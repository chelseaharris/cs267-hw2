#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"


#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

extern double size;
int numBins;
int numRows;
int numCols;
//
//  benchmarking program
//

#define FOR(i,n) for( int i=0; i<n; i++ )

struct bin_t {
	int row;
	int col;
	int num;
	int* ids;
};

void binning(bin_t* _bins, particle_t* _particles, int _n) {
	FOR (i, numBins)
    _bins[i].num = 0;
	FOR (i, _n) {
		int id = floor(_particles[i].y/cutoff) * numCols +  floor(_particles[i].x/cutoff);
		_bins[id].ids[_bins[id].num++] = i;
	}
}


void init_bins(bin_t*& _bins, int _n) {
	numRows = ceil(size/cutoff);
	numCols = ceil(size/cutoff);
	numBins = numRows * numCols;
	_bins = (bin_t*) malloc(numBins * sizeof(bin_t));
	FOR (i, numBins)
    _bins[i].ids = (int*) malloc(_n*sizeof(int));
	FOR (i, numRows) {
		FOR (j, numCols) {
			int id = i * numCols + j;
			_bins[id].row = i;
			_bins[id].col = j;
			_bins[id].num = 0;
		}
	}
}
void init_local_bins(bin_t*& _bins, int _nlocal_bins, int _n, int _offset) {
    //	numRows = ceil(size/cutoff);
	numCols = ceil(size/cutoff);
    //	numBins = numRows * numCols;
	_bins = (bin_t*) malloc(_nlocal_bins * sizeof(bin_t));
	FOR (i, _nlocal_bins) {
        //        _bins[i].ids = (int*) malloc(_n*sizeof(int));
        //        _bins[i].ids[0] = 2123;
        int id = _offset + i;
        int row = floor(id/numCols);
        int col = id - row * numCols;
        _bins[i].row = row;
        _bins[i].col = col;
        _bins[i].num = 0;
    }
    
}
void free_bins(bin_t*& _bins) {
	FOR (i, numBins)
    free(_bins[i].ids);
    
	free(_bins);
}


void simple_apply_force(particle_t& _particle, particle_t& _neighbor) {
	double dx = _neighbor.x - _particle.x;
	double dy = _neighbor.y - _particle.y;
	double r2 = dx * dx + dy * dy;
	if( r2 > cutoff*cutoff )
		return;
	
	r2 = fmax( r2, min_r*min_r );
	double r = sqrt( r2 );
	//
	//  very simple short-range repulsive force
	//
	double coef = ( 1 - cutoff / r ) / r2 / mass;
	_particle.ax += coef * dx;
	_particle.ay += coef * dy;
}

void get_statistics( particle_t &particle, particle_t &neighbor ,
					double *dmin, double *davg, int *navg) {
	double dx = neighbor.x - particle.x;
	double dy = neighbor.y - particle.y;
	double r2 = dx * dx + dy * dy;
	if( r2 > cutoff*cutoff )
		return;
	if (r2 != 0)
	{
		if (r2/(cutoff*cutoff) < *dmin * (*dmin))
			*dmin = sqrt(r2)/cutoff;
		(*davg) += sqrt(r2)/cutoff;
		(*navg) ++;
	}
}

void get_statistics_bin(particle_t* _particles, bin_t* _bins, int _id, double *dmin, double *davg, int *navg) {
	bin_t* bin = _bins + _id;
	int dx[9] = {-1, -1, -1, 0, 0, 0, 1, 1, 1};
	int dy[9] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
	FOR (i, bin->num) {
		FOR (k, 9) {
			int row = bin->row + dx[k];
			int col = bin->col + dy[k];
			if (row >= 0 && col >= 0 && row < numRows && col < numCols ) {
				bin_t* bin_neigh = _bins + row * numCols + col;
				FOR (j, bin_neigh->num) {
					get_statistics(_particles[bin->ids[i]],_particles[bin_neigh->ids[j]], dmin, davg, navg);
				}
			}
		}
	}
}

void simple_apply_force_bin(particle_t* _particles, bin_t* _bins, bin_t& bin) {
	//bin_t* bin = _bins + _id;
	int dx[9] = {-1, -1, -1, 0, 0, 0, 1, 1, 1};
	int dy[9] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
	FOR (i, bin.num) {
		FOR (k, 9) {
			int row = bin.row + dx[k];
			int col = bin.col + dy[k];
			if (row >= 0 && col >= 0 && row < numRows && col < numCols) {
				bin_t* bin_neigh = _bins + row * numCols + col;
				FOR (j, bin_neigh->num)
                simple_apply_force(_particles[bin.ids[i]],_particles[bin_neigh->ids[j]]);
			}
		}
	}
}
//
//  benchmarking program
//
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{
    
    //
    //  process command line parameters
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );
    set_size( n );
    char *savename = read_string( argc, argv, "-o", NULL );
    
    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    //
    //  allocate generic resources
    //
    //    printf("rank = %d, point 1\n", rank);
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    
    bin_t *bins;
    if (rank == 0)
        init_bins(bins, n);
    //    printf("rank = %d, point 2\n", rank);
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 7, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    MPI_Datatype BIN;
    MPI_Type_contiguous(7, MPI_DOUBLE, &BIN);
    MPI_Type_commit(&BIN);
    
    //    printf("rank = %d, point 3\n", rank);
    //
    //  set up the data partitioning across processors
    //
    int particle_per_proc = (n + n_proc - 1) / n_proc;
    int bin_per_proc = (numBins + n_proc - 1)/n_proc;
    
    int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    int *bin_offsets = (int*) malloc((n_proc + 1 ) * sizeof(int) );
    for( int i = 0; i < n_proc+1; i++ ) {
        partition_offsets[i] = min( i * particle_per_proc, n );
        bin_offsets[i] = min(i * bin_per_proc, numBins);
    }
    
    //    printf("rank = %d, point 4\n", rank);
    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    int *bin_sizes = (int*) malloc(n_proc * sizeof(int));
    for( int i = 0; i < n_proc; i++ ) {
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
        bin_sizes[i] = bin_offsets[i+1] - bin_offsets[i];
    }
    
    //    printf("rank = %d, point 5\n", rank);
    //
    //  allocate storage for local partition
    //
    int nlocal = partition_sizes[rank];
    //int nlocal_bin = bin_sizes[rank];
    //int local_offset = bin_offsets[rank];
    particle_t *local = (particle_t*) malloc( nlocal * sizeof(particle_t) );
    //bin_t *local_bins;
    //init_local_bins(local_bins, nlocal_bin, n, local_offset);
    
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    
    if( rank == 0 )
        init_particles( n, particles );
    // printf("rank = %d, point a6, testing %d\n", rank, local_bins[0].ids[0]);
    MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
    //MPI_Scatterv(bins,bin_sizes,bin_offsets,BIN,local_bins,nlocal_bin,BIN,0,MPI_COMM_WORLD);
    
    //    MPI_Allgatherv( local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );
    
    //    printf("rank = %d, point 7\n", rank);
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        //        if (rank == 0)
        //        binning(bins, particles, n);
        //            printf("rank = %d, point 7.3\n", rank);
        //
        //  collect all global data locally (not good idea to do)
        //
        //printf("rank = %d, point a7, testing %d\n", rank, local_bins[0].ids[0]);
        //printf("rank = %d, step = %d, point 7.6\n", rank, step);
        MPI_Allgatherv( local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );
        //            printf("rank = %d, point 7.6\n", rank);
        //printf("rank = %d, point a8, testing %d\n", rank, local_bins[0].ids[0]);
        //printf("rank = %d, step = %d, point 8\n", rank, step);
        //        local_binning(local_bins, particles, nlocal_bin, n, local_offset);
        //                binning(bins, particles,  n);
        //printf("rank = %d, point a9, testing %d\n",rank, local_bins[0].ids[0]);
        //
        //  save current step if necessary (slightly different semantics than in other codes)
        
        //
        //
        //  compute all forces
        //
        //        printf("rank = %d, step =%d, point 9\n", rank, step);
        for( int i = 0; i < nlocal; i++ )
        {
            local[i].ax = local[i].ay = 0;
            for (int j = 0; j < n; j++ )
                simple_apply_force( local[i], particles[j]);
        }
        //printf("rank = %d, point a10, testing %d\n",local_bins[0].ids[0]);
        if (rank == 0) {
            for (int i = 0; i < numBins; i++) {
                //            int id = local_offset + i;
                //            simple_apply_force_bin(particles, bins, bins[i]);
            }
        }
        //        printf("rank = %d, step = %d, point 10\n", rank, step);
        //printf("rank = %d, point a11, testing %d\n", rank, local_bins[0].ids[0]);
        //
        //  move particles
        //
        for( int i = 0; i < nlocal; i++ )
            move( local[i] );
        
        //        printf("rank = %d, step = %d, point 11\n", rank, step);
    }
    simulation_time = read_timer( ) - simulation_time;
    
    if (rank == 0) {
        printf( "n = %d, simulation time = %g seconds", n, simulation_time);
        
        printf("\n");
        
    }
    
    //
    //  release resources
    free( partition_offsets );
    free( partition_sizes );
    free( local );
    free( particles );
    free_bins(bins);
    //free_bins(local_bins);
    //free(bin_offsets);
    //free(bin_sizes);
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
