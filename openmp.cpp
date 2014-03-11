#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"


#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

extern double size; 
int numBins; 
int numRows;
int numCols; 
int dx[9] = {-1, -1, -1, 0, 0, 0, 1, 1, 1};
int dy[9] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};

#define FOR(i,n) for( int i=0; i<n; i++ )

// added function
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
		int id = floor(_particles[i].x/cutoff) * numCols +  floor(_particles[i].y/cutoff);
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
	//  very simple short-range repulsive force
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

void simple_apply_force_bin(particle_t* _particles, bin_t* _bins, int _id) {
	bin_t* bin = _bins + _id; 

	FOR (i, bin->num) {
		_particles[bin->ids[i]].ax = 0;
		_particles[bin->ids[i]].ay = 0;
		FOR (k, 9) {
			int row = bin->row + dx[k];
			int col = bin->col + dy[k];
			if (row >= 0 && col >= 0 && row < numRows && col < numCols) {
				bin_t* bin_neigh = _bins + row * numCols + col;
				FOR (j, bin_neigh->num)
					simple_apply_force(_particles[bin->ids[i]],_particles[bin_neigh->ids[j]]);
			}
		}
		//move(_particles[bin
	}
}



//
//  benchmarking program
//
int main( int argc, char **argv )
{   
    int navg,nabsavg=0,numthreads; 
    double dmin, absmin=1.0,davg,absavg=0.0;
	
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" ); 
        printf( "-no turns off all correctness checks and particle output\n");   
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
	
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;      

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
	bin_t* bins; 
	init_bins(bins, n);
	//binning(bins, particles, n);
	int p = read_int( argc, argv, "-p", 1 );
	omp_set_num_threads(p);
	bool isNo = ( find_option( argc, argv, "-no" ) == -1 );
	omp_lock_t* locks = (omp_lock_t*) malloc(numBins * sizeof(omp_lock_t) );
	FOR (i, numBins)
		omp_init_lock(locks+i);
	//omp_init_lock(locks);
    //
    //  simulate a number of time steps
    //
	double simulation_time = read_timer();

    #pragma omp parallel
    {
    numthreads = omp_get_num_threads();
    for( int step = 0; step < 1000; step++ )
    {
        navg = 0;
        davg = 0.0;
		dmin = 1.0;
        //
        //  compute all forces
		#pragma omp for
		FOR (i, numRows) {
			FOR (j, numCols) {
				int id = i * numCols + j; 
				bins[id].num = 0; 
			}
		}
		#pragma omp for
		FOR (i, n) {
			int id = floor(particles[i].y/cutoff) * numCols +  floor(particles[i].x/cutoff);
			omp_set_lock(locks+id);			
			bins[id].ids[bins[id].num++] = i;
			omp_unset_lock(locks+id);		
		}


        #pragma omp for
		FOR (i, numBins)   // include ax = 0; ay = 0; 
			simple_apply_force_bin(particles, bins, i);

        //
        //  move particles
        //
        #pragma omp for
		FOR (i, n)
            move( particles[i] );

		
		if (isNo) {
			#pragma omp master
			FOR (i, numBins)
				get_statistics_bin(particles, bins, i, &dmin, &davg, &navg);
		}
		//#pragma omp master
		//binning(bins, particles, n);
	


        if( find_option( argc, argv, "-no" ) == -1 ) 
        {
          //
          //  compute statistical data
          //
          #pragma omp master
          if (navg) { 
            absavg += davg/navg;
            nabsavg++;
          }

          #pragma omp critical
	  if (dmin < absmin) absmin = dmin; 
		
          //
          //  save if necessary
          //
          #pragma omp master
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
}
   	simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");
    
    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );

    free( particles );
    if( fsave )
        fclose( fsave );
    free_bins(bins);
	FOR (i, numBins) 
		omp_destroy_lock(locks+i);
    return 0;
}
