#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"

double size;
int binSize; 
int numBins; 
int* globalIds; 
int* globalNums; 
//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define cutoff2 (cutoff*cutoff)
#define min_r   (cutoff/100)
#define min_r2  (min_r * min_r)
#define dt      0.0005

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n ) {
	size = sqrt( density * n );
	binSize = (int)ceil(size / cutoff);
	numBins = binSize * binSize; 
	globalIds =  (int*) malloc(n * sizeof(int));
	globalNums = (int*)malloc(numBins * sizeof(int));
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p )
{
    srand48( time( NULL ) );
        
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);
		
        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
		p[i].ax = 0; 
		p[i].ay = 0;
		p[i].id = i;
    }
    free( shuffle );
}

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor)
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff2 )
        return;

		
	//printf("%f\n",min_r);
	
    r2 = fmax( r2, double(min_r2) );

    double r = sqrt(r2);
    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}


//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }
}


void move_and_update( particle_t &p)
{
	//
	//  slightly simplified Velocity Verlet integration
	//  conserves energy better than explicit Euler method
	//
	//printf("x = %f, y = %f, vx = %f, vy = %f, ax = %f, ay = %f\n",  p.x, p.y, p.vx, p.vy, p.ax, p.ay);
	p.vx += p.ax * dt;
	p.vy += p.ay * dt;
	p.x  += p.vx * dt;
	p.y  += p.vy * dt;

	//
	//  bounce from walls
	//
	while( p.x < 0 || p.x > size )
	{
		p.x  = p.x < 0 ? -p.x : 2*size-p.x;
		p.vx = -p.vx;
	}
	while( p.y < 0 || p.y > size )
	{
		p.y  = p.y < 0 ? -p.y : 2*size-p.y;
		p.vy = -p.vy;
	}

	p.ax = 0; 
	p.ay = 0;
	//int id = ;
	globalIds[p.id] = (int)(floor(p.x / cutoff) * binSize 
		+ floor(p.y / cutoff));
	//if (globalIds[p.id] < 0 || globalIds[p.id] >= numBins) {
		//printf("id = %d, bin id = %d x = %f, y = %f\n", p.id, globalIds[p.id], p.x, p.y);
		//exit(-1);
	//}
	
}


//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}


void binning(particle_t* _particles, bin_t* _bins, int _num) {
	#pragma omp parallel for 
	FOR (i, numBins) {
		_bins[i].num_particles = 0;
		globalNums[i] = 0; 
	}
	
	
	FOR (i, _num)
		globalNums[globalIds[i]]++;

	#pragma omp parallel for 
	FOR (i, numBins) {
		/*if (_bins[i].particle_ids)
			free(_bins[i].particle_ids);*/
		if (globalNums[i] > 0)
			_bins[i].particle_ids = (particle_t**)malloc(globalNums[i] * sizeof(particle_t*));
	}
	#pragma omp parallel for 
	FOR (i, _num) {
		int id = globalIds[i];
		_bins[id].particle_ids[_bins[id].num_particles] = _particles+i;
		_bins[id].num_particles++;
	}

}



void apply_force_bin(particle_t* _particles, bin_t& _bin, int _binId) {
	//bin_t* bin = _bins + _binId;

	FOR (i, _bin.num_particles) {
		particle_t& p = *_bin.particle_ids[i];
		FOR (k, _bin.num_neigh) {
			bin_t* new_bin = _bin.neighbors_ids[k]; 
			for(int j = 0; j < new_bin->num_particles; j++)
				apply_force(p, *new_bin->particle_ids[j]);
		}
		//move_and_update(p);
	}
}


void get_statistics( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg ) {
	double dx = neighbor.x - particle.x;
	double dy = neighbor.y - particle.y;
	double r2 = dx * dx + dy * dy;
	if( r2 > cutoff2 )
		return;

	if (r2 != 0) {
		if (r2/(cutoff2) < *dmin * (*dmin))
			*dmin = sqrt(r2)/cutoff;

		(*davg) += sqrt(r2)/cutoff;
		(*navg) ++;
	}
}

void get_statistics_bin( particle_t* _particles, bin_t& _bin, int _binId, double *dmin, double *davg, int *navg ) {
	FOR (i, _bin.num_particles) {
		FOR (k, _bin.num_neigh) {
			bin_t* new_bin = _bin.neighbors_ids[k]; 
			for(int j = 0; j < new_bin->num_particles; j++)
				get_statistics(*_bin.particle_ids[i], *new_bin->particle_ids[j], dmin, davg, navg);
		}
	}
}

void init_bins( bin_t* _bins ) {
	int dx[] = {-1, -1, -1, 0, 0, 0, 1, 1, 1};
	int dy[] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
	FOR (i, numBins) {
		_bins[i].particle_ids = 0; 
		_bins[i].num_neigh = 0; 
		_bins[i].neighbors_ids = (bin_t**) malloc(9 * sizeof(bin_t*));
		int x = i % binSize; 
		int y = (i - x) / binSize; 
		FOR (k, 9) {
			int new_x = x + dx[k]; 
			int new_y = y + dy[k];
			if (new_x >= 0 && new_y >= 0 && new_x < binSize && new_y < binSize) {
				int new_id = new_x + new_y * binSize;				
				_bins[i].neighbors_ids[_bins[i].num_neigh] = _bins + new_id; 
				_bins[i].num_neigh++;
			}
		}

	}
}



void apply_force_bin_batch( particle_t* _particles, bin_t* _bins, int _binId, int _batchSize, int _num) {
	int start = _binId * _batchSize; 
	int end = min((_binId+1) * _batchSize, _num); 
	for (int i = start; i < end; i++ )
		apply_force_bin(_particles, _bins[i], i);
}

void move_and_update_batch( particle_t* _particles, int _binId, int _batchSize, int _num ) {
	int start = _binId * _batchSize; 
	int end = min((_binId+1) * _batchSize, _num); 
	for (int i = start; i < end; i++ )
		move_and_update(_particles[i]);
}
