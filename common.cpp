#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <vector>
#include "common.h"

double size;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
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
void set_size( int n )
{
    size = sqrt( density * n );
}

double get_size()
{
  return size;
}

double get_cutoff()
{
  return cutoff;
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
    }
    free( shuffle );
}

//
// Binning functions
//

// function to bin particles (square bins)
bins_t bin_particles( particle_t* particles, const int n ) 
{

  bins_t particle_bins;

  double bin_wid = get_cutoff(); // ideal bin width
  double grid_size = get_size();
  int num_bins_side = floor( grid_size/bin_wid ); // number of bins in one direction (some grid unbinned)
  bin_wid = grid_size/num_bins_side; // adjust bin width so that bins fill the grid
  int num_bins = num_bins_side*num_bins_side; // total number of bins

  particle_bins.num_bins = num_bins;
  particle_bins.bin_wid = bin_wid;

  // offsets to make a 3x3 box of bins
  // goes from lower left corner to upper right corner, row-wise
  int shifts[9] = { -num_bins_side-1, -num_bins_side, -num_bins_side+1, -1, 0, 1, num_bins_side-1, num_bins_side, num_bins_side+1 };

  particle_bins.shiftlist = shifts;

  /* //for testing the binning:
  printf("Grid size: %e\n",grid_size);
  printf("Cutoff: %e\n",get_cutoff());
  printf("Bin width: %e\n",bin_wid);
  */

  BinnedParticles binned_particles(num_bins);

  // bin the particles
  particle_t* end_ptr = particles + n;
  particle_t* p_ptr = particles;
  do {
    int bin_x = (int)floor( (*p_ptr).x / bin_wid ); // bin number in x direction
    int bin_y = (int)floor( (*p_ptr).y / bin_wid ); // bin number in y direction
    int bin_i = num_bins_side*bin_y + bin_x; //linear bin index;
    binned_particles[bin_i].push_back(p_ptr) ; 
    p_ptr++;
  } while (p_ptr != end_ptr);

  particle_bins.binned_parts = binned_particles;

  return particle_bins;
}


//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor )//, double *dmin, double *davg, int *navg)
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
        return;
    /* 
     * Move this part to get_stats()
    if (r2 != 0)
        {
	   if (r2/(cutoff*cutoff) < *dmin * (*dmin))
	      *dmin = sqrt(r2)/cutoff;
           (*davg) += sqrt(r2)/cutoff;
           (*navg) ++;
        }
    */
		
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );
 	
    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

// apply force to the particles in a bin; only look at 
// adjacent bins
/*void apply_force_in_bin( bins_t &part_bins, int i_bin )
{
    

}*/

// separate out the part that gets statistics, since these are not always wanted
void get_stats( particle_t &particle, particle_t &neighbor, double *dmin, double *davg, int *navg )
{
  double dx = neighbor.x - particle.x;
  double dy = neighbor.y - particle.y;
  double r  = sqrt( dx*dx + dy*dy );
  double frac_r = r/cutoff;
  if (r!=0)
    {
      if (frac_r < *dmin) 
	*dmin = frac_r;
      (*davg) += frac_r;
      (*navg) ++;
    }
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


