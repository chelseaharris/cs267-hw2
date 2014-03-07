#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__
#define FOR(i,n) for( int i=0; i<n; i++ )
inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;

//
// particle data structure
//
typedef struct 
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
} particle_t;

typedef struct{
	int num_particles;
	int num_neigh; 
	int* neighbors_ids; 
	int* particle_ids;
} bin_t;




//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
void init_particles( int n, particle_t *p );
void apply_force( particle_t &particle, particle_t &neighbor);
void move( particle_t &p );


// added functions
void move_and_update( particle_t &p, int _id);
void init_bins(bin_t* _bins);
void binning(particle_t* _particles, bin_t* _bins, int _num);
void apply_force_bin(particle_t* _particles,  bin_t* _bins, int _binId);
//void apply_force_bin_openmp(particle_t* _particles,  bin_t* _bins, int _binId);
//void apply_force_bin_batch(particle_t* _particles,  bin_t* _bins, int _binId, int _batchSize, int _num);
void get_statistics( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void get_statistics_bin(particle_t* _particles,  bin_t* _bins, int _binId, double *dmin, double *davg, int *navg);
//void move_and_binning( particle_t &p);


//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );




//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
