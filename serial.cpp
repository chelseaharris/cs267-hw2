#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"
extern int numBins; 
//
//  benchmarking program
//
int main( int argc, char **argv ) {    
	int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;
	int numthreads = 1; 
    //
    // argument parsing
    // 
    if( find_option( argc, argv, "-h" ) >= 0 ) {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
	bool isCheck = (find_option( argc, argv, "-no" ) == -1);
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
	
	set_size(n);
	bin_t* bins = (bin_t*) malloc(numBins * sizeof(bin_t));
	init_bins(bins);
    init_particles(n, particles );
	
	for (int i = 0; i < n; i++)
		move_and_update(particles[i]);

	binning(particles, bins, n);
    //
    //  simulate a number of time steps
    //
	int batchSz = ceil(numBins / (float)numthreads);
	int batchSzParticles = ceil(n / (float)numthreads);
    double simulation_time = read_timer();
	FOR (step, NSTEPS) {
		//#pragma omp parallel for 
		FOR (i, n) {
			particles[i].ax = 0; 
			particles[i].ay = 0;
		}
	
		//#pragma omp parallel for 
		FOR (i, numthreads)
			apply_force_bin_batch(particles, bins, i, batchSz, numBins);

		if (isCheck) {
			navg = 0;
			davg = 0.0;
			dmin = 1.0;

			FOR (i, numBins)
				get_statistics_bin(particles, bins[i], i, &dmin, &davg, &navg);
		}
    
		//#pragma omp parallel for 
		FOR (i, numthreads)
			move_and_update_batch(particles, i, batchSzParticles, n);

		binning(particles, bins, n);

        if (isCheck)
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }

          if (dmin < absmin) 
			  absmin = dmin;
		
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }

	
	simulation_time = read_timer( ) - simulation_time;

	    
  printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if (isCheck)
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
    
	if (absmin < 0.4) 
		printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    
	if (absavg < 0.8) 
		printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }

    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );    
    free( particles );
    if( fsave )
        fclose( fsave );
    free(bins);
    return 0;

}
