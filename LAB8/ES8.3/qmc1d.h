#include <iostream>
#include <fstream>
#include <cmath>
#include <TRandom3.h>


#ifndef __qmc1d__
#define __qmc1d__


class QMC1D {

	public:
	//  Definition of the physical constants. The constants have been put to 1 in order to work with more comfortable numbers
	double hbar = 1;
	double boltzmann = 1;
	double mass = 1;

	//********** GLOBAL VARIABLES ************/
	double sigma = 0.61;
	double mu = 0.77;
	int func_type = 0;
	double lambda;			// lambda is hbar*hbar/2m, a constant widely used in Path Integral Theory.
	double dtau, alpha;		// dtau is the timestep, the "small imaginary-time" by which the total propagation time is divided by the Path Integral.
	int PIGS;				// PIGS = flag variable used to determine whether you are running a zero temperature or a finite temperature Path Integral.
	TRandom3* generator;

	// The following declarations are the variables used by QMC1D
	int timeslices, brownianBridgeReconstructions, brownianBridgeAttempts, brownianMotionReconstructions;
	int MCSTEPS, equilibration, blocks, histogram_bins;
	int timeslices_averages_start, timeslices_averages_end;
	double temperature, imaginaryTimePropagation, delta_variational, delta_translation;
	double histogram_start, histogram_end;
	int acceptedTranslations, acceptedVariational, acceptedBB, acceptedBM;
	int totalTranslations, totalVariational, totalBB, totalBM;

    double* positions;
    double* potential_energy;
    double* potential_energy_accumulator;
    double* potential_energy_square_accumulator;                                                                                                             
    double* kinetic_energy;
    double* kinetic_energy_accumulator;
	double* kinetic_energy_square_accumulator;                                                                                                            
	double* positions_histogram;
	double* positions_histogram_accumulator;
	double* positions_histogram_square_accumulator;

	// Constructor
	QMC1D();
	// Destructor
	~QMC1D();
	// Methods
	void readInput(int argc, char ** argv);		// reads input from the file "input.dat"
	void deleteMemory();	// handles the dynamic allocation of memory
	void initialize();		// initializes the variables
	void consoleOutput();	// writes the output on the screen                                                                                                       
	double potential_density_matrix(double val, double val_next); // potential_density_matrix returns only the potential part of the correlation between two adjacent timeslices.
	double u_prime(double x, int m);
	double u_sec(double x, int m);

	double external_potential(double);			// External potential definition
	double external_potential_prime(double);	// First derivative
	double external_potential_second(double);	// Second derivative 
                                                                                                                               
	void translation();								// performs a rigid translation
	void brownianBridge();							// reconstructs a segment of the polymer with a free particle propagation. 
	void brownianMotion(int);						// reconstructs a segment at the extremities of the polymer with a free particle propagation.                                                                                                                 
	double variationalWaveFunction(double);			// variationalWaveFunction is the variational wave function that is projected in a PIGS simulation.
	double variationalWaveFunction_second(double);	// As for the potential, you have to specify its first and second derivative for the evaluation of the kinetic local energy.
	double variationalLocalEnergy(double val);

	int index_mask(int);		// index_mask is just a compatibility function that takes into account whether the polymer is open (PIGS) or closed in periodic boundary contitions (PIMC-ring polymer).
	void upgradeAverages();		// at every MCSTEP accumulates the estimators values.
	void upgradeHistogram();	// fills the histogram of positions foreach MCSTEP
	void endBlock();			// finalizes the averages at the end of each block
	double kineticEstimator(double,double);  // evaluates the kinetic energy along the polymer

	// Called at the end of the simulation, basically they average over each block and evaluate the error on the block average.
	void finalizePotentialEstimator();
	void finalizeKineticEstimator();
	void finalizeHistogram();
};

#endif
