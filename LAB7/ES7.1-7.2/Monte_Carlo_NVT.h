#ifndef __NVT__
#define __NVT__

#include "random.h"

class mcNVT {

    private:
        Random rnd; // Random numbers

    public:
    	// Parameters, observables
    	int n_props, iv, iw, igofr;
    	double vtail, ptail, bin_size, nbins, sd;
		double *walker;

		// Averages
		double blk_norm, accepted, attempted;
		double *blk_av, *glob_av, *glob_av2;
		double stima_pot, stima_pres, err_pot, err_press, err_gdir;

		// Configuration
		const int m_part = 108;
		double *x, *y, *z;

		// Thermodynamical state
		int npart;
		double beta, temp, vol, rho, box, rcut;

		// Simulation
		int nstep, nblk;
		double delta;

		// Pigreco
		const double pi = 3.1415927;

		// Constructors
		mcNVT();
		// Destructor
		~mcNVT();
		// Methods
		void Input(int, char **);
		void Reset(int);
		void Accumulate(void);
		void Averages(int);
		void Move(void);
		void ConfFinal(void);
		void ConfXYZ(int);
		void Measure(void);
		double Boltzmann(double, double, double, int);
		double Pbc(double);
		double Error(double,double,int);
};

#endif