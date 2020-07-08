#ifndef __NVT__
#define __NVT__

#include "random.h"

class VCM {

    private:
        Random rnd; // Random numbers

    public:

		// Variabili
		double mu = 0.76;
		double sigma = 0.6;
		double delta = 2.2;

		int accepted = 0, attempted = 0;
		int blk_norm = 0;
		double blk_avg = 0., walker = 0.;
		int n_blok = 100, n_step = 10000;
		double glob_av = 0., glob_av2 = 0.;
		double stima_histo = 0., err_histo = 0.;
	    double stima_h = 0., err_h = 0.;

		int n_bins = 100;
		int *histogram, *blk_histo;
		double *glob_histo, *glob_histo2;

		// Limits for the x in the histogram
		double inf = -3, sup = 3;
		double bin_length;

		// Number of total points in the histogram
		double tot = 0;

		// Starting position
		double x = 0.;

		// Constructors
		VCM();
		// Destructor
		~VCM();
		// Methods
		void Input(int, char **);
		void Reset(int);
		void Accumulate();
		void Averages(int);
		void Move();
		void Measure();
		double Error(double,double,int);
		double Wave_function(double);
		double Potential(double);
		double Second_derivative(double);
};

#endif