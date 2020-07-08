#ifndef __ISING__
#define __ISING__

#include "random.h"

class Ising {

    private:
        Random rnd; // Random numbers

    public:
        // Parameters, observables
        int n_props, iu, ic, im, ix, ig;
        double nbins;
        double *walker;

        // Averages
        double blk_norm, accepted, attempted;
        double *blk_av, *glob_av, *glob_av2;
        double stima_u, stima_c, stima_m, stima_x, stima_g;
        double err_u, err_c, err_m, err_x, err_g;

        // Configuration
        const int m_spin = 50;
        double *s;

        // Thermodynamical state
        int nspin;
        double beta, temp, J, h;

        // Simulation
        int nstep, nblk, metro;

        // Constructors
        Ising();
        // Destructor
        ~Ising();
        // Methods
        void Input(void);
        void Reset(int);
        void Accumulate(void);
        void Averages(int);
        void Move();
        void ConfFinal(void);
        void Measure(void);
        double Boltzmann(int, int);
        int Pbc(int);
        double Error(double, double, int);

};

#endif // __Ising__