#ifndef __MolDyn__
#define __MolDyn__

class MolDyn{

private:
	// Parameters, observables
	double stima_pot, stima_kin, stima_etot, stima_temp;

	// Configuration
	int m_part=108;
	double *x, *y, *z, *xold, *yold, *zold, *vx, *vy, *vz;
	double *ave_pot, *ave_kin, *ave_etot, *ave_temp;
	double *sum_pot, *sum_kin, *sum_etot, *sum_temp;
	double *sum2_pot, *sum2_kin, *sum2_etot, *sum2_temp;
	double *err_pot, *err_kin, *err_etot, *err_temp;
	int counter=0;

protected:

public:
	// Simulation
	int nstep, nblock=100, iprint, seed;
	int nbins = 100;
	int n_props = nbins;
	double bin_size;
	int attempted, accepted;
	int stepinblock;
	double delta;
	int restart;
	int numsim;
	double *walker;

	// Thermodynamical state
	int npart;
	double energy,temp,vol,rho,box,rcut;

	// Constructors
	MolDyn();
	// Destructor
	~MolDyn();

	// Methods
	void Input(int, char **);
	void InputConfig(int);
	void ChangeInput0();
	void Move();
	void ConfFinal();
	void ConfFinalOld();
	void ConfXYZ(int);
	void Measure(int);
	double Force(int, int);
	double Pbc(double);
};


#endif // __MolDyn__
