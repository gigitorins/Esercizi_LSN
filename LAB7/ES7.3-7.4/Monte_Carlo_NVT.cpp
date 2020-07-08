#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_NVT.h"

using namespace std;


mcNVT :: mcNVT(){

}


mcNVT :: ~mcNVT(){

	delete [] walker;
	delete [] blk_av;
	delete [] glob_av;	
	delete [] glob_av2;
	delete [] x;
	delete [] y;
	delete [] z;
}

void mcNVT :: Input(int argc, char ** argv){

	ifstream ReadInput, ReadConf;

	rnd.Set();	// Read seed for random numbers

	// Read input informations
    if(argc == 1)	ReadInput.open("input.dat");
    else			ReadInput.open(argv[1]);
	
	ReadInput >> temp;
	ReadInput >> npart;
	ReadInput >> rho;
	ReadInput >> rcut;
	ReadInput >> delta;
	ReadInput >> nblk;
	ReadInput >> nstep;
	ReadInput.close();

	// Read nblk & nstep
	if(argc > 2)	nblk = atoi(argv[2]);
    if(argc > 3)	nstep = atoi(argv[3]);
	
	beta = 1.0/temp;
	vol = (double)npart/rho;
	box = pow(vol,1.0/3.0);

	// Tail corrections for potential energy and pressure
	vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
	ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));

	cout << "Number of blocks = " << nblk << endl;
	cout << "Number of steps in one block = " << nstep << endl << endl;

	// Prepare arrays for measurements
	iv = 0;			// Potential energy
	iw = 1;			// Virial
	n_props = 2;	// Number of observables

	// Measurement of g(r)
	igofr = 2;
	nbins = 100;
	n_props += nbins;
	bin_size = (box/2.0)/(double)nbins;

	// Creazione vettori
	walker = new double[n_props];
	blk_av = new double[n_props];
	glob_av = new double[n_props];
	glob_av2 = new double[n_props];
	x = new double[npart];
	y = new double[npart];
	z = new double[npart];

	// Read initial configuration
	cout << "Read initial configuration from file config.0 " << endl << endl;
	ReadConf.open("config.0");
	for(int i=0; i<npart; ++i){
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = Pbc( x[i] * box );
		y[i] = Pbc( y[i] * box );
		z[i] = Pbc( z[i] * box );
	}
	ReadConf.close();
  
	// Evaluate potential energy and virial of the initial configuration
	Measure();

}


void mcNVT :: Move(void){
	int o;
	double p, energy_old, energy_new;
	double xold, yold, zold, xnew, ynew, znew;

	for(int i=0; i<npart; ++i){
		// Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
		o = (int)(rnd.Rannyu()*npart);

		// Old
		xold = x[o], yold = y[o], zold = z[o];
		energy_old = Boltzmann(xold,yold,zold,o);

		// New
		xnew = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
		ynew = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
		znew = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );
		energy_new = Boltzmann(xnew,ynew,znew,o);

		// Metropolis test
		p = exp(beta*(energy_old-energy_new));
		if(p >= rnd.Rannyu()){
			// Update
			x[o] = xnew, y[o] = ynew, z[o] = znew;
			accepted += 1.0;
		}
		attempted += 1.0;
	}
}

double mcNVT :: Boltzmann(double xx, double yy, double zz, int ip){
	double ene = 0.0;
	double dx, dy, dz, dr;

	for (int i=0; i<npart; ++i){
		if(i != ip){
			// Distance ip-i in pbc
			dx = Pbc(xx - x[i]);
			dy = Pbc(yy - y[i]);
			dz = Pbc(zz - z[i]);
			dr = sqrt(dx*dx + dy*dy + dz*dz);

			if(dr < rcut){
				ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
			}
		}
	}

	return 4.0*ene;
}

void mcNVT :: Measure(){

	double v = 0.0, w = 0.0, vij, wij;
	double dx, dy, dz, dr;
    double inf = 0, sup = 0; // g(r) hystogram

	// Reset the hystogram of g(r)
	for(int k=igofr; k<igofr+nbins; ++k)	walker[k] = 0.0;

	// Cycle over pairs of particles
	for(int i=0; i<npart-1; ++i){
		for (int j=i+1; j<npart; ++j){
			// Distance i-j in pbc
			dx = Pbc(x[i] - x[j]), dy = Pbc(y[i] - y[j]), dz = Pbc(z[i] - z[j]);
			dr = sqrt(dx*dx + dy*dy + dz*dz);

			// Update of the histogram of g(r)
            for(int k=igofr; k<igofr+nbins; ++k){
				inf = (k-igofr)*bin_size;
				sup = (k-igofr+1)*bin_size;
				attempted++;
                if(dr >= inf && dr < sup){
					walker[k] += 2;		// Check if dr is in the k-th interval
					accepted++; 
            	}
			}

			// Contribution to energy and virial
			if(dr < rcut){
				vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
				wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
				v += vij;	
				w += wij;
			}
		}          
	}

	walker[iv] = 4.0 * v;
	walker[iw] = 48.0 * w / 3.0;
}


void mcNVT :: Reset(int iblk){	// Reset block averages
   
	if(iblk == 1){
		for(int i=0; i<n_props; ++i)	glob_av[i] = 0, glob_av2[i] = 0;
	}

	for(int i=0; i<n_props; ++i)	blk_av[i] = 0;
	blk_norm = 0;
	attempted = 0;
	accepted = 0;
}


void mcNVT :: Accumulate(void){	// Update block averages

   for(int i=0; i<n_props; ++i)		blk_av[i] += walker[i];
   blk_norm += 1.0;
}


void mcNVT :: Averages(int iblk){	// Print results for current block
    
	double r, gdir, DeltaV;
	ofstream Gofr, Gave, Epot, Pres;
	const int wd=12;
    
	// cout << "Block number " << iblk << endl;
	// cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
	Epot.open("output.epot.0",ios::app);
	Pres.open("output.pres.0",ios::app);
	Gofr.open("output.gofr.0",ios::app);
	Gave.open("output.gave.0");
    
    stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot = Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_press = Error(glob_av[iw],glob_av2[iw],iblk);

	// Potential energy per particle
	Epot << iblk << " " << stima_pot << " " << glob_av[iv]/(double)iblk << " " << err_pot << endl;
	// Pressure
	Pres << iblk << " " <<  stima_pres << " " << glob_av[iw]/(double)iblk << " " << err_press << endl;

	// g(r)
    for(int k=igofr; k<igofr+nbins; k++){

        r = bin_size*(k-igofr);
        DeltaV = 4./3 * M_PI * r*r*r;

		if(DeltaV == 0)	gdir = 0;
		else			gdir = blk_av[k]/rho/(double)npart/DeltaV/nbins;
		
        glob_av[k] += gdir;
        glob_av2[k] += gdir*gdir;
        err_gdir = Error(glob_av[k], glob_av2[k], iblk);
        
        Gofr << " " << r << " " << gdir << " " << endl;
        Gave << " " << r << " " << glob_av[k]/(double)iblk << " " << err_gdir << endl;
    }
	
	// cout << "----------------------------" << endl << endl;

	Epot.close();
	Pres.close();
	Gofr.close();
	Gave.close();
}


void mcNVT :: ConfFinal(void){

	ofstream WriteConf("config.final");
	cout << "Print final configuration to file config.final " << endl << endl;
	for (int i=0; i<npart; ++i)		WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
	WriteConf.close();
	rnd.SaveSeed();
}


void mcNVT :: ConfXYZ(int nconf){ // Write configuration in .xyz format

  ofstream WriteXYZ("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  for(int i=0; i<npart; ++i)	WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  WriteXYZ.close();
}


double mcNVT :: Pbc(double r){	// Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}


double mcNVT :: Error(double sum, double sum2, int iblk){
    if(iblk == 1)	return 0.0;
    else			return	sqrt( abs( (sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1) ));
}