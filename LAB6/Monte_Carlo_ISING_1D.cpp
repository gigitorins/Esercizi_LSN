#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;


Ising :: Ising(){


}


Ising :: ~Ising(){

	delete [] walker;
	delete [] blk_av;
	delete [] glob_av;	
	delete [] glob_av2;
	delete [] s;
}


void Ising :: Input(void){

	ifstream ReadInput;
	rnd.Set();	// Read seed for random numbers
  
	// Read input informations
	ReadInput.open("input.dat");
	ReadInput >> temp;
	beta = 1.0/temp;
	ReadInput >> nspin;
	ReadInput >> J;
	ReadInput >> h;
	ReadInput >> metro; // if=1 Metropolis else Gibbs
	ReadInput >> nblk;
	ReadInput >> nstep;

	ReadInput.close();

	// Prepare arrays for measurements
	iu = 0;			// Energy
	ic = 1;			// Heat capacity
	im = 2;			// Magnetization
	ix = 3;			// Magnetic susceptibility
	n_props = 4;	// Number of observables

	walker = new double[n_props];
	blk_av = new double[n_props];
	glob_av = new double[n_props];
	glob_av2 = new double[n_props];
	s = new double[m_spin];

	// Initial configuration
	ifstream ReadConf("Files/config.final");	// If config.final already exists it is loaded as initial configuration
	if(ReadConf.is_open()){
			for(int i=0; i<nspin; ++i)	ReadConf >> s[i];
	}else	{
			for(int i=0; i<nspin; ++i){
				if(rnd.Rannyu() >= 0.5) s[i] = 1;
				else s[i] = -1;
			}
	}
  
	ReadConf.close();
	Measure();	// Evaluate energy etc. of the initial configuration
}


void Ising :: Move(){
  
	int o;
	double p, energy_old, energy_new, sm, energydiff;
	double s_new;

	for(int i=0; i<nspin; ++i){
		// Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
		o = (int)(rnd.Rannyu()*nspin);
		sm = s[o];

		if(metro==1){ 	// Metropolis
			energy_old = Boltzmann(sm,o);
			energy_new = Boltzmann(-sm,o);
			energydiff = energy_new - energy_old;
			p = exp(-beta*energydiff);
			if( rnd.Rannyu() < min(1., p) ){
				s[o] = -sm;
				accepted++;
			}
			attempted++;
		}else	{	// Gibbs
            attempted++;
            // Select a random spin
            if(rnd.Rannyu() >= 0.5) s_new = 1;
            else s_new = -1;
            
            p = exp(-beta*Boltzmann(s_new, o))/(exp(-beta*Boltzmann(1,o)) + exp(-beta*Boltzmann(-1,o)));
            
            if(rnd.Rannyu() < p){
                s[o] = s_new;
                accepted++;
            }
        }
	}

}


double Ising :: Boltzmann(int sm, int ip){
	return -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
}


void Ising :: Measure(){
	int bin;
	double u = 0.0, m = 0.0;

	// Cycle over spins
	for (int i=0; i<nspin; ++i){
		u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
		m += s[i];
	}

	walker[iu] = u;
	walker[ic] = u*u;
	walker[im] = m;
	walker[ix] = m*m;
}


void Ising :: Reset(int iblk){		// Reset block averages
   
	if(iblk == 1){
		for(int i=0; i<n_props; ++i){	glob_av[i] = 0, glob_av2[i] = 0;	}
	}

	for(int i=0; i<n_props; ++i)   blk_av[i] = 0;
	blk_norm = 0;
	attempted = 0;
	accepted = 0;
}


void Ising :: Accumulate(void){		// Update block averages

	for(int i=0; i<n_props; ++i)   blk_av[i] += walker[i];
	blk_norm += 1.0;

}


void Ising :: Averages(int iblk){	// Print results for current block
    
	ofstream Ene, Heat, Mag, Chi;
    
	// Energy
	Ene.open("Files/output.ene.0",ios::app);
	stima_u = blk_av[iu]/blk_norm/(double)nspin;
	glob_av[iu]  += stima_u;
	glob_av2[iu] += stima_u*stima_u;
	err_u = Error(glob_av[iu],glob_av2[iu],iblk);
	Ene << iblk <<  " " << stima_u << " " << glob_av[iu]/(double)iblk << " " << err_u << endl;
	Ene.close();

	// Heat
	Heat.open("Files/output.heat.0",ios::app);
	stima_c = beta*beta*(blk_av[ic]/blk_norm - (blk_av[iu]/blk_norm)*(blk_av[iu]/blk_norm) )/(double)nspin; // Heat capacity
	glob_av[ic]  += stima_c;
	glob_av2[ic] += stima_c*stima_c;
	err_c = Error(glob_av[ic],glob_av2[ic],iblk);
	Heat << iblk <<  " " << stima_c << " " << glob_av[ic]/(double)iblk << " " << err_c << endl;
	Heat.close();

	// Magnetization
	Mag.open("Files/output.mag.0",ios::app);
	stima_m = blk_av[im]/blk_norm/(double)nspin;
	glob_av[im]  += stima_m;
	glob_av2[im] += stima_m*stima_m;
	err_m = Error(glob_av[im],glob_av2[im],iblk);
	Mag << iblk <<  " " << stima_m << " " << glob_av[im]/(double)iblk << " " << err_m << endl;
	Mag.close();

	// Susceptibility
	Chi.open("Files/output.chi.0",ios::app);
	stima_x = beta*blk_av[ix]/blk_norm/(double)nspin; 
	glob_av[ix]  += stima_x;
	glob_av2[ix] += stima_x*stima_x;
	err_x = Error(glob_av[ix],glob_av2[ix],iblk);
	Chi << iblk <<  " " << stima_x << " " << glob_av[ix]/(double)iblk << " " << err_x << endl;
	Chi.close();

}


void Ising :: ConfFinal(void){

	ofstream WriteConf;
	WriteConf.open("Files/config.final");
	for(int i=0; i<nspin; ++i)   WriteConf << s[i] << endl;
	WriteConf.close();
	rnd.SaveSeed();

}


int Ising :: Pbc(int i){	// Algorithm for periodic boundary conditions

	if(i >= nspin) i = i - nspin;
	else if(i < 0) i = i + nspin;

	return i;
}


double Ising :: Error(double sum, double sum2, int iblk){	// Function for statistical uncertainty estimation

	if(iblk==1) return 0;
	else		return sqrt( abs( (sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1) ) );
}
