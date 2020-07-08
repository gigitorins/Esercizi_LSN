#include <stdlib.h>   
#include <iostream>     
#include <fstream>      
#include <string>
#include <sstream>
#include <cmath>        
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(){ 

	Ising is1d;
	is1d.Input();	// Inizialization
	int numblocks = is1d.nblk;
	int numsteps = is1d.nstep;

	for(int iblk=1; iblk <= numblocks; ++iblk){	// Simulation

		is1d.Reset(iblk);	// Reset block averages

		for(int istep=1; istep <= numsteps; ++istep){
			is1d.Move();
			is1d.Measure();
			is1d.Accumulate();	// Update block averages
		}

		is1d.Averages(iblk);	// Print results for current block
	}

	is1d.ConfFinal();	// Write final configuration

	return 0;
}
