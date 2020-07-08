#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <string>
#include <sstream>
#include <cmath>        
#include "MolDyn_NVE.h"

using namespace std;

int main(int argc, char ** argv){

		MolDyn mdyn;
		mdyn.Input(argc, argv);			//Inizialization
		int M = mdyn.nstep;		// Number of time-steps

		for(int istep=1; istep <= M; ++istep){
			mdyn.Move();		 		// Move particles with Verlet algorithm
			if(istep%mdyn.iprint == 0) cout << "Number of time-steps: " << istep << endl;
			if(istep%10 == 0){
				mdyn.Measure(istep);		// Properties measurement
			}	
			if(istep == (M-1)){
				mdyn.ConfFinalOld();		// Write actual configuration in XYZ format
			}
		}

		mdyn.ConfFinal();			// Write final configuration to restart

	return 0;
}
