#include <stdlib.h>   
#include <iostream>     
#include <fstream>      
#include <string>
#include <sstream>
#include <cmath>
#include <ostream>
#include <iomanip>
#include "Monte_Carlo_NVT.h"

using namespace std;

int main(int argc, char ** argv){

    mcNVT mcnvt; 
    mcnvt.Input(argc, argv);  // Inizialization
	int numblocks = mcnvt.nblk;
	int numsteps = mcnvt.nstep;

    for(int iblk = 1; iblk <= numblocks; ++iblk){  // Simulation

        mcnvt.Reset(iblk);  // Reset block averages
        for(int istep=1; istep <= numsteps; ++istep){
            mcnvt.Move();
            mcnvt.Measure();
            mcnvt.Accumulate(); // Update block averages
            cout << "-- Progress: " << int(double(istep+(iblk-1)*numsteps)/double(numsteps*numblocks)*10000)/(double)100<< "% \r";
        }

        mcnvt.Averages(iblk);   // Print results for current block
    }

    mcnvt.ConfFinal(); // Write final configuration

  return 0;
}