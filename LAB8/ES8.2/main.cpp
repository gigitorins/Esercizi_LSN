#include <stdlib.h>   
#include <iostream>     
#include <fstream>      
#include <string>
#include <sstream>
#include <cmath>
#include <ostream>
#include <iomanip>
#include "VCM.h"

using namespace std;

int main(int argc, char ** argv){

    VCM vcm1; 
    vcm1.Input(argc, argv);  // Inizialization
	int numblocks = vcm1.n_blok;
	int numsteps = vcm1.n_step;

    for(int iblk = 1; iblk <= numblocks; ++iblk){  // Simulation

        vcm1.Reset(iblk);  // Reset block averages
        for(int istep=1; istep <= numsteps; ++istep){
            vcm1.Move();
            vcm1.Measure();
            vcm1.Accumulate(); // Update block averages
            cout << "-- Progress: " << int(double(istep+(iblk-1)*numsteps)/double(numsteps*numblocks)*10000)/(double)100<< "% \r";
        }

        vcm1.Averages(iblk);   // Print results for current block
    }

  return 0;
}