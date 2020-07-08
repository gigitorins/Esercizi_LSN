#include <cstdlib>   
#include <iostream>     
#include <fstream>      
#include <string>
#include <sstream>
#include <cmath>
#include <ostream>
#include <iomanip>
#include "genetic.h"

using namespace std;

int main(int argc, char ** argv){

    genetic gen1;
	gen1.Initialize(argc, argv);
    int ngens = gen1.stepnumber, wd = 12;

    cout << endl << "Best walk iniziale: " << gen1.votominimo << endl;
    ofstream out("walks.dat");
    double beta;
    int tempsteps = gen1.tsteps;
    
    for(int i=1; i<=tempsteps; i++){
        beta = 1.0*i;

        for(int step=0; step<ngens; step++){
            gen1.Move(beta);
            cout << "-- Progresso: " << (int)( ( (double)(step+i*ngens)/(ngens*tempsteps) )*100.) << "% \r";
        }

        out << setw(wd) << beta << setw(wd) << gen1.votiPop << endl;
    }

    cout << endl << setw(40) << "--- RISULTATI FINALI ---" << setw(40) << endl;
    cout << setw(wd) << "Temperatura" << setw(wd) << "Cammino" << endl;
    cout << setw(wd) << 1./beta << setw(wd) << gen1.votominimo << endl;
    cout << "Rate di accettazione: " << (double)gen1.accepted/gen1.attempted << endl;
    out.close();

	gen1.PrintConfig();
	gen1.DeleteMemory();

	return 0;
}
