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
    int ngens = gen1.numgenerations;
    
    cout << endl << "Best walk iniziale: " << gen1.votominimoassoluto << endl;

    for(int i=0; i<ngens; i++){
        gen1.Generazione(i);
        cout << "-- Progress: " << int(double(i)/double(ngens)*10000)/(double)100<< "% \r";
    }
    
    cout << "Best walk trovato: " << gen1.votominimoassoluto << endl;

	gen1.PrintConfig();
	gen1.DeleteMemory();

	return 0;
}