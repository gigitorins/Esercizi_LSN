#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include "random.h"

using namespace std;


int main (int argc, char *argv[]){

	Random rnd;		// Inizializzazione classe random
	rnd.Set();

	int M = 1E6;			// Numero di lanci
	int N[4] = {1,2,10,100};		
	double sum_u = 0, sum_e = 0, sum_l = 0;
	ofstream out1;

	for(int l=0; l<4; l++){

		out1.open("file" + to_string(N[l]) + ".txt");

		for(int i=0; i<M; i++){

			sum_u=0, sum_e=0, sum_l=0;

			for(int j=0; j<N[l]; j++){
				sum_u += rnd.Rannyu();
				sum_e += rnd.Exponential(1.0);
				sum_l += rnd.Lorentz(0.0, 1.0);
			}

			sum_u /= N[l];
			sum_e /= N[l];
			sum_l /= N[l];

			out1 << sum_u << " " << sum_e << " " << sum_l << endl;
		}

		out1.close();
	}

  return 0;
}
