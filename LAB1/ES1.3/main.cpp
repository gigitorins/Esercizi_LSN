#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include "random.h"

using namespace std;


int main (int argc, char *argv[]){

	Random rnd;
	rnd.Set();

	int M = int(1E5);		// Numero totale di lanci
	int N = 100;			// Numero di blocchi
	int L = int(M/N);		// Numero di lanci in ogni blocco (USARE PER M UN MULTIPLO DI N)

	if(argc > 1)	M = atoi(argv[1]);
	if(argc > 2)	M = atoi(argv[1]), N = atoi(argv[2]);
	cout << endl << "M = " << M << "   N = " << N << endl << endl;
	L = int(M/N);

	double x, theta, y;			// Coordinata x del centro dell'ago, angolo tra la verticale e l'ago
	double l = 10., d = 30.;	// Lunghezza dell'ago, distanza tra due linee
	double p = 0;				// Conteggi positivi

	ofstream out1("pi.txt");
	ofstream out2("pos.dat");
	double * ave = new double[N];
	double * ave2 = new double[N];
	double * sum_prog = new double[N];
	double * su2_prog = new double[N];
	double * err_prog = new double[N];
	
	for(int i=0; i<N; i++){

		if((i+1)%100 == 0) cout << endl << "Blocco " << i+1;
		p = 0;
		for(int j=0; j<L; j++){
			y = rnd.Rannyu(0.,d);
			theta = rnd.Rannyu(0., 3.14); 
			x = rnd.Rannyu(0.,d);
			if(x-(l/2)*sin(theta)<=0 || x+(l/2)*sin(theta)>=d)		p++;
			if((i+1)%100 == 0){
				out2 << x << " " << y << " " << theta << endl;
			}
		}
		ave[i] = (2*l*L)/(p*d);
		ave2[i] = pow(ave[i],2);

		for(int j=0; j<(i+1); j++){
			sum_prog[i]+= ave[j];				// SUM_{j=0,i} r_j
			su2_prog[i]+= ave2[j];				// SUM_{j=0,i} (r_j)^2
		}  	
		sum_prog[i]/=(i+1);									// Cumulative average
		su2_prog[i]/=(i+1);									// Cumulative square average
		err_prog[i] = rnd.error(sum_prog,su2_prog,i); 		// Statistical uncertainty
		out1 << sum_prog[i] << " " << err_prog[i] << endl;
	}

	out1.close();
	out2.close();
	cout << endl;

// ------------------------------- Memory cleaning -------------------------------- //

	delete [] ave;
	delete [] ave2;
	delete [] sum_prog;
	delete [] su2_prog;
	delete [] err_prog;

  return 0;
}
