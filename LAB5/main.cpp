#include <sstream> 
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include "random.h"

#define bohr0 0.0529E-9

using namespace std;


int main (int argc, char *argv[]){

	Random rnd;
	rnd.Set();

	int M = 1000000;	// Total number of throws
	int N = 100;		// Number of blocks
	int L = int(M/N);	// Number of throws in each block, please use for M a multiple of N

	double * ave = new double[N];
	double * sum_prog = new double[N];
	double * su2_prog = new double[N];
	double * err_prog = new double[N];
	double * ave2 = new double[N];
	double * sum_prog2 = new double[N];
	double * su2_prog2 = new double[N];
	double * err_prog2 = new double[N];

	double xprimo[3] = {bohr0*1., bohr0*1., bohr0*1.};
	double xenne[3] = {bohr0*0.5, bohr0*0.5, bohr0*0.5};
	double xprimo2[3] = {bohr0*1., bohr0*1., bohr0*1.};
	double xenne2[3] = {bohr0*0., bohr0*0., bohr0*2.};
	double probfrac = 0;
	double probfrac2 = 0;
	
	cout << "\nInizio programma\n\n";

// ---------------------------------- ( Parte 1 ) ---------------------------------- //

	cout << "\n\n// -------------------- ( Parte 1 ) -------------------- //\n\n";

	// Stream
	ofstream out1("file1.txt");
	ofstream out2("file2.txt");
	ofstream outrv("fileR1.txt");
	ofstream outrv2("fileR2.txt");

	// Azzero variabili
	for(int i=0; i<N; i++){
		ave[i] = 0, sum_prog[i] = 0;
		su2_prog[i] = 0, err_prog[i] = 0;
		ave2[i] = 0, sum_prog2[i] = 0;
		su2_prog2[i] = 0, err_prog2[i] = 0;
	}

	// Equilibrazione
	for(int i=0; i<1000; i++){

		// Genero x' casuale
		for(int s=0; s<3; s++){
			xprimo[s] = xenne[s] + bohr0*rnd.Rannyu(-1.22,1.22);
			xprimo2[s] = xenne2[s] + bohr0*rnd.Rannyu(-2.98,2.98);
		}
		probfrac = pow(rnd.psi110(xprimo),2) / pow(rnd.psi110(xenne),2);
		probfrac2 = pow(rnd.psi210(xprimo2),2) / pow(rnd.psi210(xenne2),2);

		// Accetto o rifiuto la mossa
		if( rnd.Rannyu() <= min(1., probfrac) ){  
			for(int s=0; s<3; s++){ xenne[s] = xprimo[s]; }
		}

		if( rnd.Rannyu() <= min(1., probfrac2) ){  
			for(int s=0; s<3; s++){ xenne2[s] = xprimo2[s]; }
		}
	}


	for(int i=0; i<N; i++){

		for(int j=0; j<L; j++){	
			// Genero x' casuale
			for(int s=0; s<3; s++){
				xprimo[s] = xenne[s] + bohr0*rnd.Rannyu(-1.22,1.22);
				xprimo2[s] = xenne2[s] + bohr0*rnd.Rannyu(-2.98,2.98);
			}
			probfrac = pow(rnd.psi110(xprimo),2) / pow(rnd.psi110(xenne),2);
			probfrac2 = pow(rnd.psi210(xprimo2),2) / pow(rnd.psi210(xenne2),2);

			// Accetto o rifiuto la mossa
			if( rnd.Rannyu() <= min(1., probfrac) ){  
				for(int s=0; s<3; s++){ xenne[s] = xprimo[s]; }
				outrv << xenne[0] << " " << xenne[1] << " " << xenne[2] << endl;
			}

			if( rnd.Rannyu() <= min(1., probfrac2) ){  
				for(int s=0; s<3; s++){ xenne2[s] = xprimo2[s]; }
				outrv2 << xenne2[0] << " " << xenne2[1] << " " << xenne2[2] << endl;
			}
			
			ave[i] += sqrt( pow(xenne[0],2) + pow(xenne[1],2) + pow(xenne[2],2) );
			ave2[i] += sqrt( pow(xenne2[0],2) + pow(xenne2[1],2) + pow(xenne2[2],2) );
		}

		// Media dei raggi in blocco
		ave[i] /= L;  
		ave2[i] /= L; 
		
		if( (i+1)%50 == 0){ cout << "Calcolo medie " << i+1 << endl; }
	
		// Block averaging
		for(int k=0; k<(i+1); k++){
  			sum_prog[i]+= ave[k];				// SUM_{j=0,i} r_j
  			sum_prog2[i]+= ave2[k];				
  			su2_prog[i]+= pow(ave[k],2);			// SUM_{j=0,i} (r_j)^2
			su2_prog2[i]+= pow(ave2[k],2);			
		}  	
		sum_prog[i]/=(i+1);					// Cumulative average
  		su2_prog[i]/=(i+1);					
		sum_prog2[i]/=(i+1);					// Cumulative average
  		su2_prog2[i]/=(i+1);					
  		err_prog[i] = rnd.error(sum_prog,su2_prog,i); 		// Statistical uncertainty
  		out1 << sum_prog[i] << " " << err_prog[i] << endl;
  		err_prog2[i] = rnd.error(sum_prog2,su2_prog2,i); 	// Statistical uncertainty
  		out2 << sum_prog2[i] << " " << err_prog2[i] << endl;
	}

	out1.close();
	out2.close();
	outrv.close();
	outrv2.close();


// ---------------------------------- ( Parte 2 ) ---------------------------------- //

	cout << "\n\n// -------------------- ( Parte 2 ) -------------------- //\n\n";

	// Stream
	out1.open("file3.txt");
	out2.open("file4.txt");
	outrv.open("fileR3.txt");
	outrv2.open("fileR4.txt");

	// Nuove starting positions
	xenne[0] = bohr0*0.5, xenne[1] = bohr0*0.5, xenne[2] = bohr0*0.5;
	xenne2[0] = bohr0*0., xenne2[1] = bohr0*0., xenne2[2] = bohr0*2.;

	// Azzero variabili
	for(int i=0; i<N; i++){
		ave[i] = 0, sum_prog[i] = 0;
		su2_prog[i] = 0, err_prog[i] = 0;
		ave2[i] = 0, sum_prog2[i] = 0;
		su2_prog2[i] = 0, err_prog2[i] = 0;
	}

	// Equilibrazione
	for(int i=0; i<1000; i++){

		// Genero x' casuale
		for(int s=0; s<3; s++){
			xprimo[s] = rnd.Gauss(xenne[s], bohr0*0.76);
			xprimo2[s] = rnd.Gauss(xenne2[s], bohr0*1.875);
		}
		probfrac = pow(rnd.psi110(xprimo),2) / pow(rnd.psi110(xenne),2);
		probfrac2 = pow(rnd.psi210(xprimo2),2) / pow(rnd.psi210(xenne2),2);

		// Accetto o rifiuto la mossa
		if( rnd.Rannyu() <= min(1., probfrac) ){  
			for(int s=0; s<3; s++){ xenne[s] = xprimo[s]; }
		}

		if( rnd.Rannyu() <= min(1., probfrac2) ){  
			for(int s=0; s<3; s++){ xenne2[s] = xprimo2[s]; }
		}
	}


	for(int i=0; i<N; i++){

		for(int j=0; j<L; j++){
			// Genero x' casuale
			for(int s=0; s<3; s++){
				xprimo[s] = rnd.Gauss(xenne[s], bohr0*0.76);
				xprimo2[s] = rnd.Gauss(xenne2[s], bohr0*1.875);
			}
			probfrac = pow(rnd.psi110(xprimo),2) / pow(rnd.psi110(xenne),2);
			probfrac2 = pow(rnd.psi210(xprimo2),2) / pow(rnd.psi210(xenne2),2);

			// Accetto o rifiuto la mossa
			if( rnd.Rannyu() <= min(1., probfrac) ){  
				for(int s=0; s<3; s++){ xenne[s] = xprimo[s]; }
				outrv << xenne[0] << " " << xenne[1] << " " << xenne[2] << endl;
			}

			if( rnd.Rannyu() <= min(1., probfrac2) ){  
				for(int s=0; s<3; s++){ xenne2[s] = xprimo2[s]; }
				outrv2 << xenne2[0] << " " << xenne2[1] << " " << xenne2[2] << endl;
			}
			
			ave[i] += sqrt( pow(xenne[0],2) + pow(xenne[1],2) + pow(xenne[2],2) );
			ave2[i] += sqrt( pow(xenne2[0],2) + pow(xenne2[1],2) + pow(xenne2[2],2) );
		}

		// Media dei raggi in blocco
		ave[i] /= L;  
		ave2[i] /= L; 
		
		if( (i+1)%50 == 0){ cout << "Calcolo medie " << i+1 << endl; }
	
		// Block averaging
		for(int k=0; k<(i+1); k++){
  			sum_prog[i]+= ave[k];				// SUM_{j=0,i} r_j
  			sum_prog2[i]+= ave2[k];				
  			su2_prog[i]+= pow(ave[k],2);			// SUM_{j=0,i} (r_j)^2
			su2_prog2[i]+= pow(ave2[k],2);			
		}  	
		sum_prog[i]/=(i+1);					// Cumulative average
  		su2_prog[i]/=(i+1);					
		sum_prog2[i]/=(i+1);					// Cumulative average
  		su2_prog2[i]/=(i+1);					
  		err_prog[i] = rnd.error(sum_prog,su2_prog,i); 		// Statistical uncertainty
  		out1 << sum_prog[i] << " " << err_prog[i] << endl;
  		err_prog2[i] = rnd.error(sum_prog2,su2_prog2,i); 	// Statistical uncertainty
  		out2 << sum_prog2[i] << " " << err_prog2[i] << endl;
	}

	out1.close();
	out2.close();
	outrv.close();
	outrv2.close();

// ------------------------------- Memory cleaning -------------------------------- //

	delete [] ave;
	delete [] sum_prog;
	delete [] su2_prog;
	delete [] err_prog;
	delete [] ave2;
	delete [] sum_prog2;
	delete [] su2_prog2;
	delete [] err_prog2;	

	return 0;
}
