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

  int M = int(1E6);		// Total number of throws
  int N = 100;			// Number of blocks
  int L = int(M/N);		// Number of throws in each block, please use for M a multiple of N

  double * r = new double[M];
  rnd.random_unif(r,M);

  double * ave = new double[N];
  double * sum_prog = new double[N];
  double * su2_prog = new double[N];
  double * err_prog = new double[N];
  int k;

// ---------------------------------- ( ES 1.1 ) ---------------------------------- //

  ofstream out1("file1.txt");
  for(int i=0; i<N; i++){

  	for(int j=0; j<L; j++){
  		k = j+i*L;
		ave[i] += r[k]/L;
	}

  	for(int k=0; k<(i+1); k++){
  		sum_prog[i]+= ave[k];				// SUM_{j=0,i} r_j
  		su2_prog[i]+= pow(ave[k],2);		// SUM_{j=0,i} (r_j)^2
	}  	
	sum_prog[i]/=(i+1);					// Cumulative average
  	su2_prog[i]/=(i+1);					// Cumulative square average
  	err_prog[i] = rnd.error(sum_prog,su2_prog,i); 		// Statistical uncertainty
  	out1 << sum_prog[i] << " " << err_prog[i] << endl;
  }

  out1.close();

// ---------------------------------- ( ES 1.2 ) ---------------------------------- //

  ofstream out2("file2.txt");
  for(int i=0; i<N; i++){
	ave[i] = 0;
	sum_prog[i] = 0;
  	su2_prog[i] = 0;
	err_prog[i] = 0;
  }

  for(int i=0; i<N; i++){

  	for(int j=0; j<L; j++){
  		k = j+i*L;
  		ave[i] += (pow((r[k]-0.5),2))/L;
	}

  	for(int j=0; j<(i+1); j++){
  		sum_prog[i]+= ave[j];			// SUM_{j=0,i} r_j
  		su2_prog[i]+= pow(ave[j],2);		// SUM_{j=0,i} (r_j)^2
	}  	
	sum_prog[i]/=(i+1);				// Cumulative average
  	su2_prog[i]/=(i+1);				// Cumulative square average
  	err_prog[i] = rnd.error(sum_prog,su2_prog,i); 	// Statistical uncertainty
  	out2 << sum_prog[i] << " " << err_prog[i] << endl;
  }

  out2.close();

// ---------------------------------- ( ES 1.3 ) ---------------------------------- //

	M = 100;			// Numero di sottointervalli in cui dividiamo [0,1]
	int n = 1E6; 		// Numero di numeri pseudo-random generati
	int n_throw = 1E4; 	// Numero di numeri pseudo-random usati per calcolare il chi^2 una sola volta
	int j = 0;
	int *count = new int[M]; // Vettore che conta il numero di eventi in ogni sottointervallo
	double *chi_quadro = new double[n/n_throw];
	double * r2 = new double[n];
	
	rnd.random_unif(r2,n);
	ofstream out3("file3.txt");
	
	for(int counter=0; counter<n/n_throw; counter++){
		for(int i=0;i<M;i++){		// Pongo a zero il conteggio del numero di eventi in ogni sottointervallo
			count[i] = 0;
		}

		for(int i=counter*n_throw; i<counter*n_throw + n_throw; i++){
			j = 0;
			while(r2[i] < double(j)/M || r2[i] >= double(j+1)/M){	// Trovo l'indice del sottointervallo in cui è caduto il mio numero random
				j++;					     
			}						    
			count[j] += 1;
		}

		// Calcolo chi^2
		chi_quadro[counter] = 0;
		for(int i=0; i<M; i++){
			chi_quadro[counter] += pow((count[i]-n_throw/M),2)/(n_throw/M);	
		}
		out3 << chi_quadro[counter] << endl;
	}

  out3.close();

// ------------------------------- Memory cleaning -------------------------------- //

  delete [] r;
  delete [] r2;
  delete [] count;
  delete [] chi_quadro;
  delete [] ave;
  delete [] sum_prog;
  delete [] su2_prog;
  delete [] err_prog;

  return 0;
}
