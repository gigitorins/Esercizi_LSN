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

	int M = 100000;			// Total number of throws
	int N = 250;			// Number of blocks
	int L = int(M/N);		// Number of throws in each block, please use for M a multiple of N
	ofstream out1("file1.txt");	// Stream
	ofstream out2("file2.txt");	
	double * r = new double[M];
	rnd.random_unif(r,M);

	double * ave = new double[N];
	double * ave2 = new double[N];
	double * sum_prog = new double[N];
	double * su2_prog = new double[N];
	double * err_prog = new double[N];
	double sum;
	int k;

// ---------------------------------- ( ES 2.1 parte 1 ) ---------------------------------- //

	for(int i=0; i<N; i++){
		sum = 0;
		for(int j=0; j<L; j++){
			k = j+i*L;
			sum+=(M_PI/2)*cos((M_PI*r[k])/2);
		}
		ave[i] = sum/L;
		ave2[i] = pow(ave[i],2);
	}
  
	for(int i=0; i<N; i++){
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

// ---------------------------------- ( ES 2.1 parte 2 ) ---------------------------------- //

	double x;

	for(int i=0; i<N; i++){
		sum = 0;
		for(int j=0; j<L; j++){
			x = 1 - sqrt(1-rnd.Rannyu());
			sum += (M_PI/2)*cos(x*M_PI/2)/(2*(1-x));
		}
		ave[i] = sum/L;
		ave2[i] = pow(ave[i],2);
	}
  

	for(int i=0; i<N; i++){
		sum_prog[i] = 0;
		su2_prog[i] = 0;	
		err_prog[i] = 0;

		for(int j=0; j<(i+1); j++){
			sum_prog[i]+= ave[j];				// SUM_{j=0,i} r_j
			su2_prog[i]+= ave2[j];				// SUM_{j=0,i} (r_j)^2
		}  	
		sum_prog[i]/=(i+1);					// Cumulative average
		su2_prog[i]/=(i+1);					// Cumulative square average
		err_prog[i] = rnd.error(sum_prog,su2_prog,i); 		// Statistical uncertainty
		out2 << sum_prog[i] << " " << err_prog[i] << endl;
	}

	out2.close();

// ------------------------------- Memory cleaning -------------------------------- //

  delete [] r;
  delete [] ave;
  delete [] ave2;
  delete [] sum_prog;
  delete [] su2_prog;
  delete [] err_prog;

  return 0;
}
