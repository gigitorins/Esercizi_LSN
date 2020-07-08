#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <cstdlib>
#include "random.h"

using namespace std;

int main (int argc, char *argv[]){

	Random rnd;
	rnd.Set();

	int M = 1E5;				// Total number of throws
	int N = 250;				// Number of blocks
	int L = int(M/N);			// Number of throws in each block, please use for M a multiple of N
	ofstream outC("fileC.txt");
	ofstream outP("fileP.txt");	

	double * aveC = new double[N];
	double * aveP = new double[N];
	double * ave2C = new double[N];
	double * ave2P = new double[N];
	double * sum_progC = new double[N];
	double * sum_progP = new double[N];
	double * su2_progC = new double[N];
	double * su2_progP = new double[N];
	double * err_progC = new double[N];
	double * err_progP = new double[N];
	double sumC, sumP;

	double S0 = 100;	// Asset price at t=0
	double T = 1.;		// Delivery time
	double K = 100;		// Strike price 
	double r = 0.1;		// Risk-free interest rate
	double sigma = 0.25;	// Volatility
	double S_cur = 0.;	// Dummy variable


// ---------------------------------- ( ES 2.1 parte 1 ) ---------------------------------- //

	for(int i=0; i<N; i++){
		sumC = 0, sumP = 0;
		for(int j=0; j<L; j++){
 			S_cur = S0 * exp(T*(r-0.5*sigma*sigma)+sqrt(T)*sigma*rnd.Gauss(0.,1.));		// Spot price via the Brownian motion final distribution
			sumC+= std::max(S_cur - K, 0.0), sumP+= std::max(K - S_cur, 0.0);		// Sum of the option pay-offs
		}
		aveC[i] = (sumC*exp(-r*T))/L, aveP[i] = (sumP*exp(-r*T))/L;				// Average of the pay-off sums, discounted the risk-free rate from the price			
		ave2C[i] = pow(aveC[i],2), ave2P[i] = pow(aveP[i],2);
	}
  

	for(int i=0; i<N; i++){
		for(int j=0; j<(i+1); j++){
			sum_progC[i]+= aveC[j], sum_progP[i]+= aveP[j];			
			su2_progC[i]+= ave2C[j], su2_progP[i]+= ave2P[j];	
		}  
		sum_progC[i]/=(i+1), sum_progP[i]/=(i+1);				// Cumulative average					
		su2_progC[i]/=(i+1), su2_progP[i]/=(i+1);				// Cumulative square average
		err_progC[i] = rnd.error(sum_progC,su2_progC,i);		// Statistical uncertainty
		err_progP[i] = rnd.error(sum_progP,su2_progP,i); 	
		outC << sum_progC[i] << " " << err_progC[i] << endl;
		outP << sum_progP[i] << " " << err_progP[i] << endl;
	}

	outC.close(), outP.close();

// ---------------------------------- ( ES 2.1 parte 2 ) ---------------------------------- //

	ofstream outC2("fileC2.txt"), outP2("fileP2.txt");	// Stream

	// Azzero variabili per medie, somme cumulative ed errore
	for(int i=0; i<N; i++){
		aveC[i] = 0., aveP[i] = 0., ave2C[i] = 0., ave2P[i] = 0.;
		sum_progC[i] = 0., sum_progP[i] = 0., su2_progC[i] = 0., su2_progP[i] = 0.;
		err_progC[i] = 0., err_progP[i] = 0.;
	}
	
	for(int i=0; i<N; i++){
		sumC = 0, sumP = 0;
		for(int j=0; j<L; j++){
			S_cur = S0;
			for(int l=1; l<=100; l++){  S_cur*= exp((r-0.5*sigma*sigma)/100 + sqrt(1./100)*sigma*rnd.Gauss(0.,1.));  }	// Spot price via the Brownian motion final distribution
			sumC+= std::max(S_cur - K, 0.0),  sumP+= std::max(K - S_cur, 0.0);						// Sum of the option pay-offs
		}
		aveC[i] = (sumC*exp(-r*T))/L, aveP[i] = (sumP*exp(-r*T))/L;				// Average of the pay-off sums, discounted the risk-free rate from the price			
		ave2C[i] = pow(aveC[i],2), ave2P[i] = pow(aveP[i],2);
	}
  
	for(int i=0; i<N; i++){
		for(int j=0; j<(i+1); j++){
			sum_progC[i]+= aveC[j], sum_progP[i]+= aveP[j];			// SUM_{j=0,i} r_j				
			su2_progC[i]+= ave2C[j], su2_progP[i]+= ave2P[j];		// SUM_{j=0,i} (r_j)^2
		}  	
		sum_progC[i]/=(i+1), sum_progP[i]/=(i+1);		// Cumulative average					
		su2_progC[i]/=(i+1), su2_progP[i]/=(i+1);		// Cumulative square average
		err_progC[i] = rnd.error(sum_progC,su2_progC,i); 	// Statistical uncertainty
		err_progP[i] = rnd.error(sum_progP,su2_progP,i); 	
		outC2 << sum_progC[i] << " " << err_progC[i] << endl;
		outP2 << sum_progP[i] << " " << err_progP[i] << endl;
	}

	outC2.close(), outP2.close();

// ------------------------------- Memory cleaning -------------------------------- //

  delete [] aveC;
  delete [] aveP;
  delete [] ave2C;
  delete [] ave2P;
  delete [] sum_progC;
  delete [] sum_progP;
  delete [] su2_progC;
  delete [] su2_progP;
  delete [] err_progC;
  delete [] err_progP;

  return 0;
}
