/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"

#define bohr0 0.0529E-9

using namespace std;



Random :: Random(){}

Random :: ~Random(){}

void Random :: Set(){
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
		input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	SaveSeed();
}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Lorentz(double mean, double width){
	double y=Rannyu();
	return width*(tan(M_PI*(y-0.5)))+mean;
}

double Random :: Exponential(double lambda){
	double y=Rannyu();
	return ((-1.)/lambda)*log(1-y);
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

void Random :: random_unif(double * r, int n){
  
  for(int i=0; i<n; i++){
  	r[i] = Rannyu();
  }
  return;
}

void Random :: random_exp(double * r, int n, double lambda){
  
  for(int i=0; i<n; i++){
  	r[i] = Exponential(lambda);
  }
  return;
}

void Random :: random_lorentz(double * r, int n, double mean, double width){
  
  for(int i=0; i<n; i++){
  	r[i] = Lorentz(mean,width);
  }
  return;
}

void Random :: random_gauss(double * r, int n, double mean, double sigma){
  
  for(int i=0; i<n; i++){
  	r[i] = Gauss(mean,sigma);
  }
  return;
}

double Random :: error(double *AV, double *AV2, int n){  //Function for statistical uncertainty estimation
	if(n==0){
	return 0;
	}else	{
	return sqrt((AV2[n] - pow(AV[n],2))/n);
 	}
}

void Random :: angolo_rand(double ** omega, int dim){

	for(int i=0; i<dim; i++){
		omega[i][1] = acos(1-2*Rannyu());
		omega[i][2] = Rannyu(0., 2.*M_PI);
	}

	return;
}

double Random :: psi110(double * vettore){
	
	double raggio = sqrt( pow(vettore[0],2) + pow(vettore[1],2) + pow(vettore[2],2) );
	return ( sqrt(1./(M_PI * pow(bohr0,3))) * exp(-raggio/bohr0) );

}
	
double Random :: psi210(double * vettore){

	double raggio = sqrt( pow(vettore[0],2) + pow(vettore[1],2) + pow(vettore[2],2) );
	return ( (1./8) * sqrt(2./(M_PI * pow(bohr0,5))) * vettore[2] * exp(-raggio/(2*bohr0)) );

}

