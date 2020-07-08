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

	int M = 1E5;				// Numero totale di lanci
	int N = 100;				// Numero di blocchi
	int L = int(M/N);			// Numero di lanci in ogni blocco (USARE PER M UN MULTIPLO DI N)
	int sort;				// Variabile d'appoggio per la prob. lungo x/y/z
	double a = 1.;				// Passo del RW
	int n_passi = 101;			// Numero di passi del RW

	double ** ave = new double * [N];	// Creo matrice per medie delle posizioni
	for(int i=0; i<N; i++){ ave[i] = new double[n_passi]; }
	double ** ave2 = new double * [N];	// Creo matrice per medie del quadrato delle posizioni
	for (int i = 0; i < N; i++){ ave2[i] = new double[n_passi]; }
	double * sum_prog = new double[n_passi];
	double * su2_prog = new double[n_passi];
	double * err_prog = new double[n_passi];
	double pos[3] = {0.,0.,0.};
	double x = 0, y = 0, z = 0;		// Posizione iniziale (origine)
	double * dist = new double[n_passi];	// Vettore delle distanze

	ofstream out1("rw3dD.txt");		// Streams
	ofstream out2("rw3dC.txt");
	ofstream out3("rw3dpD.txt");
	ofstream out4("rw3dpC.txt");

// ------------------------------------ Random Walk discreto ------------------------------------ //

	// Produco le medie su ogni blocco per i passi del RW 
	for(int g=0; g<N; g++){									// Ciclo sul numero di blocchi
		for(int ik=0; ik<n_passi; ik++){ dist[ik]=0; }		// Azzero vettore distanze per ogni blocco

		for(int j=0; j<L; j++){					// **************** INIZIO BLOCCO: ciclo sul numero di lanci L per ogni blocco

			for (int i = 0; i < n_passi; i++){	// ***** INIZIO RW: simulo 1 random walk, salvo le distanze percorse su dist[i]
				sort = int(rnd.Rannyu(0.,3.));							// Generandoli di seguito ottengo sempre RW diversi
				if(i==0){
					for(int ggg=0; ggg<3; ggg++){	pos[ggg] = 0;	}								// Ogni RW parte da 0
					dist[i]+= pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2];						// Contatore delle distanze a 0
					if(j==0 && g==0){ out1 << pos[0] << " " << pos[1] << " " << pos[2] << endl; }	// Salvo su file posizione del primo passo (esempio)
				}else	{
					pos[sort] += (2 * rnd.Bernoulli(0.5) - 1) * a;
					dist[i]+= pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2];						// Contatore delle distanze al passo i-esimo
					if(j==0 && g==0){ out1 << pos[0] << " " << pos[1] << " " << pos[2] << endl; }  	// Salvo su file posizione del primo RW (esempio)
				}
			}									// ***** FINE RW:
		}										// **************** FINE BLOCCO: ho le somme di posizioni per ogni passo i del RW sul g-esimo blocco

		for(int kk=0; kk<n_passi; kk++){		// Ciclo sui passi del RW alla fine di ogni ciclo L sul g-esimo blocco
			ave[g][kk] = dist[kk]/L;				// Assegno media della posizione al passo i=kk sul g-esimo blocco
			ave2[g][kk] = pow(ave[g][kk],2);
		}
	}
 

	// Produco le medie finali per i passi del RW 
	for(int i=0; i<n_passi; i++){				// Ciclo sul numero di passi
		for(int j=0; j<N; j++){					// Ciclo sul numero di blocchi
			sum_prog[i]+= ave[j][i];			// Somme dell'i-esimo passo 
			su2_prog[i]+= ave2[j][i];			
		}  	
		sum_prog[i]/=N;							// Media dell'i-esimo passo su tutti i blocchi
		su2_prog[i]/=N;					
		err_prog[i] = sqrt((su2_prog[i] - pow(sum_prog[i],2))/(i+1));	// Errore sull'i-esimo passo
		out3 << sum_prog[i] << " " << err_prog[i] << endl;		  
	}

	out1.close();
	out3.close();


// ------------------------------------ Random Walk continuo ------------------------------------ //
	
	int dim_omega=10000000;
	x = 0, y = 0, z = 0;								// Azzero le variabili posizione
	double ** omega = new double * [dim_omega];			// Creo matrice per variabili angolari di 10^6 = 10^4 colonne * 100 righe
	for(int i = 0; i < dim_omega; i++){ omega[i] = new double[2]; }
	rnd.angolo_rand(omega, dim_omega);					// Assegno a ogni riga una coppia di angoli casuali theta e phi
	int cont = 0;										// Contatore

	// Produco le medie su ogni blocco per i passi del RW 
	for(int g=0; g<N; g++){							// Ciclo sul numero di blocchi
		for(int ik=0; ik<n_passi; ik++)	dist[ik]=0;	// Azzero vettore distanze per ogni blocco

		for(int j=0; j<L; j++){					// **************** INIZIO BLOCCO: ciclo sul numero di lanci L per ogni blocco
			
			for (int i = 0; i < n_passi; i++){				// ***** INIZIO RW: simulo 1 random walk, salvo le distanze percorse su dist[i]
				cont=i+j*n_passi+g*L;						// Contatore che si sposta su omega
				if(i==0){
					x=0, y=0, z=0;													// Ogni RW parte da 0
					dist[i]+= x*x + y*y + z*z;								// Contatore delle distanze a 0
					if(j==0 && g==0)	out2 << x << " " << y << " " << z << endl;	// Salvo su file posizione del primo passo (esempio)
				}else	{
				x+= a*(sin(omega[cont][1]))*(cos(omega[cont][2]));			// Calcolo delle posizioni ai passi successivi
				y+= a*(sin(omega[cont][1]))*(sin(omega[cont][2]));
				z+= a*(cos(omega[cont][1]));
				dist[i]+= x*x + y*y + z*z;								// Contatore delle distanze al passo i-esimo
				if(j==0 && g==0)	out2 << x << " " << y << " " << z << endl;	// Salvo su file posizione del primo RW (esempio)
				}
			}							// ***** FINE RW:
		}							// **************** FINE BLOCCO: ho le somme di posizioni per ogni passo i del RW sul g-esimo blocco

		for(int kk=0; kk<n_passi; kk++){			// Ciclo sui passi del RW alla fine di ogni ciclo L sul g-esimo blocco
			ave[g][kk] = dist[kk]/L;				// Assegno media della posizione al passo i=kk sul g-esimo blocco
			ave2[g][kk] = pow(ave[g][kk],2);
		}
	}

	// Azzero variabili per somme cumulative ed errore
	for(int i=0; i<n_passi; i++){	sum_prog[i] = 0, su2_prog[i] = 0, err_prog[i] = 0;	}

	// Produco le medie finali per i passi del RW 
	for(int i=0; i<n_passi; i++){				// Ciclo sul numero di passi
		for(int j=0; j<N; j++){					// Ciclo sul numero di blocchi
			sum_prog[i]+= ave[j][i];			// Somme dell'i-esimo passo 
			su2_prog[i]+= ave2[j][i];			
		}  	
		sum_prog[i]/=N;							// Media dell'i-esimo passo su tutti i blocchi
		su2_prog[i]/=N;					
		err_prog[i] = sqrt((su2_prog[i] -pow(sum_prog[i],2))/(i+1));	// Errore sull'i-esimo passo
		out4 << sum_prog[i] << " " << err_prog[i] << endl;		  
	}

	out2.close();
	out4.close();


// ------------------------------------ Memory cleaning ------------------------------------ //

	for(int i = 0; i < N; i++){ delete [] ave[i]; }
	delete [] ave;

	for(int i = 0; i < N; i++){ delete [] ave2[i]; }
	delete [] ave2;
	
	delete [] dist;
	delete [] sum_prog;
	delete [] su2_prog;
 	delete [] err_prog;

	for(int i = 0; i < dim_omega; i++){ delete [] omega[i]; }
	delete [] omega;

	return 0;
}
