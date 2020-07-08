#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <algorithm>    // std::random_shuffle
#include "genetic.h"


using namespace std;



genetic ::  genetic(){	}



genetic :: ~genetic(){	}



// Inizializzazione
void genetic :: Initialize(int argc, char ** argv){

	// Set seed for random numbers
	generator = new TRandom3(0);

	// Initialize vectors
	citylist = new int[citynumber];
	for(int i = 0; i < citynumber; ++i)		citylist[i] = i+1;

	// Cities position
    citypositions = new double * [citynumber];
    for(int i=0; i<citynumber; i++)  citypositions[i] = new double[2];

	// Square or circular configuration
    if(argc > 1)		func_type = atoi(argv[1]);
	if(func_type == 0){
		cout << "Square configuration" << endl << endl;
		ReadConfig(func_type);
		// RandomSquare();
	}else	{
		cout << "Circular configuration" << endl << endl;
		ReadConfig(func_type);
		// RandomCirc();
	}

	individuo = new int[citynumber];
	nuovoindividuo = new int [citynumber];
	bestind = new int[citynumber];
	
	// Creo individuo iniziale
	int nswaps = (int)(citynumber*generator->Rndm());
	int elem1 = 0, elem2 = 0, kk = 0;
	for(int i = 0; i < citynumber; ++i)		individuo[i] = citylist[i];

	while(kk<nswaps){
		elem1 = (int)(generator->Uniform(1.,citynumber));
		elem2 = (int)(generator->Uniform(1.,citynumber));

		if(elem1 != elem2){
			Swap(individuo, elem1, elem2);
			kk++;
		}
	}

	// Nuovo individuo & bestind = individuo iniziale
	for(int k = 0; k < citynumber; ++k){
		nuovoindividuo[k] = individuo[k];
		bestind[k] = individuo[k];
	}

	votiPop = L1(individuo,citynumber);
	votominimo = votiPop;

}



// ---------------------------------------------------------- PRINT INDIVIDUI ---------------------------------------------------------- //

void genetic :: PrintInd(int * individuo){
	
	for(int k=0; k<citynumber; k++){
		if(k == 0)		cout << "[ " << individuo[k] << ", ";
		else {
			if(k == citynumber-1)	cout << individuo[k] << " ]";
			else					cout << individuo[k] << ", ";
		}
	}
}


// -------------------------------------------------------- METROPOLIS -------------------------------------------------------- //

void genetic :: Move(double bb){

	attempted++;
    double beta = bb;
    double energy_old, energy_new;
    double r, A;

	for(int j=0; j<citynumber; ++j)		nuovoindividuo[j] = individuo[j];

	energy_old = L1(individuo,citynumber);
	Mutazioni();
    energy_new = L1(nuovoindividuo,citynumber);
    
    if((energy_new-energy_old) < 0){

		accepted++;
		for(int j=0; j<citynumber; ++j)		individuo[j] = nuovoindividuo[j];

    }else	{
        r = generator->Rndm();
        A = exp(-beta*(energy_new-energy_old));
        
        if(r < min(1.,A)){
			accepted++;
			for(int j=0; j<citynumber; ++j)		individuo[j] = nuovoindividuo[j];
        }
    }

	votiPop = L1(individuo,citynumber);

	if(votiPop < votominimo){
		votominimo = votiPop;
		for(int i=0; i<citynumber; ++i)		bestind[i] = individuo[i];
	}

}


// -------------------------------------------------------------- SORTING -------------------------------------------------------------- //

// Swaps
void genetic :: Swap(int *individuo, int pos1, int pos2){
	int	appo = individuo[pos1];
	individuo[pos1] = individuo[pos2];
	individuo[pos2] = appo;
}

void genetic :: Swap(double *individuo, int pos1, int pos2){
	double appo = individuo[pos1];
	individuo[pos1] = individuo[pos2];
	individuo[pos2] = appo;
}



// ----------------------------------------------------------- CONFIGURAZIONI ----------------------------------------------------------- //

// Configurazione quadrata
void genetic :: RandomSquare(){

	ofstream WriteConfig("config_s.dat");
	double appo1 = 0, appo2 = 0;

	for(int i = 0; i < citynumber; ++i){
		for(int k = 0; k < 2; ++k)	citypositions[i][k] = generator->Uniform(0.,squareside);
		appo1 = citypositions[i][0], appo2 = citypositions[i][1];
		WriteConfig << appo1 << " " << appo2 << endl;
	}

	WriteConfig.close();
}

// Configurazione circolare
void genetic :: RandomCirc(){

	ofstream WriteConfig("config_c.dat");
	double appo1 = 0, appo2 = 0;
	double theta = 0;

	for(int i = 0; i < citynumber; ++i){
		theta = generator->Uniform(0.,2.*M_PI);
		citypositions[i][0] = circleradius + circleradius*cos(theta);
		citypositions[i][1] = circleradius + circleradius*sin(theta);
		appo1 = citypositions[i][0], appo2 = citypositions[i][1];
		WriteConfig << appo1 << " " << appo2 << endl;
	}

	WriteConfig.close();
}



// -------------------------------------------------------------- FITNESS -------------------------------------------------------------- //

// Funzione L1 di una certa configurazione di #citynum città
double genetic :: L1(int * configuration, int citynum){

	double dvec[2], path = 0;
	int index = 0, indexprimo = 0;

	for(int i=0; i<citynum; ++i){
		if(i == (citynum-1)){
			index = configuration[i] - 1;
			dvec[0] = abs( citypositions[index][0] - citypositions[0][0]),	dvec[1] = abs( citypositions[index][1] - citypositions[0][1]);
			path += sqrt(dvec[0]*dvec[0] + dvec[1]*dvec[1]);
		}else	{
			index = configuration[i] - 1,	indexprimo = configuration[i+1] -1;
			dvec[0] = abs( citypositions[index][0] - citypositions[indexprimo][0]),	dvec[1] = abs( citypositions[index][1] - citypositions[indexprimo][1]);
			path += sqrt(dvec[0]*dvec[0] + dvec[1]*dvec[1]);
		}
	}

	return path;
}


// Funzione L2 di una certa configurazione di #citynum città
double genetic :: L2(int * configuration, int citynum){

	double dvec[2], path = 0;
	int index = 0, indexprimo = 0;

	for(int i=0; i<citynum; ++i){
		if(i == (citynum-1)){
			index = configuration[i] - 1;
			dvec[0] = abs( citypositions[index][0] - citypositions[0][0]),	dvec[1] = abs( citypositions[index][1] - citypositions[0][1]);
			path += dvec[0]*dvec[0] + dvec[1]*dvec[1];
		}else	{
			index = configuration[i] - 1,	indexprimo = configuration[i+1] -1;
			dvec[0] = abs( citypositions[index][0] - citypositions[indexprimo][0]),	dvec[1] = abs( citypositions[index][1] - citypositions[indexprimo][1]);
			path += dvec[0]*dvec[0] + dvec[1]*dvec[1];
		}
	}

	return path;
}



// ------------------------------------------------------------- MUTAZIONI ------------------------------------------------------------- //

// Applica mutazioni di tipo casuale alla popolazione
void genetic :: Mutazioni(){

	int mutationnum = (int)(generator->Uniform(0.,5.));

	switch(mutationnum){
		case 1:{	if(generator->Rndm() <= pm1)	Mutazione1(nuovoindividuo);
					break;
		}
		case 2:{	if(generator->Rndm() <= pm2)	Mutazione2(nuovoindividuo);
					break;
		}
		case 3:{	if(generator->Rndm() <= pm3)	Mutazione3(nuovoindividuo);
					break;
		}
		case 4:{	if(generator->Rndm() <= pm4)	Mutazione4(nuovoindividuo);
					break;
		}
	}
}


// 1) --> Swap di due parti casuali scelte nell'individuo (esclusa la prima)
void genetic :: Mutazione1(int *individuo){

	int elem1 = 0, elem2 = 0, kk = 0;

	while(kk == 0){
		elem1 = (int)(generator->Uniform(1.,citynumber));
		elem2 = (int)(generator->Uniform(1.,citynumber));

		if(elem1 != elem2){
			Swap(individuo, elem1, elem2);
			kk++;
		}
	}
}


// 2) --> Shifta di n posizioni #m città contigue
void genetic :: Mutazione2(int *individuo){

	int posinit = (int)(generator->Uniform(1.,citynumber-1));				// Posizione iniziale da cui shiftare
	int ncitt = (int)(generator->Uniform(1.,citynumber-posinit));			// Numero m di città contigue	
	int shift = (int)(generator->Uniform(1.,citynumber+1-ncitt-posinit));	// Valore dello shift

	// Copio lo shift su individuo e lo completo
	for(int i=posinit+ncitt-1; i>=posinit; i--){
		for(int j=0; j<shift; j++){
			Swap(individuo,i+j,i+j+1);
		}
	}
}


// 3) --> Permuta di #m città con altre #m città contigue
void genetic :: Mutazione3(int *individuo){
	
	int pos1 = 0;				// Posizione iniziale da cui contare m
	if(citynumber % 2 == 0)		pos1 = (int)(generator->Uniform(1.,citynumber/2));
	else						pos1 = (int)(generator->Uniform(1.,(double)citynumber/2));
		
	int m = (int)(generator->Uniform(1.,citynumber/2-pos1));
	for(int i=0; i<m; i++)		Swap(individuo, pos1+i, pos1+i+m);
}


// 4) --> Inverte l'ordine di m città contigue
void genetic :: Mutazione4(int *individuo){
	
	int iniziale = 2, finale = 1;

	while(iniziale > finale || iniziale == finale){
		iniziale = (int)(generator->Uniform(1.,citynumber-1));		// Posizione iniziale
		finale = (int)(generator->Uniform(1.,citynumber-1));		// Posizione finale	
	}

	for(int i=0; i<finale-iniziale; i++){
		if(iniziale+i == finale-i)	break;
		else						Swap(individuo,iniziale+i,finale-i);
	}
}


// -------------------------------------------------------------- ALTRO -------------------------------------------------------------- //

// Legge la configurazione iniziale
void genetic :: ReadConfig(int tipo){

	ifstream ReadConfig;
	double appo1 = 0, appo2 = 0;

	if(tipo == 0)	ReadConfig.open("config_s.dat");
	else			ReadConfig.open("config_c.dat");
	
	for(int i = 0; i < citynumber; ++i){
		ReadConfig >> appo1 >> appo2;
		citypositions[i][0] = appo1,	citypositions[i][1] = appo2;
	}

	ReadConfig.close();
}



// Stampa la configurazione
void genetic :: PrintConfig(){

	ofstream WriteConfig("solution.dat");
	for(int i = 0; i < citynumber; ++i)		WriteConfig << bestind[i] << endl;
	WriteConfig.close();
}


// Dealloca la memoria usata
void genetic ::DeleteMemory(){

	delete generator;
    delete [] citylist;
	for(int i = 0; i < citynumber; i++)	delete [] citypositions[i];
	delete [] citypositions;
	delete [] individuo;
	delete [] nuovoindividuo;
}