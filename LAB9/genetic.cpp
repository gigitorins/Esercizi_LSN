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
	bestind = new int[citynumber];
	for(int i = 0; i < citynumber; ++i){
		citylist[i] = i+1;
		bestind[i] = i+1;
	}

	// Cities position
    citypositions = new double * [citynumber];
    for(int i=0; i<citynumber; i++)  citypositions[i] = new double[2];

	// Best walks and averages
	bestwalks = new double[numgenerations+1];
	averages = new double[numgenerations+1];
	for(int i=0; i<numgenerations+1; i++){
		bestwalks[i] = 0, averages[i] = 0;
	}

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

	votiPop = new double[popnumber];
	for(int i = 0; i < popnumber; ++i)		votiPop[i] = 0;

	// Initial population
	PopIniziale();

	// Print initial population
	//PrintPop();

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

void genetic :: PrintPop(){
	for(int i=0; i<popnumber; ++i){
        cout << i+1 << " = ";
		PrintInd(population[i]);
		cout << endl;
	}
}


// ------------------------------------------------------ GENERAZIONI E SELEZIONE ------------------------------------------------------ //

// Genera popolazione iniziale
void genetic :: PopIniziale(){

	double avg = 0;
	population = new int * [popnumber];
	newpopulation = new int * [popnumber];
	
	// Creo popolazione iniziale e newpopulation (serve per Crossover e Mutazioni)
	for(int i=0; i<popnumber; ++i){
		GeneraInd(i);
		newpopulation[i] = new int[citynumber];
		for(int k = 0; k < citynumber; ++k)		newpopulation[i][k] = 0;
	}

	// Calcolo voti iniziali, li ordino in modo crescente
	for(int i=0; i<popnumber; ++i)		votiPop[i] = L1(population[i],citynumber);
	QuickSort(votiPop,0,popnumber-1);

	// Calcolo medie della migliore metà della popolazione (la prima metà)
	for(int i=0; i<popnumber/2; i++)	avg += votiPop[i];
	averages[0] = avg/((int)popnumber/2);

	votominimo = votiPop[0];
	bestwalks[0] = votominimo;
	votominimoassoluto = votominimo;
}


// Genera individuo random in posizione #numero nella popolazione iniziale
void genetic :: GeneraInd(int numero){

	int nswaps = (int)(citynumber*generator->Rndm());
	int elem1 = 0, elem2 = 0, kk = 0;
	population[numero] = new int[citynumber];
	for(int i = 0; i < citynumber; ++i)		population[numero][i] = citylist[i];

	while(kk<nswaps){
		elem1 = (int)(generator->Uniform(1.,citynumber));
		elem2 = (int)(generator->Uniform(1.,citynumber));

		if(elem1 != elem2){
			Swap(population[numero], elem1, elem2);
			kk++;
		}
	}
}


// Generazione di una popolazione a partire da quella precedente
void genetic :: Generazione(int num){

	double avg = 0;
	double copiavoti[popnumber];

	// Seleziono i genitori dalla population e li assegno a newpopulation (tra parentesi c'è il criterio)
	Selezione(2);	

	// Applico crossover ai genitori, genero i figli
	for(int i=0; i<(popnumber/2); i++)	Crossover(newpopulation[2*i], newpopulation[2*i+1]);

	// Applico mutazioni ai figli
	Mutazioni();	

	// Calcolo voti, li ordino in modo crescente
	for(int i=0; i<popnumber; ++i){
		votiPop[i] = L1(newpopulation[i],citynumber);
		copiavoti[i] = votiPop[i];
	}

	QuickSort(votiPop,0,popnumber-1);

	// Calcolo medie della migliore metà della popolazione (la prima metà)
	for(int i=0; i<popnumber/2; i++)	avg += votiPop[i];
	averages[num+1] = avg/((int)popnumber/2);

	// Calcolo voto minimo e salvo l'individuo se è il migliore trovato fino a questa generazione
	votominimo = votiPop[0];
	bestwalks[num+1] = votominimo;
	int pos = Position(votominimo, copiavoti, popnumber);

	if(votominimo < votominimoassoluto){
		votominimoassoluto = votominimo;
		for(int i=0; i<citynumber; ++i)		bestind[i] = newpopulation[pos][i];
	}

	// Vecchia popolazione == Nuova popolazione
	for(int i=0; i<popnumber; ++i){
		for(int j=0; j<citynumber; ++j) population[i][j] = newpopulation[i][j];
	}

}


// Selezione dei genitori
void genetic :: Selezione(int criterio){

	switch(criterio){
		case 1:{	// ---------------- Applico selezione Fitness Proportionate ---------------- //

			double totalevotiG = 0, randnum = 0;
			int rrr = 0;

			// Calcolo i voti come inverso di L1
			for(int i=0; i<popnumber; ++i){
				votiPop[i] = 1./(L1(population[i],citynumber));
				totalevotiG += votiPop[i];
			}

			// Normalizzo i voti e li riassegno in modo da avere una somma cumulativa progressiva
			for(int i=0; i<popnumber; ++i){
				votiPop[i] /= totalevotiG;
				if(i != 0)	votiPop[i] += votiPop[i-1];
			}

			// Estraggo gli indici dei genitori
			for(int i=0; i<popnumber; ++i){
				randnum = generator->Rndm();
				rrr = 0;
				for(int j=0; j<popnumber; ++j){
					if(votiPop[j] < randnum)	rrr++;
				}

				// Creo genitori
				for(int j=0; j<citynumber; ++j) newpopulation[i][j] = population[rrr][j];
			}
					
			break;
		}	// ----------------------------------------------------------------------------- //

		case 2:{	// ---------------- Estrazione casuale con voti ordinati ---------------- //

			int s = 0, pos = 0;
    		double random = 0;
			double copiavoti[popnumber];

			// Calcolo i voti
			for(int i=0; i<popnumber; ++i){
				votiPop[i] = (L1(population[i],citynumber));
				copiavoti[i] = votiPop[i];
			}

			QuickSort(votiPop,0,popnumber-1);

			for(int i=0; i<popnumber; ++i){
    			random = pow(generator->Rndm(),selectionexp);
    			s = int(popnumber * random);
				pos = Position(votiPop[s], copiavoti, popnumber);
				for(int j=0; j<citynumber; ++j) newpopulation[i][j] = population[pos][j];
			}

			break;
		}	// ----------------------------------------------------------------------------- //
	}
}


// -------------------------------------------------------------- SORTING -------------------------------------------------------------- //

// QuickSort: arr[] --> array to be sorted, low  --> starting index, high  --> ending index
void genetic :: QuickSort(double * arr, int low, int high){

    if(low < high){
		int pi = Partition(arr, low, high);	// pi is partitioning index, arr[p] is now at right place
        QuickSort(arr, low, pi - 1);		// Separately sort elements before partition and after partition 
        QuickSort(arr, pi + 1, high); 
    } 
} 


int genetic :: Partition(double *arr, int low, int high){
	 
    double pivot = arr[high];	// Pivot
	int i = (low-1);		// Index of smaller element
  
    for(int j = low; j <= high-1; j++){ 
        if(arr[j] <= pivot){	// If current element is smaller than or equal to pivot
            i++;    			// increment index of smaller element 
            Swap(arr, i, j); 
        } 
    } 

    Swap(arr, i+1, high);
    return (i+1); 
} 


int genetic :: Position(double valore, double * lista, int dimensione){

	int pos = 0;

	for(int i=0; i<dimensione; i++){
		if(lista[i] == valore){
			pos = i;
			break;
		}
	}

	if(pos < 0 || pos > dimensione)		cout << endl << "ERRORE IN POSITION" << endl;
	return pos;
}


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

	for(int i = 0; i < citynumber; ++i){
		for(int k = 0; k < 2; ++k)	citypositions[i][k] = generator->Uniform(0.,squareside);
	}
}

// Configurazione circolare
void genetic :: RandomCirc(){

	double theta = 0;

	for(int i = 0; i < citynumber; ++i){
		theta = generator->Uniform(0.,2.*M_PI);
		citypositions[i][0] = circleradius + circleradius*cos(theta);
		citypositions[i][1] = circleradius + circleradius*sin(theta);
	}
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



// ------------------------------------------------------------- CROSSOVER ------------------------------------------------------------- //

void genetic :: Crossover(int* genitore1, int* genitore2){

    int errore1 = 0, errore2 = 0;
	int crosspos = (int)(generator->Uniform(1.,citynumber-1));
	int copia1[citynumber], copia2[citynumber];
	
	// Copie degli individui
	for(int i=0; i<citynumber; i++)		copia1[i] = genitore1[i],	copia2[i] = genitore2[i];
	
    // CROSSOVER
	for(int i=crosspos; i<citynumber; i++){

        // Crossover individuo 1
		for(int j=1; j<citynumber; j++){
            errore1 = 0;
            for(int k=1; k<i; k++){
				if(genitore1[k] == copia2[j]){
					    errore1 = 1;
					    break;
				}
			}

            if(errore1 == 0){
                genitore1[i] = copia2[j];
                break;
            }
        }

        // Crossover individuo 2
		for(int j=1; j<citynumber; j++){
            errore2 = 0;
            for(int k=1; k<i; k++){
				if(genitore2[k] == copia1[j]){
					    errore2 = 1;
					    break;
				}
			}

            if(errore2 == 0){
                genitore2[i] = copia1[j];
                break;
            }
		}
	}

}


// ------------------------------------------------------------- MUTAZIONI ------------------------------------------------------------- //

// Applica mutazioni di tipo casuale alla popolazione
void genetic :: Mutazioni(){

	int mutationnum = 0;

	for(int i=0; i<popnumber; ++i){
	
		mutationnum = (int)(generator->Uniform(0.,5.));

		switch(mutationnum){
			case 1:{	if(generator->Rndm() <= pm1)	Mutazione1(newpopulation[i]);
						break;
			}
			case 2:{	if(generator->Rndm() <= pm2)	Mutazione2(newpopulation[i]);
						break;
			}
			case 3:{	if(generator->Rndm() <= pm3)	Mutazione3(newpopulation[i]);
						break;
			}
			case 4:{	if(generator->Rndm() <= pm4)	Mutazione4(newpopulation[i]);
						break;
			}
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

	ofstream WriteConfig("solution.dat"), WriteBest("walks.dat"), WriteAvg("averages.dat");
	
	for(int i = 0; i < citynumber; ++i)		WriteConfig << bestind[i] << endl;

	for(int i = 0; i < numgenerations+1; ++i){
		WriteBest << bestwalks[i] << endl;
		WriteAvg << averages[i] << endl;
	}

	WriteAvg.close();
	WriteBest.close();
	WriteConfig.close();
}


// Dealloca la memoria usata
void genetic ::DeleteMemory(){

	delete generator;
    delete [] citylist;
	delete [] bestind;
	delete [] bestwalks;
	delete [] averages;
	for(int i = 0; i < citynumber; i++)	delete [] citypositions[i];
	delete [] citypositions;
	for(int i = 0; i < popnumber; i++)	delete [] population[i];
	delete [] population;
	for(int i = 0; i < popnumber; i++)	delete [] newpopulation[i];
	delete [] newpopulation;
	delete [] votiPop;
}