#include <iostream>
#include <fstream>
#include <cmath>
#include <TRandom3.h>

#ifndef __genetic__
#define __genetic__

class genetic {

    public:
		// Random numbers
		TRandom3* generator;

		// Variables
		int citynumber = 32;
		int popnumber = 3000;
		int numgenerations = 200;
		int * citylist;				// Lista ordinata delle città
		double ** citypositions;	// Posizioni delle città		
		int ** population;			// Popolazione
		int ** newpopulation;		// Popolazione dopo crossover/mutazione
		double * votiPop;
		double votominimo;
		double votominimoassoluto;
		int * bestind;
		double * averages;
		double * bestwalks;

		double pm1 = 0.1;
		double pm2 = 0.1;
		double pm3 = 0.1;
		double pm4 = 0.1;
		double probcross = 0.7;
		int selectionexp = 2;

		// Square and circular configuration
		int func_type = 0;
		double squareside = 10.;
		double circleradius = squareside/2.;

		// Constructors
		genetic();
		// Destructor
		~genetic();
		// Methods
		void Initialize(int, char **);	
		void RandomSquare();
		void RandomCirc();
		void PrintInd(int *);
		void PrintPop();
		void PopIniziale();
		void GeneraInd(int);
		void Generazione(int);
		void Selezione(int);	
		void QuickSort(double *, int, int);
		int Partition(double *, int, int);
		int Position(double, double *, int);
		void Swap(int*, int, int);
		void Swap(double*, int, int);
		double L1(int*, int);
		double L2(int*, int);
		void Crossover(int*, int*);
		void Mutazioni();
		void Mutazione1(int *);
		void Mutazione2(int *);
		void Mutazione3(int *);
		void Mutazione4(int *);
		void DeleteMemory();
		void ReadConfig(int);
		void PrintConfig();

};

#endif