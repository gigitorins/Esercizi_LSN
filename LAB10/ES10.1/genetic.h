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
		int stepnumber = 1000;
		int * citylist;				// Lista ordinata delle città
		double ** citypositions;	// Posizioni delle città		
		int * individuo;			// Individuo
		int * nuovoindividuo;		// Nuovo individuo
		int * bestind;
		double votiPop = 0;	
		double votominimo;
		int tsteps = 2000;
		int attempted = 0, accepted = 0;

		double pm1 = 0.2, pm2 = 0.2, pm3 = 0.2, pm4 = 0.2;

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
		void PrintInd(int *);
		void RandomSquare();
		void RandomCirc();	
		void Swap(int*, int, int);
		void Swap(double*, int, int);
		double L1(int*, int);
		double L2(int*, int);
		void Move(double);
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