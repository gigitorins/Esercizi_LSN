#include <iostream>
#include <fstream>
#include <cmath>
#include <TRandom3.h>
#include "qmc1d.h"

#define LEFT 0
#define RIGHT 1

using namespace std;

int main(int argc, char**argv){

    QMC1D qmcd1d;
    qmcd1d.readInput(argc, argv);  
	qmcd1d.initialize();  

    for(int i=0; i<qmcd1d.equilibration; i++){

        if(qmcd1d.PIGS){  // only a PIGS polymer has a start and an end. 
            qmcd1d.brownianMotion(LEFT);
            qmcd1d.brownianMotion(RIGHT);
        }
		
        qmcd1d.translation();
		
        for(int j=0;j<qmcd1d.brownianBridgeAttempts;j++)   qmcd1d.brownianBridge();
	}
	
	for(int b=0;b<qmcd1d.blocks;b++){
        for(int i=0;i<qmcd1d.MCSTEPS;i++){
            if(qmcd1d.PIGS){
                qmcd1d.brownianMotion(LEFT);
                qmcd1d.brownianMotion(RIGHT);
			}

			qmcd1d.translation();
		
			for(int j=0;j<qmcd1d.brownianBridgeAttempts;j++)   qmcd1d.brownianBridge();
			
            qmcd1d.upgradeAverages();
		}

		cout << "Completed block: " << b+1 << "/" << qmcd1d.blocks << endl;
        qmcd1d.endBlock();
	}
	
    qmcd1d.consoleOutput();
    qmcd1d.finalizePotentialEstimator();
	qmcd1d.finalizeKineticEstimator();
    qmcd1d.finalizeHistogram();

    qmcd1d.deleteMemory();  // de-allocate dynamic variables

	return 0;
}