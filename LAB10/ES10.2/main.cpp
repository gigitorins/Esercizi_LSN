#include <cstdlib>   
#include <iostream>     
#include <fstream>      
#include <string>
#include <sstream>
#include <cmath>
#include <ostream>
#include <iomanip>
#include <vector>
#include "mpi.h"
#include "genetic.h"

using namespace std;

int main(int argc, char ** argv){

    MPI::Init(argc, argv);
    int rank = MPI::COMM_WORLD.Get_rank();
    int size = MPI::COMM_WORLD.Get_size();
    MPI_Status stat;
    string r = to_string(rank);

    genetic gen1;
	gen1.Initialize(argc, argv);
    int ngens = gen1.numgenerations, cityN = gen1.citynumber;

    std::srand ( unsigned ( std::time(0) ) );
    std::vector<int> swaps;
    for(int i=0; i<size; ++i)  swaps.push_back(i);
    int appo1[cityN], appo2[cityN];
    int nmigr = 10;
    
    cout << endl << "Best walk iniziale: " << gen1.votominimoassoluto << endl;

    for(int i=0; i<ngens; i++){

        if(i%nmigr == 0){
            // We randomly shuffle the vector swaps, which we will use to swap the best individuals between the processes
			if(rank == 0) random_shuffle(swaps.begin(), swaps.end());
			MPI_Bcast(&swaps.front(), 4, MPI_INT, 0, MPI_COMM_WORLD);
			
			for(int k=0; k<2; k++){     // Swap the best individuals between the processes with rank swaps[i] and swaps[i+1] (swaps[0] <-> swaps[1], swaps[2] <-> swaps[3])

				if(rank == swaps[2*k]){
                    for(int j=0; j<cityN; ++j)  appo1[j] = gen1.bestind[j];
					MPI_Send(&appo1, cityN, MPI_INT, swaps[2*k+1], 2*k, MPI_COMM_WORLD);
					MPI_Recv(&appo2, cityN, MPI_INT, swaps[2*k+1], (2*k+1), MPI_COMM_WORLD, &stat);
                    for(int j=0; j<cityN; ++j)  gen1.population[0][j] = appo2[j];
				}
				
				if(rank == swaps[2*k+1]){
                    for(int j=0; j<cityN; ++j)  appo2[j] = gen1.bestind[j];
					MPI_Recv(&appo1, cityN, MPI_INT, swaps[2*k], 2*k, MPI_COMM_WORLD, &stat);
					MPI_Send(&appo2, cityN, MPI_INT, swaps[2*k], (2*k+1), MPI_COMM_WORLD);
                    for(int j=0; j<cityN; ++j)  gen1.population[0][j] = appo1[j];
				}
			}

		}

        gen1.Generazione(i);
        cout << "-- Progress: " << int(double(i)/double(ngens)*10000)/(double)100<< "% \r";
    }
    
    cout << "Best walk trovato: " << gen1.votominimoassoluto << endl;

	gen1.PrintConfig(r);
	gen1.DeleteMemory();

    MPI::Finalize();
    if(rank == 0)   cout << endl << "*** PROGRAM COMPLETED *** " << endl << endl;

	return 0;
}