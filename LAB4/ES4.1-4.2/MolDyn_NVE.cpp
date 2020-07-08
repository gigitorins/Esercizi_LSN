#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <string>
#include <sstream>
#include "MolDyn_NVE.h"

using namespace std;



MolDyn :: MolDyn(){
	x = new double[m_part], y = new double[m_part], z = new double[m_part];
	xold = new double[m_part], yold = new double[m_part], zold = new double[m_part];
	vx = new double[m_part], vy = new double[m_part], vz = new double[m_part];
	ave_pot = new double[nblock], ave_kin = new double[nblock], ave_etot = new double[nblock], ave_temp = new double[nblock];
	sum_pot = new double[nblock], sum_kin = new double[nblock], sum_etot = new double[nblock], sum_temp = new double[nblock];
	sum2_pot = new double[nblock], sum2_kin = new double[nblock], sum2_etot = new double[nblock], sum2_temp = new double[nblock];
	err_pot = new double[nblock], err_kin = new double[nblock], err_etot = new double[nblock], err_temp = new double[nblock];
	walker = new double[n_props];
}



MolDyn :: ~MolDyn(){
	delete [] x;
	delete [] y;
	delete [] z;	
	delete [] xold;
	delete [] yold;
	delete [] zold;
	delete [] vx;
	delete [] vy;
	delete [] vz;
	delete [] ave_pot;
	delete [] ave_kin;
	delete [] ave_etot;
	delete [] ave_temp;
	delete [] sum_pot;
	delete [] sum_kin;
	delete [] sum_etot;
	delete [] sum_temp;
	delete [] sum2_pot;
	delete [] sum2_kin;
	delete [] sum2_etot;
	delete [] sum2_temp;
	delete [] err_pot;
	delete [] err_kin;
	delete [] err_etot;
	delete [] err_temp;
	delete [] walker;
}



void MolDyn :: Input(int argc, char ** argv){		//Prepare all stuff for the simulation

	ifstream ReadInput("input.dat");
	// cout << endl << "/******** Molecular dynamics simulation in NVE ensemble ********/" << endl;

	// Simulation number
	if(argc > 1)	numsim = atoi(argv[1]);

/*	cout << endl << "Classic Lennard-Jones fluid" << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "The program uses Lennard-Jones units " << endl;
*/
	seed = 1;			//Set seed for random numbers
	srand(seed);			//Initialize random number generator
  
	ReadInput >> temp;
	ReadInput >> npart;
	ReadInput >> rho;	
	vol = (double)npart/rho;
	box = pow(vol,1.0/3.0);
	ReadInput >> rcut;
	ReadInput >> delta;
	ReadInput >> nstep;
	stepinblock = nstep/nblock;
	ReadInput >> iprint;
	ReadInput >> restart;		// Restart
	ReadInput.close();	

/*	cout << "Number of particles = " << npart << endl;
	cout << "Density of particles = " << rho << endl;
	cout << "Volume of the simulation box = " << vol << endl;
	cout << "Edge of the simulation box = " << box << endl;
	cout << "The program integrates Newton equations with the Verlet method " << endl;
	cout << "Time step = " << delta << endl;
	cout << "Number of steps = " << nstep << endl << endl;
*/
	InputConfig(restart);

	// Azzeramento vettori per le grandezze termodinamiche
	for(int i=0; i<nblock; ++i){	
		ave_pot[i] = 0, ave_kin[i] = 0, ave_etot[i] = 0, ave_temp[i] = 0;
		sum_pot[i] = 0, sum_kin[i] = 0, sum_etot[i] = 0, sum_temp[i] = 0;
		sum2_pot[i] = 0, sum2_kin[i] = 0, sum2_etot[i] = 0, sum2_temp[i] = 0;
		err_pot[i] = 0, err_kin[i] = 0, err_etot[i] = 0, err_temp[i] = 0;
	}

	for(int i=0; i<n_props; i++)	walker[i] = 0.;

	bin_size = (box/2.0)/(double)nbins;
}



void MolDyn :: InputConfig(int restart){

if(restart == 0){	/********************* RESTART == 0 *********************/

	//Read initial configuration
	ifstream ReadConf("config.0");
	cout << "Read initial configuration from file config.0 " << endl;
	for (int i=0; i<npart; ++i){
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = x[i] * box;
		y[i] = y[i] * box;
		z[i] = z[i] * box;
	}
	ReadConf.close();

	//Prepare initial velocities
	cout << "Prepare random velocities with center of mass velocity equal to zero " << endl;
	double sumv[3] = {0.0, 0.0, 0.0}, sumv2 = 0.0, fs;
	for (int i=0; i<npart; ++i){
		vx[i] = rand()/double(RAND_MAX) - 0.5;
		vy[i] = rand()/double(RAND_MAX) - 0.5;
		vz[i] = rand()/double(RAND_MAX) - 0.5;
		sumv[0] += vx[i];
		sumv[1] += vy[i];
		sumv[2] += vz[i];
	}
	for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
	for (int i=0; i<npart; ++i){
		vx[i] = vx[i] - sumv[0];
		vy[i] = vy[i] - sumv[1];
		vz[i] = vz[i] - sumv[2];
		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	}
	sumv2 /= (double)npart;
	fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
	for (int i=0; i<npart; ++i){
		vx[i] *= fs;
		vy[i] *= fs;
		vz[i] *= fs;
		xold[i] = Pbc(x[i] - vx[i] * delta);
		yold[i] = Pbc(y[i] - vy[i] * delta);
		zold[i] = Pbc(z[i] - vz[i] * delta);
	}

	ChangeInput0();	// Cambia restart da 0 a 1 in input.dat

}else	{	/********************* RESTART == 1 *********************/

	//Read initial configuration
	ifstream ReadConf("config.final");
	ifstream ReadConf0("old.final");
	cout << "Read initial configurations from file config.final and old.final " << endl << endl;

	for (int i=0; i<npart; ++i){
		ReadConf >> x[i] >> y[i] >> z[i];
		ReadConf0 >> xold[i] >> yold[i] >> zold[i];
		x[i] = x[i] * box;
		y[i] = y[i] * box;
		z[i] = z[i] * box;
		xold[i] = xold[i] * box;
		yold[i] = yold[i] * box;
		zold[i] = zold[i] * box;
	}
	ReadConf.close();
	ReadConf0.close();

	Move();		// Move once with Verlet, x(t)=x(t+dt) and xold=x(t)

	// Compute initial velocities with center of mass velocity equal to zero
	double sumv[3] = {0.0, 0.0, 0.0}, sumv2 = 0.0;
	for(int i=0; i<npart; ++i){			// Velocities at time = t+dt/2
		vx[i] = Pbc(x[i] - xold[i])/(delta);			
		vy[i] = Pbc(y[i] - yold[i])/(delta);
		vz[i] = Pbc(z[i] - zold[i])/(delta);
		sumv[0] += vx[i];
		sumv[1] += vy[i];
		sumv[2] += vz[i];
	}

	for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;

	for (int i=0; i<npart; ++i){
		vx[i] = vx[i] - sumv[0];
		vy[i] = vy[i] - sumv[1];
		vz[i] = vz[i] - sumv[2];
		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	}

	sumv2 /= (double)npart;	
	double fs = sqrt((3 * temp) / sumv2);		// Scale factor
	for (int i=0; i<npart; ++i){			// New velocities
		vx[i] *= fs;
		vy[i] *= fs;
		vz[i] *= fs;				
		xold[i] = Pbc(x[i] - vx[i] * delta);	// New r(old)
		yold[i] = Pbc(y[i] - vy[i] * delta);
		zold[i] = Pbc(z[i] - vz[i] * delta);
	}
}		/********************* END RESTART *********************/
}



void MolDyn :: ChangeInput0(){	// Cambia restart da 0 a 1 nel file input.dat

	string readout, search="0", replace="1";
	string oldname = "input.dat";
	string newname1 = "input1.dat", newname2 = "input2.dat";
	ofstream outFile(newname1);
	ifstream readFile(oldname);
	while(getline(readFile,readout)){
		if(readout == search)	outFile << replace << endl;
		else	outFile << readout << endl;
	}
	outFile.close();
	readFile.close();

	int result;
	result = rename(oldname.c_str(), newname2.c_str());
	if(result != 0)	cout << endl << "Error renaming file " << oldname << endl;
	result = rename(newname1.c_str(), oldname.c_str());
	if(result != 0)	cout << endl << "Error renaming file " << newname1 << endl;
	if(remove(newname2.c_str()) != 0)	cout << endl << "Error deleting file " << newname2 << endl;
	cout << "Input.dat successfully modified " << endl << endl;
}



void MolDyn :: Move(){		// Move particles with Verlet algorithm

	double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

	for(int i=0; i<npart; ++i){		// Force acting on particle i
		fx[i] = Force(i,0), fy[i] = Force(i,1), fz[i] = Force(i,2);
	}

	for(int i=0; i<npart; ++i){		// Verlet integration scheme

		xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );	// New positions
		ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
		znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );
		vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);			// Velocities
		vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
		vz[i] = Pbc(znew - zold[i])/(2.0 * delta);
		xold[i] = x[i];							// Old positions == actual positions
		yold[i] = y[i];
		zold[i] = z[i];
		x[i] = xnew;							// Actual position == new position
		y[i] = ynew;
		z[i] = znew;
	}
}



double MolDyn :: Force(int ip, int idir){	// Compute forces as -Grad_ip V(r)
	
	double f=0.0, dvec[3], dr;
	for(int i=0; i<npart; ++i){
		if(i != ip){
			dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
			dvec[1] = Pbc( y[ip] - y[i] );
			dvec[2] = Pbc( z[ip] - z[i] );
			dr = sqrt(dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2]);
		if(dr < rcut){	f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); }	// -Grad_ip V(r)
		}
	} 
	
	return f;
}



void MolDyn :: Measure(int step){	// Properties measurement

	// Reset observables
	int bin = 0, sstep = 0;
	attempted = 0, accepted = 0;
    double inf = 0, sup = 0;	// Boundaries for blocks in g(r) hystogram
	double v = 0.0, t = 0.0, vij, dx, dy, dz, dr;	
	for(int k=0; k<0+nbins; ++k)	walker[k] = 0.0;

	ofstream Epot("output_epot.dat",ios::app);
	ofstream Ekin("output_ekin.dat",ios::app);
	ofstream Temp("output_temp.dat",ios::app);
	ofstream Etot("output_etot.dat",ios::app);
	ofstream aveEpot("ave_epot.out",ios::app);
	ofstream aveEkin("ave_ekin.out",ios::app);
	ofstream aveTemp("ave_temp.out",ios::app);
	ofstream aveEtot("ave_etot.out",ios::app);

	// Potential energy (cycle over pairs of particles)
	for(int i=0; i<npart-1; ++i){
		for (int j=i+1; j<npart; ++j){
			dx = Pbc( x[i] - x[j] );
			dy = Pbc( y[i] - y[j] );
			dz = Pbc( z[i] - z[j] );
			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);

			// Update of the histogram of g(r)
            for (int k=0; k<0+nbins; ++k){
				sstep = k-0;
				inf = (sstep)*bin_size;
				sup = (sstep+1)*bin_size;      
				bin_size = (box/2.0)/(double)nbins;
                attempted++;
                // I see if the dr is in the k-th interval (the size is divided into 100 intervals).
                if(dr >= inf &&  dr < sup){
                    walker[k] += 2;
                    accepted++;   
                }
                
            }

			if(dr < rcut){
			vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
			v += vij;	
			}
		}          
	}

	// Kinetic energy
	for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
	stima_pot = v/(double)npart;			// Potential energy per particle
	stima_kin = t/(double)npart;			// Kinetic energy per particle
	stima_temp = (2.0 / 3.0) * t/(double)npart;	// Temperature
	stima_etot = (t+v)/(double)npart;		// Total energy per particle

	Epot << stima_pot  << endl;
	Ekin << stima_kin  << endl;
	Temp << stima_temp << endl;
	Etot << stima_etot << endl;
	Epot.close();
	Ekin.close();
	Temp.close();
	Etot.close();

	/************************************ Block averaging ************************************/

	ave_pot[counter] += stima_pot*10./stepinblock;	// Calcolo le medie per ogni blocco, counter corre su nblock == numero dei blocchi
	ave_kin[counter] += stima_kin*10./stepinblock;	 
	ave_etot[counter] += stima_etot*10./stepinblock;
	ave_temp[counter] += stima_temp*10./stepinblock;

	if(step%stepinblock == 0){		// Numero di step == multiplo di numero di step per blocco

		for(int k=0; k<(counter+1); k++){
	  		sum_pot[counter]+= ave_pot[k], sum2_pot[counter]+= pow(ave_pot[k],2);		// SUM & SUM2 {k=0,counter}
	  		sum_kin[counter]+= ave_kin[k], sum2_kin[counter]+= pow(ave_kin[k],2);		
	  		sum_etot[counter]+= ave_etot[k], sum2_etot[counter]+= pow(ave_etot[k],2);		
	  		sum_temp[counter]+= ave_temp[k], sum2_temp[counter]+= pow(ave_temp[k],2);		
		}  	
		
		sum_pot[counter]/=(counter+1), sum2_pot[counter]/=(counter+1);		// Cumulative average & cumulative square average
  		sum_kin[counter]/=(counter+1), sum2_kin[counter]/=(counter+1);		// Cumulative average & cumulative square average
		sum_etot[counter]/=(counter+1), sum2_etot[counter]/=(counter+1);	// Cumulative average & cumulative square average
		sum_temp[counter]/=(counter+1), sum2_temp[counter]/=(counter+1);	// Cumulative average & cumulative square average	

		if(counter==0){			// Statistical uncertainty		
			err_pot[counter] = 0;
			err_kin[counter] = 0;
			err_etot[counter] = 0;
			err_temp[counter] = 0;
		}else	{
			err_pot[counter] = sqrt((sum2_pot[counter] - pow(sum_pot[counter],2))/counter);
			err_kin[counter] = sqrt((sum2_kin[counter] - pow(sum_kin[counter],2))/counter);
			err_etot[counter] = sqrt((sum2_etot[counter] - pow(sum_etot[counter],2))/counter);
			err_temp[counter] = sqrt((sum2_temp[counter] - pow(sum_temp[counter],2))/counter);
		}
			
  		aveEpot << sum_pot[counter] << " " << err_pot[counter] << endl;		// Stampo medie & errori su file
  		aveEkin << sum_kin[counter] << " " << err_kin[counter] << endl;
  		aveEtot << sum_etot[counter] << " " << err_etot[counter] << endl;
  		aveTemp << sum_temp[counter] << " " << err_temp[counter] << endl;

		counter += 1;			// Passo al blocco successivo
	}

	aveEpot.close();
	aveEkin.close();
	aveTemp.close();
	aveEtot.close();

	/*****************************************************************************************/

	return;
}



void MolDyn :: ConfFinal(){		// Write final configuration

	ofstream WriteConf;
	cout << "Print final configuration to file config.final " << endl << endl;
	WriteConf.open("config.final");

	for(int i=0; i<npart; ++i)	WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
	WriteConf.close();
}



void MolDyn :: ConfFinalOld(){		// Write old.final configuration

	ofstream WriteConf;
	cout << "Print configuration to file old.final " << endl;
	WriteConf.open("old.final", std::ofstream::trunc);

	for (int i=0; i<npart; ++i)	WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
	WriteConf.close();
}



void MolDyn :: ConfXYZ(int nconf){	// Write configuration in .xyz format

	ofstream WriteXYZ("frames/config_" + to_string(nconf) + ".xyz");
	WriteXYZ << npart << endl;
	WriteXYZ << "This is only a comment!" << endl;
	for (int i=0; i<npart; ++i){
		WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
	}
	WriteXYZ.close();
}



double MolDyn :: Pbc(double r){		// Algorithm for periodic boundary conditions with side L=box
	return r - box * rint(r/box);
}



