#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <string>
#include <sstream>
#include "MolDyn_NVE.h"

using namespace std;



MolDyn :: MolDyn(){

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
	delete [] ave_gofr;
	delete [] glob_av;
	delete [] glob_av2;
}



void MolDyn :: Input(int argc, char ** argv){		//Prepare all stuff for the simulation

	ifstream ReadInput("input.dat");

	// Simulation number
	if(argc > 1)	numsim = atoi(argv[1]);
	
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

	if(argc > 2)	nblock = atoi(argv[2]);
    if(argc > 3)	nstep = atoi(argv[3]);

	if(argc > 2)	stepinblock = nstep;
	else	stepinblock = nstep/nblock;

	ReadInput >> iprint;
	ReadInput >> restart;		// Restart
	ReadInput.close();	

	// Creazione vettori
	x = new double[npart], y = new double[npart], z = new double[npart];
	xold = new double[npart], yold = new double[npart], zold = new double[npart];
	vx = new double[npart], vy = new double[npart], vz = new double[npart];
	ave_pot = new double[nblock], ave_kin = new double[nblock], ave_etot = new double[nblock], ave_temp = new double[nblock];
	sum_pot = new double[nblock], sum_kin = new double[nblock], sum_etot = new double[nblock], sum_temp = new double[nblock];
	sum2_pot = new double[nblock], sum2_kin = new double[nblock], sum2_etot = new double[nblock], sum2_temp = new double[nblock];
	err_pot = new double[nblock], err_kin = new double[nblock], err_etot = new double[nblock], err_temp = new double[nblock];
	walker = new double[nbins];
	ave_gofr = new double[nbins], glob_av = new double[nbins], glob_av2 = new double[nbins];

	InputConfig(restart);

	// Azzeramento vettori per le grandezze termodinamiche
	for(int i=0; i<nblock; ++i){	
		ave_pot[i] = 0, ave_kin[i] = 0, ave_etot[i] = 0, ave_temp[i] = 0;
		sum_pot[i] = 0, sum_kin[i] = 0, sum_etot[i] = 0, sum_temp[i] = 0;
		sum2_pot[i] = 0, sum2_kin[i] = 0, sum2_etot[i] = 0, sum2_temp[i] = 0;
		err_pot[i] = 0, err_kin[i] = 0, err_etot[i] = 0, err_temp[i] = 0;
	}

	// Azzero walker
	for(int i=0; i<nbins; i++){
		walker[i] = 0.;
		ave_gofr[i] = glob_av[i] = glob_av2[i] = 0.;
	}

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
	// cout << "Prepare random velocities with center of mass velocity equal to zero " << endl;
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
	// cout << "Input.dat successfully modified " << endl << endl;
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
			if(dr < rcut)	f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); 	// -Grad_ip V(r)
		}
	} 
	
	return f;
}



void MolDyn :: Measure(int step){	// Properties measurement

	// Reset observables
	attempted = 0, accepted = 0;
    double inf = 0, sup = 0;	// Boundaries for blocks in g(r) hystogram
	double r, gdir, DeltaV;
	double v = 0.0, t = 0.0, vij, dx, dy, dz, dr;	

	// Reset the hystogram of g(r)
	for(int i=0; i<nbins; i++)	walker[i] = 0.;

	ofstream Gave("output.gave.0",ios::app);

	// Potential energy (cycle over pairs of particles)
	for(int i=0; i<npart-1; ++i){
		for (int j=i+1; j<npart; ++j){
			dx = Pbc( x[i] - x[j] );
			dy = Pbc( y[i] - y[j] );
			dz = Pbc( z[i] - z[j] );
			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);

			// Update of the histogram of g(r)
            for (int k=0; k<nbins; ++k){
				inf = k*bin_size;
				sup = (k+1)*bin_size;      
                attempted++;

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


	for(int k=0; k<nbins; ++k)	ave_gofr[k] += walker[k]*10.;

	if(step%stepinblock == 0){		// Numero di step == multiplo di numero di step per blocco

		for(int k=0; k<nbins; k++){
			r = bin_size*k;
			DeltaV = 4./3 * M_PI * r*r*r;
			if(DeltaV == 0)	gdir = 0;
			else			gdir = ave_gofr[k]/rho/(double)npart/DeltaV/nbins;

			glob_av[k] += gdir;
			glob_av2[k] += gdir*gdir;
			err_gdir = Error(glob_av[k], glob_av2[k], counter);
			Gave << " " << r << " " << glob_av[k]/(double)counter << " " << err_gdir << endl;
			ave_gofr[k] = 0.;
		}


		counter += 1;			// Passo al blocco successivo
	}

	Gave.close();

	/*****************************************************************************************/

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
	//cout << "Print configuration to file old.final " << endl;
	WriteConf.open("old.final", std::ofstream::trunc);

	for (int i=0; i<npart; ++i)	WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
	WriteConf.close();
}



void MolDyn :: ConfXYZ(int nconf){	// Write configuration in .xyz format

	ofstream WriteXYZ("frames/config_" + to_string(nconf) + ".xyz");
	WriteXYZ << npart << endl;
	for (int i=0; i<npart; ++i){
		WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
	}
	WriteXYZ.close();
}



double MolDyn :: Pbc(double r){		// Algorithm for periodic boundary conditions with side L=box
	return r - box * rint(r/box);
}


double MolDyn :: Error(double sum, double sum2, int iblk){
    if(iblk == 1)	return 0.0;
    else			return	sqrt( abs( (sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1) ));
}