#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "VCM.h"

using namespace std;


VCM :: VCM(){

}


VCM :: ~VCM(){
	delete [] histogram;
	delete [] blk_histo;
	delete [] glob_histo;
	delete [] glob_histo2;
}


void VCM :: Input(int argc, char ** argv){

	rnd.Set();	// Read seed for random numbers
    if(argc != 1){
		sigma = atof(argv[1]);
		mu = atof(argv[2]);
    }
    cout << "Sigma: " << sigma << "   mu: " << mu << endl;
    
    bin_length = (sup - inf)/(double)n_bins;

	histogram = new int[n_bins];
	blk_histo = new int[n_bins];
	glob_histo = new double[n_bins];
	glob_histo2 = new double[n_bins];

}


void VCM :: Move(void){
    
	double p_old = pow(Wave_function(x), 2);
	double x_new = x + delta*(rnd.Rannyu(-1.,1.));
    double p_new = pow(Wave_function(x_new), 2);
    	
	// Metropolis test
	double test = min(p_new/p_old, 1.);
        
	if(rnd.Rannyu() <= test){
		x = x_new;
		accepted++;
	}
	
	attempted++;
}


double VCM :: Wave_function(double x){
    return exp(-pow((x-mu),2)/(2.*sigma*sigma)) + exp(-pow((x+mu),2)/(2.*sigma*sigma));
}


double VCM :: Second_derivative(double x){

    double part1 = exp(-pow((x-mu),2)/(2.*sigma*sigma)) * ( pow((x-mu),2)/(pow(sigma,4)) - 1./(sigma*sigma));
	double part2 = exp(-pow((x+mu),2)/(2.*sigma*sigma)) * ( pow((x+mu),2)/(pow(sigma,4)) - 1./(sigma*sigma));
    
    return part1 + part2;
}


double VCM :: Potential(double x){
    
    return pow(x,4.)-2.5*pow(x,2.);
}


void VCM :: Measure(){

    walker = (-0.5 * Second_derivative(x) + Potential(x) * Wave_function(x)) / Wave_function(x);

    if(x > inf && x < sup){
        int position = int( (x - inf)/bin_length );
        histogram[position]++;
    }
}


void VCM :: Reset(int iblk){	// Reset block averages
   
	if(iblk == 1){
		for(int i=0; i<n_bins; ++i)	glob_histo[i] = 0, glob_histo2[i] = 0;
	}

    blk_avg = 0.;
    blk_norm = 0;
    attempted = 0;
    accepted = 0;
    tot = 0;

	for(int i=0; i<n_bins; i++){   
		histogram[i] = 0;
		blk_histo[i] = 0;
    }
}


void VCM :: Accumulate(void){	// Update block averages

    blk_avg += walker;
    blk_norm++;
    
    for(int i=0; i<n_bins; i++){
        blk_histo[i] += histogram[i];
        tot += histogram[i];
    }
        
    return;
}


void VCM :: Averages(int i){

    ofstream Ham("output.hamilton.0", ios::app);
	ofstream Histo("output.histo.0");
    
    stima_h = blk_avg/((double)blk_norm);
    glob_av += stima_h;
    glob_av2 += stima_h*stima_h;
    err_h = Error(glob_av, glob_av2, i);
    Ham << " " << i << " " << stima_h << " " << glob_av/(double)i << " " << err_h << endl;
    
    for(int k=0; k<n_bins; k++){
        stima_histo = blk_histo[k]/(bin_length * tot);
        glob_histo[k] += stima_histo;
        glob_histo2[k] += stima_histo*stima_histo;
        err_histo = Error(glob_histo[k], glob_histo2[k], i);
        Histo << " " << (k*bin_length + inf) << " " << stima_histo << " " << glob_histo[k]/(double)i << " " << err_histo << endl;
    }

	Ham.close();
    Histo.close();
}


double VCM :: Error(double sum, double sum2, int iblk){
    if(iblk == 1)	return 0.0;
    else			return	sqrt( abs( (sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk) ));
}