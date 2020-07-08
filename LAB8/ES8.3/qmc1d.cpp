/*********************** QMC1D **************************
************** PATH INTEGRAL GROUND STATE ***************
************** PATH INTEGRAL MONTE CARLO ****************
************ APPLIED TO A SINGLE PARTICLE ***************
************** IN AN EXTERNAL POTENTIAL ****************/

#include "qmc1d.h"

#define LEFT 0
#define RIGHT 1

using namespace std;


QMC1D :: QMC1D(){

}


QMC1D :: ~QMC1D(){

}


double QMC1D :: potential_density_matrix(double val, double val_next){	// This is the primitive approximation without the kinetic correlation

	return (-dtau/2)*(external_potential(val)+external_potential(val_next));
}


void QMC1D :: initialize(){		// Initialization of the variables and allocation of the memory

	lambda = hbar*hbar/(2*mass);
	if(temperature==0)	PIGS=1;
	else				PIGS=0;
	
	if(PIGS)	dtau = imaginaryTimePropagation/(timeslices-1);
	else		dtau = hbar/(boltzmann*temperature*timeslices);
	
	acceptedTranslations = 0,	acceptedVariational = 0;
	acceptedBB = 0,				acceptedBM = 0;
	totalTranslations = 0,		totalVariational = 0;
	totalBB = 0,				totalBM = 0;
	
	generator = new TRandom3();

	positions = new double[timeslices];
	potential_energy = new double[timeslices];
	potential_energy_accumulator = new double[timeslices];
	potential_energy_square_accumulator = new double[timeslices];
                                                                                                                 
	kinetic_energy = new double[timeslices];
	kinetic_energy_accumulator = new double[timeslices];
	kinetic_energy_square_accumulator=new double[timeslices];
                                                                                                                
	positions_histogram=new double[histogram_bins];
	positions_histogram_accumulator=new double[histogram_bins];
	positions_histogram_square_accumulator=new double[histogram_bins];
	
	for(int i=0;i<timeslices;i++)	positions[i]=0.0;
	
	for(int i=0;i<timeslices;i++){

		potential_energy[i] = 0;
		potential_energy_accumulator[i] = 0;
		potential_energy_square_accumulator[i] = 0;

		kinetic_energy[i] = 0;
		kinetic_energy_accumulator[i] = 0;
		kinetic_energy_square_accumulator[i] = 0;
	}
	
	for(int i=0;i<histogram_bins;i++){
		positions_histogram[i] = 0;
		positions_histogram_accumulator[i] = 0;
		positions_histogram_square_accumulator[i] = 0;
	}

	alpha=0;
}


// External potential
double QMC1D :: external_potential(double val){
    return pow(val,4.) - 5./2. * pow(val,2.);
}

double QMC1D :: external_potential_prime(double val){
    return 4*pow(val,3.) - 5.*val;
}

double QMC1D :: external_potential_second(double val){
    return 12*pow(val,2.) - 5.;
}


// Variational Wave Function
double QMC1D :: variationalWaveFunction(double v){
    if(func_type == 0)	return 1.0;
    else				return exp(-0.5*pow((v-mu),2)/(sigma*sigma)) + exp(-0.5*pow((v+mu),2)/(sigma*sigma));
}

double QMC1D :: variationalWaveFunction_second(double v){

    if(func_type == 0)	return 0.0;
    else{
			double first_term = exp(-pow((v-mu),2)/(2.*sigma*sigma)) * (pow((v-mu),2)/pow(sigma,4) - 1./(sigma*sigma));
			double second_term = exp(-pow((v+mu),2)/(2.*sigma*sigma)) * (pow((v+mu),2)/pow(sigma,4) - 1./(sigma*sigma));
			return first_term + second_term;
    }
}


void QMC1D :: translation(){

	totalTranslations++;
	double delta = generator->Uniform(-delta_translation,delta_translation);
	double acc_density_matrix_difference=0;
	int last = timeslices;
	if(PIGS)	last=timeslices-1;
	
	int inext;
	double newcorr,oldcorr;

	for(int i=0;i<last;i++){
		inext = index_mask(i+1);
		newcorr = potential_density_matrix(positions[i]+delta,positions[inext]+delta);
		oldcorr = potential_density_matrix(positions[i],positions[inext]);
		acc_density_matrix_difference += oldcorr-newcorr;
	}

	// metropolis: PIGS contains also the statistical weight of the variational Wave Function.
	double acceptance_probability = exp(-acc_density_matrix_difference);
	
	if(PIGS)
		acceptance_probability *= variationalWaveFunction(positions[0]+delta)*variationalWaveFunction(positions[timeslices-1]+delta)/
									(variationalWaveFunction(positions[0])*variationalWaveFunction(positions[timeslices-1]));
	
	if(generator->Rndm()<acceptance_probability){
		for(int i=0;i<timeslices;i++)	positions[i]+=delta;
		acceptedTranslations++;
	}
}


/* BB removes a segment of the polymer, in this case from "starting_point+1" to "endpoint-1"
and replaces it with a free particle propagation. The free particle propagation is achieved
with the gaussian sampling of the kinetic part of the density matrix.
The function index_mask handles the compatibility between PIGS and PIMC: in PIGS the polymer is 
open, so you can't have a starting index greater than an ending index. In PIMC, instead, you
have a ring polymer so when you reach the end you can continue from the beginning. 
The compatibility solution that has been chosen consists in viewing the ring polymer as an open
polymer that has been closed on periodic boundary contitions. index_mask takes this into account. */
void QMC1D :: brownianBridge(){

	totalBB++;
	int available_starting_points = timeslices-brownianBridgeReconstructions-1; // for PIGS simulation
	if(!PIGS)	available_starting_points = timeslices-1;

	int starting_point = (int)(generator->Rndm()*available_starting_points);
	int endpoint = index_mask(starting_point + brownianBridgeReconstructions + 1);
	
	double starting_coord = positions[starting_point];
	double ending_coord = positions[endpoint];
	double new_segment[brownianBridgeReconstructions+2];
	new_segment[0] = starting_coord;
	new_segment[brownianBridgeReconstructions+1] = ending_coord;
	double previous_position = starting_coord;

	for(int i=0;i<brownianBridgeReconstructions;i++){

		int left_reco = brownianBridgeReconstructions-i;
		// gaussian sampling of the free particle propagator
		double average_position = previous_position + (ending_coord-previous_position)/(left_reco+1);
		double variance = 2*lambda*dtau*left_reco/(left_reco+1);
		double newcoordinate = generator->Gaus(average_position,sqrt(variance));
		new_segment[i+1] = newcoordinate;
		previous_position = newcoordinate;
	}
	
	// Metropolis --> The kinetic part has been sampled exactely, only the potential part of the density matrix determines the acceptance probability of the move
	double acc_density_matrix_difference = 0;
	int i_old, i_next_old;
	double newcorr,oldcorr;

	for(int i=0;i<brownianBridgeReconstructions+1;i++){
		i_old = index_mask(starting_point+i);
		i_next_old = index_mask(starting_point+i+1);
		newcorr = potential_density_matrix(new_segment[i],new_segment[i+1]);
		oldcorr = potential_density_matrix(positions[i_old],positions[i_next_old]);
		acc_density_matrix_difference += oldcorr-newcorr;
	}
	
	double acceptance_probability = exp(-acc_density_matrix_difference);
	if(generator->Rndm()<acceptance_probability){
		for(int i=1;i<brownianBridgeReconstructions+1;i++){
			i_old = index_mask(starting_point+i);
			positions[i_old]=new_segment[i];
		}
		acceptedBB++;
	}
}


/* BM removes a segment at one of the two ends of the polymer, and replaces it with a free particle propagation using a Brownian Bridge after the sampling of the
starting (left move) or final (right move) position. The free particle propagation is achieved with the gaussian sampling of the kinetic part of the density matrix. */
void QMC1D :: brownianMotion(int which){	// BM is called only for PIGS simulations

	int starting_point, endpoint, left_reco;
	double starting_coord, ending_coord, average_position, variance, newposition, oldposition;

	totalBM++;

	if(which==LEFT){
		starting_point = 0;
		endpoint = brownianMotionReconstructions+1;
		ending_coord = positions[endpoint];
		average_position = ending_coord;
		variance = 2*lambda*dtau*(brownianMotionReconstructions+1);
		starting_coord = generator->Gaus(average_position,sqrt(variance));
		oldposition = positions[starting_point];
		newposition = starting_coord;
	}else		{
		starting_point = timeslices-2-brownianMotionReconstructions;
		endpoint = timeslices-1;
		starting_coord = positions[starting_point];
		average_position = starting_coord;
		variance = 2*lambda*dtau*(brownianMotionReconstructions+1);
		ending_coord = generator->Gaus(average_position,sqrt(variance));
		oldposition = positions[endpoint];
		newposition = ending_coord;
	}

	double new_segment[brownianMotionReconstructions+2];
	new_segment[0] = starting_coord;
	new_segment[brownianMotionReconstructions+1] = ending_coord;
	double previous_position = starting_coord;
	double newcoordinate;

	for(int i=0; i<brownianMotionReconstructions; i++){

		left_reco = brownianMotionReconstructions-i;
		// gaussian sampling of the free particle propagator
		average_position = previous_position + (ending_coord-previous_position)/(left_reco+1);
		variance = 2*lambda*dtau*left_reco/(left_reco+1);
		newcoordinate = generator->Gaus(average_position,sqrt(variance));
		new_segment[i+1] = newcoordinate;
		previous_position = newcoordinate;
	}

	// Metropolis --> The kinetic part has been sampled exactely, only the potential part of the density matrix determines the acceptance probability of the move
	double acc_density_matrix_difference = 0;
	double newcorr,oldcorr;

	for(int i=0;i<brownianMotionReconstructions+1;i++){

		newcorr = potential_density_matrix(new_segment[i],new_segment[i+1]);
		oldcorr = potential_density_matrix(positions[starting_point+i],positions[starting_point+i+1]);
		acc_density_matrix_difference += oldcorr-newcorr;
	}

	double acceptance_probability = exp(-acc_density_matrix_difference)*variationalWaveFunction(newposition)/variationalWaveFunction(oldposition);
	if(generator->Rndm()<acceptance_probability){

		for(int i=0;i<brownianMotionReconstructions+2;i++)		positions[starting_point+i]=new_segment[i];
		acceptedBM++;
	}
}


int QMC1D :: index_mask(int ind){

	if(PIGS)	return ind;  // no pbc over indices
	else	{
		int new_ind = ind;
		while(new_ind>=timeslices)	new_ind-=timeslices;	// pbc over indices
		return new_ind;
	}
}


void QMC1D :: consoleOutput(){
	cout << "Acceptances: " << endl;
	if(PIGS)	cout << "BM: " << ((double)acceptedBM)/totalBM << endl;
	cout << "Transl: " << ((double)acceptedTranslations)/totalTranslations << endl;
	cout << "BB: " << ((double)acceptedBB)/totalBB << endl;
}


/* This function accumulates the expectation values in their respective variables. At the end of the block, these variables are divided by the MCSTEPS
	value and the block average and its squared value are accumulated in apposite variables */
void QMC1D :: upgradeAverages(){

	for(int i=0;i<timeslices;i++)	potential_energy[i]+=external_potential(positions[i]);
	
	int flag = 0;
	if(PIGS)	flag=1;
		
	for(int i=flag;i<timeslices-flag;i++){
		int i_mod = index_mask(i);
		int i_mod_next = index_mask(i+1);
		kinetic_energy[i_mod] += kineticEstimator(positions[i_mod],positions[i_mod_next]);
	}
	if(flag){
		kinetic_energy[0] += variationalLocalEnergy(positions[0]);
		kinetic_energy[timeslices-1] += variationalLocalEnergy(positions[timeslices-1]);
	}
	
	upgradeHistogram();
}


// Fill in the values of an histogram. 
void QMC1D :: upgradeHistogram(){

	double delta_pos = (histogram_end-histogram_start)/histogram_bins;
	for(int i=timeslices_averages_start; i<=timeslices_averages_end; i+=1){
		int k=1;
		while((histogram_start + k*delta_pos)<positions[i])
			k++;
		positions_histogram[k-1]+=1;
	}
}


void QMC1D :: endBlock(){  // Calculating and accumulating block averages

	for(int i=0;i<timeslices;i++){
		potential_energy[i] /= MCSTEPS;
		potential_energy_accumulator[i] += potential_energy[i];
		potential_energy_square_accumulator[i] += potential_energy[i]*potential_energy[i];
		potential_energy[i] = 0;
		kinetic_energy[i] /= MCSTEPS;
		kinetic_energy_accumulator[i] += kinetic_energy[i];
		kinetic_energy_square_accumulator[i] += kinetic_energy[i]*kinetic_energy[i];
		kinetic_energy[i] = 0 ;
	}
	
	for(int i=0;i<histogram_bins;i++){
		positions_histogram[i] /= MCSTEPS;
		positions_histogram_accumulator[i] += positions_histogram[i];
		positions_histogram_square_accumulator[i] += positions_histogram[i]*positions_histogram[i];
		positions_histogram[i] = 0;
	}
}

/* finalize**** functions average the accumulators of the block averages and their square values.
Then the error is calculated by the usual formula err(A)=sqrt(|<A><A>-<A*A>|/Nblocks) */
void QMC1D :: finalizePotentialEstimator(){
	
	ofstream out("potential.dat");
	double potential_energy_average, potential_energy_square_avg, p_error;

	for(int i=0;i<timeslices;i++){
		potential_energy_average = potential_energy_accumulator[i]/blocks;
		potential_energy_square_avg = potential_energy_square_accumulator[i]/blocks;
		p_error = sqrt( abs( pow(potential_energy_average,2) - potential_energy_square_avg )/blocks );
		out << i << " " << potential_energy_average << " " << p_error << endl;
	}
	out.close();
}

void QMC1D :: finalizeKineticEstimator(){

	ofstream out("kinetic.dat");
	double kinetic_energy_average, kinetic_energy_square_avg, k_error;

	for(int i=0;i<timeslices;i++){
		kinetic_energy_average = kinetic_energy_accumulator[i]/blocks;
		kinetic_energy_square_avg = kinetic_energy_square_accumulator[i]/blocks;
		k_error = sqrt( abs( pow(kinetic_energy_average,2) - kinetic_energy_square_avg )/blocks );
		out << i << " " << kinetic_energy_average << " " << k_error <<endl;
	}
	out.close();
}

void QMC1D :: finalizeHistogram(){

	ofstream out("probability.dat");
	double current_position, hist_average, hist_square_avg, error;
	double delta_pos = (histogram_end-histogram_start)/histogram_bins;
	double norma = 0.0;
	for(int i=0; i<histogram_bins; i++)		norma += positions_histogram_accumulator[i]/blocks;

	norma *= delta_pos;

	for(int i=0; i<histogram_bins; i++){
		current_position = histogram_start + (i+0.5)*delta_pos;
		hist_average = positions_histogram_accumulator[i]/blocks;
		hist_square_avg = positions_histogram_square_accumulator[i]/blocks;
		error = sqrt(abs(hist_average*hist_average-hist_square_avg)/blocks);
		out << current_position << " " << hist_average/norma << " " << error/norma << endl;
                
        }
	out.close();
}

// (-hbar*hbar/2m)d^2/dx^2G(x,x',dtau)
double QMC1D :: kineticEstimator(double value,double next_value){

	double term_1 = (dtau/2)*external_potential_prime(value) + (value-next_value)/(2*lambda*dtau);
	double term_2 = (dtau/2)*external_potential_second(value) + 1./(2*lambda*dtau);
	return -(hbar*hbar/(2*mass))*(term_1*term_1 - term_2);
}

// (-hbar*hbar/2m)(d^2/dx^2G(x,x',dtau))/G(x,x',dtau)
double QMC1D :: variationalLocalEnergy(double val){

	return -(hbar*hbar/(2*mass))*variationalWaveFunction_second(val)/variationalWaveFunction(val);
}


void QMC1D :: readInput(int argc, char ** argv){

	ifstream input_file;

    if(argc == 1)	input_file.open("input.dat");
    else	{
		if(atoi(argv[1]) == 0){
            input_file.open("input.pigs");
            cout << endl << "---- PIGS ----" << endl << endl;
        }else	if(atoi(argv[1]) == 1){
					input_file.open("input.pimc");
					cout << endl << "---- PIMC ----" << endl << endl;
        }
        
        if(argc >= 3){
			func_type = atoi(argv[2]);
			if(func_type == 0)	cout << "Constant wavefunction" << endl << endl;
			else				cout << "Trial wavefunction" << endl << endl;		
        }

        if(argc >= 5){
			sigma = atof(argv[3]);
			mu = atof(argv[4]);
			cout << "Mu = " << mu << "   Sigma = " << sigma << endl << endl;
        }
    }

	char* string_away = new char[60];

	input_file >> string_away >> timeslices;
	input_file >> string_away >> temperature;
	input_file >> string_away >> imaginaryTimePropagation;
	input_file >> string_away >> brownianMotionReconstructions;
	input_file >> string_away >> delta_translation;
	input_file >> string_away >> brownianBridgeReconstructions;
	input_file >> string_away >> brownianBridgeAttempts;
	input_file >> string_away >> MCSTEPS;
	input_file >> string_away >> equilibration;
	input_file >> string_away >> blocks;
	input_file >> string_away >> histogram_bins;
	input_file >> string_away >> histogram_start;
	input_file >> string_away >> histogram_end;
	input_file >> string_away >> timeslices_averages_start>>timeslices_averages_end;
	input_file.close();

	if(argc >= 6){
		imaginaryTimePropagation = atof(argv[5]);
		cout << "Imaginary time = " << imaginaryTimePropagation << endl << endl;
	}

	delete [] string_away;
}


void QMC1D :: deleteMemory(){

	delete [] positions;
	delete [] potential_energy;
	delete [] potential_energy_accumulator;
	delete [] potential_energy_square_accumulator;
                                                                                                                 
	delete [] kinetic_energy;
	delete [] kinetic_energy_accumulator;
	delete [] kinetic_energy_square_accumulator;
                                                                                                                 
	delete [] positions_histogram;
	delete [] positions_histogram_accumulator;
	delete [] positions_histogram_square_accumulator;

	delete generator;
}
