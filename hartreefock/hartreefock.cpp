/********************************************************************
	This file calculates the final eigen energies for an input
		single particle file and two body interaction, V, file.
		The input is assumed to be COUPLED with an optional
		conversion to UNCOUPLED format.

	KEY: 	c/unc 	== 	coupled/uncoupled
			sp 		== 	single particle
			tb		== 	two body
*********************************************************************/

#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<string>
#include<sstream>
#include "armadillo"
#include "time.h"
#include "memory.h"
#include "hartreefock.h"

#define TOLERANCE 1e-8
#define maxITERATIONS 30

int main(int argc, char* argv[]){
	//Declare basic variables
	char* spfilename, *tbfilename;
	int i, j, k, l, 
		tz, size, 
		iterations=0;
	bool tb;
	clock_t start, finish;
	double time, rho, energy;

	//Verify proper input
	if(argc<5){
		std::cout << "Bad usage. Enter also:	infilename_sp infilename_tb tz bool_tb" << std::endl;
		exit(1);
	}
	else{
		spfilename = argv[1];
		tbfilename = argv[2];
		tz=atoi(argv[3]);
		tb = atoi(argv[4]);
				
		if(tz==0){
			size = 80;
		}
		else if(abs(tz)==1){
			size = 40;
		}
	}

	/********************************************************************
		Declare required matrices and vectors

			mat: 	states 		--> 	lists state information (n,l,2j,etc)
					H 			--> 	hamiltonian
					cgcoeff		-->		Clebsch-Gordan coefficient

			vec:	E 			-->		final eigen energies
					mask		-->		translates "i" to state label
											from file
	********************************************************************/
	arma::mat 	states_c = arma::zeros<arma::mat>(80,5),
				states_unc = states_c,
				H_c = arma::zeros<arma::mat>(size,size),
				H0 = H_c,
				H_unc = H0,
				cgcoeff = H_unc,
				C = arma::eye<arma::mat>(size,size);
	arma::vec	E_c = arma::zeros<arma::vec>(size),
				E_unc = E_c,
				mask = E_unc;
	matrix4D<double> V(80,80,80,80,0.0);

	/********************************************************************
		Extract data from files
		Create state number mask for use if abs(tz)==1
		Initialize coupled single particle H0
	********************************************************************/
	extract_from_spfile(spfilename, states_unc);
	if(tb==true){
		extract_from_tbfile(tbfilename,V);
	}
	create_mask(states_unc,tz,mask);
	sp_energies(states_unc,tz,H0);

	/********************************************************************
		Run solver
		Prints iteration count, time, and diff.max() < TOLERANCE to screen
		H0, V, tz, mask are inputs. Returns H_unc, E_unc
	********************************************************************/
	solve_iterations(H0,V,tz,mask,H_unc,E_unc);

	/*************************************
		Print eigen energies to screen
	*************************************/
	for(i=0;i<size;i++){
		std::cout << std::setw(10) << E_unc(i);
		if((i+1)%10==0){
			std::cout << std::endl;
		}
	}
	std::cout << std::endl;

	/********************************************************************
		Calculates uncoupled states 
	********************************************************************/
/*
	states_uncoupled_from_coupled(states_c,states_unc);
	CG_coeff(states_c,states_unc,tz,cgcoeff);
	H_uncoupled_from_coupled(H_c,cgcoeff,tz,H_unc);
	arma::eig_sym(E_unc,C,H_unc);
	E_unc = arma::sort(E_unc);
*/
	return 0;
}

