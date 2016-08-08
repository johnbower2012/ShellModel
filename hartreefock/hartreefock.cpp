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

#define TOLERANCE 1e-8
#define maxITERATIONS 30

void extract_from_spfile(char*,arma::mat&);
void extract_from_tbfile(char*,matrix4D<double>&);
void create_uncoupled_from_coupled(arma::mat&,arma::mat&);
void create_mask(arma::mat&,int,arma::vec&);
void sp_energies(arma::mat&,int,arma::mat&);
void compute_densityMatrix(arma::mat&,int,arma::mat&);
void solve_iterations(arma::mat&,matrix4D<double>&,int,arma::vec&,arma::mat&,arma::vec&);
bool triangle_rule(int,int,int);
int factorial(int);
int factorial_half(int);
double cg_coeff(int,int,int,int,int,int);
void states_uncoupled_from_coupled(arma::mat&,arma::mat&);
void CG_coeff(arma::mat&,arma::mat&,int,arma::mat&);
void H_uncoupled_from_coupled(arma::mat&,arma::mat&,int,arma::mat&);

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
	extract_from_spfile(spfilename, states_c);
	if(tb==true){
		extract_from_tbfile(tbfilename,V);
	}
	create_mask(states_c,tz,mask);
	sp_energies(states_c,tz,H0);


	/********************************************************************
		Run solver
		Prints iteration count, time, and diff.max() < TOLERANCE to screen
		H0, V, tz, mask are inputs. Returns H_c, E_c
	********************************************************************/
	solve_iterations(H0,V,tz,mask,H_c,E_c);

	/********************************************************************
		Print eigen energies to screen
	********************************************************************/
	for(i=0;i<size;i++){
		std::cout << std::setw(10) << E_c(i);
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

/**************************************************************************************************

	BEGIN FUNCTION DEFINITIONS

**************************************************************************************************/

void extract_from_spfile(char* filename, arma::mat& states){
	int i, num;
	std::string empty, words;
	std::ifstream spfile;
	spfile.open(filename);
	for(i=0;i<6;i++){
		spfile >> empty;
	}
	for(i=0;i<80;i++){
		spfile >> empty; spfile >> words;
		spfile >> num;
		spfile >> states(i,0);
		spfile >> states(i,1);
		spfile >> states(i,2);
		spfile >> states(i,3);
		spfile >> states(i,4);
	}
	spfile.close();
}
void extract_from_tbfile(char* filename, matrix4D<double>& V){
	int a,b,c,d;
	std::ifstream tbfile;
	tbfile.open(filename);
	do{
		tbfile >> a; tbfile >> b;
		tbfile >> c; tbfile >> d;
		tbfile >> V.memory[a-1][b-1][c-1][d-1];
	} while (!tbfile.eof());
	tbfile.close();
}
void create_uncoupled_from_coupled(arma::mat& uncoupled, arma::mat& coupled){
	int i, j2, marker=0, prevj2=-1, tz = 1;

	for(i=0;i<80;i++){
		uncoupled(i,0) = coupled(i,0);
		uncoupled(i,1) = j2 = 2*coupled(i,1);
		if(prevj2==j2){
			if(i%2==0){
				marker += 2;
				if(marker>j2*2+1){
					tz *= -1;
					marker = 0;
				}
			}
		}
		else{
			marker = 0;
			tz *= -1;
		}
		uncoupled(i,2) = j2 - marker;
		uncoupled(i,3) = 1 - 2*(i%2);
		uncoupled(i,4) = tz;
		prevj2 = j2;
	}
}
void sp_energies(arma::mat& input, int tz, arma::mat& H0){
	int i, j;
	if(tz==0){
		for(i=0;i<80;i++){
			H0(i,i) = 10*(2*input(i,0) + input(i,1) + 1.5);
		}
	}
	else if(tz==1||tz==-1){
		j=0;
		for(i=0;i<80;i++){
			if(input(i,4)==tz){
				H0(j,j) = 10*(2*input(i,0) + input(i,1) + 1.5);
				j++;
			}
		}
	}
	else{
		std::cout << "Improper use of 'sp_energies':	tz must be 1 (n), 0 (b), -1 (p)" << std::endl;
		exit(1);
	}
}
void create_mask(arma::mat& states, int tz, arma::vec& mask){
	int i, j = 0;
	for(i=0;i<80;i++){
		if(states(i,4)==tz){
			mask[j] = i;
			j++;
		}
	}
}
void compute_densityMatrix(arma::mat& C, int size, arma::mat& densityMatrix){
	int i, j, k, part;
	double rho;	
	part = size/5;
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			rho = 0;
			for(k=0;k<part;k++){
				rho += C(i,k)*C(j,k);
			}
			densityMatrix(i,j) = rho;
		}
	}
}
void solve_iterations(arma::mat& H0, matrix4D<double>& V, int tz, arma::vec& mask, arma::mat& H, arma::vec& E){	
	//Run solver for tz==+-1 or 0
	//Note only the n/p only case uses the mask
	//There is no need for it in the total case
	int i, j, k, l,
		I, J, K, L,
		iterations = 0,
		size = 80 - fabs(tz)*40;
	double energy, time;

	arma::mat 	C = arma::eye<arma::mat>(size,size),
				densityMatrix = C;
	arma::vec	prevE = arma::zeros<arma::vec>(size),
				diff = prevE;
	
	clock_t start, finish;
	start = clock();
	if(abs(tz)==1){		
		while(iterations < maxITERATIONS){
			H = arma::zeros<arma::mat>(size,size);
			compute_densityMatrix(C,size,densityMatrix);
			for(i=0;i<size;i++){
				I = mask(i);
				for(j=i;j<size;j++){
					J = mask(j);
					energy = 0;
					for(k=0;k<size;k++){
						K = mask(k);
						for(l=0;l<size;l++){
							L = mask(l);
							energy += densityMatrix(k,l)*V.memory[I][K][J][L];
						}
					}
					H(i,j) = H(j,i) = energy + H0(i,j);
				}
			}
			arma::eig_sym(E,C,H);
			diff = prevE - E;
			prevE = E;
			if(fabs(diff.max()) < TOLERANCE){
				break;
			}
			iterations++;
		}
	}
	else if(tz==0){
		while(iterations < maxITERATIONS){
			H = arma::zeros<arma::mat>(size,size);
			diff = arma::zeros<arma::vec>(size);
			compute_densityMatrix(C,size,densityMatrix);
			for(i=0;i<size;i++){
				for(j=i;j<size;j++){
					energy = 0;
					for(k=0;k<size;k++){
						for(l=0;l<size;l++){
							energy += densityMatrix(k,l)*V.memory[i][k][j][l];
						}
					}
					H(i,j) = H(j,i) = energy + H0(i,j);
				}
			}
			arma::eig_sym(E,C,H);
			diff = prevE - E;
			prevE = E;
			if(fabs(diff.max()) < TOLERANCE){
				break;
			}
			iterations++;
		}
	}
	finish = clock();
	time = (finish - start)/(double) CLOCKS_PER_SEC;
	std::cout << std::endl;
	std::cout << "Convergence in " << time << " seconds after " << iterations;
	std::cout << " iterations with " << fabs(diff.max()) << " < " << TOLERANCE << std::endl;
	std::cout << std::endl;
	
	E = arma::sort(E);
}
bool triangle_rule(int j1, int j2, int j3){
	bool test = true;
	if((j1+j2+j3)%2==1){
		test = false;
	}
	else if((j1+j2)<j3){
		test = false;
	}
	else if((j1+j3)<j2){
		test = false;
	}
	else if((j2+j3)<j1){
		test = false;
	}
	return test;
}

int factorial(int n){
	int i, fac = 1;
	if(n<0){
		fac = 0;
	}
	else{
		for(i=1;i<n+1;i++){
			fac *= i;
		}
	}
	return fac;
}
int factorial_half(int n){
	int i, fac = 1;
	n /= 2;
	if(n<0){
		fac = 0;
	}
	else{
		for(i=1;i<n+1;i++){
			fac *= i;
		}
	}
	return fac;
}
double cg_coeff(int j1, int j2, int m1, int m2, int j, int m){
	int k, max, min;
	double fac, fac_k, fac1, fac2, fac3, fac4, fac5, fac6,
			coeff = 0.0;
	if(triangle_rule(j1,j2,j)&&(m1+m2==m)){
		max = j1 + j2 - j;
		min = 0;
		if( (j1-m1) < max){
			max = j1 - m1;
		}
		if( (j2+m2) < max){
			max = j2 + m2;
		}
		if( (-j + j2 - m1) > min){
			min = -j + j2 - m1;
		}
		if( (-j + j1 + m2) > min){
			min = -j + j1 + m2;
		}
		for(k=min;k<max+1;k+=2){
			fac_k = factorial_half(k);
			fac1 = factorial_half(j1+j2-j-k);
			fac2 = factorial_half(j1-m1-k);
			fac3 = factorial_half(j2+m2-k);
			fac4 = factorial_half(j-j2+m1+k);
			fac5 = factorial_half(j-j1-m2+k);
			fac = fac_k*fac1*fac2*fac3*fac4*fac5;
			coeff += pow(-1,k/2)/fac;
		}

		fac1 = factorial_half(j+m); fac2 = factorial_half(j-m);
		fac3 = factorial_half(j1-m1); fac4 = factorial_half(j1+m1);
		fac5 = factorial_half(j2-m2); fac6 = factorial_half(j2+m2);
		fac = fac1*fac2*fac3*fac4*fac5*fac6;		
		coeff *= sqrt(fac);
		
		fac1 = factorial_half(j+j1-j2); fac2 = factorial_half(j-j1+j2);
		fac3 = factorial_half(j1+j2-j); fac4 = factorial_half(j1+j2+j+2);
		fac = fac1*fac2*fac3/fac4;
		coeff *= sqrt(((double) j+1.0)*fac);
		
	}

	return coeff;
}
void states_uncoupled_from_coupled(arma::mat& coupled, arma::mat& uncoupled){
	int i, j2, marker=0, prevj2=-1, tz = 1;

	for(i=0;i<80;i++){
		uncoupled(i,0) = coupled(i,0);
		uncoupled(i,1) = j2 = 2*coupled(i,1);
		if(prevj2==j2){
			if(i%2==0){
				marker += 2;
				if(marker>j2*2+1){
					tz *= -1;
					marker = 0;
				}
			}
		}
		else{
			marker = 0;
			tz *= -1;
		}
		uncoupled(i,2) = j2 - marker;
		uncoupled(i,3) = 1 - 2*(i%2);
		uncoupled(i,4) = tz;
		prevj2 = j2;
	}
}
void CG_coeff(arma::mat& c_states ,arma::mat& unc_states, int tz, arma::mat& cgcoeff){
	int i, k, nc, nunc, jc, junc,
		I, K,
		m1, m2, J, M,
		size = 40;
	arma::vec	c_mask = arma::zeros<arma::vec>(size),
				unc_mask = c_mask;
	if(fabs(tz)==1){
		create_mask(c_states,tz,c_mask);
		create_mask(unc_states,tz,unc_mask);
		for(i=0;i<size;i++){
			I = c_mask(i);
			for(k=0;k<size;k++){
				K = unc_mask(k);
				nc = c_states(I,0);
				nunc = unc_states(K,0);
				if(nc==nunc){
					jc = 2*c_states(I,1);
					junc = unc_states(K,1);
					if(jc==junc){
						J = c_states(I,2);
						M = c_states(I,3);
						m1 = unc_states(K,2);
						m2 = unc_states(K,3);
						cgcoeff(i,k) = cg_coeff(jc,1,m1,m2,J,M);
					}
				}
			}
		}
	}
	else if(tz==0){
		for(tz=-1;tz<2;tz+=2){
			create_mask(c_states,tz,c_mask);
			create_mask(unc_states,tz,unc_mask);
			for(i=0;i<size;i++){
				I = c_mask(i);
				for(k=0;k<size;k++){
					K = unc_mask(k);
					nc = c_states(I,0);
					nunc = unc_states(K,0);
					if(nc==nunc){
						jc = 2*c_states(I,1);
						junc = unc_states(K,1);
						if(jc==junc){
							J = c_states(I,2);
							M = c_states(I,3);
							m1 = unc_states(K,2);
							m2 = unc_states(K,3);
							cgcoeff(I,K) = cg_coeff(jc,1,m1,m2,J,M);
						}
					}
				}
			}
		}
	}
}
void H_uncoupled_from_coupled(arma::mat& H_c,arma::mat& cgcoeff,int tz,arma::mat& H_unc){
	int i, j, k, l,
		size = 80 - fabs(tz)*40;
	double energy;
	H_unc = arma::zeros<arma::mat>(size,size);
	for(i=0;i<size;i++){
		for(j=i;j<size;j++){
			energy = 0.0;
			for(k=0;k<size;k++){
				for(l=0;l<size;l++){
					energy += cgcoeff(k,i)*cgcoeff(l,j)*H_c(k,l);
				}
			}
			H_unc(i,j) = H_unc(j,i) = energy;
		}
	}
}
