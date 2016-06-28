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
void create_mask(arma::mat&,int,arma::vec&);
void sp_energies(arma::mat&,int,arma::vec&);
void compute_densityMatrix(arma::mat&,int,arma::mat&);

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

	//Declare required matrices and vectors
	arma::mat 	states = arma::zeros<arma::mat>(80,5),
				H = arma::zeros<arma::mat>(size,size),
				C = arma::eye<arma::mat>(size,size),
				densityMatrix = C;
	arma::vec	spE = arma::zeros<arma::vec>(size),
				E = spE,
				prevE = E,
				mask = prevE,
				diff = arma::ones<arma::vec>(size);
	matrix4D<double> V(80,80,80,80,0.0);

	//Extract data from files
	//Set up array with sp energies
	extract_from_spfile(spfilename, states);
	if(tb==true){
		extract_from_tbfile(tbfilename,V);
	}
	create_mask(states,tz,mask);
	sp_energies(states,tz,spE);
	prevE = spE;
 

	//Run solver for tz==+-1 or 0
	//Note only the n/p only case uses the mask
	//There is no need for it in the total case
	start = clock();
	if(abs(tz)==1){
		while(iterations < maxITERATIONS){
			H = arma::zeros<arma::mat>(size,size);
			compute_densityMatrix(C,size,densityMatrix);
			for(i=0;i<size;i++){
				for(j=i;j<size;j++){
					energy = 0;
					for(k=0;k<size;k++){
						for(l=0;l<size;l++){
							energy += densityMatrix(k,l)*V.memory[(int)mask(i)][(int)mask(k)][(int)mask(j)][(int)mask(l)];
						}
					}
					H(i,j) = energy;
				}
				H(i,i) += spE(i);
			}
			arma::eig_sym(E,C,H);
			diff = prevE - E;
			prevE = E;
			iterations++;
			if(fabs(diff.max()) < TOLERANCE){
				break;
			}
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
					H(i,j) = H(j,i) = energy;
				}
				H(i,i) += spE(i);
			}
			arma::eig_sym(E,C,H);
			diff = prevE - E;
			prevE = E;
			iterations++;
			if(fabs(diff.max()) < TOLERANCE){
				break;
			}
		}
	}
	finish = clock();
	time = (finish - start)/(double) CLOCKS_PER_SEC;
	
	E = arma::sort(E);

	std::cout << std::endl;
	std::cout << "Convergence in " << time << " seconds after " << iterations << " iterations with " << fabs(diff.max()) << " < " << TOLERANCE << std::endl;
	std::cout << std::endl;
	for(i=0;i<size;i++){
		std::cout << std::setw(10) << E(i);
		if((i+1)%10==0){
			std::cout << std::endl;
		}
	}
	std::cout << std::endl;
	
	return 0;
}

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
void sp_energies(arma::mat& input, int tz, arma::vec& spE){
	int i, j;

	if(tz==0){
		for(i=0;i<80;i++){
			spE(i) = 10*(2*input(i,0) + input(i,1) + 1.5);
		}
	}
	else if(tz==1||tz==-1){
		j=0;
		for(i=0;i<80;i++){
			if(input(i,4)==tz){
				spE(j) = 10*(2*input(i,0) + input(i,1) + 1.5);
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
//	arma::mat C_dagger = arma::trans(C);
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

