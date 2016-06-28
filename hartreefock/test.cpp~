#include<iostream>
#include<iomanip>
#include<cmath>
#include "armadillo"
#include "memory.h"

int** basis_uncoupled(int,int);
int** basis_coupled(int,int);
void basis_print(int**,int);
bool triangle_rule(int,int,int);
int factorial(int);
double cg_coeff(int,int,int,int,int,int);
void extract_from_spfile(char* filename, arma::mat& states);
void create_uncoupled_from_coupled(arma::mat&, arma::mat&);

int main(int argc, char* argv[]){
	int i, j, j2, prevj2=-1, tz=-1, marker=0, size;
	char* infile, *outfile;
	
	if(argc<3){
		std::cout << "Bad usage. Enter also: 	coupledinfile uncoupledoutfile" << std::endl;
		exit(1);
	}
	else{
		infile = argv[1];
		outfile = argv[2];
	}

	arma::mat 	unc_states = arma::zeros<arma::mat>(80,5),
				c_states = unc_states;

	extract_from_spfile(infile,c_states);
	create_uncoupled_from_coupled(unc_states,c_states);
	for(i=0;i<80;i++){
		for(j=0;j<5;j++){
			std::cout << std::setw(5) << unc_states(i,j);
		}
		std::cout << std::endl;
	}

/*
	for(i=0;i<80;i++){
		unc_states(i,0) = c_states(i,0);
		unc_states(i,1) = j2 = 2*c_states(i,1);
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
		unc_states(i,2) = j2 - marker;
		unc_states(i,3) = 1 - 2*(i%2);
		unc_states(i,4) = tz;
		prevj2 = j2;
	}
*/
	
	return 0;
}

//////////////////////////////////////////////////////////////////////////////////

int** basis_uncoupled(int j1, int j2){
	int i, mj1, mj2, J2, size,
		**uncoupled;

	mj1 = j1; J2 = j2 + 1; 
	size = (j1+1)*(j2+1);
	uncoupled = new int*[size];

	for(i=0;i<size;i++){
		mj2 = j2 - (i%J2)*2;
		uncoupled[i] = new int[4];
		uncoupled[i][0] = j1;
		uncoupled[i][1] = j2;
		uncoupled[i][2] = mj1;
		uncoupled[i][3] = mj2;
		if(mj2==-j2){
			mj1 -= 2;
		}
	}
	return uncoupled;
}

int** basis_coupled(int j1, int j2){
	int i, J, M, size,
		**coupled;
 
	size = (j1+1)*(j2+1);
	coupled = new int*[size];

	J = j1 + j2; M = J;
	for(i=0;i<size;i++){
		coupled[i] = new int[4];
		coupled[i][0] = j1;
		coupled[i][1] = j2;
		coupled[i][2] = J;
		coupled[i][3] = M;
		M -= 2;
		if(M==-(J+2)){
			J -= 2;
			M = J;
		}
	}
	return coupled;
}

void basis_print(int**memory, int size){
	int i, j;
	for(i=0;i<size;i++){
		for(j=0;j<4;j++){
			std::cout << std::setw(5) << memory[i][j];
		}
		std::cout << std::endl;
	}
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




