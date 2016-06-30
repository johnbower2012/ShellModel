#include<iostream>
#include<iomanip>
#include<cmath>
#include "armadillo"
#include "time.h"
#include "memory.h"

int** basis_uncoupled(int,int);
int** basis_coupled(int,int);
void basis_print(int**,int);
bool triangle_rule(int,int,int);
int factorial(int);
int factorial_half(int);
double cg_coeff(int,int,int,int,int,int);
void extract_from_spfile(char* filename, arma::mat& states);
void extract_from_tbfile(char*,matrix4D<double>&);
void sp_energies(arma::mat& input,int,arma::mat&);
void create_mask(arma::mat&,int,arma::vec&);
void states_uncoupled_from_coupled(arma::mat&,arma::mat&);
void V_uncoupled_from_coupled(matrix4D<double>&,arma::mat&,int,matrix4D<double>&);
void CG_coeff(arma::mat&,arma::mat&,int,arma::mat&);
void H_uncoupled_from_coupled(arma::mat&,arma::mat&,int,arma::mat&);

int main(int argc, char* argv[]){
	int a, b, c, d,
		i, j, k, l,
		A, B, C, D,
		I, J, K, L, 
		tz, size;
	double energy, V, time;
	char* spinfile, *tbinfile;
	clock_t start, finish;
	
	if(argc<4){
		std::cout << "Bad usage. Enter also: 	sbfilename tbfilename tz" << std::endl;
		exit(1);
	}
	else{
		spinfile = argv[1];
		tbinfile = argv[2];
		tz = atoi(argv[3]);
		size = 80 - fabs(tz)*40;
	}

	arma::mat 	unc_states = arma::zeros<arma::mat>(80,5),
				c_states = unc_states,
				cgcoeff = arma::zeros<arma::mat>(size,size),
				H_unc = cgcoeff,
				H_c = H_unc,
				C_ = arma::eye<arma::mat>(size,size);
	arma::vec	unc_n_mask = arma::zeros<arma::vec>(size),
				unc_p_mask = unc_n_mask,
				c_n_mask = unc_p_mask,
				c_p_mask = c_n_mask,
				spE = arma::zeros<arma::vec>(size),
				E = spE;
	matrix4D<double> 	V_c(80,80,80,80,0.0),
						V_unc(80,80,80,80,0.0);
	
	extract_from_spfile(spinfile,c_states);
	extract_from_tbfile(tbinfile,V_c);
	sp_energies(c_states,tz,H_c);
	states_uncoupled_from_coupled(unc_states,c_states);

/*
	if(fabs(tz)==1){
		create_mask(c_states,tz,c_n_mask);
		create_mask(unc_states,tz,unc_n_mask);
		for(a=0;a<size;a++){
			A = unc_n_mask(a);
			for(b=0;b<size;b++){
				B = unc_n_mask(b);
				for(c=0;c<size;c++){
					C = unc_n_mask(c);
					for(d=0;d<size;d++){
						D = unc_n_mask(d);
						V = 0.0;
	start = clock();
						for(i=0;i<size;i++){
							I = c_n_mask(i);
							for(j=0;j<size;j++){
								J = c_n_mask(j);
								for(k=0;k<size;k++){
									K = c_n_mask(k);
									for(l=0;l<size;l++){
										L = c_n_mask(l);
										V += cgcoeff(i,a)*cgcoeff(j,b)*cgcoeff(k,c)*cgcoeff(l,d)*V_c.memory[I][J][K][L];
									}
								}
							}
						}
						V_unc.memory[A][B][C][D] = V;
	finish = clock();
	time = (finish - start)/(double) CLOCKS_PER_SEC;
	std::cout << time << " seconds" << std::endl;
					}
				}
			}
		}
	}
	else if(tz==0){
		for(a=0;a<size;a++){
			for(b=0;b<size;c++){
				for(c=0;c<size;c++){
					for(d=0;d<size;d++){
						V = 0.0;
						for(i=0;i<size;i++){
							for(j=0;j<size;j++){
								for(k=0;k<size;k++){
									for(l=0;l<size;l++){
										V += cgcoeff(i,a)*cgcoeff(j,b)*cgcoeff(k,c)*cgcoeff(l,d)*V_c.memory[i][j][k][l];
									}
								}
							}
						}
						V_unc.memory[a][b][c][d] = V;
					}
				}
			}
		}
	}
*/

	CG_coeff(c_states,unc_states,tz,cgcoeff);
	H_uncoupled_from_coupled(H_c,cgcoeff,tz,H_unc);
	arma::eig_sym(E,C_,H_unc);

	return 0;
}

/**************************************************************************************************

	BEGIN FUNCTION DEFINITIONS

**************************************************************************************************/

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
void states_uncoupled_from_coupled(arma::mat& uncoupled, arma::mat& coupled){
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
void V_uncoupled_from_coupled(matrix4D<double>&,arma::mat&,int,matrix4D<double>&){}
void create_mask(arma::mat& states, int tz, arma::vec& mask){
	int i, j = 0;
	for(i=0;i<80;i++){
		if(states(i,4)==tz){
			mask[j] = i;
			j++;
		}
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
		for(j=0;j<size;j++){
			energy = 0.0;
			for(k=0;k<size;k++){
				for(l=0;l<size;l++){
					energy += cgcoeff(k,i)*cgcoeff(l,j)*H_c(k,l);
				}
			}
			H_unc(i,j) = energy;
		}
	}
}
	


