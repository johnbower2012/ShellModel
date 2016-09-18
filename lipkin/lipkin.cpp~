#include<iomanip>
#include<iostream>
#include "count.h"
#include "memory.h"
#include "stateset_lipkin.h"
#include "matrix_solver.h"

int main(int argc, char*argv[]){
	int i, j, k, l,
		p, s, n, m, omega,
		coun, selection,
		Jpm, Jmp, Jpp, Jmm, Jz, N;
	count number;
	operators oper;

	/***  DECLARE CONSTANTS  ***/
	double 	eps = 2.0,
			V = -1.0/3.0,
			W = -1.0/4.0,
			tolerance = 1e-8;

	if(argc<4){
		std::cout << "Enter also 'p n m'." << std::endl;
		exit(1);
	}
	else{
		p = atoi(argv[1]);
		n = atoi(argv[2]);
		m = atoi(argv[3]);
		number.choose(n,m);
		s = pow(p,m);
//		s = number.value();

		omega = m;
	}

	array<int> 		list(p,0);
	matrix<double> 	hamil(s,s,0), 
					vec(s,s,0), 
					sort(s,2,0);
	for(i=0;i<s;i++){
		vec.memory[i][i] = 1.0;
	}
	for(i=0;i<p;i++){
		list.memory[i] = m;
	}

	stateset family(p,list.memory,n,m,s);

	oper.reset(n,0,family.states.memory[0]);
	for(i=0;i<s;i++){
		std::cout << i << std::setw(15) << family.mask.memory[i] << std::setw(15);
		Jz = 0.0; Jpm = 0.0;
		for(j=0;j<m;j++){
			for(k=0;k<p;k++){
				oper.update(i,family.states.memory[i]);
				Jz += oper.J_z(family.states.memory[i],j,k,p,omega);
			}
			for(k=0;k<m;k++){
				oper.update(i,family.states.memory[i]);
				Jpm += oper.J_plus_J_minus(family.states.memory[i],j,k,omega);
			}
		}
		Jz /= 2.0;
		std::cout << Jpm + Jz*(Jz - 1.0) << std::setw(15) << family.zedangmom(i)/2.0 << std::setw(15);
		for(j=0;j<n;j++){
			family.print(i,j);
		}
		std::cout << std::endl;
	}


	oper.reset(n,0,family.states.memory[0]);
	for(i=0;i<s;i++){
		for(j=0;j<s;j++){
			Jz = 0.0; 
			Jpm = 0.0; Jmp = 0.0;
			Jpp = 0.0; Jmm = 0.0;
			for(k=0;k<m;k++){
				oper.update(j,family.states.memory[j]);
				N = oper.N_ac(family.states.memory[i]);
				for(l=0;l<p;l++){
					oper.update(j,family.states.memory[j]);
					Jz += oper.J_z(family.states.memory[i],k,l,p,omega);
				}
				for(l=0;l<m;l++){
					oper.update(j,family.states.memory[j]);
					Jpm += oper.J_plus_J_minus(family.states.memory[i],k,l,omega);
				}
				for(l=0;l<m;l++){
					oper.update(j,family.states.memory[j]);
					Jmp += oper.J_minus_J_plus(family.states.memory[i],k,l,omega);
				}
				for(l=0;l<m;l++){
					oper.update(j,family.states.memory[j]);
					Jpp += oper.J_plus_J_plus(family.states.memory[i],k,l,omega);
				}
				for(l=0;l<m;l++){
					oper.update(j,family.states.memory[j]);
					Jmm += oper.J_minus_J_minus(family.states.memory[i],k,l,omega);
				}
			}
			hamil.memory[i][j] = eps*Jz/2.0 + 0.5*V*(Jpp + Jmm) + 0.5*W*(Jpm + Jmp - N);
		}
	}
	for(i=0;i<s;i++){
		for(j=0;j<s;j++){
			std::cout << std::setw(10) << hamil.memory[i][j];
		}
		std::cout << std::endl;
	}

	jacobi_simrot_eigen_solver(hamil.memory,vec.memory,s,tolerance,coun);
	matrix_diag_sort(hamil.memory,sort.memory,s);

	for(i=0;i<s;i++){
		std::cout << std::setw(10) << family.zedangmom(i)/2.0 << std::setw(15) << sort.memory[i][0] << std::endl;
	}

 	return 0;
}
