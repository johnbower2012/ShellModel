#include<iostream>
#include<iomanip>
#include "stateset.h"
#include "time.h"

int main(int argc, char* argv[]){
	int i,j,k,n,m,s,omega,p;
	count z;
	int __builtin_popcount (unsigned int);
	bool copy;
	clock_t start, finish;
	double time;

	if(argc<3){
		std::cout << "Bad usage. Include 'p omega m'." << std::endl;
		exit(1);
	}
	else{
		p = atoi(argv[1]);
		omega = atoi(argv[2]);
		m = atoi(argv[3]);
		n = p*omega;
	}

	z.choose(n,m);
	s = z.value();

	stateset family(n,m,s);

	unsigned int twod;

	start = clock();
	for(i=0;i<s;i++){
		for(k=0;k<s;k++){
			twod = 0;
			for(j=0;j<n;j++){
				twod += __builtin_popcount(family.states.memory[i][j] ^ family.states.memory[k][j]);
			}
		}
	}
	finish = clock();
	time = (finish - start)/(double) CLOCKS_PER_SEC;
	std::cout << time << std::endl;

	start = clock();
	for(i=0;i<s;i++){
		for(k=0;k<s;k++){
			twod = 0;
			for(j=0;j<n;j++){
				twod += (family.states.memory[i][j] ^ family.states.memory[k][j]);
			}
		}
	}
	finish = clock();
	time = (finish - start)/(double) CLOCKS_PER_SEC;
	std::cout << time << std::endl;

	return 0;
}

