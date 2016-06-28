#include<iostream>
#include<iomanip>
#include "stateset.h"
#include "states.h"
#include "time.h"



int main(int argc, char* argv[]){
	int i,j,k,n,m,s,omega,p,q,w,e;
	count z;
	int* list, diff;
	bool fixed;
	clock_t start, finish;
	double time;

	if(argc<4){
		std::cout << "Bad usage. Include 'p m bool_fixed_omega'." << std::endl;
		exit(1);
	}
	else{
		p = atoi(argv[1]);
		m = atoi(argv[2]);
		fixed = atoi(argv[3]);
		q = atoi(argv[4]);
	}

	array<int> omegal(p);
	if(fixed==true){
		std::cout << "Enter omega: ";
		std::cin >> omega;
		for(i=0;i<p;i++){
			omegal.memory[i] = omega;
		}
		n = p*omega;
	}
	else if(fixed==false){
		n = 0;
		for(i=0;i<p;i++){
			std::cout << "level" << i << ":";
			std::cin >> omegal.memory[i];
			n += omegal.memory[i];
		}
	}

	z.choose(n,m);
	s = z.value();
	stateset family(p,omegal.memory,n,m,s);

	array<int> holes, particles, hp;
	for(i=0;i<s;i++){
		diff = family.differences(q,i);
		holes.resize(diff,0);
		particles.resize(diff,0);
		hp.resize(diff*2,0);
		holes.memory = family.hole_loc(q,i,diff);
		particles.memory = family.particle_loc(q,i,diff);
		hp.memory = family.hp_loc(diff,holes.memory,particles.memory);

		std::cout << std::setw(5) << i << std::setw(5);
		for(j=0;j<n;j++){
			std::cout << family.states.memory[i][j];
		}
		std::cout << std::setw(5) << family.zedangmom(i);
		std::cout << std::setw(5) << family.phase(q,diff,holes.memory,particles.memory) << std::setw(5);
		std::cout << std::setw(5) << family.phase_(q,diff,hp.memory) << std::setw(5);
		std::cout << std::setw(5) << diff << ":" << std::setw(5);
		for(j=0;j<diff*2;j++){
			std::cout << hp.memory[j];
		}
		std::cout << std::setw(5);
		for(j=0;j<diff;j++){
			std::cout << holes.memory[j];
		}
		std::cout << std::setw(5) << "-->" << std::setw(5);
		for(j=0;j<diff;j++){
			std::cout << particles.memory[j];
		}

		std::cout << std::endl;
	}
	
	array<int> list2;
	int count,min,max;
	max = family.zedang_max();
	min = family.zedang_min();
	for(j=min;j<max+1;j+=2){
		count = family.zedang_count(j);
		list2.resize(count,0);
		list2.memory = family.zedang_loc(j,count);
		std::cout << std::setw(5) << j << ": ";
		for(i=0;i<count;i++){
			std::cout << std::setw(5) << list2.memory[i];
		}
		std::cout << std::endl;
	}

	start = clock();
	for(i=0;i<s;i++){		
		diff = family.differences(q,i);
		holes.resize(diff,0);
		particles.resize(diff,0);
		hp.resize(diff*2,0);
		holes.memory = family.hole_loc(q,i,diff);
		particles.memory = family.particle_loc(q,i,diff);
		family.phase(q,diff,holes.memory,particles.memory);
	}
	finish = clock();
	time = (finish - start)/(double) CLOCKS_PER_SEC;

	std::cout << time << std::endl;

	start = clock();
	for(i=0;i<s;i++){		
		diff = family.differences(q,i);
		holes.resize(diff,0);
		particles.resize(diff,0);
		hp.resize(diff*2,0);
		holes.memory = family.hole_loc(q,i,diff);
		particles.memory = family.particle_loc(q,i,diff);
		hp.memory = family.hp_loc(diff,holes.memory,particles.memory);
		family.phase_(q,diff,hp.memory);
	}
	finish = clock();
	time = (finish - start)/(double) CLOCKS_PER_SEC;

	std::cout << time << std::endl;


	return 0;
}

