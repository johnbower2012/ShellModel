#ifndef STATES_H
#define STATES_H

#include<iomanip>
#include<iostream>
#include "count.h"
#include "memory.h"

class states{
	public:
		int levels, degeneracy;
		matrix<bool> state;

		states(int,int);
		states(int,int*);
};

states::states(int p, int omega) : state(p,omega){
	int i, j;	
	levels = p; degeneracy = omega;
	for(i=0;i<levels;i++){
		for(j=0;j<degeneracy;j++){
			state.memory[i][j] = false;
		}
	}
}


#endif
