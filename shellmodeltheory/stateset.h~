#ifndef STATESET_H
#define STATESET_H

#include<iomanip>
#include<iostream>
#include<cmath>
#include "memory.h"
#include "count.h"

class stateset{
	public:
		int levels,
			spstates, 
			particles, 
			mpstates;
		
		array<int> omega;
		matrix<bool> states;

		stateset();
		stateset(int,int*,int,int,int);

		inline void print(int,int);

		int differences(int,int);
		int* hole_loc(int,int,int);
		int* particle_loc(int,int,int);
		int* hp_loc(int,int*,int*);
		int phase(int,int,int*,int*);
		int phase_(int,int,int*);

		int zedangmom(int);
		int zedang_min();
		int zedang_max();
		int zedang_min(int,int*);
		int zedang_max(int,int*);
		int zedang_count(int);
		int* zedang_loc(int,int);
};

class operators{
	public:
		int length,
			reference;
		array<bool> state;
		bool phase;

		operators();
		operators(int,int,bool*);

		void print();
		void update(int,bool*);
		void reset(int,int,bool*);

	 	void annihilate(int);
		void create(int);
};



/***************************************

		STATESET CLASS DEFINITIONS

***************************************/
/***************************************
	Construct null stateset
***************************************/
stateset::stateset() : states(0,0), omega(0) {
	levels = 0; 
	spstates = 0; 
	particles = 0; 
	mpstates = 0;
}
/***************************************
	Construct state set with n single
	particle states, m particles, and
	s many particle states -- note s
	should match n choose m.
***************************************/
stateset::stateset(int p, int* list, int n, int m, int s) : states(s,n), omega(p) {
	int i, j, count, temp;
	array<int> m_list(m);

	count = 0;
	for(j=0;j<p;j++){
		count += list[j];
		omega.memory[j] = list[j];
	}
	if(count!=n){
		std::cout << "Error:	The number of single particle states from array 'list' does not match the input of 'n'" << std::endl;
		exit(1);
	}

	levels = p; 
	spstates = n; 
	particles = m;
	mpstates = s;

	i=0; count = 0;
	for(j=0;j<m;j++){
		m_list.set(j,j);
	}

	while(i<m){
		if(m_list.value(m-1-i)<(n-i)){
			if(i==0){
				for(j=0;j<m;j++){
					states.set(count,m_list.value(j),true);
				}
				temp = m_list.value(m-1) + 1;
				m_list.set(m-1,temp);
				count++;
			}
			else if(i<m){
				temp = m_list.value(m-1-i) + 1;
				m_list.set(m-1-i,temp);
				while(i>0){
					m_list.set(m-i, m_list.value(m-1-i) + 1);
					i -= 1;
				}
			}
		}
		else if(m_list.value(m-1-i)==n-i){
			i += 1;
		}
	}
}
/***************************************
	Print occupation of position y in 
	state x
***************************************/
inline void stateset::print(int x, int y){
	std::cout << states.memory[x][y];
}
/***************************************
	Calculate the number of moved 
	particles for same-sized states
***************************************/
int stateset::differences(int a, int b){
	int i, excited = 0;
	for(i=0;i<spstates;i++){
		excited += (states.memory[a][i]^states.memory[b][i]);
	}
	return excited/2;
}/***************************************
	Return list of holes from state a
	to state b
***************************************/
int* stateset::hole_loc(int a, int b, int diff){
	int i, hole, *list, marker = 0;
	list = new int[diff];
	for(i=0;i<spstates;i++){
		hole = (states.memory[a][i]^states.memory[b][i])&states.memory[a][i];
		if(hole!=0){
			list[marker] = i;
			marker++;
		}
	}
	return list;
}
/***************************************
	Return list of particles from 
	state a to state b
***************************************/
int* stateset::particle_loc(int a, int b, int diff){
	int i, part, *list, marker = 0;
	list = new int[diff];
	for(i=0;i<spstates;i++){
		part = (states.memory[a][i]^states.memory[b][i])&states.memory[b][i];
		if(part!=0){
			list[marker] = i;
			marker++;
		}
	}
	return list;
}
/***************************************
	Return ordered list of holes and
	particles
***************************************/
int* stateset::hp_loc(int diff, int* holes, int* particles){
	int i, 
		marker_holes = 0, 
		marker_particles = 0,
		length = diff*2;
	int* list = new int[length];
	for(i=0;i<length;i++){
		if(holes[marker_holes] < particles[marker_particles]){
			list[i] = holes[marker_holes];
			marker_holes++;
		}
		else{
			list[i] = particles[marker_particles];
			marker_particles++;
		}
		
		if(marker_holes==diff){
			i++;
			while(i<length){
				list[i]=particles[marker_particles];
				marker_particles++;
				i++;
			}
		}
		else if(marker_particles==diff){
			i++;
			while(i<length){
				list[i]=holes[marker_holes];
				marker_holes++;
				i++;
			}
		}
	}
	return list;
}
/***************************************
	Returns the phase between a state 
	a and b given sorted lists of the 
	holes and particles
***************************************/
int stateset::phase(int a, int n, int* holes, int* parts){
	int i, j, max, min, temp, phase;
	bool fase = false;

	array<bool> list(spstates);
	for(i=0;i<spstates;i++){
		list.memory[i] = states.memory[a][i];
	}

	max = maximum(n,holes);
	temp = maximum(n,parts);
	if(max<temp){
		max = temp;
	}
	min = minimum(n,holes);
	temp = minimum(n,parts);
	if(min>temp){
		min = temp;
	}
	
	for(i=0;i<n;i++){
		list.memory[holes[i]] = false;
		for(j=min;j<holes[i];j++){
			fase = fase^list.memory[j];
		}
	}
	for(i=n-1;i>-1;i--){
		list.memory[parts[i]] = true;
		for(j=min;j<parts[i];j++){
			fase = fase^list.memory[j];
		}
	}
	phase = pow(-1,fase);
	return phase;
}
/***************************************
	Returns the phase between a state 
	a and b given a single sorted list
	of the holes and particles
***************************************/
int stateset::phase_(int a, int diff, int* hp_loc){
	int i, j, phase, start, end, twoi;
	bool fase = false;
	for(i=0;i<diff;i++){
		twoi = i*2;
		start = hp_loc[twoi]+1;
		end = hp_loc[twoi+1];
		for(j=start;j<end;j++){
			fase = fase^states.memory[a][j];
		}
	}

	phase = pow(-1,fase);
	return phase;
}

/***************************************
	Return total zprojection of 
	angular momentum
***************************************/
int stateset::zedangmom(int a){
	int i, marker1=0, marker2=0, zproj=0;
	for(i=0;i<spstates;i++){
		if(states.memory[a][i]==true){
			zproj += (omega.memory[marker2] - 1 - 2*(marker1%omega.memory[marker2]));
		}
		marker1++;
		if(marker1==omega.memory[marker2]){
			marker2++;
			marker1=0;
		}
	}
	return zproj;
}
/***************************************
	Return maximum zprojection of
	angular momentum
***************************************/
int stateset::zedang_min(){
	int i,min;
	min = zedangmom(0);
	for(i=0;i<mpstates;i++){
		if(min>zedangmom(i)){
			min = zedangmom(i);
		}
	}
	return min;
}
/***************************************
	Return minimum zprojection of
	angular momentum
***************************************/
int stateset::zedang_max(int count, int* list){
	int i,max;
	max = zedangmom(list[0]);
	for(i=0;i<count;i++){
		if(max<zedangmom(list[i])){
			max = zedangmom(list[i]);
		}
	}
	return max;
}
/***************************************
	Return minimum zprojection of
	angular momentum from states in
	list
***************************************/
int stateset::zedang_min(int count, int* list){
	int i,min;
	min = zedangmom(list[0]);
	for(i=0;i<count;i++){
		if(min>zedangmom(list[i])){
			min = zedangmom(list[i]);
		}
	}
	return min;
}
/***************************************
	Return maximum zprojection of
	angular momentum from states in 
	list
***************************************/
int stateset::zedang_max(){
	int i,max;
	max = zedangmom(0);
	for(i=0;i<mpstates;i++){
		if(max<zedangmom(i)){
			max = zedangmom(i);
		}
	}
	return max;
}
/***************************************
	Return count of states with total 
	zprojection = a
***************************************/
int stateset::zedang_count(int a){
	int i, zproj, count = 0;
	for(i=0;i<mpstates;i++){
		zproj = zedangmom(i);
		if(zproj==a){
			count++;
		}
	}
	return count;
}
/***************************************
	Return array for list of states
	with total zprojection = a
***************************************/
int* stateset::zedang_loc(int a, int count){
	int i, counter, zproj, *list;
	counter = 0;
	list = new int[count];
	for(i=0;i<mpstates;i++){
		zproj = zedangmom(i);
		if(zproj==a){
			list[counter] = i;
			counter++;
		}
	}
	return list;
}


/***************************************

		OPERATORS CLASS DEFINITIONS

***************************************/

/***************************************
	Null constructor
***************************************/
operators::operators(){
	length = 0; reference = -1; 
	state.memory = nullptr;
	phase = false;
}
/***************************************
	Construct state of length n, index
	a in stateset, and structure of copy
***************************************/
operators::operators(int n, int a, bool* copy) : state(n) {
	int i;
	length = n;
	reference = a;
	phase = false;
	for(i=0;i<length;i++){
		state.memory[i] = copy[i];
	}
}
/***************************************
	Print state with no spaces
***************************************/
void operators::print(){
	for(int i=0; i<length; i++){ std::cout << state.memory[i];}
}
/***************************************
	New index reference and state 
	structure
	Keep length
***************************************/
void operators::update(int a, bool* copy){
	int i;
	reference = a;
	phase = false;
	for(i=0;i<length;i++){
		state.memory[i] = copy[i];
	}
}
/***************************************
	Wipe old information, replace with
	full new set
***************************************/
void operators::reset(int n, int a, bool* copy){
	int i;
	if(n!=length){
		state.resize(n,0);
	}
	length = n;
	reference = a;
	phase = false;
	for(i=0;i<length;i++){
		state.memory[i] = copy[i];
	}
}
/***************************************
	Annihilate space j
	Test for full space
		--if true, set false & 
			calculate phase
		--if false, collapse to nullptr
***************************************/
void operators::annihilate(int j){
	int i;
	if(j>(length-1)){
		std::cout << "Error:	operators::annihilate improperly used. Entry value larger than bit array length." << std::endl;
		exit(1);
	}
	if(state.memory[j]==true){
		state.memory[j] = false;
		for(i=0;i<j;i++){
			phase = phase^state.memory[i];
		}
	}
	else if(state.memory[j]==false){
		delete[] state.memory;
		state.memory = nullptr;
		length = 0;
		reference = -1;
		phase = false;
	}
}
/***************************************
	Annihilate space j
	Test for empty space
		--if false, set true & 
			calculate phase
		--if true, collapse to nullptr
***************************************/
void operators::create(int j){
	int i;	
	if((j>length-1)){
		std::cout << "Error:	operators::create improperly used. Entry value larger than bit array length." << std::endl;
		exit(1);
	}
	if(state.memory[j]==false){
		state.memory[j] = true;
		for(i=0;i<j;i++){
			phase = phase^state.memory[i];
		}
	}
	else if(state.memory[j]==true){
		delete[] state.memory;
		state.memory = nullptr;
		length = 0;
		reference = -1;
		phase = false;
	}
}

#endif
