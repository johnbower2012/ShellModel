#ifndef COUNT_H
#define COUNT_H

#include<iomanip>
#include<iostream>

class count{
		int number;
	public:
		void factorial(int);
		void choose(int,int);
		void broken_pairs(int,int,int,int);

		inline int value();
		inline void print();
};

void count::broken_pairs(int n, int omega, int m, int brokenpairs){
	int i, holder;
	if(m%2==0){
		holder = m/2 - brokenpairs;
		choose(n/omega, holder);
		for(i=0;i<2*brokenpairs;i++){
			number *= (n - 2*(i+holder));
			number /= (i+1);
		}
	}
	else if(m%2==1){
		holder = (m-1)/2 - brokenpairs;
		choose(n/omega, holder);		
		for(i=0;i<2*brokenpairs+1;i++){
			number *= (n - 2*(i+holder));
			number /= (i+1);
		}
	}
}

void count::factorial(int n){
	int i;
	number = 1;
	if(n<0){
		std::cout << "count.factorial: Bad usage. Choose a positive number." << std::endl;
		exit(1);
	}
	else if(n==0){
		number = 1;
	}
	else if(n>0){
		for(i=1;i<n+1;i++){
			number *= i;
		}
	}
}
void count::choose(int n, int m){
	int i;
	number = n;
	if(n<m||m<0){
		number = 0;
	}
	else if(m==n||m==0){
		number = 1;
	}
	else{
		if(m*2 > n){
			m = n-m;
		}
		for(i=1;i<m;i++){
			number *= (n-i);
			number /= (i+1);
		}
	}
}
inline int count::value(){
	return number;
}
inline void count::print(){
	std::cout << number;
}
int maximum(int n, int* list){
	int i, max;
	max = list[0];
	for(i=0;i<n;i++){
		if(max<list[i]){
			max = list[i];
		}
	}
	return max;
}
int minimum(int n, int* list){
	int i, min;
	min = list[0];
	for(i=0;i<n;i++){
		if(min>list[i]){
			min = list[i];
		}
	}
	return min;
}

#endif
