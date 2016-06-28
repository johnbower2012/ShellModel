#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include "project3_library.h"
#include "time.h"

using namespace std;

ofstream ofile;
ofstream ofile2;

clock_t start, finish;

int main(int argc, char* argv[]){
	bool test, loop;
	bool** stateset, **temp;

	int n, p, omega, m, count, nchoosem, counter;
	int i, j, k, q;
	int* list;

	char* outfile, *outfile_screenprint;

	double energy, g, d, tolerance, time;
	double** hamiltonian, **hamil, **vectors, **eigenvalues;
	tolerance = 1e-10;

	/************************************
		Check for proper file input
	************************************/
	if(argc<7){
		cout << "Bad usage. Enter also 'level_count degeneracy level_spacing g outfilename screenprint_filename' on same line" << endl;
		exit(1);
	}
	else{
		p = atoi(argv[1]);
		omega = atoi(argv[2]);
		n = p*omega;
		d = atof(argv[3]);
		g = atof(argv[4]);
		outfile = argv[5];
		outfile_screenprint = argv[6];
	}

	/************************************
		Print to screen constants
		and labels for coming chart
	************************************/
	ofile2.open(outfile_screenprint);

	ofile2 << endl;
	ofile2 << "p=" << p << setw(10) << "omega=" << omega<< setw(10) << "d=" << d << setw(10) << "g=" << g << endl;
	ofile2 << endl;
	ofile2 << "part" << setw(6) << "bp" << setw(10) << "states" << setw(10) << "gs en." << setw(10) << "time" << endl;
	
	cout << endl;
	cout << "p=" << p << setw(10) << "omega=" << omega<< setw(10) << "d=" << d << setw(10) << "g=" << g << endl;
	cout << endl;
	cout << "part" << setw(6) << "bp" << setw(10) << "states" << setw(10) << "gs en." << setw(10) << "time" << endl;

	/************************************
		Write labels to file
	************************************/
	ofile.open(outfile);
	ofile << "#m";
	for(i=0;i<p/2+1;i++){
		ofile << setw(8) << i << "bp";
	}
	ofile << endl;

	/************************************
		Begin looping over number of
		particles. Allocate stateset
	************************************/
	for(i=1;i<n+1;i++){
		m=i;
		choose(nchoosem, n, m);
		matrix_alloc(stateset, nchoosem, n);
		construct_stateset(stateset, count, n, m);
		ofile << m;
		ofile2 << endl << m << endl;
		cout << endl << m << endl;

		/************************************
			Begin looping over number of
			broken pairs. Allocate and
			initialize matrices.
		************************************/
		for(q=0;q<m/2+1;q++){

			/************************************
				Start clock.
				Calculate how many states will
				have q broken pairs given
					n, omega, m
				If count==0, break loop
			************************************/
			start = clock();
			count_brokenpairstates(count, n, omega, m, q);
			if(count==0){
				break;
			}

			/************************************
				**list holds the location of 
					appropriate states from 
					stateset
				**temp is created from the **list
					of states
				**hamiltonian is energy matrix
				**hamil copies hamiltonian so it
					may be printed to screen
				**vectors are the state vectors
					of hamiltonian--shows mixing
				**eigenvalues is to be an ordered
					list of energies, sm-->lg
			************************************/
			array_alloc(list,count);
			matrix_alloc(temp, count, n);
			matrix_alloc(hamiltonian, count, count);
			matrix_alloc(hamil, count, count);
			matrix_alloc(vectors, count, count);
			matrix_alloc(eigenvalues, count, 2);
			for(j=0;j<count;j++){
				vectors[j][j] = 1.0;
			}
			for(j=0;j<count;j++){
				eigenvalues[j][0] = 0;
				eigenvalues[j][1] = -1;
			}

			/************************************
				construct list from stateset of 
				those states which have q broken
				pairs
			************************************/
			construct_brokenpairs_list(stateset, list, nchoosem, n, omega, q);
			for(j=0;j<count;j++){
				for(k=0;k<n;k++){
					temp[j][k] = stateset[list[j]][k];
				}
			}
			
			/************************************
				construct and diagonalize the
				energy matrix given our pairing
				hamiltonian
			************************************/
			derive_hamiltonian_matrix_zproj0(hamiltonian, temp, count, n, omega, d, g);
			jacobi_simrot_eigen_solver(hamiltonian, vectors, count, tolerance, counter);
			matrix_diag_sort(hamiltonian, eigenvalues, count);


			/************************************
				Calculate time to find all states
				with q broken pairs, construct
				hamiltonian, and diagonalize,
				all steps necessary to find
				groundstate energy.
			************************************/
			finish = clock();
			time = (finish - start)/((double) CLOCKS_PER_SEC);

			/************************************
				Print to screen information about
				bp, states, gs en., time
			************************************/
			ofile2 << setw(10) << q << setw(10) << count << setw(10) <<  setw(10) << eigenvalues[0][0] << setw(15) << time << endl;			
			cout << setw(10) << q << setw(10) << count << setw(10) <<  setw(10) << eigenvalues[0][0] << setw(15) << time << endl;
			ofile << setw(10) << eigenvalues[0][0];

			/************************************
				Delete all but stateset to allow
				for next broken pair count
			************************************/
			array_delete(list);
			matrix_delete(temp, count);
			matrix_delete(hamiltonian, count);
			matrix_delete(vectors, count);
			matrix_delete(eigenvalues, count);
		}

		/************************************
			Next line in output file
			Delete stateset to allow for next
			next particle count
		************************************/
		ofile << endl;
		matrix_delete(stateset, nchoosem);
	}

	/************************************
		close output file
		create space for tidy screen
	************************************/
	ofile.close();
	ofile2.close();
	cout << endl;
	return 0;

}


