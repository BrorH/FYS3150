#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "time.h"

using namespace std;
using namespace arma;

mat fill_A(int);

int main(int argc, char *argv[])
{
	int n = atoi(argv[1]);
	string name = "data/";
	name += argv[2];
	mat A = fill_A(n);

	cx_vec eigval;
	cx_mat eigvec;

	clock_t start = clock();
	eig_gen(eigval, eigvec, A);
	double time = ((clock() - start)/(double)CLOCKS_PER_SEC);

	ofstream datafile;
	datafile.open(name, ios::app);
	datafile << n << "," << time << endl;
	for (int i = 0; i < n; i++)
	{
		datafile << real(eigval[i]);
		for (int j = 0; j < n; j++)
		{
			datafile << "," << real(eigvec(j, i));
		}
		datafile << endl;
	}
	return 0;
}

mat fill_A(int n)
{
	mat A(n, n, fill::zeros);
	double h = 1.0 / (double)(n + 1);
	double hh = pow(h, -2);
	for (int i = 0; i < n; i++)
	{
		A(i, i) = 2 * hh;
		if (i != n - 1)
		{
			A(i, i + 1) = -hh;
			A(i + 1, i) = -hh;
		}
	}
	return A;
}