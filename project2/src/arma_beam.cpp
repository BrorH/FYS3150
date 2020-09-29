#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>

using namespace std;
using namespace arma;

mat fill_A(int);

int main(int argc, char *argv[])
{
	int n = atoi(argv[1]);
	string name = "data/";
	name += argv[2];
	mat A = fill_A(n);

	vec eigval;
	mat eigvec;

	eig_sym(eigval, eigvec, A);

	ofstream datafile;
	datafile.open(name, ios::app);
	datafile << n << endl;
	for (int i = 0; i < n; i++)
	{
		datafile << eigval[i];
		for (int j = 0; j < n; j++)
		{
			datafile << "," << eigvec(j, i);
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