/*
*/

#include <iostream>
#include <string>
#include <fstream>
#include <armadillo>
#include "solver.h"

using namespace std;
using namespace arma;

char *garg[5];

mat diag(int, double);


int main(int argc, char *argv[])
{
    if (argc <= 4) {
        cout << "Bad usage" << endl;
        exit(1);
    }
    int n = atoi(argv[1]); // read n from commandline
    double rho_max = atof(argv[2]);
    for (int i = 0; i < argc; i++)
    {
        garg[i] = argv[i];
    }
    double tolerance = 1e-8;

    mat d = diag(n, rho_max);
    Solver problem(n, rho_max, tolerance, d);
    problem.solve();
    mat eigvals = problem.eigenvalues();
    eigvals.print();

    return 0;
}

mat diag(int n, double rho_max)//, int behaviour = 0, double rho_max = 0, double omega = 0)
{
    /*
    int n, double h, int particles = 0, double rho_max = 0, double omega = 0
    Diagonal of matrix A, for the three cases our solver will be used for.
    Arguments:
    -----------
    n: int
        number of points
    h: double
        step length
    behaviour: int (optional)
        extra terms in the diagonal of A
        0: No extra term, for buckling beam problem
        1: rho^2 term, for single electron in HO
        2: w^2p^2 + 1/p, for two electron
    rho_max: double (optional)
        for use in quantum problems
    omega: double (optional)
        for use in two-electron system
    */
    double h = rho_max / n;
    int behaviour = atoi(garg[3]);
    double omega = atof(garg[4]);

    mat d(n, 1, fill::zeros);

    for (int i = 0; i < n; i++)
    {
        double rho = i * h + h;
        d[i] = 2 * pow(h, -2);
        if (behaviour == 1)
        {
            d[i] += pow(rho, 2);
        }
        if (behaviour == 2)
        {
            d[i] += pow(omega * rho, 2);
            d[i] += pow(rho, -1);
        }
    }
    return d;
}
