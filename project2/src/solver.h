#ifndef SOLVER_H
#define SOLVER_H

#include <armadillo>
using namespace arma;

class Solver
{
private:
    double pmax; // defaults to 1
    double h;

    int counts = 0;

    double eps; // tolerance, maximum allowed non-diagonal element

    mat rho; // (n, 1) matrix with radial domain
    mat A; // (n, n) tridiagonal matrix
    mat B; // (n, n) diagonal matrix with eigenvalues of A
    mat EigVec; // (n, n) Eigenvectormatrix of A
    mat diags;
    void fillA(int, mat, double);
    mat Givens(int, int, double);
    void Jacobi_algorithm();

public:
    int n;
    Solver(int, double, double, mat); // Constructor
    void solve(); // calls Jacobi_algorithm
    void write(); // writes all to data.dat
    mat get_A(); // returns A
    mat get_rho(); // returns rho
    mat eigenvalues(); // returns diagonal of EigVal
    mat eigenvectors(); // returns EigVec
};

#endif