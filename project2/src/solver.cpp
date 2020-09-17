#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "solver.h"

using namespace std;
using namespace arma;

Solver::Solver(int _n, double rho_max, double tolerance, mat diag)
{
    n = _n;
    eps = tolerance;
    pmax = rho_max;
    h = pmax / n;

    double non_diag = -pow(h, -2);
    fillA(n, diag, non_diag);
    rho = mat(n, 1);
    for (int i = 1; i < n + 1; i++)
    {
        rho[i] = i * h;
    }
}

void Solver::fillA(int n, mat d, double a)
{
    A = mat(n, n);
    for (int i = 0; i < n; i++)
    {
        A(i, i) = d(i);
        if (i != n - 1)
        {
            A(i + 1, i) = a;
            A(i, i + 1) = a;
        }
    }
}

void Solver::solve()
{
    Jacobi_algorithm();
}

mat Solver::Givens(int i, int k, double t)
{
    mat S(n, n, fill::eye);
    double c = 1 / sqrt(1 + pow(t, 2));
    double s = t * c;
    S(i, i) = c;
    S(k, k) = c;
    S(i, k) = s;
    S(k, i) = -s;
    return S;
}

void Solver::Jacobi_algorithm()
{
    EigVal = A;
    EigVec = mat(n,n,fill::eye);
    mat Eig_prime;
    int max_idx;
    int i;
    int k;
    uvec idxs;
    double t1;
    double t2;
    double tau;
    double max_aik;

    while (true){
        Eig_prime = EigVal - diagmat(EigVal.diag());
        max_idx = abs(Eig_prime).index_max();

        idxs = ind2sub(size(A), max_idx);
        i = idxs(0);
        k = idxs(1);

        max_aik = pow(Eig_prime(i, k), 2);
        if (max_aik < eps)
        {
            break;
        }
        tau = (EigVal(i,i) - EigVal(k,k)) / (2 * EigVal(i, k));
		t1 = -tau + sqrt(1 + pow(tau, 2));
		t2 = -tau - sqrt(1 + pow(tau, 2));
        mat S = Givens(i, k, min(t1, t2));
        EigVec *= S;
        EigVal = S.t() * EigVal * S;

    }
}

mat Solver::get_A()
{
    return A;
}
mat Solver::get_rho()
{
    return rho;
}
mat Solver::eigenvalues()
{
    return EigVal.diag();
}
mat Solver::eigenvectors()
{
    return EigVec;
}
