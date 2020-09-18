#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "solver.h"

using namespace std;
using namespace arma;

double static pi = 3.1415926538;

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

mat Solver::Givens(int i, int j, double t)
{
    mat S(n, n, fill::eye);
    double c = 1 / sqrt(1 + pow(t, 2));
    double s = t * c;
    S(i, i) = c;
    S(j, j) = c;
    S(i, j) = s;
    S(j, i) = -s;
    return S;
}

void Solver::Jacobi_algorithm()
{
    B = A;
    EigVec = mat(n, n, fill::eye);
    mat Eig_prime;
    int max_idx;
    int i;
    int j;
    uvec idxs;
    double t1;
    double t2;
    double tau;
    double max_Aij;
    int count = 0;
    
    while (true)
    {
        count ++;
        max_idx = abs(B - diagmat(B.diag())).index_max(); // Finds linear index of max non-diag element

        idxs = ind2sub(size(A), max_idx); // Decomposes linear index
        i = idxs(0); // Decomposed indecies
        j = idxs(1);

        max_Aij = pow(B(max_idx), 2); // max non-diag value
        if (max_Aij < eps) break;

        tau = (B(i, i) - B(j, j)) / (2 * B(i, j));
        t1 = -tau + sqrt(1 + pow(tau, 2));
        t2 = -tau - sqrt(1 + pow(tau, 2));
        // double t = t2;
        // if (abs(t1) < 1) t = t1;
        // cout << atan(t1) * 180 / pi << "  " << atan(t2) * 180 / pi << endl;
        // cout << "sent " << atan(t) * 180 / pi << endl << endl;
        // double t;
        // if (tau > 0 ) {
        //     t = t2;
        // } else {
        //     t = t1;
        // }
        // cout << tau << endl;
        // cout << t1 << "   " << t2 << endl;
        // cout << "sent " << t << endl << endl;
        mat S = Givens(i, j, min(abs(t1), abs(t2)));
        // mat S = Givens(i, j, t);
        EigVec *= S;
        B = S.t() * B * S;
    }
    cout << n << ":  " <<count << endl;
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
    return B.diag();
}
mat Solver::eigenvectors()
{
    return EigVec;
}
