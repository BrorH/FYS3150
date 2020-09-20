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
    diags = diag;
    double non_diag = -pow(h, -2);
    fillA(n, diag, non_diag);
    rho = mat(n, 1);
    I = mat(n, n, fill::eye);
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

void Solver::Givens(int i, int j, double t)
{
    S = I;
    double c = 1 / sqrt(1 + pow(t, 2));
    double s = t * c;
    S(i, i) = c;
    S(j, j) = c;
    S(i, j) = s;
    S(j, i) = -s;
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
    cout << "Start loop" << endl;
    while (true)
    {
        counts ++;
        max_idx = abs(B - diagmat(B.diag())).index_max(); // Finds linear index of max non-diag element

        idxs = ind2sub(size(A), max_idx); // Decomposes linear index
        i = idxs(0); // Decomposed indecies
        j = idxs(1);

        max_Aij = pow(B(max_idx), 2); // max non-diag value
        if (max_Aij < eps) break;

        tau = (B(i, i) - B(j, j)) / (2 * B(i, j));
        t1 = -tau + sqrt(1 + pow(tau, 2));
        t2 = -tau - sqrt(1 + pow(tau, 2));

        Givens(i, j, min(abs(t1), abs(t2)));
        EigVec *= S;
        // B.print();
        // cout << i << "  " << j << endl;
        B = S.t() * B * S;
        // B.print();
        // cout << endl << endl;
    }
    cout << n << ":  " <<counts << endl;
}
void Solver::write(){
    ofstream datafile;
    datafile.open("data.dat", ios::app);
    datafile << n << "," <<pmax << ","<< eps << ","<< counts<< endl;
    //datafile << "diags";
    for(int i = 0; i < n; i++){
        datafile << diags[i] << ",";
    }
    // datafile << endl <<"eigvals: ";
    // for(int i = 0; i < n; i++){
    //     datafile << B.diag()[i] << ",";
    // }
    datafile << endl;
    for(int i = 0; i < n; i++){
        datafile << B.diag()[i] << ",";
        for (int j =0; j < n; j++){
            datafile << EigVec(j,i) << ",";
        }
        datafile << endl;
    }

    datafile <<"*" << endl;



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
