/*
main.cpp: Setus up the instance and the problem based on passed parameters:

Parameters:
    n: int
        - Number of points we solve for. n = N - 1, where N is the size of our system
    epsilon: int
        - Tolerance of jacobi algorithm (in -log10).
            i.e epsilon = 8 -> tolerance of 1e-8.
    rhomax: float
        - Upper bound on rho variable
    behaviour: int (0, 1, 2)
        - Which problem to solve for.
            0: Buckling beam problem
            1: Quantum system; one electron
            2: Quantum system; two electrons
    omega: float
        - Frequency of harmonic oscillator in the quantum systems
           (NB! Should be set to 0 when solving for behaviour 0!)
*/



#include "solver.h"


int n, method;
double rho_max, epsilon, omega;
string name;


void fill_d(double *d)
{
    /*
    Fills the diagonal elements of d.
    In our three cases, only the diagonals differ from case to case.
    We therefore only have to pass the diagonals to our solver.
    */
    double h = rho_max /(double) n;
    double rho = 0;

    int c; // numerator in fraction in term for two electron
    if (method == 2) c = 1;
    else c = 0;

    for (int i = 0; i < n; i++)
    {
        rho += h;
        d[i] = 2 * pow(h, -2) + pow(omega * rho, 2) + c/rho;
    }

}



int main(int argc, char *argv[])
{
    if (argc <= 6){
        cout << "Not correct amount of (proceeding) args. Expected 5, got " << argc<<endl;
        exit(1);
    }

    // read arguments from commandline
    name = argv[1];
    n = atoi(argv[2]);
    epsilon = pow(10, -atoi(argv[3]));
    rho_max = (double)atof(argv[4]);
    method = atoi(argv[5]);
    omega = (double)atof(argv[6]);

    if ((method == 0) && (omega != 0))cout << "method set to 0 but omega is " << omega << ". Should usually be 0. Proceeding." <<endl;
    if ((method == 0) && (rho_max != 1))cout << "method set to 0 but rho_max is " << rho_max << ". Should usually be 1. Proceeding." <<endl;
    if ((method == 1) && (omega != 1))cout << "method set to 1 but omega is " << omega << ". Should usually be 1. Proceeding." <<endl;

    double *d = new double[n]; // array to contain all diagonal elements

    fill_d(d);
    Solver problem(n, rho_max, epsilon, d, name); // create problem
    problem.solve(); // solve problem
    problem.write(); // write solution

    delete[] d;
    return 0;
}

