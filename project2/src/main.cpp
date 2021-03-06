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
           (NB! Should be set to 0 when solving for behaviour 0 and set to 1 when solving behaviour 1)
*/



#include "solver.h"


int n, method;
double rho_max, epsilon, omega;
string name, datafile;


void fill_d(double*);



int main(int argc, char *argv[])
{
   
    // read arguments from commandline
    name = argv[1];
    n = atoi(argv[2]);
    epsilon = pow(10, -atof(argv[3]));
    rho_max = (double)atof(argv[4]);
    method = atoi(argv[5]);
    omega = (double)atof(argv[6]);
    datafile = argv[7];
   
    double *d = new double[n]; // array to contain all diagonal elements

    fill_d(d);
    Solver problem(n, rho_max, epsilon, d, name); // create problem
    problem.solve(); // solve problem
    problem.write(datafile); // write solution

    delete[] d;
    return 0;
}

void fill_d(double *d)
{
    /*
    Fills the diagonal elements of d.
    In our three cases, only the diagonals differ from case to case.
    We therefore only have to pass the diagonals to our solver.
    */
    double h = rho_max /((double) n+1);
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
