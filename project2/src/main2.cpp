/*


Parameters:
    n: int
        - number of points we solve for. N - 1
    epsilon: float
        - tolerance of algorithm
    rhomax: float
        - length of system
    behaviour: int (0, 1, 2)
        - type of problem we solve. (Beam, one electron, two electron)
    omega: float
        - frequency of HO in two electron problem
*/

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

#include "solver2.h"

//char *garg[6];
int n, method;
double rho_max, epsilon, omega;
using namespace std;

void fill_d(double *d)
{
    /*
    Diagonal of matrix A, for the three cases our solver will be used for.
    Arguments set in commandline:
    -----------
    n: int
        - number of points
    rho_max: float
        - length of system
    behaviour: (0, 1, 2)
        - problem to solve
    omega: float
        - frequency in last problem
    */
    double h = rho_max /(double) n;
    double omega = 0;
    int c = 0; // numerator in fraction in term for two electron
    if (method == 2) c = 1;
    
    double rho = 0;
    for (int i = 0; i < n; i++)
    {
        rho += h;
        d[i] = 2 * pow(h, -2) + pow(omega * rho, 2) + c/rho;
    }
    
}


void print(double *v){
    // prints vector v
    for(int i = 0; i < n; i ++){
        cout << v[i]<<" ";
    }
    cout <<  endl;
}


int main(int argc, char *argv[])
{
    if (argc <= 4)
    {
        cout << "Bad usage" << endl;
        exit(1);
    }
    n = atoi(argv[1]); // read n from commandline
    epsilon = pow(10, -atoi(argv[2]));
    rho_max = (double)atof(argv[3]);
    method = atoi(argv[4]);
    omega = (double)atof(argv[5]);

    if ((method == 0) && (omega != 0))cout << "method set to 0 but omega is " << omega << ". Should usually be 0. Proceeding." <<endl; 
    if ((method == 0) && (rho_max != 1))cout << "method set to 0 but rho_max is " << rho_max << ". Should usually be 1. Proceeding." <<endl; 
    if ((method == 1) && (omega != 1))cout << "method set to 1 but omega is " << omega << ". Should usually be 1. Proceeding." <<endl; 
    
    double *d = new double[n];

    fill_d(d);
    Solver problem(n, rho_max, epsilon, d);
    problem.solve();
    //problem.write();
    //Solver problem(n, rho_max, tolerance, d);
    //problem.solve();
    //problem.write();
    //mat eigvals = problem.eigenvalues();

    //sort(eigvals).print();
    delete[] d;
    return 0;
}

