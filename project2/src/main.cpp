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
#include "solver.h"

char *garg[6];

mat diag(int, double);

int main(int argc, char *argv[])
{
    if (argc <= 4)
    {
        cout << "Bad usage" << endl;
        exit(1);
    }
    int n = atoi(argv[1]); // read n from commandline
    double tolerance = pow(10, -atoi(argv[2]));
    double rho_max = atof(argv[3]);
    for (int i = 0; i < argc; i++)
    {
        garg[i] = argv[i];
    }

    mat d = diag(n, rho_max);
    Solver problem(n, rho_max, tolerance, d);
    problem.solve();
<<<<<<< HEAD
    problem.write();
    //mat eigvals = problem.eigenvalues();

    //sort(eigvals).print();

=======
    mat eigvals = problem.eigenvalues();
    sort(eigvals).print();
>>>>>>> 945b25ee5e3e8ca90612dd1a9831991cd2df2f80
    return 0;
}

mat diag(int n, double rho_max)
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
    double h = rho_max / n;
    int behaviour = atoi(garg[4]);
    double omega = 0;
    int c = 0; // numerator in fraction in term for two electron
    if (behaviour == 1)
    {
        omega = 1;
    }
    else if (behaviour == 2)
    {
        omega = atof(garg[5]);
        c = 1;
    }

    mat d(n, 1, fill::zeros);

    double rho = 0;
    for (int i = 0; i < n; i++)
    {
        rho += h;
        d[i] = 2 * pow(h, -2);
        d[i] += pow(omega * rho, 2);
        d[i] += c / rho;
    }
    return d;
}
