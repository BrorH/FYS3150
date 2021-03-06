#ifndef SOLVER_H
#define SOLVER_H
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include "time.h"

using namespace std;

class Solver
{
private:
    double **A_init;
    double time;
    string name;

    int maxidx[2];

    void fillA(double **);
    void findMaxIdx();
    void rotateA(double, double);
    void print(double **);
    bool unitTest();
    void checkOrthog(double tol);
    void checkSymmetry(double tol);
    void checkMaxIdx();


public:
    double rho_max, h, eps, *diag, **A, **V;
    int n;
    int transforms = 0;
    Solver(int, double, double, double *, string); // Constructor
    void solve(); // calls Jacobi_algorithm
    void write(string filename); // writes all to data.dat
       ~Solver();
};

#endif