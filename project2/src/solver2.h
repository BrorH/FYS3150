#ifndef SOLVER2_H
#define SOLVER2_H


class Solver
{
private:
    double **A_init;

    int maxidx[2];

    void fillA(double **);
    void findMaxIdx();
    void rotateA(int, int, double, double);
    void printA();
    bool unitTest();

public:
    double rho_max, h, eps, *diag, **A; 
    int n;
    int transforms = 0;
    Solver(int, double, double, double *); // Constructor
    void solve(); // calls Jacobi_algorithm
    void write(); // writes all to data.dat
   
};

#endif