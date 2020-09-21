#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include "solver2.h"

using namespace std;



Solver::Solver(int _n, double _rho_max, double _eps, double *_diag)
{
    n = _n;
    eps = _eps;
    rho_max = _rho_max;
    h = rho_max / n;
    diag = _diag;

    A = new double *[n];
    A_init = new double* [n]; // To store the original A
    fillA(A);
    fillA(A_init);
 
  
 
}

void Solver::fillA(double **_A)
{   double off_diag = -1*pow(h,-2);
    for (int i = 0; i < n; i++){ _A[i] = new double[n];}
    for (int i = 0; i < n; i++)
    {
        _A[i][i] = diag[i];
        if (i != n - 1)
        {
            _A[i + 1][i] = off_diag;
            _A[i][i + 1] = off_diag;
        }
    }
}



bool Solver::unitTest(){
    /* First check if eigenvectors still are orthonormal
    armadillo's function 'accu' sumes over all elements in a matrix */
    //cout << "Orthog? Should be 0: " << accu(V.t()*V)-n <<endl;
    /* Then checking if V are still eigenvectors */
    //cout << "Eigvecs? Should be 0: "<< accu(armaEigvals-eig_sym(B)) << endl;
    //armaEigvals=eig_sym(B);
    /* Lastly checking is A is symmetric */
    //cout << "Symm? Should be 0: "<< accu(B-B.t()) << endl;
    cout << endl;
    return true;
}   



void Solver::findMaxIdx(){
    // burde optimaliseres
    double maxval = 0;
    maxidx[0] = 0; maxidx[1] = 0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(i==j) continue;
            if (abs(A[i][j]) > maxval){
                maxidx[0] = i;
                maxidx[1] = j;
                maxval = abs(A[i][j]);
            }
        }
    }
    
}

void Solver::rotateA(int k, int l, double c, double s){
    // instead of passing k and l they can just be retrieved from the global maxidx array
    double A_kk, A_ll, A_ik, A_il;
    A_kk = A[k][k];
    A_ll = A[l][l];
    A[k][k]=A_kk*pow(c,2) - 2*A[k][l]*c*s + A_ll*pow(s,2);
    A[l][l]=A_ll*pow(c,2) + 2*A[k][l]*c*s + A_kk*pow(s,2);
    A[k][l] = 0;
    A[l][k] = 0;
    for(int i = 0; i < n ; i ++){
        if((i != k) && (i != l)){
            A_ik = A[i][k];
            A_il = A[i][l];
            A[i][k] = A_ik*c - A_il*s;
            A[k][i] = A[i][k];
            A[i][l] = A_il*c + A_ik*s;
            A[l][i] = A[i][l];
        }
    }
   

}

void Solver::printA(){
    for(int i = 0; i < n; i ++){
        for(int j = 0; j < n; j ++){
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
}

void Solver::solve()
{

    int k, l;
    double t,c,s;
    double tau;
    transforms = 0;
    double max_Aij= 2*eps; // make it bigger for now 
    while (max_Aij > eps)
    {
        transforms ++;
        findMaxIdx();
        k = maxidx[0]; l = maxidx[1];
     
        max_Aij = pow(A[k][l], 2); // max non-diag value squared
       
        tau = (A[l][l] - A[k][k])/(2*A[k][l]);
        if(tau > 0){
            t = 1/(tau+sqrt(1+pow(tau,2)));
        } else {
            t = 1/(tau-sqrt(1+pow(tau,2)));
        }
        c = 1/sqrt(1+pow(t,2));
        s = c*t;
        rotateA(k,l,c,s);
        
        
        
    }
    cout << n << ":  " <<transforms << endl;
}
void Solver::write(){
    ofstream datafile;
    datafile.open("data.dat", ios::app);
    datafile << n << "," <<rho_max << ","<< eps << ","<< transforms<< endl;
    //datafile << "diags";
    for(int i = 0; i < n; i++){
        datafile << diag[i] << ",";
    }
    // datafile << endl <<"eigvals: ";
    // for(int i = 0; i < n; i++){
    //     datafile << B.diag()[i] << ",";
    // }
    datafile << endl;
    for(int i = 0; i < n; i++){
        //datafile << B.diag()[i] << ",";
        for (int j =0; j < n; j++){
            datafile << A[i][j] << ", ";
        }
        datafile << endl;
    }
    
    datafile <<"*" << endl;



}
