#include "solver.h"


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
    V = new double *[n]; // To store eigenvectors as columns
    for(int i = 0; i < n; i++){V[i] = new double[n];}
    for(int i = 0; i < n; i++){
        V[i][i] = 1; // V is initially identity matrix
    }
   
}
Solver::~Solver(){
    // Memory handling
    delete[] V, A, A_init;
}

void Solver::fillA(double **_A)
{   /*
    Fills the passed double pointer, square matrix _A with the correct values
    */
    double off_diag = -1*pow(h,-2);
    for (int i = 0; i < n; i++){ _A[i] = new double[n];}
    for (int i = 0; i < n; i++){
        _A[i][i] = diag[i];
        if (i != n - 1){
            _A[i + 1][i] = off_diag;
            _A[i][i + 1] = off_diag;
        }
    }
}

void Solver::checkOrthog(double tol = 1e-5){
    // Matrix multiplication of V^T * V and summing together
    double resVTV = 0;
    for(int i = 0; i < n; i ++){
        for(int j = 0; j < n; j ++){
            for(int k = 0; k < n; k ++){
                resVTV += V[k][i]*V[k][j]; // calculate cummulative sum of VTV
            }
        }
    }
    if (abs(resVTV-n) > tol){
        cout << "Not orthogonal! Got sum " << resVTV << endl;
    }

}

void Solver::checkSymmetry(double tol = 1e-5){
    // checks the symmetry of A
    double diff;
    for(int i = 0; i < n; i ++){
        for(int j = 0; j < n; j ++){
            diff = abs(A[i][j]-A[j][i]);
            if( diff > tol){
                cout << "Not symmetric! Got diff of " << diff << endl;
            }
        }
    }
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

void Solver::rotateA(double c, double s){
    /*
    Performs the unitary rotation on A and interates the eigenmatrix V.
    args c and s are cos and sin from the Jacobi analysis.
    */
    int k = maxidx[0]; int l = maxidx[1]; // retrieve indecies from global array
    double A_kk, A_ll, A_ik, A_il, V_ik, V_il; // define variables for the indecies as the changes are static
    
    A_kk = A[k][k]; A_ll = A[l][l];

    A[k][k]=A_kk*pow(c,2) - 2*A[k][l]*c*s + A_ll*pow(s,2); // blindly follow Jacobi theory
    A[l][l]=A_ll*pow(c,2) + 2*A[k][l]*c*s + A_kk*pow(s,2); // --""--
    A[k][l] = 0; // manually set to 0
    A[l][k] = 0; // --""--
    for(int i = 0; i < n ; i ++){

        if((i != k) && (i != l)){ // dont change the already changed values

            A_ik = A[i][k]; A_il = A[i][l];
            // blindly follow Jacobi theory
            A[i][k] = A_ik*c - A_il*s;
            A[k][i] = A[i][k];
            A[i][l] = A_il*c + A_ik*s;
            A[l][i] = A[i][l];
        }
        // then update the eigenvectors accordingly
        V_ik = V[i][k];
        V_il = V[i][l];
        V[i][k] = V_ik*c - V_il*s;
        V[i][l] = V_il*c + V_ik*s;
    }
}


void Solver::solve()
{
    /*
    Solves the system

    MORE HERE

    */
    int k, l;
    double t, c, s, tau, max_A;

    transforms = 0;
    max_A = 2*eps; // Set a higher start value

    while (max_A > eps) // loop while the maximum non-diagonal value in A is below threshold
    {
        findMaxIdx(); // find max non-idagonal index and update dynamic class array "maxidx"
        k = maxidx[0]; l = maxidx[1]; // retrieve new max indecies
        max_A = pow( A[k][l], 2); // max non-diag value squared
       
        tau = (A[l][l] - A[k][k])/(2*A[k][l]);

        if(tau > 0){ // pick the smallest and most fitting theta
            t = 1/(tau+sqrt(1+pow(tau,2)));
        } else {
            t = 1/(tau-sqrt(1+pow(tau,2)));
        }
        c = 1/sqrt(1+pow(t,2));
        s = c*t;

        rotateA(c,s); // perform the unitary transform

        if(transforms%10 == 0){
            // Do unit tests every so often.
            checkOrthog();
            checkSymmetry();
            // ONE MORE!
        }
        transforms ++; // update counter
    }
    //cout << "n: " << n<< " transforms: " <<transforms << endl;
}
void Solver::write(){
    // Simple data-writing system. Raw data not meant to be read by humans. 
    
    ofstream datafile;
    datafile.open("data.dat", ios::app);
    datafile << n << "," <<rho_max << ","<< eps << ","<< transforms<< endl;
    for(int i = 0; i < n; i++){datafile << diag[i] << ",";}
    datafile << endl;
    for(int i = 0; i < n; i++){
        datafile << A[i][i] << ",";
        for (int j =0; j < n; j++){
            datafile << V[j][i] << ",";
        }
        datafile << endl;
    }
    
    datafile <<"*" << endl; // separator



}
