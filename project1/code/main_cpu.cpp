#include "time.h"
#include <fstream>
#include <iostream>
#include <math.h>

// OPTIMIZED FOR CPU USAGE

using namespace std;
static double e = 2.718281828459045;

void LUdcmp(double **A, int n, double **L, double **U) {
  // this is the optimized LU decomposition that accounts for the
  // tri-diagonal properties of A.
  // U[i][k] = (u_ik if k=i), (-1 if k=i+1), (0 else) [u_lk def. in book]
  // L[i][k] = (1 if k=i), (l_ik if k=i-1), (0 else) [l_ik def. in book]
  // this optimization is just a de-foorlooping of the complete decomposition.
  // To achieve the best optimizitaion the algorithm may have to be reworked
  // completely from a mathematical standpoint.

  // Study of the patters yield (PROVE THIS!!):
  // u_ik (where k=i) = (i+2)/(i+1) for i >= 0.
  // l_ik (where k=i-1) = -i/(i+1) for i > 0.
  double sum = 0;
  for (int i = 0; i < n; i++) {
    U[i][i] = (double)(i + 2) / ((double)(i + 1));
    U[i][i + 1] = -1;
    L[i][i] = 1;
    if (i > 0) {
      L[i][i - 1] = -(double)i / ((double)(i + 1));
    }
  }
}

double *fwdsub(double **L, double *b, int n) {
  // The optimized version of forward substitution
  // due to the definition of L (see comment in LUdcmp) many elements of the
  // sum in y will be 0. In fact only the L[i][i-1] and L[i][i]-terms will be
  // non-zero. And since the sum goes from j € [0, i-1], only one term in the
  // sum is nonzero, i.e y[i] = (b[i]-L[i][i-1]*y[i-1])/L[i][i]
  double *y = new double[n];
  y[0] = b[0] / L[0][0];
  for (int i = 1; i < n; i++) {
    y[i] = (b[i] - L[i][i - 1] * y[i - 1]) / L[i][i];
  }
  return y;
}

double *bwdsub(double **U, double *b, int n) {
  // optimized version of bwdsub. Like in fdwsub some of the elements in the
  // sum will be 0, as only U[i][i] and U[i][i+1] will be nonzero. Unlike in
  // fwdsubOpt both of these terms will be included in the sum, and must be
  // accounted for. However U[i][i+1] = -1, so no need to actually look up its
  // value when summing
  //
  double *v = new double[n];
  v[n - 1] = b[n - 1] / U[n - 1][n - 1];
  double sum = 0;
  for (int i = n - 2; i >= 0; i--) {
    sum += U[i][i] * v[i];
    sum += -1 * v[i + 1];
    v[i] = (b[i] - sum) / U[i][i];
    sum = 0;
  }
  return v;
}

void construct_A(double **A, int n) {
  for (int i = 0; i < n; i++) {
    A[i][i] = 2;
    if (i != n - 1) {
      A[i + 1][i] = -1;
      A[i][i + 1] = -1;
    }
  }
}

void write(double *v, double *x, int n) {

  ofstream file("num.dat");
  file << n << endl;
  for (int i = 0; i < n; i++) {
    file << v[i] << " " << x[i] << endl;
  }
  file.close();
}

void print2d(double **A, int n) {

  cout << endl;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      cout << A[i][j] << " ";
    }
    cout << endl;
  }
}
void print1d(double *a, int n) {
  for (int i = 0; i < n; i++) {
    cout << a[i] << " ";
  }
  cout << endl;
}

double f(double x) { return 100 * pow(e, -10 * x); }

void printTime(clock_t start) {

  clock_t now = clock();
  double elapsed = ((now - start) / (double)CLOCKS_PER_SEC);
  cout << "elapsed time: " << elapsed << endl;
}

int main() {
  // the process that takes the longest time is the LUdcmp
  // how to make faster:
  //  -Do not define L and U in outter scope (unnecesary mem)
  //  -skip the matrix elements that are 0 anyway
  clock_t start, finish;
  start = clock();
  int n = (int)1e5;
  double h = 1 / (double)n;
  double **A, **L, **U, *b, *x;
  A = new double *[n];
  L = new double *[n];
  U = new double *[n];
  b = new double[n];
  x = new double[n];
  for (int i = 0; i < n; i++) {
    A[i] = new double[n];
    L[i] = new double[n];
    U[i] = new double[n];
    x[i] = i * h;
    b[i] = pow(h, 2) * f(x[i]);
  }
  printTime(start);
  construct_A(A, n);
  LUdcmpOpt(A, n, L, U);
  printTime(start);
  double *y = fwdsub(L, b, n);
  printTime(start);

  double *v = bwdsub(U, y, n);
  printTime(start);
  // write(v, x, n);

  delete[] A;
  delete[] L;
  delete[] U;
  delete[] b;
  delete[] x;
  
  return 0;
}
