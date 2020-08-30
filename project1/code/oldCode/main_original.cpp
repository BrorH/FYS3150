#include "time.h"
#include <fstream>
#include <iostream>
#include <math.h>

using namespace std;
static double e = 2.718281828459045;

void LUdcmp(double **A, int n, double **L, double **U) {
  // this is general LU decomposition.
  double sum = 0;
  for (int i = 0; i < n; i++) {
    for (int k = i; k < n; k++) {
      for (int j = 0; j < i; j++) {

        sum += L[i][j] * U[j][k];
      }
      U[i][k] = A[i][k] - sum;
      sum = 0;
    }
    for (int k = i; k < n; k++) {
      if (i == k) {
        L[i][i] = 1;
      }
      for (int j = 0; j < i; j++) {

        sum += L[k][j] * U[j][i];
      }
      L[k][i] = (A[k][i] - sum) / U[i][i];
      sum = 0;
    }
  }
}

double *fwdsub(double **L, double *b, int n) {
  double *y = new double[n];
  y[0] = b[0] / L[0][0];
  double sum = 0;
  for (int i = 1; i < n; i++) {
    for (int j = 0; j < i; j++) {
      sum += L[i][j] * y[j];
    }
    y[i] = (b[i] - sum) / L[i][i];
    sum = 0;
  }
  return y;
}

double *bwdsub(double **U, double *b, int n) {
  double *v = new double[n];
  v[n - 1] = b[n - 1] / U[n - 1][n - 1];
  double sum = 0;
  for (int i = n - 2; i >= 0; i--) {
    for (int j = i; j < n; j++) {
      sum += U[i][j] * v[j];
    }
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

void write_all(double *v, double *x, float time, int n) {
  ofstream file("num.dat");
  file << n << "," << time << endl;
  for (int i = 0; i < n; i++) {
    file << v[i] << ",";
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

float printTime(clock_t start, string msg) {

  clock_t now = clock();
  double elapsed = ((now - start) / (double)CLOCKS_PER_SEC);
  cout << msg << ": elapsed time: " << elapsed << endl;
  return elapsed;
}

int main(int argc, char *argv[]) {
  clock_t start, finish;
  start = clock();
  // int n = (int)1e1;
  int n = pow(10, atoi(argv[1]));
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
  printTime(start, "init");
  construct_A(A, n);
  LUdcmp(A, n, L, U);

  printTime(start, "lu");
  double *y = fwdsub(L, b, n);
  printTime(start, "fwdsub");
  double *v = bwdsub(U, y, n);

  float done = printTime(start, "done");
  write_all(v, x, done, n);

  delete[] A;
  delete[] L;
  delete[] U;
  delete[] b;
  delete[] x;

  return 0;
}
