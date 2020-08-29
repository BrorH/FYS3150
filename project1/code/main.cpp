#include "time.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
using namespace std;

// arguments: n, method, (int)write_solution

int n;
double h;
static double e = 2.718281828459045;

double f(double x) { return 100 * pow(e, -10 * x); }

// ======================================================== //

void LUdcmp_cpu_ram(double *L, double *U) {
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

  // to ram optimize we need to do some tricks with the indexes. The conversions
  // are as follows:
  // L[i][i] -> L[0] = 1
  // L[i][i-1] -> L[i+1]
  // U[i][i] -> U[i+1]
  // U[i][i+1] -> U[0] = -1

  double sum = 0;
  for (int i = 0; i < n; i++) {
    U[i + 1] = (double)(i + 2) / ((double)(i + 1));
    if (i > 0) {
      L[i] = -(double)i / ((double)(i + 1));
    }
  }
}

void LUdcmp_cpu(double **A, double **L, double **U) {
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

void LUdcmp(double **A, double **L, double **U) {
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

// ======================================================== //

double *fwdsub_cpu_ram(double *L, double *b) {
  // The optimized version of forward substitution
  // due to the definition of L (see comment in LUdcmp) many elements of the
  // sum in y will be 0. In fact only the L[i][i-1] and L[i][i]-terms will be
  // non-zero. And since the sum goes from j € [0, i-1], only one term in the
  // sum is nonzero, i.e y[i] = (b[i]-L[i][i-1]*y[i-1])/L[i][i]

  // see LUdcmp for the index trickery regarding ram optimization

  double *y = new double[n];
  y[0] = b[0];
  for (int i = 1; i < n; i++) {
    y[i] = (b[i] -
            L[i + 1] * y[i - 1]); // originally also included /L[i][i]  = / 1
  }
  return y;
}

double *fwdsub_cpu(double **L, double *b) {
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

double *fwdsub(double **L, double *b) {
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

// ======================================================== //

double *bwdsub_cpu_ram(double *U, double *b) {
  // optimized version of bwdsub. Like in fdwsub some of the elements in the
  // sum will be 0, as only U[i][i] and U[i][i+1] will be nonzero. Unlike in
  // fwdsub both of these terms will be included in the sum, and must be
  // accounted for. However U[i][i+1] = -1, so no need to actually look up its
  // value when summing
  //
  // see LUdcmp for the index trickery regarding ram optimization
  double *v = new double[n];
  v[n - 1] = b[n - 1] / U[n];
  double sum = 0;
  for (int i = n - 2; i >= 0; i--) {
    sum += U[i + 1] * v[i];
    sum += -1 * v[i + 1];
    v[i] = (b[i] - sum) / U[i + 1];
    sum = 0;
  }
  return v;
}

double *bwdsub_cpu(double **U, double *b) {
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

double *bwdsub(double **U, double *b) {
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

// ======================================================== //

void construct_A(double **A) {
  // this is the same for non-optimized and cpu-optimized
  for (int i = 0; i < n; i++) {
    A[i][i] = 2;
    if (i != n - 1) {
      A[i + 1][i] = -1;
      A[i][i + 1] = -1;
    }
  }
}

double max_rel_err(double *x, double *v) {
  double exact;
  double epsilon;
  double max = 0;
  for (int i = 1; i < n - 1; i++) {
    exact = 1 - (1 - (double)pow(e, -10)) * x[i] - (double)pow(e, -10 * x[i]);
    epsilon = log10(abs((v[i] - exact) / exact));
    if (epsilon > max)
      max = epsilon;
  }
  cout << "max err = " << max << endl;
  return max;
}

// ======================================================== //

void write_solution(double *v) {
  // writese just the solution points to the solutions.dat file
  ofstream file("solution.dat");
  for (int i = 0; i < n; i++) {
    file << v[i] << ",";
  }

  file.close();
}

void write_specs(float time, double max_err) {
  // writes specifics about the simulations, n, time and max_err
  ofstream file("specs.dat");
  file << n << "," << time << "," << max_err << endl;
  file.close();
}

float printTime(clock_t start, string msg) {
  clock_t now = clock();
  double elapsed = ((now - start) / (double)CLOCKS_PER_SEC);
  cout << msg << ": elapsed time: " << elapsed << endl;
  return elapsed;
}

// ======================================================== //

void solve_original(bool write_sol = false) {
  clock_t start = clock();
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
  construct_A(A);
  LUdcmp(A, L, U);
  double *y = fwdsub(L, b);
  double *v = bwdsub(U, y);

  float done = printTime(start, "done");
  write_specs(done, max_rel_err(x, v));
  if (write_sol) {
    write_solution(v);
  }
  delete[] A;
  delete[] L;
  delete[] U;
  delete[] b;
  delete[] x;
}

void solve_cpu(bool write_sol = false) {
  clock_t start = clock();
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
  construct_A(A);
  LUdcmp_cpu(A, L, U);
  double *y = fwdsub_cpu(L, b);
  double *v = bwdsub_cpu(U, y);

  float done = printTime(start, "done");
  write_specs(done, max_rel_err(x, v));
  if (write_sol) {
    write_solution(v);
  }
  delete[] A;
  delete[] L;
  delete[] U;
  delete[] b;
  delete[] x;
}

void solve_cpu_ram(bool write_sol = false) {
  clock_t start = clock();
  double *L, *U, *b, *x;
  L = new double[n];
  L[0] = 1;
  ;

  // fill U
  U = new double[n + 1];
  U[0] = -1;

  // fill b and x
  x = new double[n];
  b = new double[n];
  for (int i = 0; i < n; i++) {
    x[i] = i * h;
    b[i] = pow(h, 2) * f(x[i]);
  }

  LUdcmp_cpu_ram(L, U);
  double *y = fwdsub_cpu_ram(L, b);
  double *v = bwdsub_cpu_ram(U, y);

  float done = printTime(start, "done");
  double err = max_rel_err(x, v);
  cout << "error: " << err << endl;
  write_specs(done, err);
  if (write_sol) {
    write_solution(v);
  }

  delete[] L;
  delete[] U;
  delete[] b;
  delete[] x;
  delete[] y;
  delete[] v;
}

int main(int argc, char *argv[]) {

  n = pow(10, atoi(argv[1]));
  h = 1 / (double)n;
  string method = argv[2];

  bool write_sol = atoi(argv[3]); // if true  solved data points will be written

  if (method == "original") {
    solve_original(write_sol);
  } else if (method == "cpu") {
    solve_cpu(write_sol);
  } else if (method == "cpu_ram") {
    solve_cpu_ram(write_sol);
  }
  return 0;
}
