#include "time.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
using namespace std;

/*
Program arguments:
 - n (int): size of matrix (in log10)
 - method ("original" or "optimized"): method of solution
 - write_solution (boolean int): write all solution points to solutions.dat. Takes long time for high n

 Written by
 - Bror Hjemgaard, bahjemga@uio.no
 - Håkon Olav Torvik, haakooto@uio.no

 In most loops, the individual FLOPS are noted. Condition statements are omitted.
*/

int n;
long double h;
static long double e = 2.718281828459045;

double f(double x) { return 100 * pow(e, -10 * x); } // 3 FLOPS

// ======================================================== //

void LUdcmp_optimized(double *L, double *U)
{
  /*
  Optimized LU decomposition that accounts for the tri-diagonal properties of A.
  The analytic expressions for L and U are in the paper. Efficiency O(n)
  */

  for (double i = 0; i < n; i++) // n iterations: 6n FLOPS
  {
    U[(int)i + 1] = (i + 2) / (i + 1); // 3 FLOPS
    if (i > 0)
    {
      L[(int)i] = -i / (i + 1); // 3 FLOPS
    }
  }
}

void LUdcmp(double **A, double **L, double **U)
{
  /*
  General LU-decomposition for matrix A.
  This is methodically identical to the mathematical decomposition,
  and according to Wikipedia using this algorithm requires 2/3 n^3 FLOPS: O(n^3).
  (https://en.wikipedia.org/wiki/LU_decomposition#Using_Gaussian_elimination)
  */

  double sum = 0;
  int count = 0;  // FLOP-counter
  for (int i = 0; i < n; i++)
  {
    for (int k = i; k < n; k++)
    {
      for (int j = 0; j < i; j++)
      {
        sum += L[i][j] * U[j][k];
        count += 2;
      }
      count ++;
      U[i][k] = A[i][k] - sum;
      sum = 0;
    }
    for (int k = i; k < n; k++)
    {
      if (i == k)
      {
        L[i][i] = 1;
      }
      for (int j = 0; j < i; j++)
      {
        sum += L[k][j] * U[j][i];
        count += 2;
      }
      count += 2;
      L[k][i] = (A[k][i] - sum) / U[i][i];
      sum = 0;
    }
  }
cout << "Count for FLOPS in LU for n = " << n << " is " << count << endl;
}

// ======================================================== //

double *fwdsub_optimized(double *L, double *b)
{
  /*
  Optimized forward substitution. Efficiency O(n)
  */

  double *y = new double[n];
  y[0] = b[0];
  for (int i = 1; i < n; i++) // n-1 iterations: 2n-2 FLOPS
  {
    y[i] = b[i] - L[i] * y[i - 1]; //2 FLOPS
  }
  return y;
}
double *fwdsub(double **L, double *b)
{
  /*
  General forward substitution. Efficiency O(n^2)
  */

  double *y = new double[n];
  y[0] = b[0] / L[0][0]; // 1 FLOP
  double sum = 0;
  for (int i = 1; i < n; i++) // n iterations: n^2+3n FLOPS

  {
    for (int j = 0; j < i; j++) // (n^2+n)/2 iterations: (n^2+n)/2 * 2 FLOPS
    {
      sum += L[i][j] * y[j]; // 2 FLOPS
    }
    y[i] = (b[i] - sum) / L[i][i]; // 2 FLOPS
    sum = 0;
  }
  return y;
}

// ======================================================== //

double *bwdsub_optimized(double *U, double *b)
{
  /*
  Optimized backwards substitution. Efficiency O(n)
  */

  double *v = new double[n];
  v[n - 1] = b[n - 1] / U[n];      // 1 FLOP
  for (int i = n - 2; i >= 0; i--) // n-2 iterations: 2n - 4 FLOPS
  {
    v[i] = (b[i] + v[i + 1]) / U[i + 1]; // 2 FLOPS
  }
  return v;
}

double *bwdsub(double **U, double *b)
{
  /*
  General backwards substitution. Efficiency O(n^2)
  */

  double *v = new double[n];
  v[n - 1] = b[n - 1] / U[n - 1][n - 1]; // 1 FLOP
  double sum = 0;
  for (int i = n - 2; i >= 0; i--) // n-2 iterations: n^2+3n FLOPS
  {
    for (int j = i; j < n; j++) //  (n^2+n)/2 *2 FLOPS
    {
      sum += U[i][j] * v[j]; // 2 FLOPS
    }
    v[i] = (b[i] - sum) / U[i][i]; // 2 FLOPS
    sum = 0;
  }
  return v;
}

// ======================================================== //

void construct_A(double **A)
{
  // Fills A with its matrix elements
  for (int i = 0; i < n; i++)
  {
    A[i][i] = 2;
    if (i != n - 1)
    {
      A[i + 1][i] = -1;
      A[i][i + 1] = -1;
    }
  }
}

double max_rel_err(double *x, double *v)
{
  /*
  Computes the max relative error
  */
  double exact, epsilon, max;
  max = 0;
  for (int i = 1; i < n - 1; i++)
  {
    exact = (1 - (1 - (double)pow(e, -10)) * x[i] - (double)pow(e, -10 * x[i])) / pow(h, 2);
    epsilon = abs((v[i] - exact) / exact);
    if (epsilon > max)
      max = epsilon;
  }
  return abs(log10(max));
}

// ======================================================== //

void write_solution(double *v)
{
  // Writes the contents of v to solutions.dat. Takes some time for large n
  ofstream file("solution.dat");
  for (int i = 0; i < n; i++)
  {
    file << v[i] << ",";
  }
  file.close();
}

void write_specs(float time, double max_err)
{
  // Writes n, elapsed time and max_err to specs.dat
  ofstream file("specs.dat");
  file << n << "," << time << "," << max_err << endl;
  file.close();
}

float getTime(clock_t start)
{
  // returns the passed time from given start event
  clock_t now = clock();
  double elapsed = ((now - start) / (double)CLOCKS_PER_SEC);
  return elapsed;
}

// ======================================================== //

void solve_original(bool write_sol = false)
{
  /*
  Solves the system with the general, non-optimized method
  */
  double **A, **L, **U, *b, *x;
  A = new double *[n];
  L = new double *[n];
  U = new double *[n];
  b = new double[n];
  x = new double[n];
  for (int i = 0; i < n; i++)
  {
    A[i] = new double[n];
    L[i] = new double[n];
    U[i] = new double[n];
    x[i] = i * h;
    b[i] = pow(h, 2) * f(x[i]);
  }
  construct_A(A); // Fill A

  clock_t start = clock(); // start timer
  LUdcmp(A, L, U);         // Solve general LU decomposition
  double *y = fwdsub(L, b);
  double *v = bwdsub(U, y);
  float done = getTime(start);

  write_specs(done, max_rel_err(x, v));
  if (write_sol)
  {
    write_solution(v);
  }
  for (int i = 0; i < n; i++)
  {
    delete[] A[i];
    delete[] L[i];
    delete[] U[i];
  }
  // To avoid memory leak
  delete[] A;
  delete[] L;
  delete[] U;
  delete[] y;
  delete[] v;
  delete[] b;
  delete[] x;
}

void solve_optimized(bool write_sol = false)
{
  /*
  Solves the system for the special diagonal matrix with an optimized algorithm.
  */
  double *L, *U, *b, *x;
  L = new double[n];
  L[0] = 1; // The first element of L is its constant diagonal element 1

  // fill U
  U = new double[n + 1];
  U[0] = -1; // The first element of U is its constant upper diagonal element -1

  // fill b and x
  x = new double[n];
  b = new double[n];
  for (int i = 0; i < n; i++)
  {
    x[i] = i * h;
    b[i] = pow(h, 2) * f(x[i]);
  }

  LUdcmp_optimized(L, U);  // Optimized LU-decomposition
  clock_t start = clock(); // Start time at start of numerical computation
  double *y = fwdsub_optimized(L, b);
  double *v = bwdsub_optimized(U, y);
  float done = getTime(start);

  double err = max_rel_err(x, v);
  write_specs(done, err);
  if (write_sol)
  {
    write_solution(v);
  }
  delete[] L;
  delete[] U;
  delete[] b;
  delete[] x;
  delete[] y;
  delete[] v;
}

int main(int argc, char *argv[])
{
  n = pow(10, atoi(argv[1])); // Read n from command line
  h = 1 / (long double)n;
  string method = argv[2]; // Solve either optimized or general
  cout << "Solving " << method << " for log(n)= " << argv[1] << endl;

  bool write_sol = atoi(argv[3]); // if true solved data points will be written to file

  if (method == "original")
  {
    solve_original(write_sol);
  }
  else if (method == "optimized")
  {
    solve_optimized(write_sol);
  }
  return 0;
}
