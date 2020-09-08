#include "time.h"
#include <algorithm>
#include <armadillo>
#include <fstream>
#include <iostream>

using namespace std;
using namespace arma;

double h;
static double e = 2.718281828459045;

double f(double x) { return 100 * pow(e, -10 * x); }

void construct_b(vec *b, mat *xn, int n)
{
  // makes and fills the array b
  for (int i = 0; i < n; i++)
  {
    (*b)(i) = pow(h, 2) * f((*xn)(i));
  }
}

void construct_A(mat *A, int n)
{
  // makes and fills the matrix A
  // note that it is very easy to convert this to a funciton that
  // allows for any values along, and next to, the diagonal!

  // A->set_size(n, n);
  A->fill(0);
  for (int i = 0; i < n; i++)
  {
    // fill A with the diagonal elements
    (*A)(i, i) = 2;
    if (i != n - 1)
    {
      (*A)(i + 1, i) = -1;
      (*A)(i, i + 1) = -1;
    }
  }
}
void write_specs(int n, double avg, double err)
{
  ofstream file("arma.dat");
  file << n << "," << avg << "," << err << endl;
  file.close();
}

void printTime(clock_t start, string msg)
{

  clock_t now = clock();
  double elapsed = ((now - start) / (double)CLOCKS_PER_SEC);
  cout << msg << ": elapsed time: " << elapsed << endl;
}

double max_rel_err(vec *x, vec *v, int n)
{
  double exact;
  double epsilon;
  double max = 0;
  for (int i = 1; i < n - 1; i++)
  {
    exact = 1 - (1 - (double)pow(e, -10)) * x->operator[](i) -
            (double)pow(e, -10 * x->operator[](i));
    exact = exact / h / h;
    epsilon = log10(abs((v->operator[](i) - exact) / exact));
    epsilon = abs(epsilon);
    if (epsilon > max)
      max = epsilon;
  }
  // cout << "max err = " << max << endl;
  return max;
}
void runs(int n, int m = 5)
{
  // does m timed runs of scale n and averages their time
  clock_t start;
  double avg;
  mat A(n, n);
  vec x(n);
  vec b(n);
  vec v;

  for (int i = 0; i < m; i++)
  {
    start = clock();
    for (double i = 0; i < n; i++)
    {
      x(i) = h * i;
    }

    construct_A(&A, n);

    construct_b(&b, &x, n);
    v = solve(A, b);

    // stop = clock();
    avg += ((clock() - start) / (double)CLOCKS_PER_SEC);
  }
  avg = avg / (double)m;
  double err = max_rel_err(&x, &v, n);
  cout << "n = " << n << ", t = " << avg << " s. err = " << err << endl;
  write_specs(n, avg, err);
}

int main(int argc, char *argv[])
{
  int n = pow(10, atoi(argv[1]));
  h = 1 / ((double)n);
  // double res;
  int m = 5;
  runs(n, m);
  //
  // clock_t start, finish;
  // start = clock();
  // cout << "n = " << n << endl;
  // mat A(n, n);
  // vec x(n);
  // vec b(n);
  // for (double i = 0; i < n; i++) {
  //   x(i) = h * i;
  // }
  // printTime(start, "setup");
  //
  // construct_A(&A);
  // printTime(start, "A done");
  //
  // construct_b(&b, &x);
  // printTime(start, "b done");
  //
  // vec v = solve(A, b);
  // printTime(start, "solved");
  //
  // write(&v);
  // printTime(start, "writen");

  return 0;
}
