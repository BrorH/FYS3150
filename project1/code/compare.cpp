#include "time.h"
#include <algorithm>
#include <armadillo>
#include <fstream>
#include <iostream>

using namespace std;
using namespace arma;

int n = (int)2e4;
double h = 1 / ((double)n);
static double e = 2.718281828459045;

double f(double x) { return 100 * pow(e, -10 * x); }

void construct_b(vec *b, mat *xn) {
  // makes and fills the array b
  for (int i = 0; i < n; i++) {
    (*b)(i) = pow(h, 2) * f((*xn)(i));
  }
}

void construct_A(mat *A) {
  // makes and fills the matrix A
  // note that it is very easy to convert this to a funciton that
  // allows for any values along, and next to, the diagonal!

  // A->set_size(n, n);
  A->fill(0);
  for (int i = 0; i < n; i++) {
    // fill A with the diagonal elements
    (*A)(i, i) = 2;
    if (i != n - 1) {
      (*A)(i + 1, i) = -1;
      (*A)(i, i + 1) = -1;
    }
  }
}
void write(vec *v) {

  ofstream file("num.dat");
  file << n << endl;
  for (int i = 0; i < n; i++) {
    file << (*v)(i) << ",";
  }

  file.close();
}

void printTime(clock_t start, string msg) {

  clock_t now = clock();
  double elapsed = ((now - start) / (double)CLOCKS_PER_SEC);
  cout << msg << ": elapsed time: " << elapsed << endl;
}

int main() {
  clock_t start, finish;
  start = clock();
  cout << "n = " << n << endl;
  mat A(n, n);
  vec x(n);
  vec b(n);
  for (double i = 0; i < n; i++) {
    x(i) = h * i;
  }
  printTime(start, "setup");

  construct_A(&A);
  printTime(start, "A done");

  construct_b(&b, &x);
  printTime(start, "b done");

  vec v = solve(A, b);
  printTime(start, "solved");

  write(&v);
  printTime(start, "writen");

  return 0;
}
