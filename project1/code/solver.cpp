#include <algorithm>
#include <armadillo>
#include <cmath>
#include <iostream>
using namespace std;
using namespace arma;

static double e = 2.71828;
double f(double x) { return 100 * pow(e, -10 * x); }

Mat<double> make_b(int n, double h, double *x_i, double (*f)(double)) {
  Mat<double> b(n, 1);
  for (int i = 0; i < n; i++) {
    b(i) = pow(h, 2) * f(x_i[i]);
  }
  return b;
}

int main() {
  Mat<double> A, v;
  static const int n = 5;      // no. of steps (matrix width)
  static const double h = 0.1; // funcion step length

  // Sete the proper sizes of the objects
  A.set_size(n, n);
  v.set_size(n, 1);

  // make x_i and fill it with properly spaced elements on (0,1)
  double x_i[n];
  for (int i = 0; i < n; i++) {
    x_i[i] = (1 / (double)n) * (double)i;
  }

  // get the b_ values
  Mat<double> b_ = make_b(n, h, x_i, *f);
  // the below dat arrays can be filled with any datapoints preffered
  double adat[n - 1], bdat[n], cdat[n - 1];
  fill_n(adat, n - 1, -1);
  fill_n(bdat, n, 2);
  fill_n(cdat, n - 1, -1);
  A.fill(0);
  for (int i = 0; i < n; i++) {
    // fill A with the diagonal elements (and the ones above and below)
    A(i, i) = bdat[i];
    if (i != n - 1) {
      A(i + 1, i) = adat[i];
      A(i, i + 1) = cdat[i];
    }
  }
  cout << A << endl;

  Mat<double> L, U; // L and U in LU decomposition
  lu(L, U, A);      // do the decomp
  L.print();
  U.print();
  return 0;
}
