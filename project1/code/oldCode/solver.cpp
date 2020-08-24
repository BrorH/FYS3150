#include "LUsolverPoisson.h"
#include "common.h"
static double e = 2.718281828459045;

int n = 10;
double h = 1 / ((double)n);

double f(double x) { return 100 * pow(e, -10 * x); }

Mat<double> fwdsub(Mat<double> L, Mat<double> b) {
  // performs forward subsititution on the lower triangular matrix L in the
  // sytem Lx=b
  int size = (int)b.size();
  Mat<double> x(size, 1);
  x(0) = b(0) / L(0, 0);
  double sum;
  for (int m = 1; m < size; m++) {
    sum = 0;
    for (int i = 0; i < m; i++) {
      sum += L(m, i) * x(i);
    }
    x(m) = (b(m) - sum) / L(m, m);
  }
  return x;
}
Mat<double> bwdsub(Mat<double> U, Mat<double> b) {
  // performs backwards subsititution on the upper triangular matrix L in the
  // sytem Ux=b
  int size = (int)b.size();
  Mat<double> x(size, 1);
  x(size - 1) = b(size - 1) / U(size - 1, size - 1);
  double sum;
  for (int m = size - 2; m >= 0; m--) {
    sum = 0;
    for (int i = m; i < size; i++) {
      sum += U(m, i) * x(i);
    }
    x(m) = (b(m) - sum) / U(m, m);
  }
  return x;
}

Mat<double> make_b(double *x_i, double (*f)(double)) {
  Mat<double> b(n, 1);
  for (int i = 0; i < n; i++) {
    b(i) = pow(h, 2) * f(x_i[i]);
  }
  return b;
}

int main() {
  Mat<double> A;
  // static const int n = 1000;     // no. of steps (matrix width)
  // static const double h = 0.001; // funcion step length
  // Sete the proper sizes of the objects
  A.set_size(n, n);

  // make x_i and fill it with properly spaced elements on (0,1)
  double x_i[n];
  for (int i = 0; i < n; i++) {
    x_i[i] = (1 / (double)n) * (double)i;
  }

  // get the b_ values
  Mat<double> b_ = make_b(x_i, *f);
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
  solve(A, b_).print();
  // Mat<double> L, U; // L and U in LU decomposition
  // lu(L, U, A);      // do the decomp
  // Mat<double> y = fwdsub(L, b_);
  // Mat<double> v = bwdsub(U, y);
  // // v.print();
  // // for (int i = 0; i < n; i++) {
  // //   cout << x_i[i] << " ";
  // // }
  // // cout << endl;
  // Mat<double> xw;
  // xw.set_size(n, 1);
  // for (int i = 0; i < n; i++) {
  //   xw[i] = x_i[i];
  // }
  // ofstream file1("x.dat");
  // file1 << h << " " << n << endl;
  // xw.save(file1, csv_ascii);
  // file1.close();
  // ofstream file("v.dat");
  // file << h << " " << n << endl;
  // v.save(file, csv_ascii);
  // file.close();
  // // U.print();
  return 0;
}
