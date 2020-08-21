#include <algorithm>
#include <armadillo>
#include <fstream>
#include <iostream>
using namespace std;
using namespace arma;

class LUsolverPoisson {
  Mat<double> A, b, xn;
  double h, x_min, x_max, u0, u1;
  int n;
  double (*f)(double);

public:
  LUsolverPoisson(double (*_f)(double), double _u0,
                  double _u1) { // Mat<double> _A, Mat<double> _b) {
    f = _f;
    u0 = _u0;
    u1 = _u1;
  }
  void setup(double _h, double _n, Mat<double> _xn) {
    h = _h;
    n = _n;
    xn = _xn;
    // construct_A();
    // construct_b();
    // if a vector of x_n values are provided, the yare chosen

    // set the parameters of the system
  }
  void setup(double _h, int _n, double _x_min, double _x_max) {
    // if a vector of x_n values are NOT provided, x is evenly distributed on
    // range
    double delta = _x_max - _x_min;
    Mat<double> _xn;
    _xn.set_size(_n, 1);
    for (int i = 0; i < _n; i++) {
      _xn(i) = (delta / (double)_n) * (double)i;
    }
    setup(_h, _n, _xn);
  }

  Mat<double> compute(bool write = true) {
    // runs the numerical integration and writes to the files v.dat and x.dat
    Mat<double> L, U;
    lu(L, U, A);
    Mat<double> y = fwdsub(L, b);
    Mat<double> v = bwdsub(U, y);
    if (write) {
      ofstream xfile("x.dat");
      xfile << h << " " << n << endl;
      xn.save(xfile, csv_ascii);
      xfile.close();
      ofstream vfile("v.dat");
      vfile << h << " " << n << endl;
      v.save(vfile, csv_ascii);
      vfile.close();
    }
    return v;
  }

private:
  void construct_A() {
    // makes and fills the matrix A
    // note that it is very easy to convert this to a funciton that
    // allows for any values along, and next to, the diagonal!
    A.set_size(n, n);
    A.fill(0);
    for (int i = 0; i < n; i++) {
      // fill A with the diagonal elements
      A(i, i) = 2;
      if (i != n - 1) {
        A(i + 1, i) = -1;
        A(i, i + 1) = -1;
      }
    }
  }
  void construct_b() {
    // makes and fills the array b
    b.set_size(n, 1);
    for (int i = 0; i < n; i++) {
      b(i) = pow(h, 2) * f(xn(i));
    }
  }

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
};
