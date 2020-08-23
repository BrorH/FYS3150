#include <algorithm>
#include <armadillo>
#include <fstream>
#include <iostream>
using namespace std;
using namespace arma;

class LUsolverPoisson
{
  mat A, b, xn, v;
  double h, x_min, x_max, u0, u1;
  int n;
  double (*f)(double);

public:
  LUsolverPoisson(double (*_f)(double), mat _xn)
  {
    f = _f;
    xn = _xn;
    n = xn.n_rows;
    x_min = xn[0];
    x_max = xn[n - 1];

    // construct_A();
    // construct_b();
    // if a vector of x_n values are provided, they aare chosen

    // set the parameters of the system
  }
  LUsolverPoisson(double (*_f)(double), int _n, double _x_min, double _x_max)
  {
    // if a vector of x_n values are NOT provided, x is evenly distributed on
    // range
    f = _f;
    n = _n;
    h = (_x_max - _x_min) / (double)n;
    xn.set_size(n, 1);
    for (double i = 0; i < n; i++)
    {
      xn(i) = h * i;
    }
  }
  void set_boundary(double _u0, double _u1)
  {
    u0 = _u0;
    u1 = _u1;
    construct_A();
    construct_b();
  }

  void compute(bool write = true)
  {
    // runs the numerical integration and writes to the files v.dat and x.dat
    mat L, U;
    lu(L, U, A);
    mat y = fwdsub(L, b);
    v = bwdsub(U, y);
  }
  void write(string filename)
  {
    ofstream xfile("x.dat");
    xfile << h << " " << n << endl;
    xn.save(xfile, csv_ascii);
    xfile.close();
    ofstream vfile("v.dat");
    vfile << h << " " << n << endl;
    v.save(vfile, csv_ascii);
    vfile.close();
  }

private:
  void construct_A()
  {
    // makes and fills the matrix A
    // note that it is very easy to convert this to a funciton that
    // allows for any values along, and next to, the diagonal!
    A.set_size(n, n);
    A.fill(0);
    for (int i = 0; i < n; i++)
    {
      // fill A with the diagonal elements
      A(i, i) = 2;
      if (i != n - 1)
      {
        A(i + 1, i) = -1;
        A(i, i + 1) = -1;
      }
    }
  }
  void construct_b()
  {
    // makes and fills the array b
    b.set_size(n, 1);
    for (int i = 0; i < n; i++)
    {
      b(i) = pow(h, 2) * f(xn(i));
    }
  }

  mat fwdsub(mat L, mat b)
  {
    // performs forward subsititution on the lower triangular matrix L in the
    // sytem Lx=b
    int size = b.n_rows;
    mat x(size, 1);
    x(0) = b(0) / L(0, 0);
    double sum;
    for (int i = 1; i < size; i++)
    {
      sum = 0;
      for (int j = 0; j < i; j++)
      {
        sum += L(i, j) * x(j);
      }
      x(i) = (b(i) - sum) / L(i, i);
    }
    return x;
  }
  mat bwdsub(mat U, mat b)
  {
    // performs backwards subsititution on the upper triangular matrix L in the
    // sytem Ux=b
    int size = b.n_rows;
    mat x(size, 1);
    x(size - 1) = b(size - 1) / U(size - 1, size - 1);
    double sum;
    for (int i = size - 2; i >= 0; i--)
    {
      sum = 0;
      for (int j = i; j < size; j++)
      {
        sum += U(i, j) * x(j);
      }
      x(i) = (b(i) - sum) / U(i, i);
    }
    return x;
  }
};
