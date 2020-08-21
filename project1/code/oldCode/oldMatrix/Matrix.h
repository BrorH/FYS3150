#include "common.h"

class Matrix {
  // create and work with matrices of any shape.

public:
  vector<vector<double>> vals;

  int n, m;
  Matrix(int _n, int _m) {
    // M, N are size of vector
    n = _n;
    m = _m;
    vals.resize(n);
    for (int i = 0; i < n; ++i)
      vals[i].resize(_m);
  }
  ~Matrix() { delete[] vals }
  double operator()(int i, int j) {
    // return value at indexes i,j
    return vals[i][j];
  }

  void set(int i, int j, double val) {
    // sets array element at i,j to val
    vals[i][j] = val;
  }

  void print() {
    for (int i = 0; i < n; i++) {
      printf("|");
      for (int j = 0; j < m; j++) {
        printf("%.2f", vals[i][j]);
        if (j != m - 1) {
          printf(" ");
        }
        // cout << round(this->operator()(i, j)) << " ";
      }
      cout << "|" << endl;
    }
  }
};

Matrix operator*(Matrix M1, Matrix M2) {
  // calculates the product M1*M2
  assert(M1.m == M2.n); // must be right shape
  double element;
  Matrix mulRes(M1.n, M2.m);
  for (int i = 0; i < M1.n; i++) {
    for (int j = 0; j < M2.m; j++) {
      element = 0;
      for (int k = 0; k < M1.m; k++) {
        element += M1(i, k) * M2(k, j);
      }
      mulRes.set(i, j, element);
    }
  }
  delete[] element;
  return mulRes;
}

Matrix operator*(Matrix M, double a) {
  // calculates product between matrix and number a
  Matrix mulRes(M.n, M.m);
  for (int i = 0; i < M.n; i++) {
    for (int j = 0; j < M.m; j++) {
      mulRes.set(i, j, a * (M.vals[i][j]));
    }
  }
  return mulRes;
}

Matrix operator*(double a, Matrix M) {
  // makes sure you can do M*a as well as a*M
  return M * a;
}
