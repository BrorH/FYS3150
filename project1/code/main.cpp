#include "LUsolverPoisson.h"
static double e = 2.718281828459045;

double f(double x) { return 100 * pow(e, -10 * x); }

int main() {
  double h = 0.001;
  int n = 1000;

  LUsolverPoisson solver(*f, 0, 0);
  solver.setup(h, n, 0, 1);
  // solver.compute();
  return 0;
}
