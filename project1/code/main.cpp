#include "LUsolverPoisson.h"
#include <cmath>
static double e = 2.718281828459045;

double f(double x) { return 100 * pow(e, -10 * x); }

int main(int argc, char *argv[])
{
  int exp = atoi(argv[1]);
  int n = pow(10, exp);

  LUsolverPoisson solver(*f, n, 0, 1);
  solver.set_boundary(0, 0);
  solver.compute();
  solver.write("data.txt");
  return 0;
}
