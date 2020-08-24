import numpy as np
import matplotlib.pyplot as plt


def read(filename):
    with open(filename, "r+") as file:
        lines = file.readlines()
    n = int(lines[0].strip())

    X = np.zeros((2, n - 1))
    for i, line in enumerate(lines[1:-1]):
        foo = [float(a) for a in line.split(" ")]
        X[:, i] = foo
    return *X, n

f = lambda y: 1 - (1 - np.exp(-10)) * y - np.exp(-10 * y)

y, x, n = read("num.dat")
print("n = %G" % n)
plt.plot(x, y, label="num")
plt.plot(x, f(x), label="ana")
plt.grid()
plt.legend()
plt.show()
