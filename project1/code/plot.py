import numpy as np
import matplotlib.pyplot as plt


def read(filename):
    with open(filename, "r+") as file:
        lines = file.readlines()
    n = int(lines[0].strip())
    v = np.asarray([float(i) for i in lines[1].strip().split(",")[:-1]])
    return v, n

f = lambda y: 1 - (1 - np.exp(-10)) * y - np.exp(-10 * y)

y, n = read("num.dat")
x = np.linspace(0, 1, n)
assert(len(x) == len(y))
print("n = %G" % n)
plt.plot(x, y, label="num")
plt.plot(x, f(x), label="ana")
plt.grid()
plt.legend()
plt.show()
