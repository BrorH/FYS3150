import numpy as np
import matplotlib.pyplot as plt


def read(filename):
    x = []
    y = []
    with open(filename, "r+") as file:
        lines = file.read().split("\n")
    n = [float(i) for i in lines[0].split("\n")[0].split(" ")][0]
    for line in lines[1:-1]:
        foo = [float(a) for a in line.split(" ")]
        y.append(foo[0])
        x.append(foo[1])
    return np.array(x), np.array(y), n
    # dat = [float(i) for i in lines[1:-1]]
    # return np.array(dat), conf


f = lambda y: 1 - (1 - np.exp(-10)) * y - np.exp(-10 * y)

x, y, n = read("num.dat")
plt.plot(x, y, label="num")
plt.plot(x, f(x), label="ana")
plt.grid()
plt.legend()
plt.show()
