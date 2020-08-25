import numpy as np
import matplotlib.pyplot as plt
import subprocess
import sys


def read_big_data(filename):
    with open(filename, "r+") as file:
        lines = file.readlines()
    n, time = [float(i) for i in lines[0].strip().split(",")]
    n = int(n)
    v = np.asarray([float(i) for i in lines[1].strip().split(",")[:-1]])
    return v, n, time


def read_small_data(filename):
    with open(filename, "r+") as file:
        line = file.readline()
    n, time, max_err = [float(i) for i in line.strip().split(",")]
    n = int(n)
    return n, time, max_err


exact = lambda y: 1 - (1 - np.exp(-10)) * y - np.exp(-10 * y)


def compute_max_err(v):
    n = len(v)
    x = np.linspace(0, 1, n)
    u = exact(x)
    eps = np.log10((abs(v[1:-1] - u[1:-1]) / u[1:-1]))
    return max(eps)


def multiple_n(n):
    max_log_n = n
    ns = np.zeros(max_log_n)
    errs = np.zeros(max_log_n)
    times = np.zeros(max_log_n)
    for ni in np.arange(max_log_n):
        print(f"Solving with log10(n) = {ni + 1}")
        subprocess.run(f"./main.out {ni + 1}".split(" "))
        print()

        n, time, max_err = read_small_data("small.dat")

        ns[ni] = np.log10(n)
        times[ni] = time
        errs[ni] = np.log10(max_err)
    plt.scatter(ns, errs)
    plt.xlabel("Number of points [log10]")
    plt.ylabel("Maximum relative error [log10]")
    plt.show()


def single_n(n):
    subprocess.run(f"./main.out {n}".split(" "))
    v, n, time = read_big_data("num.dat")

    x = np.linspace(0, 1, n)
    u = exact(x)
    assert len(x) == len(v)
    max_err = compute_max_err(v)
    print("max err = ", max_err)

    print("n = %G" % n)
    plt.plot(x, v, label="numerical")
    plt.plot(x, u, label="analytic")
    plt.title(f"Maximum relative error log10 = {max_err}")
    plt.grid()
    plt.legend()
    plt.show()


try:
    n = int(sys.argv[1])
except:
    print("Wrong usage. First command line arg must be int, maximum log(n)")
else:
    if "single" in sys.argv:
        single_n(n)
    else:
        multiple_n(n)
