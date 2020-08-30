import numpy as np
import matplotlib.pyplot as plt
import subprocess
import sys

exact = lambda y: 1 - (1 - np.exp(-10)) * y - np.exp(-10 * y)


def read_solution():
    with open("solution.dat", "r+") as file:
        lines = file.readlines()
    v = np.asarray([float(i) for i in lines[0].strip().split(",")[:-1]])
    return v


def read_specs():
    with open("specs.dat", "r+") as file:
        line = file.readlines()[0].strip().split(",")
    n = int(line[0])
    time = float(line[1])
    max_err = float(line[2])
    return n, time, max_err


# ================================== #


def plot_comparisons_original(_n):
    for n in range(1, _n + 1):
        subprocess.run(f"./main.out {n} original 1".split(" "))
        v = read_solution()
        n, time, max_err = read_specs()
        x = np.linspace(0, 1, n)
        u = exact(x)
        assert len(x) == len(v)
        plt.plot(x, v, label="$\log (n) = %i$, $\epsilon_{max}$ = %.2e" % (int(np.log10(n)), max_err))

    plt.plot(x, u, "--", label="Analytical solution")
    plt.title(f"General solution algorithm, comparison")
    plt.xlabel("$x$")
    plt.ylabel("$u(x)$")
    plt.grid()
    plt.legend()
    plt.savefig("figures/generalAlgComparison.png")
    plt.show()


plot_comparisons_original(3)
