import numpy as np
import matplotlib.pyplot as plt
import subprocess
import pandas as pd
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


def comparisons(_n, type="original", solPlot=False, table=False, errPlot=False, timePlot=False):
    # compares the runs of the original solver (unoptimized) for _n values.
    # Is able to print a latex table of the specs.dat file
    print(f"comparing. type: {type}, n: {_n}, plot: {solPlot}, table: {table}")
    specs = []
    if solPlot:
        fig, ax = plt.subplots(_n, 1, sharex=True, sharey=True)
    for n in range(1, _n + 1):
        subprocess.run(f"./main.out {n} {type} {int(solPlot)}".split(" "))
        N, time, max_err = read_specs()
        specs.append([N, time, max_err])
        if solPlot:# and n < 4:
            v = read_solution()
            x = np.linspace(0, 1, N)
            ax_ = ax[n - 1]
            assert len(x) == len(v)
            ax_.plot(x, v, label="$\log (n) = %i$" % int(np.log10(N)), c="r")
            ax_.plot(x, exact(x), "--", label="Analytical", c="k")
            ax_.grid(True)
            ax_.legend()
    if solPlot:
        # ax[0].set_title(f"General solution algorithm, comparison for {n} steps")
        ax[1].set_ylabel("$u(x)$")
        ax[2].set_xlabel("$x$")
        plt.subplots_adjust(hspace=0.1)
        plt.savefig(f"f../igures/sol.{type}.{_n}.png")

        plt.show()

    if table:
        specdict = {}
        specdict[r"$\log (n)$"] = [int(np.log10(a[0])) for a in specs]
        specdict[r"$t$ [s]"] = [a[1] for a in specs]
        specdict[r"$\epsilon_{max}$"] = [a[2] for a in specs]
        df = pd.DataFrame(specdict)
        print("\n" + df.to_latex(index=False, float_format="%.2e", label=f"tab:{type}", caption=type, escape=False, column_format="c" * _n))
    if errPlot:
        n = np.log10([a[0] for a in specs])
        x = np.linspace(n[0], n[-1], 1000)
        y = np.log10([a[2] for a in specs])

        a, b = tuple(np.polyfit(n, y, deg=1))

        plt.plot(x, a * x + b, "--", c="k", label="fitted line, a = %.2e" % a)

        plt.plot(n, y, "o", c="r", label="$\epsilon_{max}$")

        plt.xlabel("$\log (n)$")
        plt.ylabel("$\log (\epsilon_{max})$")
        plt.grid()
        plt.legend()
        plt.xticks(range(1, _n + 1))
        plt.savefig(f"../figures/err.{type}.{_n}.png")
        plt.show()
    if timePlot:
        y = np.log10([a[1] for a in specs])
        n = np.log10([a[0] for a in specs])
        x = np.linspace(n[0], n[-1], 1000)
        a, b = tuple(np.polyfit(n, y, deg=1))

        plt.plot(x, a * x + b, "--", c="k", label="fitted line, a = %.2e" % a)

        plt.plot(n, y, "o", c="r", label=f"$t$")

        plt.xlabel("$\log (n)$")
        plt.ylabel(f"$\log (t)$")
        plt.grid()
        plt.legend()
        plt.xticks(range(1, _n + 1))
        plt.savefig(f"../figures/time.{type}.{_n}.png")
        plt.show()


params = {"solPlot": False, "errPlot": False, "timePlot": False, "table": False}
args = sys.argv[1:]
try:
    n = int(args[0])
except:
    print("error, first argument must be log(n), aka type int")
    sys.exit(1)
type = args[1]
assert args[1] in ["original", "optimized"], f'type must be "original" or "optimized", not {type}'
for arg in args[2:]:
    try:
        params[arg] = True
    except:
        pass

comparisons(n, type, solPlot=params["solPlot"], errPlot=params["errPlot"], timePlot=params["timePlot"], table=params["table"])
