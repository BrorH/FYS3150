import numpy as np
import matplotlib.pyplot as plt
import subprocess
import pandas as pd

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


def comparisons(_n, type="original", solPlot=False, table=True, errPlot=False, timePlot=False):
    # compares the runs of the original solver (unoptimized) for _n values.
    # Is able to print a latex table of the specs.dat file
    assert type in ["original", "cpu", "cpu_ram"], f'type must be "original", "cpu" or "cpu_ram", not {type}'
    print(f"comparing. type: {type}, n: {_n}, plot: {solPlot}, table: {table}")
    if solPlot:
        plot = 1
    else:
        plot = 0
    specs = []
    if solPlot:
        fig, ax = plt.subplots(3, 1, sharex=True, sharey=True)
    for n in range(1, _n + 1):
        subprocess.run(f"./main.out {n} {type} {plot}".split(" "))
        v = read_solution()
        N, time, max_err = read_specs()
        specs.append([N, time, max_err])
        if solPlot:
            x = np.linspace(0, 1, N)
            ax_ = ax[n - 1]
            assert len(x) == len(v)
            ax_.plot(x, v, label="$\log (n) = %i$" % int(np.log10(N)), c="r")
            ax_.plot(x, exact(x), "--", label="Analytical", c="k")
            ax_.grid(True)
            ax_.legend()
    if solPlot:
        ax[0].set_title(f"General solution algorithm, comparison for $n$ steps")
        ax[1].set_ylabel("$u(x)$")
        ax[2].set_xlabel("$x$")
        plt.subplots_adjust(hspace=0.1)
        plt.savefig(f"figures/sol.{type}.{_n}.png")

        plt.show()

    if table:
        specdict = {}
        specdict[r"$\log (n)$"] = [int(np.log10(a[0])) for a in specs]
        specdict[r"$t$ [s]"] = [a[1] for a in specs]
        specdict[r"$\epsilon_{max}$"] = [a[2] for a in specs]
        df = pd.DataFrame(specdict)
        print(df.to_latex(index=False, float_format="%.2e", label="LABEL HERE", caption="CAPTION HERE", escape=False, column_format="c" * _n))
    if errPlot or timePlot:
        # since the plotting of bot the errors and the cputime is similar, they can be handled at the same time
        # this means that only one can be done at a time, but meh
        if errPlot:
            y = np.log10([a[2] for a in specs])
            label = "error"
            fileExt = "err"
        else:
            y = np.log10([a[1] for a in specs])
            label = "CPU time"
            fileExt = "time"

        n = np.log10([a[0] for a in specs])
        plt.plot(n, y, "o", label=label)
        a, b = tuple(np.polyfit(n, y, deg=1))
        x = np.linspace(n[0], n[-1], 1000)

        plt.plot(x, a * x + b, label="fitted pol, a = %.2e" % a)
        plt.xlabel
        plt.grid()
        plt.legend()
        plt.savefig(f"figures/{fileExt}.{type}.{_n}.png")
        plt.show()


# comparisons(3, "original", producePlot=True, produceTable=True)
comparisons(3, "original", solPlot=True)
# comparisons(5, "cpu", errPlot=True)
# comparisons(5, "cpu", timePlot=True)
# plot_comparisons_original(3, True)
# table_comparisons_cpu(6)
# table_comparisons_cpu_ram(6)
