import numpy as np
import matplotlib.pyplot as plt
import subprocess
import pandas as pd
import sys
import os

font = {"family": "DejaVu Sans", "weight": "bold", "size": 22}
plt.rc("font", **font)


def pushFile(filename):
    path = os.path.abspath(filename)
    print(path)
    subprocess.run(f"git add {path}".split())
    subprocess.run(f'git commit -m "automatic_figure_update" --quiet'.split())
    subprocess.run(f"git push --quiet".split())
    print("done")


push = "push" in sys.argv
N = int(sys.argv[1])


n = np.zeros(N)
t = np.zeros(N)
e = np.zeros(N)

for i in range(N):
    subprocess.run(f"./arma.out {i + 1}".split())

    with open("arma.dat", "r") as file:
        _n, _t, _e = file.readline().strip().split(",")
    n[i] = np.log10(int(_n))
    t[i] = np.log10(float(_t))
    e[i] = np.log10(float(_e))

x = np.linspace(n[0], n[-1], 1000)
at, bt = np.polyfit(n, t, 1)

plt.plot(x, at * x + bt, "--", c="k", lw=3, label="fitted line, a = %.2e" % at)
plt.plot(n, t, "ro", ms=15, label="$t$")
plt.xlabel("$\log (n)$")
plt.ylabel(f"$\log (t)$")
plt.grid()
plt.legend()
plt.xticks(range(1, N + 1))
mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())
fig = plt.gcf()
fig.set_size_inches((9, 11), forward=False)
fig.savefig(f"../figures/time.armadillo.{N}.png")
if push:
    pushFile(f"../figures/time.armadillo.{N}.png")
plt.show()


ae, be = np.polyfit(n, e, 1)
plt.plot(x, ae * x + be, "--", c="k", lw=3, label="fitted line, a = %.2e" % ae)
plt.plot(n, e, "ro", ms=15, label="$\epsilon_{max}$")
plt.xlabel("$\log (n)$")
plt.ylabel("$\log (\epsilon_{max})$")
plt.grid()
plt.legend()
plt.xticks(range(1, N + 1))
mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())
fig = plt.gcf()
fig.set_size_inches((14, 11), forward=False)
fig.savefig(f"../figures/err.armadillo.{N}.png")
if push:
    pushFile(f"../figures/err.armadillo.{N}.png")
plt.show()
