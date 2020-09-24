import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from datareader import read_data
from bunch import Bunch
import os
import subprocess


matplotlib.use(plt.get_backend())
font = {"family": "DejaVu Sans", "weight": "bold", "size": 22}
plt.rc("font", **font)


class Plotter:
    def __init__(self, kwargs):
        self.kwargs = Bunch(kwargs)
        self.data = read_data("data.dat")
        self.rho = lambda n: np.linspace(0, self.kwargs.rho_max, n + 1, endpoint=False)[1:]

        self.args = ("vec", "vecs", "counts", "error")

        todo = self.kwargs.plot
        if todo is True:
            self.plot_vec()
            self.plot_vecs()
            self.plot_counts()
        else:
            todo = todo.split(",")
            for do in todo:
                assert do in self.args
                eval(f"self.plot_{do}()")

    def show(self, stump="tmp", xlab="x", ylab="y", size=(16, 11)):
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.grid()
        # plt.legend()

        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        fig = plt.gcf()
        fig.set_size_inches(size, forward=False)

        name = f"../figures/{stump}.png"
        if self.kwargs.savefigs:
            fig.savefig(name)
        if self.kwargs.push:
            pushFile(name)
        if self.kwargs.noshow:
            plt.clf()
        else:
            plt.show()

    def plot_vec(self):

        colors = plt.cm.viridis(np.linspace(0, 1, len(self.kwargs.names)))
        for i, run in enumerate(self.kwargs.names):
            dat = self.data[run]
            idx = np.argsort(dat.eigvals)
            idx = idx[0]
            val = dat.eigvals[idx]
            vec = dat.eigvecs[:, idx]

            plt.plot(
                self.rho(dat.n),
                vec,# * val,
                color=colors[i],
                label=fr"n = {dat.n}, $\lambda$ = {val}",
            )
            if self.kwargs.behav == 0:
                self.plot_analytical_vec(dat, colors[i])
        self.show()

    def plot_analytical_vec(self, dat, col):
        N = dat.n + 1
        # rho_max = dat.pmax
        # hh = (N / rho_max) ** 2
        # d = 2 * hh
        # a = -hh
        # eigval = d + 2 * a * np.cos(np.pi / N)
        eigvec = np.asarray([np.sin(i * np.pi / N) for i in range(1, N)])
        eigvec /= np.linalg.norm(eigvec)

        plt.plot(
            self.rho(dat.n), eigvec, "r--", color=col, label=f"Analytical eigenvector"
        )

    def plot_vecs(self):
        dat = self.data[self.kwargs.names[-1]]
        for i in range(dat.n):
            val = dat.eigvals[i]
            vec = dat.eigvecs[:, i]
            plt.plot(
                self.rho(dat.n),
                vec,
                label=fr"Eigenvector with $\lambda$ = {round(val, 4)}",
            )
        self.show()

    def plot_counts(self):
        if len(self.kwargs.n) < 2:
            return 0
        ns = [self.data[run].n for run in self.kwargs.names]
        transforms = [self.data[run].counts for run in self.kwargs.names]
        a, b, c = tuple(np.polyfit(ns, transforms, deg=2))
        x = np.linspace(ns[0], ns[-1], len(ns))
        plt.plot(ns, transforms, "ro")
        plt.plot(x, a * x ** 2 + b * x + c, "k--", label=f"2nd degree polynomial fit")
        self.show()

    def plot_error(self):
        err = []
        n = []
        for name in self.kwargs.names:
            dat = self.data[name]
            N = dat.n + 1
            h = dat.pmax / N
            hh = h ** -2
            avals = [2 * hh * (1 - np.cos(j * np.pi / N)) for j in range(1, N)]

            val = dat.eigvals

            err.append(max(abs(val - avals) / avals))
            n.append(h)

        h = np.log10(n)
        err = np.log10(err)

        x = np.linspace(min(h), max(h), len(h))
        a, b = np.polyfit(h, err, deg=1)

        plt.plot(x, x * a + b, "k--", label=f"fit. a = {a}")
        plt.plot(h, err, "ro", label="max err")
        plt.xlabel("log10(h)")
        plt.ylabel("log10(err)")
        plt.legend()
        plt.show()


def pushFile(filename):
    path = os.path.abspath(filename)
    subprocess.run(f"git add {path}".split())
    subprocess.run(f'git commit -m "automatic_figure_update" --quiet'.split())
    subprocess.run(f"git push --quiet".split())
