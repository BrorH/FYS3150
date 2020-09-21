import numpy as np
import matplotlib.pyplot as plt
from datareader import read_data
from bunch import Bunch


class Plotter:
    def __init__(self, kwargs):
        self.kwargs = Bunch(kwargs)
        self.data = read_data("data.dat")
        self.rho = lambda n: np.linspace(0, self.kwargs.rho_max, n + 1, endpoint=False)[1:]

        self.args = ("vec", "vecs", "counts")

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

    def show(self):
        plt.legend()
        plt.show()

    def plot_vec(self):

        colors = plt.cm.viridis(np.linspace(0, 1, len(self.kwargs.names)))
        for i, run in enumerate(self.kwargs.names):
            dat = self.data[run]
            idx = np.argmin(dat.eigvals)
            val = dat.eigvals[idx]
            vec = dat.eigvecs[:, idx]

            plt.plot(
                self.rho(dat.n),
                vec,
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
