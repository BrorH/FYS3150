import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import subprocess
from datareader import read_data
import sys


matplotlib.use(plt.get_backend())
font = {"family": "DejaVu Sans", "weight": "normal", "size": 18}
plt.rc("font", **font)

arma_name = "comparison.out"  # compiled armadillo c++ file
arma_file = "comparison.dat"  # datafile for armadillo
Jacobi_file = "jacobi.dat"  # datafile for jacobi method


def mat_sort_by_array(mat, vec, zero=1e-12):
    idx = np.argsort(vec)

    vec2 = np.sort(vec)
    mat2 = mat[:, idx[0]]
    if mat2[0] < 0:
        mat2 *= -1
    for i in idx[1:]:
        tmp = mat[:, i]
        if tmp[0] < 0:
            tmp *= -1
        mat2 = np.vstack((mat2, tmp))
    mat2 = np.where(abs(mat2) < zero, 0, mat2)
    return mat2.T, vec2


def read_arma():
    vals = []
    vecs = []
    with open("data/" + arma_file, "r") as file:
        lines = file.readlines()
        n, time = lines[0].split(",")
        for i in range(len(lines) - 1):
            V = [float(v) for v in lines[i + 1].split(",")]
            vals.append(V[0])
            vecs.append(V[1:])

    return np.asarray(vals), np.asarray(vecs), float(time)


def compile():
    subprocess.run(f"g++ -o {arma_name} arma_beam.cpp -larmadillo".split())
    subprocess.run(f"g++ -o main.out -O3 main.cpp solver.cpp -larmadillo".split())


def analytical(n):
    N = n + 1
    L = [2 * N ** 2 * (1 - np.cos(j * np.pi / N)) for j in range(1, N)]
    V = np.zeros((n, n))
    for j in range(1, N):
        vec = np.asarray([np.sin(i * j * np.pi / N) for i in range(1, N)])
        V[:, j - 1] = vec / np.linalg.norm(vec)
    return mat_sort_by_array(V, L)


def error(Va, Vn):
    """
    Computes the total relative difference between two solutions, using equation (20)

    Va, analytical solution
    Vn, numerical solution
    """
    error = abs(Vn - Va) / Va
    error = np.where(error == np.inf, 0, error)
    return np.sum(error)


def main(N, tol):
    arma_Jac_err = []
    analy_Jac_err = []
    analy_arma_err = []
    time_arma = []
    time_Jac = []
    for n in N:
        print("n = ", n)
        # Find solution using armadillo
        open("data/" + arma_file, "w").close()
        subprocess.run(f"./{arma_name} {n} {arma_file}".split())
        # sort eigenvectormatrix according to increasing eigenvalue
        arma_vals, arma_vecs, time = read_arma()
        arma_vecs, arma_vals = mat_sort_by_array(arma_vecs.T, arma_vals)
        time_arma.append(time)

        # Find numerical solution
        open("data/" + Jacobi_file, "w").close()
        subprocess.run(f"./main.out tmp {n} {tol} 1 0 0 {Jacobi_file}".split())
        # sort eigenvectormatrix according to increasing eigenvalue
        Jacobi = read_data(Jacobi_file)["tmp"]
        Jac_vecs, Jac_vals = mat_sort_by_array(Jacobi.eigvecs, Jacobi.eigvals)
        time_Jac.append(Jacobi.time)

        analy_vec, analy_val = analytical(n)

        # Compute total relative error in each eigenvector
        arma_Jac_err.append(error(arma_vecs, Jac_vecs))
        analy_Jac_err.append(error(analy_vec, Jac_vecs))
        analy_arma_err.append(error(analy_vec, arma_vecs))

    with plt.style.context("seaborn-darkgrid"):
        f, ax = plt.subplots(1,1, dpi=100,frameon=True)
        ax.plot(
            N,
            arma_Jac_err,
            "o",
            color="firebrick",
            label=r"$\delta_{rel}$ armadillo-Jacobi",
        )
        ax.plot(
            N,
            analy_Jac_err,
            "o",
            color="forestgreen",
            label=r"$\delta_{rel}$ analytical-Jacobi",
        )
        ax.plot(
            N,
            analy_arma_err,
            "o",
            color="mediumblue",
            label=r"$\delta_{rel}$ analytical-armadillo",
        )

        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        fig = plt.gcf()
        fig.set_size_inches((16, 11), forward=False)

        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tick_params(top='on', bottom='on', left='on', right='on', labelleft='on', labelbottom='on')
        legend = ax.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1, frameon=True, fontsize="x-small")
        frame = legend.get_frame()
        frame.set_facecolor('white')

        plt.title("Total relative difference between each method")
        plt.xlabel("N")
        plt.ylabel("relative difference")
        plt.legend()
        plt.show()

    with plt.style.context("seaborn-darkgrid"):
        f, ax = plt.subplots(1, 1, dpi=100, frameon=True)
        ax.plot(N, time_arma, "o", color="firebrick", label="armadillo")
        ax.plot(N, time_Jac, "o", color="forestgreen", label="Jacobi method")

        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tick_params(top='on', bottom='on', left='on', right='on', labelleft='on', labelbottom='on')
        legend = ax.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1, frameon=True, fontsize="x-small")
        frame = legend.get_frame()
        frame.set_facecolor('white')
        plt.title("Time of solving")
        plt.xlabel("N")
        plt.ylabel("Time, [s]")
        plt.legend()
        plt.show()


if __name__ == "__main__":
    if "compile" in sys.argv:
        compile()
    n = np.arange(10, 201, 5)
    tolerance = 12
    main(n, tolerance)
