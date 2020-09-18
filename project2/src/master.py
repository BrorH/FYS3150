# import numpy as np
import matplotlib.pyplot as plt
import subprocess
# import pandas as pd
import sys
import os
import matplotlib
from bunch import Bunch
from datareader import read_data

matplotlib.use(plt.get_backend())


def compile(compiled_name):
    subprocess.run(f"g++ -o {compiled_name} main.cpp solver.cpp -larmadillo".split())


def pushFile(filename):
    path = os.path.abspath(filename)
    subprocess.run(f"git add {path}".split())
    subprocess.run(f'git commit -m "automatic_figure_update" --quiet'.split())
    subprocess.run(f"git push --quiet".split())


def plot(xlabel="x", ylabel="y", savename="plot.pdf", size=(16, 11)):
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid()
    plt.legend()
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    fig = plt.gcf()
    fig.set_size_inches((14, 11), forward=False)
    if specs["save"]:
        fig.savefig(f"../figures/{savename}")
    if specs["push"]:
        pushFile(f"../figures/{savename}")
    if specs["noshow"]:
        plt.clf()
    else:
        plt.show()


def solve_problem(**kwargs):
    if specs.compile:
        compile(specs.compiled_name)
    print(f"solving. Problem: {specs.behaviour}, n: {specs.n}")
    open(specs.datafile, "w").close()
    subprocess.run(
        f"./{specs.compiled_name} {specs.n} {specs.tolerance} {specs.rho_max} {specs.behaviour} {specs.omega}".split()
    )
    Solved = read_data(specs.datafile)
    

# Settings
specs = Bunch(
    noshow=False,
    push=False,
    save=False,
    compile=False,
    compiled_name="main.out",
    datafile="data.dat",
    n=0,
    tolerance=8,
    behaviour="dummy",
    rho_max=0,
    omega=0,
)


def main():
    args = sys.argv
    # assert len(sys.argv) > 3
    specs.behaviour = args[1]
    assert specs.behaviour in [
        "beam",
        "q1",
        "q2",
    ], f"Type of problem must be specified. \n 'beam', 'q1', 'q2'"
    try:
        specs.n = int(args[2])
        specs.rho_max = float(args[3])
        if specs.behaviour == "q2":
            specs.omega = args[4]
    except:
        print(
            "wrong usage of file. Needs arguments n (int), rhomax (float) and omega* (float)\n*only for behaviour = 'q2'."
        )
        sys.exit(1)

    for arg in args[5:]:
        try:
            specs[arg] = True
        except KeyError:
            continue

    solve_problem()


if __name__ == "__main__":
    main()
