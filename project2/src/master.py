import numpy as np
import matplotlib.pyplot as plt
import sys
import time
from bunch import Bunch
import subprocess
from plotter import Plotter


def help():
    msg = """
    \rThis program solves an eigenvalue problem given the following command line inputs

    \rRequired inputs:
    behaviour:   Must be first in input. Determines problem to solve. Either 'beam', 'q1' or 'q2'
    n:           Must be second in input. Size of matrices. Can be single integer,
                 list split by comma , [NO PARENTHESIS OR SPACE!], or inputs to np.arange(), split by :
                 Examples:
                 $ py master.py beam 10
                                     10,20,30,40
                                     10:101:10

    \rOptional inputs. May require definiton, like 'name=some_name'.
    compile:     recompiles the c++ files. May be put anywhere in commandline input. Does not need to be defined
    name:        name of the run. Must be defined. Defaults to generic problem-descriptive name. Not recomended
    rho_max:     parameter of the system. Must be defined. Default to 1
    omega:       parameter of the system. Must be defined. Defaults to 0
    tolerance:   parameter of the solver. Must be defined. Defaults to 8 (meaning 10^-8)
    plot:        will plot eigenvector corrosponding to the smallest eigenvalue
                 of the solution, together with the analytical eigenvector.
                 If defined, will plot the defined graphs, if not, will plot all possible
                 Plottables:  vec:   smallest eigenvector for all runs
                              vecs:  all eigenvectors for last run
                              count: plot number of transformations
                              error: error in eigenvalues, only for 'beam'
    clear:      clear data.dat of all previous runs. Need not be defined.
    push:       will push saved figs to github
    savefigs:   will save figures  !!Warning, name of saved figures not unique, will overwrite. Useless atm
    noshow:     clear figure without showing

    \rExample usage of the file:
    $ python3 master.py q1 20,40,60 name=q1_20_test rho_max=3.14 tolerance=5 compile clear plot=vec,counts

    \nKnown problems:
    Not expanded to do more problems in one run (only one of beam, q1, q2)
    Can only plot data from current run. Can easily be expanded to do all data in data.dat
    One more. Forgot.
    """
    print(msg)


def compile():
    subprocess.run("g++ -o main.out main.cpp solver.cpp -Wall -larmadillo -O3".split())


def named_runs(kwargs):
    """names runs. If name not given, assigns generic name on form
        n_problem_datehourminute
    """
    if kwargs["name"] is not None:
        names = [str(n) + kwargs["name"] for n in kwargs["n"]]
    else:
        from datetime import datetime

        today = datetime.today()
        name = "_" + rev_prob[kwargs["behav"]]
        name += "_" + today.strftime("%d%H%M")
        names = [str(n) + name for n in kwargs["n"]]

    return names


def solve(kwargs):
    if kwargs["clear"]:
        open("data.dat", "w").close()
    for i, n in enumerate(kwargs["n"]):
        subprocess.run(
            f'./main.out {kwargs["names"][i]} {n} {kwargs["tolerance"]} {kwargs["rho_max"]} {kwargs["behav"]} {kwargs["omega"]}'.split()
        )
        print(f"Method: {rev_prob[kwargs['behav']]}: n = {n}: {round(100 * (i + 1) / len(kwargs['n']))} %")


default = Bunch(
    rho_max=1,
    omega=0,
    tolerance=8,  # tolerance to achieve before transformations halt
    name=None,
    plot=False,
    clear=False,
    noshow=False,
    savefigs=False,
    push=False,
)
probs = {"beam": 0, "q1": 1, "q2": 2}
rev_prob = {v: k for k, v in probs.items()}


def main(behaviour, n, *args):
    assert behaviour in probs
    if ":" in n:
        a, b, c = n.split(":")
        n = np.arange(int(a), int(b), int(c))
    else:
        n = n.split(",")
    kwargs = {"n": n, "behav": probs[behaviour]}

    for arg in args:
        if "=" not in arg:
            kwargs[arg] = True
        else:
            key, val = arg.split("=")
            kwargs[key] = val
    kwargs = {**default, **kwargs}
    kwargs["names"] = named_runs(kwargs)

    solve(kwargs)
    # print(kwargs["plot"].split(","))
    # sys.exit()
    if kwargs["plot"] is not False:
        Plotter(kwargs)
    # error(kwargs)


if __name__ == "__main__":
    args = sys.argv
    if "help" in args:
        help()
        sys.exit()
    elif "compile" in args:
        compile()
        args.remove("compile")
        if len(args) == 1:
            sys.exit()

    if len(args) < 3:
        raise SyntaxError("Insufficent arguments. Use 'help' as argument to see usage")
    elif len(args) == 3:
        main(args[1], args[2])
    else:
        main(args[1], args[2], *args[3:])
