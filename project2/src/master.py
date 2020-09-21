import numpy as np
import matplotlib.pyplot as plt
import sys
import time
from bunch import Bunch
import subprocess
from datareader import read_data


def help():
    msg = """
    \rThis program solves an eigenvalue problem given the following command line inputs

    \rRequired inputs:
    behaviour:   Must be first in input. Determines problem to solve. Either 'beam', 'q1' or 'q2'
    n:           Must be second in input. Size of matrices.

    \rOptional inputs. May require definiton, like 'name=some_name'.
    compile:     recompiles the c++ files. May be put anywhere in commandline input. Does not need to be defined
    name:        name of the run. Must be defined. Defaults to generic datetime name
    rho_max:     parameter of the system. Must be defined. Default to 1
    omega:       parameter of the system. Must be defined. Defaults to 0
    tolerance:   parameter of the solver. Must be defined. Defaults to 8 (meaning 10^-8)
    plot:        will plot eigenvector corrosponding to the smallest eigenvalue
                 of the solution, together with the analytical eigenvector. Need not be defined
    clear:      clear data.dat of all previous runs. Need not be defined.

    \rExample usage of the file:
    $ python3 master.py q1 20 name=q1_20_test rho_max=3.14 tolerance=5 compile
    """
    print(msg)


def plot(kwargs):
    if kwargs["behav"] == 0:
        N = kwargs["n"] + 1
        rhomax = float(kwargs["rho_max"])
        h = rhomax / N
        d = 2 / h ** 2
        a = -1 / h
        # a_eigval = d + 2 * a * np.cos(np.pi / N)
        a_eigvec = np.asarray([np.sin(i * np.pi / N) for i in range(1, N)])
        a_eigvec /= np.linalg.norm(a_eigvec)

    data = read_data("data.dat")[kwargs["name"]]
    min_eigval = np.argmin(data.eigvals)
    min_eigvec = data.eigvecs[:, min_eigval]

    plt.plot(min_eigvec, label="calculated eigenvector")
    plt.plot(a_eigvec, label="Analytical eigenvector")
    plt.legend()
    plt.show()


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
        name = "_" + str(kwargs["behav"])
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


default = Bunch(
    rho_max=1,
    omega=0,
    tolerance=8,  # tolerance to achieve before transformations halt
    name=None,
    plot=False,
    clear=False,
    # noshow=False,
    # savefigs=False,
    # pushfigs=False,
)


def main(behaviour, n, *args):
    # def main(behaviour, *args):
    problems = {"beam": 0, "q1": 1, "q2": 2}
    assert behaviour in problems
    try:
        n = eval(n)
    except:
        n = [int(n)]
    else:
        assert isinstance(n, (list, tuple))
    finally:
        kwargs = {"n": n, "behav": problems[behaviour]}

    for arg in args:
        if "=" not in arg:
            kwargs[arg] = True
        else:
            key, val = arg.split("=")
            kwargs[key] = val
    kwargs = {**default, **kwargs}
    kwargs["names"] = named_runs(kwargs)

    solve(kwargs)
    if kwargs["plot"]:
        plot(kwargs)


if __name__ == "__main__":
    args = sys.argv
    if "help" in args:
        help()
        sys.exit()
    elif "compile" in args:
        compile()
        args.remove("compile")

    if len(args) < 3:
        raise SyntaxError("Insufficent arguments. Use 'help' as argument to see usage")
    elif len(args) == 3:
        main(args[1], args[2])
    else:
        main(args[1], args[2], *args[3:])
