import numpy as np
import matplotlib.pyplot as plt
import sys
import time
from bunch import Bunch
import subprocess
from plotter import Plotter

defaults = {"n":[20], "rhomax":[1], "eps":[12], "omega":[0], "method":[0], "datafile":"data.dat"}

def help():
    msg = """
    \rThis program solves an eigenvalue problem given the following command line inputs
    5 parameters for the system.
        -n: matrix size. DEFAUlT: 50
        -eps: epsilon, tolerance of system in negative log10 (i.e 12 -> 1e-12). DEFAULT: 12
        -rhomax: upper value of rho. DEFAULT: 1
        -method: int, 0,1 or 2. 0= buckling, 1 =q1, 2=q1. DEFAULT: 0
        -omega: omega for q2 system. set to 1 for q1. DEFUALT: 0
        -datafile: filename in data/ folder to store the datafile. DEFAULT: data.dat
    ALL of the above arguments (except datafile) can be passed in groups as follows:
        -single values (i.e n=2),
        -as lists (i.e n=2,3,8,100), (NO WS BETWEEN COMMA SEPARATOR)
        -as ranged lists (i.e n=10:100:2) (NO WS BETWEEN COLON SEPARATOR)

    some options (passed as standalone arguments):
        -clear: clears the passed datafile before filling it (even if it is the defaulted one)
        -compile: compiles the program
        -debug: prints debug messages
        -help: prints this message
    """
    print(msg)


def compile():
    subprocess.run("g++ -o main.out main.cpp solver.cpp -Wall -larmadillo -O3".split())


def splitrange(string):
    # turns a string-formated range, i.e 10:101:10 into a python list list(range(10,101, 10))
    # can be passed as 10:101 (interval as 1) or 10:101:23 (interval as 23)
    splitted = string.split(":")
    if len(splitted) == 2:
        return splitrange(string+":1")

    return list(range(int(splitted[0]), int(splitted[1]), int(splitted[2])))


def splitlist(string):
    # splits a string of comma separated values into python list
    return [float(obj) for obj in string.split(",")]

def main(*args, **kwargs):
    #assert behaviour in probs
    final_args = {}
    for arg in args[1:]:
        for argname in defaults.keys():
            eqIdx = arg.index("=")
            if arg[:eqIdx] == argname:
                if "," in arg[eqIdx+1:]:
                    final_args[argname] = splitlist(arg[eqIdx+1:])
                elif ":" in arg[eqIdx+1:]:
                    final_args[argname] = splitrange(arg[eqIdx+1:])
                else:
                    final_args[argname] = [arg[eqIdx+1:]]
    for undefined in list(set(defaults.keys())- set(final_args.keys())):
        print(f"{undefined} defaulted to {defaults[undefined][0]}")
        final_args[undefined] = defaults[undefined]

    if kwargs["clearfile"]:
        open(f"data/{final_args['datafile']}","w").close()
    start = time.time()
    for method in final_args["method"]:
        for n in final_args["n"]:
            for eps in final_args["eps"]:
                for rhomax in final_args["rhomax"]:
                    for omega in final_args["omega"]:
                        name = f"{method}.{n}.{eps}.{rhomax}.{omega}"
                        #print(f"./main.out {name} {n} {eps} {rhomax} {int(method)} {omega} {final_args['datafile']}")
                        subprocess.run(f"./main.out {name} {n} {eps} {rhomax} {int(method)} {omega} {final_args['datafile'][0]}".split())
                        if kwargs["debug"]:
                            print(f"n:{n}, eps:{eps}, rhomax:{rhomax}, method:{int(method)}, omega:{omega}")
    print(f"Done in {round(time.time()- start, 3)} s")


if __name__ == "__main__":
    debug = False
    clearfile = False

    args = sys.argv
    if "help" in args:
        help()
        sys.exit()
    if "compile" in args:
        compile()
        args.remove("compile")
    if "clear" in args:
        clearfile = True
        args.remove("clear")
    if "debug" in args:
        debug = True
        args.remove("debug")

    main(*args, clearfile=clearfile, debug=debug)
