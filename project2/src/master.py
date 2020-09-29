import numpy as np
from argHandler import argHandler
import sys, time, subprocess

defaults = {"n":20, "rhomax":1, "eps":12, "omega":0, "method":0, "datafile":"data.dat"}

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


def main(args, **kwargs):
    if kwargs["clearfile"]:
        open(f"data/{args['datafile'][0]}","w").close()
    start = time.time()
    for method in args["method"]:
        for n in args["n"]:
            for eps in args["eps"]:
                for rhomax in args["rhomax"]:
                    for omega in args["omega"]:
                        name = f"{method}.{n}.{eps}.{rhomax}.{omega}"
                        subprocess.run(f"./main.out {name} {n} {eps} {rhomax} {int(method)} {omega} {args['datafile'][0]}".split())
                        if kwargs["debug"]:
                            print(f"n:{n}, eps:{eps}, rhomax:{rhomax}, method:{int(method)}, omega:{omega}")
    print(f"Done in {round(time.time()- start, 3)} s")


if __name__ == "__main__":
    args = sys.argv

    debug = False
    clearfile = False
    if "help" in args:
        help()
        sys.exit()
    if "compile" in args:
        print("Compiling")
        subprocess.run("g++ -o main.out main.cpp solver.cpp -Wall -larmadillo -O3".split())
        args.remove("compile")
    if "clear" in args:
        clearfile = True
        print("Clearing file")
        args.remove("clear")
    if "debug" in args:
        debug = True
        args.remove("debug")

    args = argHandler(defaults).parse(args)

    main(args, clearfile=clearfile, debug=debug)
