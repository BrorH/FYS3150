import subprocess, time, sys, inspect
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
from datareader import read_data, read_arma
from argHandler import argHandler
from colour import Color
from master import run
from scipy import optimize
mpl.use(plt.get_backend())
plotHelp = """

plotter.py, main plotting file.

each function below plots a figure.
The functions have default arguments (set in the function definition), and running the functions without any arguments will produce figures similar to those in the paper.'
In order to pass custom arguments, use same command line notation as with master.py after firstly passing the function name
i.e the function "transforms" below may be called as
 $ python3 plotter.py transforms // produces standard graph (as long as the file "data.dat" contains the right data to be plotted)
 $ python3 plotter.py transforms n=1:100:5 eps=20 // n is in range(1,100,5) and eps = 5
 $ python3 plotter.py transforms n=1:100:5 datafile=transdata.dat // same n range as above, but this time the data to be plotted is retreived from data/transdata.dat
 $ python3 plotter.py transforms eps=50 sim=True datafile=transform.dat // eps=50, and in order to plot the system it runs a simulation (sim=True) which is stored into and read from data/transform.dat

The calling of functions of custom parameters comes from a custom class in "argHandler.py" made specifically for this project.
All (some exeptions) parameters in the function definitions below can be called like this from the commandline.
Some parameters are common between (allmost) all functions:
    -n: matrix size, integration points (N-1)
    -rhomax: upper value of rho
    -eps: epsilon, tolerance
    -method: 0, 1 or 2. 0 for buckling beam problem, 1 for single-electron, and 2 for double-electron
    -datafile: Specify which file the data to be plotted is located in (must be in subdir data/)
    -sim: bool. if a dataset has not yet been generated/simulated setting sim=True will run the simulation for teh sepcified vales of n,rhomax,eps,method and omega.

generally you should not change the "method" parameter as each graphing function is custom to a method.
Any extraordinary parameters (different from the 6 mentioned above) is noted below each function.

the different functions are:
- transforms: Plots the counted transforms of a simulated system as a function of n.

- timer: Plots the time spent solving the three problems as function of n.

- optimalRhomaxSingleElectron: Plots eigval rel. error as a function of rhomax. Prints lowest relative error and corresponding rhomax

- eigvalsSingleElectron: Plots the numerical eigenvalues of the single electron system against either n or rhomax.
    special parameters:
        -n, rhomax: only ONE of these must be a list of values to be plotted against.
        -eigcount: how many eigvals are to ble plotted
- eigvalAccuracySingleElectron: Plots eigval accuracy against n and attemps a function fit of shape a*log10(n)+b
    special parameters:
        - fitIdxStart: The index of n where the function fit begins
"""

"""
MAIN PROGRAM FOR PLOTTING AND CREATING THE FIGURES IN THE PAPER

# COLOR must be installed before running
# must make a data/ folder

Consists of several functions that each create a single plot.
only one plot can be made each call of this program.
parameters:
    -plottype, must br first. is just the name of the function (see list below)
options:
    -params must follow argHandler formating (see argHandler.py)
"""

font = {"family": "DejaVu Sans", "weight": "normal", "size": 18}
plt.rc("font", **font)

def transforms(n=range(1,10), rhomax = [1,10,20], eps = 12, omega = [0,1,5], method = [0,1,2], sim=False, datafile="data.dat"):
    # plots num of transformations against N on axis ax, as well as comparing to n**2 and n**2*ln(n)
    if sim:
        print("Solving ...")

        open(f"data/{datafile}", "w+").close() # clear file
        for method_ in method:
            start = time.time()
            for N in n:
                subprocess.run(f"./main.out {N}{eps}{rhomax[method_]}{method_}{omega[method_]} {N} {eps} {rhomax[method_]} {method_} {omega[method_]} {datafile}".split())
                print(f"{round(100*(N-n[0])/(n[-1]-n[0]), 2)} %, N: {N}, eps: {eps}, method: {method_}, rhomax: {rhomax[method_]}, omega: {omega[method_]}")
            print(f"Done in {round(time.time() -start,3)} s.")
        #if datafile != "data.dat":
            # move content to designated file
            #subprocess.run(f"cp data.dat {datafile}".split())

    sols = read_data(datafile) #all solutions
    for method_ in method: # all methods are plotted
        trans = np.array([sols[f"{N}{eps}{rhomax[method_]}{method_}{omega[method_]}"].transformations for N in n]) # array of counted transformations

        color = {0:"red", 1:"blue", 2:"green"}[method_] #every method gets unique color
        label = {0:"Buckling beam", 1:r"Quantum 1, $\rho_{max} = $"+str(rhomax[method_]), 2:r"Quantum 2, $\rho_{max} = $"+f"{rhomax[method_]}, $\omega_r = {omega[method_]}$"}[method_]

        ax.plot(n, trans, "-o", markersize = 2.5, color=color, alpha=0.8, lw=2.8, label=label)

        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tick_params(top='on', bottom='on', left='on', right='on', labelleft='on', labelbottom='on')
        legend = ax.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1, frameon = True, fontsize="x-small")
        frame = legend.get_frame()
        frame.set_facecolor('white')
        ax.set_xlabel("Matrix size, $n$")
        ax.set_ylabel("Transformations, $m$")
        plt.title(r"Transformations; $m$" )

        print(f" eps: {eps}, method: {method_}, rhomax: {rhomax[method_]}, omega: {omega[method_]}:")


def timer(n=range(1, 10), rhomax=[1, 10, 20], eps=12, omega=[0, 1, [0.01, 5]], method=[0, 1, 2], sim=False, datafile="data.dat"):
    # plots elapsed time of solving against N on axis ax, as well as comparing to n**2 and n**2*ln(n)
    if sim:
        print("Solving ...")

        open(f"data/{datafile}", "w+").close() # clear file
        for method_ in method:
            start = time.time()
            if not isinstance(omega[method_], (list, np.ndarray)):
                w = [omega[method_]]
            else:
                w = omega[method_]
            for wi in w:
                for N in n:
                    subprocess.run(f"./main.out {N}{eps}{rhomax[method_]}{method_}{wi} {N} {eps} {rhomax[method_]} {method_} {wi} {datafile}".split())
                    print(f"{round(100*(N-n[0])/(n[-1]-n[0]), 2)} %, N: {N}, eps: {eps}, method: {method_}, rhomax: {rhomax[method_]}, omega: {wi}")
            print(f"Finished in {round(time.time() -start,3)} s.")
        #if datafile != "data.dat":
            # move content to designated file
            #subprocess.run(f"cp data.dat {datafile}".split())

    sols = read_data(datafile) #all solutions
    for method_ in method: # all methods are plotted
        if "__iter__" not in dir(omega[method_]):
            w = [omega[method_]]
        else:
            w = omega[method_]
        for i, wi in enumerate(w): # all omegas are plotted
            tim = np.array([sols[f"{N}{eps}{rhomax[method_]}{method_}{wi}"].time for N in n]) # array of counted transformations

            color = {0:["red"], 1:["blue"], 2:["lime", "mediumseagreen", "seagreen", "darkgreen"]}[method_][i] #every method gets unique color
            label = {0:"Buckling beam", 1:r"Quantum 1, $\rho_{max} = $"+str(rhomax[method_]), 2:r"Quantum 2, $\rho_{max} = $"+f"{rhomax[method_]}, $\omega_r = {wi}$"}[method_]

            ax.plot(n, tim, "-o", markersize=12, color=color, alpha=0.8, lw=5, label=label)

            plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
            plt.tick_params(top='on', bottom='on', left='on', right='on', labelleft='on', labelbottom='on')
            legend = ax.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1, frameon=True, fontsize="x-small")
            frame = legend.get_frame()
            frame.set_facecolor('white')
            ax.set_xlabel("Matrix size, $n$")
            ax.set_ylabel("Time, $s$")
            plt.title(r"Time of solving")#, varying $\omega_r$")

        print(f" eps: {eps}, method: {method_}, rhomax: {rhomax[method_]}, omega: {wi}:")


def time_armadillo(n=np.arange(10, 151, 10), rho_max=1, eps=12, omega=0, method=0, datafile="data.dat"):

    arma = []
    jacobi = []
    print("Solving ...")
    start = time.time()
    for _n in n:
        open(f"data/{datafile}", "w").close() # clear file
        subprocess.run(f"./arma.out {n} {datafile}".split())
        _, _, t = read_arma(datafile)
        arma.append(t)

        open(f"data/{datafile}", "w").close() # clear file1
        subprocess.run(f"./main.out {n}{eps}{rho_max}{method}{omega} {n} {eps} {rho_max} {method} {omega} {datafile}".split())
        t = read_data(datafile).time
        jacobi.append(t)

        print(f"ρ: {round(_n,3)}/{n[-1]}, {round(100*(_n-n[0])/(n[-1]-n[0]), 2)} %")
    print(f"Done in {round(time.time()- start,3)} s")

    print(arma)
    print(jacobi)

def optimalRhomaxSingleElectron(n= 100,rhomax = np.linspace(3.2,5.7, 1000), eps = 12, omega=1, method=1, sim = False, datafile = "data.dat" ):

    # given a dataset with rhomax (rhomax) finds which rhomax makes each of the first 4 eigvals of quantum 1 the most accurate
    # produces plot as well

    if sim:
        print("Solving ...")
        open(f"data/{datafile}", "w").close() # clear file
        start = time.time()
        for _rhomax in rhomax:
            subprocess.run(f"./main.out {n}{_rhomax} {n} {eps} {_rhomax} {method} {omega} {datafile}".split())
            print(f"ρ: {round(_rhomax,3)}/{rhomax[-1]}, {round(100*(_rhomax-rhomax[0])/(rhomax[-1]-rhomax[0]), 2)} %")
        print(f"Done in {round(time.time()- start,3)} s")


    sols = read_data(datafile) # read solutions


    colors = list(Color("red").range_to(Color("blue"),4)) # prettyyyy colorrssss
    colors = [str(color) for color in colors]


    analytical = np.array([4*a+3 for a in range(0,4)]) # analytical eigenvalues

    erel = np.zeros((len(rhomax),4)) #error_rel, to be filled
    for i,_rhomax in enumerate(rhomax):
        erel[i] = np.abs(np.sort(sols[f"{n}{_rhomax}"].eigvals)[:4]-analytical)/analytical # fill with all relative errors

    minidx = erel.argmin(axis=0) # find where the graph has its lowest point, a.k.a where the error is the lowest

    # plotting code
    for i in range(4):
        ax.plot(rhomax, erel[:,i], color = colors[i])
    for i,idx in enumerate(minidx):
        print(f"eig_{i+1} Best rhomax: {round(rhomax[idx],5)}, err_rel: {round(erel[idx][i],6)}") # print the results
        ax.plot(rhomax[idx], erel[idx][i], "-o", label = f"$\lambda_{i+1}: $"+r"$\rho_{max}="+f"{round(rhomax[idx],3)}$,"+r"$\epsilon_{rel}="+"%.2e$" %erel[idx][i],color = colors[i])
    legend = ax.legend(fancybox=True, framealpha=1, shadow=True, borderpad=0.4, frameon = True)
    legend.get_frame().set_facecolor('white')
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax.set_title(r"$\epsilon_{rel}$ for four first eigenvalues")
    ax.set_xlabel(r"$\rho_{max}$")
    ax.set_ylabel(r"$\epsilon_{rel}$")


def eigvalsSingleElectron(n = 100, rhomax = np.linspace(4,6, 100), eps = 16, omega = 1, method=1, eigcount = 10, sim = False, datafile = "data.dat"):
    if not isinstance(n, list) and not isinstance(n, np.ndarray):
        var = "rhomax" #rhomax varies
        vars = rhomax
        static = f"$n = {n}$"
    else:
        var = "n" #else n varies
        vars = n
        static = r"$\rho_{max}=$" + f"${rhomax}$"
    if sim: # run simulation if prompted
        print("Solving ...")
        open(f"data/{datafile}", "w").close() # clear file
        start = time.time()
        if var == "n":
            for n in vars:
                subprocess.run(f"./main.out {n}{rhomax} {n} {eps} {rhomax} {method} {omega} {datafile}".split())
                print(f"N: {n}/{vars[-1]}, {round(100*(n-vars[0])/(vars[-1]-vars[0]), 2)} %")
            print(f"Done in {round(time.time()- start,3)} s")
        elif var == "rhomax":
            for rhomax in vars:
                subprocess.run(f"./main.out {n}{rhomax} {n} {eps} {rhomax} 1 1 {datafile}".split())
                print(f"ρ: {round(rhomax,3)}/{vars[-1]}, {round(100*(rhomax-vars[0])/(vars[-1]-vars[0]), 2)} %")
            print(f"Done in {round(time.time()- start,3)} s")

    sols = read_data(datafile)


    colors = list(Color("red").range_to(Color("cyan"),len(vars))) #TO get a gradient effect
    colors = [str(color) for color in colors]


    analytical = np.array([4*a+3 for a in range(0,eigcount)]) # the analytical eigenvalues
    indexes = np.arange(1,eigcount+1)

    if var == "n":
        for i,n in enumerate(vars):
            ax.plot(indexes, np.sort(sols[f"{n}{rhomax}"].eigvals)[:eigcount],'-o', markersize=1, color=colors[i], lw=1.5, alpha=0.6)
        label = f"Matrix size, $n$"
    elif var == "rhomax":
        for i,rhomax in enumerate(vars):
            ax.plot(indexes, np.sort(sols[f"{n}{rhomax}"].eigvals)[:eigcount],'-o', markersize=1, color=colors[i], lw=1.5, alpha=0.6)
        label = r"$\rho_{max}$"

    ax.plot(indexes, analytical, label="Analytical", color="black", lw=2, alpha=0.8, linestyle="--")

    ax.set_xlabel("Index, $i$")
    ax.set_ylabel(r"$\lambda_{i}$")
    ax.set_ylim(top=analytical[-1]*1.2, bottom=analytical[0]*0.8)

    ax.set_xlim(1, eigcount)
    cmap = mpl.colors.ListedColormap(colors)
    norm = mpl.colors.Normalize(vmin=vars[0],vmax=vars[-1])
    ax.figure.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='horizontal', label=label, ticks=[int(a) for a in np.linspace(vars[0], vars[-1],5)] )
    legend = ax.legend(fancybox=True, framealpha=1, shadow=True, borderpad=0.4, frameon = True,loc='lower right')
    legend.get_frame().set_facecolor('white')
    ax.set_title(f"Eigenvalues for {static}")


def eigvalAccuracySingleElectron(n = list(range(10, 350)), rhomax=20, eps = 12, method=1, fitIdxStart = 10, sim = False, datafile="data.dat"):
    if sim:
        open(f"data/{datafile}", "w+").close()
        for n_ in n:
            subprocess.run(f"./main.out {n_} {n_} {eps} {rhomax} 1 1 {datafile}".split())
            print(f"N: {n_}, {round(100*(n_-n[0])/(n[-1]-n[0]), 2)} %")


    data = read_data(datafile)
    avals = np.array([3,7,11,15]) # analytic eigenvals


    vals = np.zeros(len(n)) # array to be filled
    for i,n_ in enumerate(n):
        vals[i] = -np.log10(1/4*np.sum(np.abs(np.sort(data[str(n_)].eigvals)[:4]-avals)))

    ax.plot(n, vals, color="red", lw=3, alpha=1, label="Numerical values")

    # using scipys optimize.curve_fit we find the best coefficient for the model a*log10(n)+b
    f = lambda t,a,b: a*np.log10(t)+b
    popt, pcov = optimize.curve_fit(f, n[fitIdxStart:], vals[fitIdxStart:])

    x = np.linspace(n[fitIdxStart], 4000, 1000)
    ax.plot(x,f(x,*popt), color="blue", lw=2, alpha=0.5, label=f"Fit ${round(popt[0],2)}\log"+r"_{10}"+f"(n){round(popt[1],3)}$")


    ax.set_xlabel("Matrix size, n")
    ax.set_ylabel("No. of correct digits")
    legend = ax.legend(fancybox=True, framealpha=1, shadow=True, borderpad=0.4, frameon = True,loc='lower right')
    legend.get_frame().set_facecolor('white')
    ax.set_title(f"Number of correct digits in eigenvalues")





if __name__ == "__main__":
    with plt.style.context("seaborn-darkgrid"):
        f, ax = plt.subplots(1,1, dpi=100,frameon=True)
        defaults = {"datafile":"data.dat"}

        if "help" in sys.argv:
            print(plotHelp)
            sys.exit()

        try:
            plotType = sys.argv[1]
            assert plotType in dir(), f"{plotType} not found as function. Call program with 'help' for help"

            flocals = inspect.getfullargspec(eval(plotType))
            for var,val in zip(flocals.args, flocals.defaults):
                if isinstance(val, np.ndarray):
                    val = [i for i in val ]
                defaults[var] = val

            args = argHandler(defaults).parse(sys.argv[1:])

            args["datafile"] = "'"+str(args["datafile"])+"'"
            call_args = ",".join([f'{arg}={args[arg]}' for arg in args.keys()])
        except Exception as e:
            print(e)
            print("plotter.py failed in initalization stage. passs arg 'help' for help and examples")
            sys.exit()
        exec(f"{plotType}({call_args})")


plt.show()





