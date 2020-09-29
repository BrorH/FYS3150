import subprocess, time, sys, inspect
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import cm
import matplotlib as mpl
from datareader import read_data
from argHandler import argHandler
from colour import Color
from master import run
from scipy import optimize

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


def plotTransforms(n=np.arange(1.10), sim=False, datafile="data.dat", eps = 12, rhomax = [1,10,20], omega=[0,1,5],**kwargs):
    # plots num of transformations against N on axis ax, as well as comparing to n**2 and n**2*ln(n)
    if sim:
        print("Solving ...")
        
        open(f"data/{datafile}", "w+").close() # clear file
        for method in [0,1,2]:
            start = time.time()
            for N in n:
                subprocess.run(f"./main.out {N}{eps}{rhomax[method]}{method}{omega[method]} {N} {eps} {rhomax[method]} {method} {omega[method]} {datafile}".split())
                print(f"{round(100*(N-n[0])/(n[-1]-n[0]), 2)} %, N: {N}, eps: {eps}, method: {method}, rhomax: {rhomax[method]}, omega: {omega[method]}")
            print(f"Done in {round(time.time() -start,3)} s.")
        #if datafile != "data.dat":
            # move content to designated file
            #subprocess.run(f"cp data.dat {datafile}".split())
    
    sols = read_data(datafile) #all solutions
    for method in [0,1,2]: # all methods are plotted
        trans = np.array([sols[f"{N}{eps}{rhomax[method]}{method}{omega[method]}"].transformations for N in n]) # array of counted transformations
        
        color = {0:"red", 1:"blue", 2:"green"}[method] #every method gets unique color
        label = {0:"Buckling beam", 1:r"Quantum 1, $\rho_{max} = $"+str(rhomax[method]), 2:r"Quantum 2, $\rho_{max} = $"+f"{rhomax[method]}, $\omega_r = {omega[method]}$"}[method]

        ax.plot(n, trans, "-o", markersize = 2.5, color=color, alpha=0.8, lw=2.8, label=label)
        
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tick_params(top='on', bottom='on', left='on', right='on', labelleft='on', labelbottom='on')
        legend = ax.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1, frameon = True, fontsize="x-small")
        frame = legend.get_frame()
        frame.set_facecolor('white')
        ax.set_xlabel("Matrix size, $n$")
        ax.set_ylabel("Transformations, $m$")
        plt.title(r"Transformations; $m$" )

        print(f" eps: {eps}, method: {method}, rhomax: {rhomax[method]}, omega: {omega[method]}:")


def plotOptimalRhomaxQuantum1(rhomax = np.linspace(3.2,5.7, 1000), omega=12, method=1, n = 100, sim = False, datafile = "data.dat", eps = 12 ):

    # given a dataset with rhomax (rhomax) finds which rhomax makes each of the first 4 eigvals of quantum 1 the most accurate
    # produces plot as well

    if sim:
        print("Solving ...")
        open(f"data/{datafile}", "w").close() # clear file
        start = time.time()
        for _rhomax in rhomax:
            subprocess.run(f"./main.out {n}{_rhomax} {n} {eps} {_rhomax} 1 1 {datafile}".split())
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


def plotEigvalsQuantum1(n = 100, rhomax=np.linspace(4,6, 100), eigcount = 10, sim = False, datafile = "data.dat", eps = 16, **kwargs):
    """
    In the paper we wish to plot the eigenvalues both as a function of varying n and varying rhomax. 
    The variable that varies is passed under "var" and the actually values of the variable under "vars". var can be "n" or "rhomax"
    the static variable (opposite of the varying) must be passed as a kwarg, i.e 'var = "n", vars = range(1,10), rhomax=20'
    eigcount is how many eigenvalues do you want to plot.
    ax is the axis to be plotted on
    """

    # one of rhomax and n has to be 1-dimensional, as it is the constant
    if not isinstance(n, list):
        var = "rhomax" #rhomax varies
        vars = rhomax
    else:
        var = "n" #else n varies
        vars = n
    if sim: # run simulation if prompted
        print("Solving ...")
        open(f"data/{datafile}", "w").close() # clear file
        start = time.time()
        if var == "n":
            for n in vars:
                subprocess.run(f"./main.out {n}{rhomax} {n} {eps} {rhomax} 1 1 {datafile}".split())
                print(f"N: {n}/{vars[-1]}, {round(100*(n-vars[0])/(vars[-1]-vars[0]), 2)} %")
            print(f"Done in {round(time.time()- start,3)} s")
            #subprocess.run(f"cp data.dat {datafile}".split())
        elif var == "rhomax":
            for rhomax in vars:
                subprocess.run(f"./main.out {n}{rhomax} {n} {eps} {rhomax} 1 1 {datafile}".split())
                print(f"ρ: {round(rhomax,3)}/{vars[-1]}, {round(100*(rhomax-vars[0])/(vars[-1]-vars[0]), 2)} %")
            print(f"Done in {round(time.time()- start,3)} s")
            #subprocess.run(f"cp data.dat {datafile}".split())
    
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
            ax.plot(indexes, np.sort(sols[f"{n}{rhomax}"].eigvals)[:eigcount],'-o', markersize=1, color=colors[i], lw=1.5, alpha=0.6, label="aaa")
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


def plotEigvalAccuracyDynamicN(n = list(np.arange(10, 350)), eps = 12, rhomax=20, datafile="data.dat", sim=False, fitIdxStart = 10):
    # WHAT IS fitIdxStart
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


functionDesc = "To be filled"


if __name__ == "__main__":
    with plt.style.context("seaborn-darkgrid"):
        f, ax = plt.subplots(1,1, dpi=100,frameon=True)
        defaults = {"datafile":"data.dat"}#"sim":False, "eps":12, "n":50, "datafile":"data.dat", "rhomax":12, "omega":0, "method":0, "vars":[], "eigcount":10}
        plotType = sys.argv[1]
        assert plotType in dir(), f"{plotType} not found as function. Available types are:\n {functionDesc}"
        
        flocals = inspect.getfullargspec(eval(plotType))
        for var,val in zip(flocals.args, flocals.defaults):
            if isinstance(val, np.ndarray):
                val = [i for i in val ]
            defaults[var] = val
        
        _args = argHandler(defaults)
    
        args = _args.parse(sys.argv[1:])
        
        args["datafile"] = "'"+str(args["datafile"])+"'"
        call_args = ",".join([f'{arg}={args[arg]}' for arg in args.keys()])

        exec(f"{plotType}({call_args})")


plt.show()


        


