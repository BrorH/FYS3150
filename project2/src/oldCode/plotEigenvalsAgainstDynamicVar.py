import numpy as np
import matplotlib.pyplot as plt 
import subprocess
from datareader import read_data
import time
import matplotlib as mpl
from colour import Color

"""
Plots the numerical eigenvalues (as well as the analytical ones) for quantum 1.
Has colorbars that show the evolution of the eigenvalues for different valeus of n and rhomax.
This specific usage produces 2 plots; one with varying n and one with varying rhomax

"""

def plotEigvalsQuantum1( ax, var = "n", vars = [], eigcount = 10, sim = False, datafile = "data.dat", eps = 16, **kwargs):
    """
    In the paper we wish to plot the eigenvalues both as a function of varying n and varying rhomax. 
    The variable that varies is passed under "var" and the actually values of the variable under "vars". var can be "n" or "rhomax"
    the static variable (opposite of the varying) must be passed as a kwarg, i.e 'var = "n", vars = range(1,10), rhomax=20'
    eigcount is how many eigenvalues do you want to plot.
    ax is the axis to be plotted on
    """

    if sim: # run simulation if prompted
        print("Solving ...")
        open("data.dat", "w").close() # clear file
        start = time.time()
        if var == "n":
            for n in vars:
                subprocess.run(f"./main.out {n}{kwargs['rhomax']} {n} {eps} {kwargs['rhomax']} 1 1".split())
                print(f"N: {n}/{vars[-1]}, {round(100*(n-vars[0])/(vars[-1]-vars[0]), 2)} %")
            print(f"Done in {round(time.time()- start,3)} s")
            subprocess.run(f"cp data.dat {datafile}".split())
        elif var == "rhomax":
            for rhomax in vars:
                subprocess.run(f"./main.out {kwargs['n']}{rhomax} {kwargs['n']} {eps} {rhomax} 1 1".split())
                print(f"ρ: {rhomax}/{vars[-1]}, {round(100*(rhomax-vars[0])/(vars[-1]-vars[0]), 2)} %")
            print(f"Done in {round(time.time()- start,3)} s")
            subprocess.run(f"cp data.dat {datafile}".split())
    
    sols = read_data(datafile)


    colors = list(Color("red").range_to(Color("cyan"),len(vars))) #TO get a gradient effect
    colors = [str(color) for color in colors]
    

    analytical = np.array([4*a+3 for a in range(0,eigcount)]) # the analytical eigenvalues
    indexes = np.arange(1,eigcount+1) 

    if var == "n":
        for i,n in enumerate(vars):
            ax.plot(indexes, np.sort(sols[f"{n}{kwargs['rhomax']}"].eigvals)[:eigcount],'-o', markersize=1, color=colors[i], lw=1.5, alpha=0.6)
        label = f"Matrix size, $n$"
    elif var == "rhomax":
        for i,rhomax in enumerate(vars):
            ax.plot(indexes, np.sort(sols[f"{kwargs['n']}{rhomax}"].eigvals)[:eigcount],'-o', markersize=1, color=colors[i], lw=1.5, alpha=0.6)
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




if __name__ == "__main__":
    with plt.style.context("seaborn-darkgrid"):
        # in thsi instance rhomax is static, while n is dynamic
        rhomax = 20
        f, ax = plt.subplots(1, 1, dpi=175, frameon=True)
        Nvars = np.array(range(10, 253, 3))
        plotEigvalsQuantum1(ax, var="n", vars=Nvars, sim = False, eigcount =10, datafile="dataNvar.dat", rhomax=rhomax)
        ax.set_title(r"$\lambda_i$ for $\rho_{max} = "+f"{rhomax}$")
        plt.show()


    with plt.style.context("seaborn-darkgrid"):
        # here n is static, rhomax changes
        n = 100
        f, ax = plt.subplots(1, 1, dpi=175, frameon=True)
        Rhovars = np.linspace(0.01, 15, 25)
        plotEigvalsQuantum1(ax, var="rhomax", vars=Rhovars, sim = False, eps=12, eigcount =10, datafile="dataPvar.dat", n=n)
        ax.set_title(r"$\lambda_i$ for $n = "+f"{n}$")
        plt.show()

