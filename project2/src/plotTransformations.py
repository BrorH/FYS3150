import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import cm
import subprocess
from datareader import read_data
import time


"""
Plots number of transformations required (m) against matrix size n.
This example produces two plots, one for a plot showing of curves at lover n-values, and one for higer n-values.

"""


def plotTransforms(Ns, ax, sim=False, datafile="data.dat", eps = 12, pmax = [1,10,20], omega=[0,1,5]):
    # plots num of transformations against N on axis ax, as well as comparing to n**2 and n**2*ln(n)
    if sim:
        print("Solving ...")
        open("data.dat", "w").close() # clear file
        for method in [0,1,2]:
            start = time.time()
            for N in Ns:
                subprocess.run(f"./main.out {N}{eps}{pmax[method]}{method}{omega[method]} {N} {eps} {pmax[method]} {method} {omega[method]}".split())
                print(f"{round(100*(N-Ns[0])/(Ns[-1]-Ns[0]), 2)} %, N: {N}, eps: {eps}, method: {method}, pmax: {pmax[method]}, omega: {omega[method]}")
            print(f"Done in {round(time.time() -start,3)} s.")
        if datafile != "data.dat":
            # move content to designated file
            subprocess.run(f"cp data.dat {datafile}".split())
    
    solutions = read_data(datafile)
    for method in [0,1,2]: # all methods are plotted
        trans = np.array([sols[f"{N}{eps}{pmax[method]}{method}{omega[method]}"].transformations for N in Ns]) # array of counted transformations
        
        color = {0:"red", 1:"blue", 2:"green"}[method] #every method gets unique color
        label = {0:"Buckling beam", 1:r"Quantum 1, $\rho_{max} = $"+str(pmax[method]), 2:r"Quantum 2, $\rho_{max} = $"+f"{pmax[method]}, $\omega_r = {omega[method]}$"}[method]

        ax.plot(Ns, trans, "-o", markersize = 2.5, color=color, alpha=0.8, lw=2.8, label=label)
        
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        plt.tick_params(top='on', bottom='on', left='on', right='on', labelleft='on', labelbottom='on')
        legend = ax.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1, frameon = True, fontsize="x-small")
        frame = legend.get_frame()
        frame.set_facecolor('white')
        ax.set_xlabel("Matrix size, $n$")
        ax.set_ylabel("Transformations, $m$")

        print(f" eps: {eps}, method: {method}, pmax: {pmax[method]}, omega: {omega[method]}:")
       
if __name__ == "__main__":
    with plt.style.context("seaborn-darkgrid"):
        f, ax = plt.subplots(1,1, dpi=250,frameon=True)
        Ns = np.arange(5, 50)
        plotTransforms(Ns, ax, sim=False, datafile="compareTransLow.dat")
        plt.title(r"Transformations; low $n$" )
        plt.show()

    with plt.style.context("seaborn-darkgrid"):
        f, ax = plt.subplots(1,1, dpi=250, frameon=True)

        Ns = np.arange(100, 450, 5)
        plotTransforms(Ns, ax, sim=False, datafile="compareTransHigh.dat")
    
        plt.title(r"Transformations; high $n$" )
        plt.show()
        


