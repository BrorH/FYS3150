import numpy as np
import matplotlib.pyplot as plt 
import subprocess
from datareader import read_data
import time
from colour import Color

"""
Produces figure that shows the optimal rhomax values for a certain n in quantum 1

"""


def plotOptimalRhomaxQuantum1(ax, pms = [], n = 100, sim = False, datafile = "data.dat", eps = 12 ):
    # pms is short for rhomax's (s for plural as it contains all the rhomax values to plot against)

    # given a dataset with rhomax (pms) finds which rhomax makes each of the first 4 eigvals of quantum 1 the most accurate
    # produces plot as well

    if sim:
        print("Solving ...")
        open("data.dat", "w").close() # clear file
        start = time.time()
        for rhomax in pms:
            subprocess.run(f"./main.out {n}{rhomax} {n} {eps} {rhomax} 1 1".split())
            print(f"ρ: {round(rhomax,3)}/{pms[-1]}, {round(100*(rhomax-pms[0])/(pms[-1]-pms[0]), 2)} %")
        print(f"Done in {round(time.time()- start,3)} s")
        subprocess.run(f"cp data.dat {datafile}".split())

    sols = read_data(datafile) # read solutions


    colors = list(Color("red").range_to(Color("blue"),4)) # prettyyyy colorrssss
    colors = [str(color) for color in colors]
    

    analytical = np.array([4*a+3 for a in range(0,4)]) # analytical eigenvalues

    erel = np.zeros((len(pms),4)) #error_rel, to be filled
    for i,rhomax in enumerate(pms):
        erel[i] = np.abs(np.sort(sols[f"{n}{rhomax}"].eigvals)[:4]-analytical)/analytical # fill with all relative errors
    
    minidx = erel.argmin(axis=0) # find where the graph has its lowest point, a.k.a where the error is the lowest

    # plotting code
    for i in range(4):
        ax.plot(pms, erel[:,i], color = colors[i])
    for i,idx in enumerate(minidx):
        print(f"eig_{i+1} Best rhomax: {round(pms[idx],5)}, err_rel: {round(erel[idx][i],6)}") # print the results
        ax.plot(pms[idx], erel[idx][i], "-o", label = f"$\lambda_{i+1}: $"+r"$\rho_{max}="+f"{round(pms[idx],3)}$,"+r"$\epsilon_{rel}="+"%.2e$" %erel[idx][i],color = colors[i])
    legend = ax.legend(fancybox=True, framealpha=1, shadow=True, borderpad=0.4, frameon = True)
    legend.get_frame().set_facecolor('white')
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    
    
    
if __name__ == "__main__":
    with plt.style.context("seaborn-darkgrid"):
        n = 100 # static
        f, ax = plt.subplots(1, 1, dpi=175, frameon=True)
        Rhovars = np.linspace(3.2,5.7, 1000) # Find which rhomax in this range gives the least error

        plotOptimalRhomaxQuantum1(ax, pms=Rhovars, sim = False, eps=12,  datafile="dataPvar.dat", n=n)
        ax.set_title(r"$\epsilon_{rel}$ for four first eigenvalues")
        ax.set_xlabel(r"$\rho_{max}$")
        ax.set_ylabel(r"$\epsilon_{rel}$")
        plt.show()


