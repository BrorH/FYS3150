import numpy as np
import matplotlib.pyplot as plt 
import subprocess
from datareader import read_data
import time
import matplotlib as mpl
from scipy import optimize

"""
Plots the igenval accuracy agains a dynamic N and estimates a curve fit in order to estimate neccesairy n for 4 digits accuracy

"""

def plotEigvalAccuracyDynamicN(ax,Ns, eps = 12, pmax=20, datafile="data.dat", sim=False):
    
    if sim:
        open("data.dat", "w+").close()
        for n in Ns:
            subprocess.run(f"./main.out {n} {n} {eps} {pmax} 1 1".split())
            print(f"N: {n}, {round(100*(n-Ns[0])/(Ns[-1]-Ns[0]), 2)} %")
        subprocess.run("cp data.dat dataFindAcc.dat".split())

    data = read_data(datafile)
    avals = np.array([3,7,11,15]) # analytic eigenvals

    
    vals = np.zeros(len(Ns)) # array to be filled
    for i,n in enumerate(Ns):
        vals[i] = -np.log10(1/4*np.sum(np.abs(np.sort(data[str(n)].eigvals)[:4]-avals)))
   
    ax.plot(Ns, vals, color="red", lw=3, alpha=1, label="Numerical values")

    # using scipys optimize.curve_fit we find the best coefficient for the model a*log10(n)+b
    f = lambda t,a,b: a*np.log10(t)+b
    popt, pcov = optimize.curve_fit(f, Ns[100:], vals[100:])

    x = np.linspace(Ns[100], 4000, 1000)
    ax.plot(x,f(x,*popt), color="blue", lw=2, alpha=0.5, label=f"Fit ${round(popt[0],2)}\log"+r"_{10}"+f"(n){round(popt[1],3)}$")
    

    ax.set_xlabel("Matrix size, n")
    ax.set_ylabel("No. of correct digits")
    legend = ax.legend(fancybox=True, framealpha=1, shadow=True, borderpad=0.4, frameon = True,loc='lower right')
    legend.get_frame().set_facecolor('white')
    ax.set_title(f"Number of correct digits in eigenvalues")


if __name__ == "__main__":
    with plt.style.context("seaborn-darkgrid"):
        f, ax = plt.subplots(1, 1, dpi=200, frameon=True)
        Ns = np.arange(10, 350)
        plotEigvalAccuracyDynamicN(ax,Ns, eps = 12, pmax=20, datafile="dataFindAcc copy 2.dat", sim=False)
        
    plt.show()