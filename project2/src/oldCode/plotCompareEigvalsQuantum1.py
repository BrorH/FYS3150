import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import subprocess
from datareader import read_data
import time
import matplotlib as mpl
from colour import Color
import matplotlib as mpl
from scipy import optimize


pmax = 20   
eps = 12
newsim = False
Nmax = 500
Nmin = 10
interval = 10
Nindexes = np.array(range(Nmin, Nmax, interval))
if newsim:
    open("data.dat", "w+").close()
    for n in Nindexes:
        subprocess.run(f"./main.out {n} {n} {eps} {pmax} 1 1".split())
        print(f"N: {n}, {round(100*(n-Nmin)/(Nmax-Nmin), 2)} %")
    subprocess.run("cp data.dat dataFindAcc.dat".split())

data = read_data("dataFindAccBIG.dat")
avals = np.array([3,7,11])
print(Nindexes)
with plt.style.context("seaborn-darkgrid"):
    f, ax = plt.subplots(1, 1, dpi=100, frameon=True)
    
    ax.set_xlabel("Matrix size, N")
    ax.set_ylabel("No. of correct digits")
    vals = np.zeros(len(Nindexes))
    for i,n in enumerate(Nindexes):
        vals[i] = -np.log10(1/3*np.sum(np.abs(np.sort(data[str(n)].eigvals)[:3]-avals)))
    f = lambda t,a,b: a*np.log10(t/b)
    popt, pcov = optimize.curve_fit(f, Nindexes[10:], vals[10:])


    x = np.linspace(Nindexes[10], 3000, 2000)
    
    ax.plot(x, f(x, *popt), color="blue", lw=3.5, alpha=0.6,label=f"fit ${round(popt[0],2)}"+r"\log_{10}(N) +"+f"{round(popt[1],2)}$")
       
    #ax.plot(Nindexes, Nanalytical, label="Analytical", color="black", lw=2, alpha=0.8, linestyle="--")
    
    ax.plot(Nindexes, vals, color="red", lw=3, alpha=0.6, label="Numerical values")
    
    #ax.plot((x[0], f(x,popt[0], popt[1])[0]), markersize=10, color="blue")
    #ax.set_ylim(top=Nanalytical[-1]*1.1, bottom=0)
    #ax.legend()
    #ax.set_title(r"$\lambda_i$ for $\rho_{max} = "+f"{rhomax}$")
    #ax.set_xlim(1, NmaxEigIdx)
    #ax.set_xticks(list(range(1,NmaxEigIdx+1)))
    #Ncmap= mpl.colors.ListedColormap(Ncolors)
    #Nnorm= mpl.colors.Normalize(vmin=Nmin,vmax=Nmax)
    #f.colorbar(mpl.cm.ScalarMappable(norm=Nnorm, cmap=Ncmap), ax=ax, orientation='horizontal', label='$N$', ticks=[Nmin, (Nmax+Nmin)//2, Nmax])
    legend = ax.legend(fancybox=True, framealpha=1, shadow=True, borderpad=0.4, frameon = True,loc='lower right')
    frame = legend.get_frame()
    frame.set_facecolor('white')
#plt.title(r"Single-electron: Num o. correct digits in first 3 eigenvalues. $\rho_{max}=$"+str(pmax) +r",$\epsilon=10^{-"+str(eps)+r"}$")
plt.show()

