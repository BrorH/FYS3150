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


pmax = 0.5
eps = 12
newsim = False
newsimP = False
omegas = np.array([0.01, 0.5, 1, 5])
Nres = 6
Ns = np.arange(10, 200, Nres)



if newsim:
    open("data/dataQuantum2.dat", "w+").close()
    for N in Ns:
        for omega in omegas:
            subprocess.run(f"./main.out {N}{omega} {N} {eps} {pmax} 2 {omega} {'dataQuantum2.dat'}".split())
            print(f"N: {N}/{Ns[-1]} ω: {round(omega,2)}, {round(100*omega/omegas[-1],2)}%")
   

pmaxs = np.linspace(0.5, 10, 50)
N = 50
if newsimP:
     open("data/dataQuantumP.dat", "w+").close()
     for rhomax in pmaxs:
        for omega in omegas:
            subprocess.run(f"./main.out {rhomax}{omega} {N} {eps} {rhomax} 2 {omega} {'dataQuantumP.dat'}".split())
            print(f"rho: {rhomax}/{pmaxs[-1]} ω: {round(omega,2)}, {round(100*omega/omegas[-1],2)}%")
   



sols = read_data("dataQuantum2.dat")
solsp = read_data("dataQuantumP.dat")
Ncolors = list(Color("red").range_to(Color("cyan"),len(omegas)))
Ncolors = [str(color) for color in Ncolors]

with plt.style.context("seaborn-darkgrid"):
    f, ax = plt.subplots(1, 1, dpi=175, frameon=True)
    ax.set_xlabel("$N$")
    ax.set_ylabel("Eigval")
    
    for N in Ns: 
        #eigvals = np.zeros(len(omegas))
        for i, omega in enumerate(omegas):
            eigval =  np.sort(sols[str(N)+str(omega)].eigvals)[0]
            #print(f"eig: {round(eigval,2)}, N: {N}, ω:{omega}")
            ax.plot(N, eigval, 'o',color = Ncolors[i])
  
    
    legend = ax.legend([f"$\omega_r={omega}$" for omega in omegas],fancybox=True, framealpha=0, shadow=True, borderpad=0.4, frameon = True,loc='upper right')
    #legend.get_frame().set_facecolor('white')
    
plt.title(r"$\rho_{max}=$"+str(pmax) +r",$\epsilon=10^{-"+str(eps)+r"}$"+f"$N  \in [{Ns[0]}, {Ns[-1]}]$")
plt.show()


with plt.style.context("seaborn-darkgrid"):
    f, ax = plt.subplots(1, 1, dpi=175, frameon=True)
    ax.set_xlabel(r"$\rho_{max}$")
    ax.set_ylabel("Eigval")
    
    for rhomax in pmaxs: 
        #eigvals = np.zeros(len(omegas))
        for i, omega in enumerate(omegas):
            eigval =  np.sort(solsp[str(rhomax)+str(omega)].eigvals)[0]
            #print(f"eig: {round(eigval,2)}, N: {N}, ω:{omega}")
            ax.plot(rhomax, eigval, 'o',color = Ncolors[i])

    
    legend = ax.legend([f"$\omega_r={omega}$" for omega in omegas],fancybox=True, framealpha=0, shadow=True, borderpad=0.4, frameon = True,loc='upper right')
    #legend.get_frame().set_facecolor('white')
    
plt.title(r"$\rho_{max}\in ["+f"{pmaxs[0]}, { pmaxs[-1]}" +r"] ,\epsilon=10^{-"+str(eps)+r"}"+f"N={N}$")
plt.show()

