import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import subprocess
from datareader import read_data
import time
import matplotlib as mpl
from scipy import optimize
N = 40 # rememebr to re-run sim when changing N
newsim = False

start = time.time()
sol = []


rhomax  = [1, 5, 10]
omega = [0, 1, 0.1]
epsilon = 1e-12
if newsim:
    print("Solving ...")
    open("data.dat", "w").close() # clear file
    for method in [0, 1, 2]:
        subtime = time.time()
        #for epsilon in np.linspace(8, 16, N):
        for n in np.linspace(60, 100, N):
            subprocess.run(f"./main.out {n} {epsilon} {rhomax[method]} {method} {omega[method]}".split())
            print(f"Method: {method}: {round(100*(n-60)/40, 2)} %")
        print(f"{method} done in {round(time.time() -subtime,3)} s.")
    print(f"All done in {round(time.time() -start,3)} s.")

    #     solutions = read_data()
    #     trfs = np.array([a.transformations for a in solutions])
    #     n = np.array([a.n for a in solutions])
    #     eps = np.array([a.eps for a in solutions])
    #     sol.append([trfs, n, eps])
    # print(sol)

unsorted = read_data()

sol1 = unsorted[0:N]
sol2 = unsorted[N:2*N]
sol3 = unsorted[2*N:3*N]

for sol_ in [sol1 , sol2, sol3]:
    trfs = np.array([a.transformations for a in sol_])
    n = np.array([a.n for a in sol_])
    #eps = np.array([a.eps for a in sol_])
    sol.append([trfs, n, epsilon])

names = [f"Buckling beam", r"Quantum 1, $\rho_{max} = $"+f"{rhomax[1]}" ,r"Quantum 2, $\rho_{max} = $"+f"{rhomax[2]}, $\omega = {omega[2]}$"]

with plt.style.context("seaborn-darkgrid"):
    f, ax1 = plt.subplots(1, 1, dpi=200, frameon=True)
    ax1.set_xlabel("Matrix size, $n$")
    ax1.set_ylabel("Transformations, $m$")

    pf = [1]
    colors = ['r', 'b', 'g']
    pfcolors = ['firebrick', 'mediumblue', "darkgreen"]
    for i in range(3):
        n = np.array(sol[i][1])#[N*epsnum:N*(epsnum+1)]))
        trfs = np.array(sol[i][0])#[N*epsnum:N*(epsnum+1)]))
        #pf = np.polyfit(n, trfs, 1)
        nfit = np.linspace(n[0], n[-1], 1000)
        #print(pf)
        fitf = lambda t,a,b:  a*t**2+b*t**2*np.log(t)
        fit = optimize.curve_fit(fitf,  n,  trfs)
        a,b = tuple(fit[0])
        ax1.plot(nfit, fitf(nfit, a, b), color=pfcolors[i], linestyle='--', alpha=1,lw=1.5)
        ax1.plot(n, trfs,'-o', markersize=1.2, color=colors[i], label=names[i], alpha=0.4, lw=2)
       
    legend = ax1.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1, frameon = True)
    frame = legend.get_frame()
    frame.set_facecolor('white')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.title(r"Transformations for $\epsilon = $"+f"{epsilon}" )
plt.tick_params(top='on', bottom='on', left='on', right='on', labelleft='on', labelbottom='on')
plt.show()
