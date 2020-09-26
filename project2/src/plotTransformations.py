import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import subprocess
from datareader import read_data
import time
import matplotlib as mpl
from scipy import optimize, stats
import pandas as pd
N = 20 # rememebr to re-run sim when changing N
newsim = False
start = time.time()
sol = []


rhomax  = [1, 5, 10]
omega = [0, 1, 0.1]
epsilon = 1e-8
if newsim:
    print("Solving ...")
    open("data.dat", "w").close() # clear file
    for method in [0, 1, 2]:
        subtime = time.time()
        #for epsilon in np.linspace(8, 16, N):
        for n in np.linspace(100, 450, N):
            subprocess.run(f"./main.out {method}{n} {n} {epsilon} {rhomax[method]} {method} {omega[method]}".split())
            #print(f"Method: {method}: {round(100*(n-100)/250, 2)} %")
        #print(f"{method} done in {round(time.time() -subtime,3)} s.")
    print(f"All done in {round(time.time() -start,3)} s.")

    #     solutions = read_data()
    #     trfs = np.array([a.transformations for a in solutions])
    #     n = np.array([a.n for a in solutions])
    #     eps = np.array([a.eps for a in solutions])
    #     sol.append([trfs, n, eps])
    # print(sol)

unsorted = list(read_data().values())
#print(unsorted)
sol1 = unsorted[0:N]
sol2 = unsorted[N:2*N]
sol3 = unsorted[2*N:3*N]

for sol_ in [sol1 , sol2, sol3]:
    trfs = np.array([a.transformations for a in sol_])
    n = np.array([a.n for a in sol_])
    #eps = np.array([a.eps for a in sol_])
    sol.append([trfs, n, epsilon])

names = [f"Buckling beam", r"Quantum 1, $\rho_{max} = $"+f"{rhomax[1]}" ,r"Quantum 2, $\rho_{max} = $"+f"{rhomax[2]}, $\omega = {omega[2]}$"]
data = {"Slope":[], "R":[], "err":[]}
with plt.style.context("seaborn-darkgrid"):
    f, ax1 = plt.subplots(1, 1, dpi=200, frameon=True)
    ax1.set_xlabel("Matrix size, $n$")
    ax1.set_ylabel("Transformations, $m$")

    pf = [1]
    colors = ['r', 'b', 'g']
    pfcolors = ['firebrick', 'mediumblue', "darkgreen"]

    for i in range(3):
        subdict = {}
        n = np.array(sol[i][1])#[N*epsnum:N*(epsnum+1)]))
        trfs = np.array(sol[i][0])#[N*epsnum:N*(epsnum+1)]))
        #pf = np.polyfit(n, trfs, 1)
        n2 = n**2
        nfit = np.linspace(n2[0], n2[-1], 1000)
       
        slope, intercept, r_value, p_value, std_err = stats.linregress(n2,trfs)
        data["Slope"].append(slope)
        data["R"].append(1-r_value)
        data["err"].append(std_err)
        print(f"a: {round(slope,3)}, b: {round(intercept,3)}, R: {r_value}, err: {std_err}")
        
        nnew = np.sqrt((trfs-intercept)/slope)
        ax1.plot(n2, slope*n2+intercept, color=pfcolors[i], linestyle='--', alpha=1,lw=1.5)
        ax1.plot(n2, trfs,'-o', markersize=1.2, color=colors[i], label=names[i], alpha=0.4, lw=2)
       
    legend = ax1.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1, frameon = True)
    frame = legend.get_frame()
    frame.set_facecolor('white')
df = pd.DataFrame(data)
print(
    "\n"
    + df.to_latex(
        index=False,
        float_format="%.2e",
        label=f"test",
        #caption=type,
        escape=False,
        #column_format="c" * _n,
    )
)
print(pf)

plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.title(r"Transformations for $\epsilon = $"+f"{epsilon}" )
plt.tick_params(top='on', bottom='on', left='on', right='on', labelleft='on', labelbottom='on')
plt.show()
