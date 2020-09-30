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


unsorted = list(read_data("datacomp.dat").values())
sol1 = unsorted[0:N]
sol2 = unsorted[N:2*N]
sol3 = unsorted[2*N:3*N]

for sol_ in [sol1 , sol2, sol3]:
    trfs = np.array([a.transformations for a in sol_])
    n = np.array([a.n for a in sol_])
    sol.append([trfs, n, epsilon])

names = [f"Buckling beam", r"Quantum 1, $\rho_{max} = $"+f"{rhomax[1]}" ,r"Quantum 2, $\rho_{max} = $"+f"{rhomax[2]}, $\omega = {omega[2]}$"]

def getshape(shape):
    data = {"Slope":[], "R":[], "err":[]}

    for i in range(3):
        subdict = {}
        n = np.array(sol[i][1])#[N*epsnum:N*(epsnum+1)]))
        trfs = np.array(sol[i][0])#[N*epsnum:N*(epsnum+1)]))
        #pf = np.polyfit(n, trfs, 1)
        n2 = eval(shape)#n**2
        nfit = np.linspace(n2[0], n2[-1], 1000)
       
        slope, intercept, r_value, p_value, std_err = stats.linregress(n2,trfs)
        data["Slope"].append(slope)
        data["R"].append(1-r_value)
        data["err"].append(std_err)
        print(f"a: {round(slope,3)}, b: {round(intercept,3)}, R: {r_value}, err: {std_err}")
    return data
s1 = getshape("n**2")
s2 = getshape("n**2*np.log(n)")
s3 = getshape("n**3")
datas = [s1, s2 ,s3]
systems = [{},{},{}]

for i,S_ in enumerate(datas):
    S = {}
    S["Slope"] = [s["Slope"][i] for s in datas]
    S["R"] = [s["R"][i] for s in datas]
    S["err"] = [s["err"][i] for s in datas]
  
    df = pd.DataFrame(S)
    print(
        "\n"
        + df.to_latex(
            index=True,
            float_format="%.2e",
            label=f"test",
            caption=type,
            escape=False,
            column_format="c"* 4,
        )
    )

