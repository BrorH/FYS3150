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
from colour import Color


n = 80
epsilon = 7
open("data.dat", "w+").close()
subprocess.run(f"./main.out {n} {epsilon} 1 0 0".split())
N = n+1
dat = read_data()[0]

def aev(j):
    # aev = analytic eigenvector
    res = np.zeros(n)
    for i in range(n):
        res[i] = np.sin(j*(i+1)*np.pi/(N))
    return res
avec = aev(1)
val, vec = dat.getLowestEigs()
print(avec)

red = Color("yellow")
colors_ = list(red.range_to(Color("blue"),n))
x = np.linspace(0,1, 1000)

with plt.style.context("seaborn-darkgrid"):
    f, ax1 = plt.subplots(1, 1, dpi=100, frameon=True)
    ax1.set_xlabel("Numerical eigenvecor coordinate")
    ax1.set_ylabel("Analytical eigenvector coordinate")

    pf = [1]
    for i,(v, av) in enumerate(zip(vec, avec)):
        print(v,av)
        print(colors_[i])
        ax1.plot(v,av, 'o',color=str(colors_[i]))
    #ax1.plot(x,x,lw=1, linestyle="--", color="black")
    #legend = ax1.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1, frameon = True)
    #frame = legend.get_frame()
    #frame.set_facecolor('white')
#plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.title(r"Analytical and numerical eigenvector components for $\epsilon = $"+f"{epsilon}" )
plt.tick_params(top='on', bottom='on', left='on', right='on', labelleft='on', labelbottom='on')
plt.show()
