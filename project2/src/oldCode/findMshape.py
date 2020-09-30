import subprocess, time, sys, inspect
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import cm
import matplotlib as mpl
from datareader import read_data
from argHandler import argHandler
from colour import Color
from master import run
from scipy import optimize

def transforms(n=np.arange(20,350,3), rhomax = [1,10,20], eps = 12, omega = [0,1,5], method = [0,1,2], sim=False, datafile="POLYdata.dat"):
    # plots num of transformations against N on axis ax, as well as comparing to n**2 and n**2*ln(n)
    if sim:
        print("Solving ...")
        
        open(f"data/{datafile}", "w+").close() # clear file
        for method_ in method:
            start = time.time()
            for N in n:
                subprocess.run(f"./main.out {N}{eps}{rhomax[method_]}{method_}{omega[method_]} {N} {eps} {rhomax[method_]} {method_} {omega[method_]} {datafile}".split())
                print(f"{round(100*(N-n[0])/(n[-1]-n[0]), 2)} %, N: {N}, eps: {eps}, method: {method_}, rhomax: {rhomax[method_]}, omega: {omega[method_]}")
            print(f"Done in {round(time.time() -start,3)} s.")
        #if datafile != "data.dat":
            # move content to designated file
            #subprocess.run(f"cp data.dat {datafile}".split())
    
    sols = read_data(datafile) #all solutions
    for method_ in [0,1,2]: # all methods are plotted
        trans = np.array([sols[f"{N}{eps}{rhomax[method_]}{method_}{omega[method_]}"].transformations for N in n]) # array of counted transformations
        print(f" eps: {eps}, method: {method_}, rhomax: {rhomax[method_]}, omega: {omega[method_]}:")
        coeff, residual, rank, sv, cond = np.polyfit(np.log10(n), np.log10(trans), deg =1, full=True)
        print(residual, rank, sv, cond)
        print(f"a= {coeff[0]} +- {residual[0]}")
        plt.plot(np.log10(n), np.log10(trans), label=str(method_))
    plt.show()
transforms(sim=False, method=[0,1])

