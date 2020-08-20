import numpy as np
import matplotlib.pyplot as plt

def read(filename):
    with open(filename, "r+") as file:
        lines = file.read().split("\n")
    conf = [float(i) for i in lines[0].split('\n')[0].split(' ')]
    dat = [float(i) for i in lines[1:-1]]
    return np.array(dat), conf



f = lambda y: 1-(1-np.exp(-10))*y-np.exp(-10*y)

x = read("x.dat")[0]
v = read("v.dat")[0]
#print(v)
plt.plot(x,v, label='num')
plt.plot(x,f(x), label='ana')
plt.grid()
plt.legend()
plt.show()
