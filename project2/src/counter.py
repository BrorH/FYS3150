import subprocess

# for n in range(1,50):
#     subprocess.run(f"./main.out {n} 1 0 0".split())
#     print(f"done {n}")

import numpy as np
import matplotlib.pyplot as plt

m = np.array([16,18,13,10,42,148,218,287,403,582,704,887,1006,1166,1334,2117,2188,2667,2775,3133,3492,3580,4161,4613,5306,5503,5842,6606,6994,7407,7852,8453,9700,9955,10369,11312,11927,12428,14172,14370,14379,15547])
n = np.array(range(1,len(m)+1))
plt.plot(n,m, "-o")
plt.show()
import numpy as np
import matplotlib.pyplot as plt

#k = np.array([88,515, 1504, 2853, 7140, 22456])
#n = np.array([5,10,15,20,30,50])

# a,b,c = np.polyfit(n,k,2)
# print(a,b,c)

# predict = lambda x: a*x**2+b*x**1+c
# # print(predict(100))
# plt.plot(n,k, "-o", label="dat")
# n_ = np.linspace(0, 70, 100)
# plt.plot(n_, predict(n_), label="polyfit")
# plt.legend()
# plt.grid()
# plt.show()
import numpy as np
import matplotlib.pyplot as plt 

from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split

X = n.reshape(len(n), 1)
y = m#X**4 + X**3 + X + 1

x_train, x_test, y_train, y_test = train_test_split(X, y, test_size=0.3)

rmses = []
degrees = np.arange(1, 10)
min_rmse, min_deg = 1e10, 0

for deg in degrees:

    # Train features
    poly_features = PolynomialFeatures(degree=deg, include_bias=False)
    x_poly_train = poly_features.fit_transform(x_train)

    # Linear regression
    poly_reg = LinearRegression()
    poly_reg.fit(x_poly_train, y_train)

    # Compare with test data
    x_poly_test = poly_features.fit_transform(x_test)
    poly_predict = poly_reg.predict(x_poly_test)
    poly_mse = mean_squared_error(y_test, poly_predict)
    poly_rmse = np.sqrt(poly_mse)
    rmses.append(poly_rmse)
    
    # Cross-validation of degree
    if min_rmse > poly_rmse:
        min_rmse = poly_rmse
        min_deg = deg

# Plot and present results
print('Best degree {} with RMSE {}'.format(min_deg, min_rmse))
        
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(degrees, rmses)
ax.set_yscale('log')
ax.set_xlabel('Degree')
ax.set_ylabel('RMSE')
plt.show()