import numpy as np
import pandas as pd 
from datareader import read_data
import subprocess

n = 100
eps = 12
open("data.dat", "w").close()
rhoms = [3.613, 4.048, 4.456, 4.827]
for rhom in rhoms:
    subprocess.run(f"./main.out {rhom} {n} {eps} {rhom} 1 1".split())
sols = read_data()
avals = np.array([4*i +3 for i in range(4)])
data = {r"$\rho_{max}$":rhoms, "$\lambda_1$":[], "$\lambda_2$":[], "$\lambda_3$":[], "$\lambda_4$":[], r"$E_{rel}$":[]}
for i in range(4):
    numvals = np.sort(sols[str(rhoms[i])].eigvals)
    data[r"$E_{rel}$"].append(np.sum(np.abs(numvals[:4]-avals)/avals))
    for j in range(4):
        data[f"$\lambda_{j+1}$"].append(numvals[j])


df = pd.DataFrame(data)
print(
        "\n"
        + df.to_latex(
            index=False,
            float_format="%.6e",
            label=f"test",
            escape=False,
        )
    )