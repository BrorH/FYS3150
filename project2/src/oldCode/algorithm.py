import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import time


def Jacobi(M, tol, V):
    # Mprime = (M - np.diag(np.diag(M))) ** 2
    n = M.shape[0]
    # idx = np.argmax(abs(Mprime))
    # i = idx // n
    # k = idx % n
    # emax = Mprime[i, k]
    # if emax < tol:
    # return M, V
    # else:
    # while emax > tol:
    te = 0
    while True:
        Mprime = (M - np.diag(np.diag(M))) ** 2
        idx = np.argmax(abs(Mprime))
        i = idx // n
        k = idx % n
        emax = Mprime[i, k]
        if emax < tol:
            print(M)
            return np.diag(M), V
        tau = (M[i, i] - M[k, k]) / 2 / M[i, k]
        t1 = np.sqrt(1 + tau ** 2) - tau
        t2 = -np.sqrt(1 + tau ** 2) - tau
        S = makeS(i, k, min(t1, t2), n)
        V = S @ V
        M = S.T @ M @ S
    # return
    # return Jacobi(B, tol, V)


def makeS(i, k, t, n):
    S = np.identity(n)
    c = 1 / np.sqrt(1 + t ** 2)
    s = t * c
    S[i, i] = c
    S[k, k] = c
    S[i, k] = -s
    S[k, i] = s
    return S


def main():
    n = 50
    pmax = 1
    h = pmax / n
    hhneg = h ** -2

    A = np.zeros((n, n))
    d = 2 * hhneg
    a = -hhneg
    for i in range(n):
        if i != 0:
            A[i, i - 1] = a
        if i != n - 1:
            A[i, i + 1] = a
        A[i, i] = d
    e = 1e-4

    # D, V = Jacobi(A, e, np.identity(n))
    # print(V)
    # print(D)
    np.set_printoptions(precision=3)
    val, vec = np.linalg.eig(A)
    val = sorted(val)
    for i in range(n):
        print(val[i])
    # for i in range(n):
    #     for k in range(n):
    #         print(f"{str(round(vec[i, k] ,4))}", end=" ")
    #     print()
    # print("Made it thourgh!")
    sys.exit()
    analytic_formula = [d + 2 * a * np.cos((i + 1) * np.pi / n) for i in range(n)]
    vector_formula = np.zeros((n, n))
    for j in range(1, n + 1):
        vector_formula[j - 1] = [np.sin(k * j * np.pi / n) for k in range(1, n + 1)]
    print(analytic_formula)
    print(vector_formula)
    # S(4, 2, np.pi / 6, n)
    for i in range(n):
        print(vec[:, 0] @ vector_formula[i])

    # print()
    # print(
    # (A @ V[:, 0]) / np.linalg.norm(A @ V[:, 0]), V[:, 0] / np.linalg.norm(V[:, 0])
    # )
    # c = ["r", "b", "g", "k", "c"]
    # for i in range(n):
    #     plt.plot(V[i], f"{c[i]}--", label=rf"$\lambda_{i}$")
    #     plt.plot(vec[i], c[i])
    # plt.legend()
    # plt.show()


if __name__ == "__main__":
    main()
