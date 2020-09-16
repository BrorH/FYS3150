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
    while True:
        Mprime = (M - np.diag(np.diag(M))) ** 2
        idx = np.argmax(abs(Mprime))
        i = idx // n
        k = idx % n
        emax = Mprime[i, k]
        if emax < tol:
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
    S[i, k] = s
    S[k, i] = -s
    return S


def main():
    n = 5
    pmax = 1
    h = pmax / n
    hhneg = h ** -2

    A = np.zeros((n, n))
    for i in range(n):
        if i != 0:
            A[i, i - 1] = -hhneg
        if i != n - 1:
            A[i, i + 1] = -hhneg
        A[i, i] = 2 * hhneg
    e = 1e-8

    D, V = Jacobi(A, e, np.identity(n))
    print("Made it thourgh!")
    val, vec = np.linalg.eig(A)
    print(val)
    print(np.diag(D))
    # S(4, 2, np.pi / 6, n)
    print()
    print(
        (A @ V[:, 0]) / np.linalg.norm(A @ V[:, 0]), V[:, 0] / np.linalg.norm(V[:, 0])
    )


if __name__ == "__main__":
    main()
