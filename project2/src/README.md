# Project 2

## Files in project
### main.cpp
Main C++ file. Contains function for making diagonal of matrix A, and main, which controlls Solver-class. Usecase desribed below.
### solver.h
Contains the outline of Solver-class, a general solver for tridiagonal eigenvalueproblems.
### solver.cpp
Flesh of the Solver-class. Contains the actual Jacobi-algorithm for solving eigenvalue problem.
### master.py
Python-script to enslave the C++-files. Not yet written.
### plot.py
Plotting. Not yet written.

## Usecase
### master.py
When written, will take inputparameters from commandline:
- n, size of system
- problem, which problem to solve. Determines "behavior" given to main.cpp
- (maybe) max_rho, maximum value of rho-domain
Will also take outputparameters from commandline. How to handle the data.
- printing
- plotting
- saving
- killing
- yeeting
- publishing

### main.cpp
Takes 4 inputparameters:
- n, size of system
- behaviour (0, 1, 2), which problem to solve. See "diag()" in file.
- rho_max, systemparameter
- omega, systemparameter (only for behaviour = 1, 2)

