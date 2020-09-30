# Project 2
Requirements:
install python module "Colour". pip install Colour
requires user to create a dir data/ in the working folder of main.cpp
## Files in project


### plotter.py
Contains all the code for all the plots and graphs used in the paper.
run
    $ python3 plotter.py help
for details on use. 

### master.py
Runs simulations of custom ranges and values directly from commandline.
run 
    $ python3 master.py help
for details on use
### main.cpp
Main C++ file. Contains function for making diagonal of matrix A, and main, which controlls Solver-class. Not used directly
### solver.h
Contains the outline of Solver-class, a general solver for tridiagonal eigenvalueproblems.
### solver.cpp
Flesh of the Solver-class. Contains the actual Jacobi-algorithm for solving eigenvalue problem.

## Usecase
### master.py
When written, will take inputparameters from commandline:
- n, size of system
- problem, which problem to solve. Determines "behavior" given to main.cpp
- (maybe) max_rho, maximum value of rho-domain



### main.cpp
Takes 4 inputparameters:
- n, size of system
- behaviour (0, 1, 2), which problem to solve. See "diag()" in file.
- rho_max, systemparameter
- omega, systemparameter (only for behaviour = 1, 2)

### Compiling
g++ main.cpp solver.cpp -larmadillo
### typical inputs for main.cpp
./a.out 5 1 0 0

PYTHON MODULE COLOUR MUST BE INSTALLED