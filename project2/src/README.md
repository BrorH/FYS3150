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


### Compiling
g++ -o main.out main.cpp solver.cpp -larmadillo -O3
g++ -o arma.out arma_beam.cpp -larmadillo -O3