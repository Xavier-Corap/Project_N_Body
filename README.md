Here is a N_body simulation of collision of galaxies with a method of cloud in cells (CIC).

You have here some fortran90 and python script. Here is a list of a little description for all scripts:

- init.f90: Create the speed and position initial conditions of two galaxies with a plummer model; 
- const.f90: Where we define all constants of our N body problem;
- var.f90: Where we define all variables of our N body problem;
- grad.f90: 


To make this simulation, you can compile all the needed script with the makefile and execute the main file "main.f90". 

You can check the result (with a plummer model) by executing the file "initial.py". 
