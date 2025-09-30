Here is a N_body simulation of collision of galaxies with a method of cloud in cells (CIC).

You have here some fortran90 and python script. Here is a list of a little description for all scripts:

-const.f90 : module comportant toutes les constantes de départ de notre projet

-var.f90: module containing the variables of our project (position, velocity, gravitational field, etc.).

-init.f90: module for setting up the initial conditions of our project.

-solv.f90: module containing the Jacobi calculation function and its residual for calculating the gravitational potential.

-grad.f90: module containing the gradient calculation function used to compute the gravitational field by deriving gravitational potential.

-pos.f90: module containing the function that updates the particle positions at each time step (using implicit Euler).

-CIC.f90: module implementing the CIC method for our problem.

-interp.f90: inverse CIC.

-main.f90: module with the main time loop. It creates the speeds and positions initial conditions in files "pos_0.dat" and "vit_0.dat" respectively. Then speed, density and position at each time step.

-main.py: module creating gifs and images of density, position and speed by using others python modules like "images.py", "vitesses.py", "rhobis.py" and "separation.py". It will give gifs of the evolution of all these variables over time.

