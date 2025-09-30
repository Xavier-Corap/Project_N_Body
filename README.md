Here is a N_body simulation of collision of galaxies with a method of cloud in cells (CIC).

You have here some fortran90 and python script. Here is a list of a little description for all scripts:
const.f90 : module comportant toutes les constantes de départ de notre projet

var.f90 : module comportant les variables de notre projet (position, vitesse, champ gravitationnel, etc..)

init.f90 : module permettant de mettre en place les conditions initiales de notre projet.

grad.f90 : module comportant la fonction de calcul de gradient permettant de calculer le champ gravitationnel.

solv.f90 : module comportant la fonction de calcul de Jacobi et de son residus.

pos.f90 : module comportant la fonction mettant à jour la position des particules lors d'un pas de temps (avec Euler implicite).

CIC.f90 : module de la méthode CIC pour notre problème.

interp.f90 : CIC inverse

main.f90 : module avec la boucle principale temporelle


To make this simulation, you can compile all the needed script with the makefile and execute the main file "main.f90". 

You can check the result (with a plummer model) by executing the file "initial.py". 
