MODULE CIC

   USE const  ! Module contenant les constantes
   USE var    ! Module contenant les variables

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE calculate_cell_density()
   
    INTEGER(xp) :: i, j, k, l, ip, jp, kp, ip1, jp1, kp1
    REAL(kind = xp) :: cell_size, fx, fy, fz

    ! Initialisation de la grille (rho) à zéro
    rho = 0

    ! Calcul de la taille de chaque cellule
    cell_size = dx
    
    ! Remplissage de la grille (cloud in cell) en 3D
    
    !$OMP PARALLEL PRIVATE(i,j,k,ip,jp,kp,ip1,jp1,kp1,fx,fy,fz,l)
    !$OMP DO SCHEDULE(DYNAMIC)

    DO l = 1, Np
      


      i = floor((x(l, 1)- cell_size/2)/ cell_size) + 1
      j = floor((x(l, 2)  - cell_size/2)/ cell_size) + 1
      k = floor((x(l, 3) - cell_size/2)/ cell_size)  + 1
      
      
      ip = MODULO(i-1, N_grille)+1
      jp = MODULO(j-1, N_grille)+1
      kp = MODULO(k-1, N_grille)+1
      
      ip1 = MODULO(i, N_grille)+1
      jp1 = MODULO(j, N_grille)+1
      kp1 = MODULO(k, N_grille)+1   

      fx = (x(l, 1) - ((i-1)*cell_size+cell_size/2))/cell_size
      fy = (x(l, 2) - ((j-1)*cell_size+cell_size/2))/cell_size
      fz = (x(l, 3) - ((k-1)*cell_size+cell_size/2))/cell_size
      
      ! Ajout de la masse de la particule à chaque cellule voisine pondérée par sa distance
      rho(ip, jp, kp) =       rho(ip, jp, kp)  + (1-fx)*(1-fy)*(1-fz) *M(l)
      rho(ip, jp, kp1) =     rho(ip, jp, kp1)  + (1-fx)*(1-fy)*fz*M(l)
      rho(ip, jp1, kp1) =   rho(ip, jp1, kp1)  + (1-fx)*(fy)  *fz*M(l)
      rho(ip, jp1, kp) =     rho(ip, jp1, kp)    + (1-fx)*(fy)  *(1-fz)*M(l)
      rho(ip1, jp1, kp) =   rho(ip1, jp1, kp)  + (fx)  *(fy)  *(1-fz)*M(l)
      rho(ip1, jp, kp) =     rho(ip1, jp, kp)  + (fx)  *(1-fy)*(1-fz) *M(l)
      rho(ip1, jp, kp1) =   rho(ip1, jp, kp1)  + (fx)  *(1-fy)*fz  *M(l)
      rho(ip1, jp1, kp1) = rho(ip1, jp1, kp1)  + (fx)  *(fy)  *fz  *M(l)
      
      
      
    END DO

    !$OMP END DO
    !$OMP END PARALLEL
    
    ! Normalisation de la densité par le volume de la cellule
    rho = rho / cell_size**3
        
  END SUBROUTINE calculate_cell_density

END MODULE CIC

