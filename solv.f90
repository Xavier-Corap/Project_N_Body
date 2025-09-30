MODULE SOLV

   USE var, only : rho, phi, M, x
   USE const

  IMPLICIT NONE

CONTAINS    ! solve_jac(CB), solv_gradconj(CB) et leurs subroutines associées

!=============================================================================================================
!=============================================================================================================
!                                       Subroutines principales
!=============================================================================================================
!=============================================================================================================

!-------------------------------------------------------------------------------------------------------------
!                                              Jacobi
!-------------------------------------------------------------------------------------------------------------
SUBROUTINE solve_jac()
  !!! Calcule le potentiel en fonction de la densité
  ! Utilise une zone buffer <=> grille de dimension (N_grille+2)**3 <=> grille voulue + coque externe (conditions aux bords)
  ! Le centre de cette grille (zone utile) est le potentiel que l'on recherche
  ! Les bords de cette grille (zone fantôme) sont des copies des bords opposés de la zone utile <=> conditions au bord périodiques

  IMPLICIT NONE
  INTEGER :: i, j, k, ib, jb, kb, n_iter, n_iter_max = 3000
  REAL(KIND=xp) :: flag, res_it, res
  REAL(KIND=xp), DIMENSION(:, :, :), allocatable:: phi_buff
  
  n_iter = 0
  res = eps+1

  allocate( phi_buff(N_grille+2, N_grille+2, N_grille+2) )

  DO WHILE  ( (n_iter < n_iter_max) .AND. res > eps )

    res_it = res
    n_iter = n_iter+1
    
    ! ---------- Définition de la zone buffer ---------------
    ! Zone utile:
    phi_buff(2:N_grille+1, 2:N_grille+1, 2:N_grille+1) = phi
    ! Zone fantome:
    ! Définition des surfaces exterieures du cube
    phi_buff(1, 2:N_grille+1, 2:N_grille+1) = phi(N_grille,:,:)
    phi_buff(N_grille+2, 2:N_grille+1, 2:N_grille+1) = phi(1,:,:)
    phi_buff(2:N_grille+1, 1, 2:N_grille+1) = phi(:,N_grille,:)
    phi_buff(2:N_grille+1, N_grille+2, 2:N_grille+1) = phi(:,1,:)
    phi_buff(2:N_grille+1, 2:N_grille+1, 1) = phi(:,:,N_grille)
    phi_buff(2:N_grille+1, 2:N_grille+1, N_grille+2) = phi(:,:,1)
    ! Définitions des coins du cube
    phi_buff(1, 1, 1) = phi(N_grille, N_grille, N_grille)
    phi_buff(1, N_grille+2, 1) = phi(N_grille, 1, N_grille)
    phi_buff(1, 1, N_grille+2) = phi(N_grille, N_grille, 1)
    phi_buff(1, N_grille+2, N_grille+2) = phi(N_grille, 1, 1)
    phi_buff(N_grille+2, 1, 1) = phi(1, N_grille, N_grille)
    phi_buff(N_grille+2, N_grille+2, 1) = phi(1, 1, N_grille)
    phi_buff(N_grille+2, 1, N_grille+2) = phi(1, N_grille, 1)
    phi_buff(N_grille+2, N_grille+2, N_grille+2) = phi (1, 1, 1)
    ! Calcul du potentiel phi en chaque point

    !On applique les calculs de phi
    phi = (phi_buff(3:N_grille+2, 2:N_grille+1, 2:N_grille+1) + phi_buff(2:N_grille+1, 3:N_grille+2, 2:N_grille+1)&
                  + phi_buff(2:N_grille+1, 2:N_grille+1, 3:N_grille+2) &
                  + phi_buff(1:N_grille, 2:N_grille+1, 2:N_grille+1) + phi_buff(2:N_grille+1, 1:N_grille, 2:N_grille+1)&
                   + phi_buff(2:N_grille+1, 2:N_grille+1, 1:N_grille) &
                  - 4._xp*pi*G*rho(:, :, :)*dx**2) * 1._xp / 6._xp 

    ! Une fois le premier pas de temps passé, les itérations étaient suffisamment nombreuses pour considérer que le résidu a bien convergé vers une valeur qu'on prendra maintenant comme
    !référence pour eps.
    if( n_iter == n_iter_max ) then

      eps = res
      n_iter_max = 50

    end if

    ! --------------------------------- Calcul du résidu ---------------------------------
    CALL residu_jac(phi_buff, res, flag)

    ! PRINT *, flag

    flag = norm2(phi - phi_buff(2:N_grille+1, 2:N_grille+1, 2:N_grille+1))/N_grille

  END DO

  
  ! Données finales
  PRINT *, 'Potentiel via Jacobi'
  PRINT *, 'Nb iterations:', n_iter, 'sur', n_iter_max, 'maximum'
  PRINT *, 'Résidu final:', res
  PRINT *, 'flag:', flag, eps

  

  
END SUBROUTINE solve_jac

!=============================================================================================================
!                                              Résidu
!=============================================================================================================

SUBROUTINE residu_jac(phi_buff, res, flag) 
  ! Calcule le résidu de la subroutine phi_jac
  ! Entrées:
  ! phi_buff (DIMENSION (N_grille+2)**3)    => non modifié
  ! res (DIMENSION N_grille**3)             => modifié

  IMPLICIT NONE
  INTEGER :: i, j, k, ib, jb, kb
  REAL(KIND=xp), DIMENSION(N_grille+2, N_grille+2, N_grille+2), INTENT(IN) :: phi_buff
  REAL(KIND=xp), INTENT(INOUT) :: res, flag
  REAL(KIND=xp), DIMENSION(:, :, :), allocatable :: lap_phi

  allocate( lap_phi(N_grille,N_grille,N_grille) )


  !parallélisation qui n'a pas réellement amélioré nos performances.
  !$OMP PARALLEL PRIVATE(i,j,k,ib,jb,kb)
  !$OMP DO SCHEDULE(DYNAMIC)
  DO i = 1, N_grille
    ib = i+1
    DO j = 1, N_grille
      jb = j+1
      DO k = 1, N_grille
        kb = k+1

        !Calculs du Laplacien de Phi.
        lap_phi(i,j,k) = (phi_buff(ib-1, jb, kb) + phi_buff(ib+1, jb, kb) - 2*phi_buff(ib, jb, kb)) / dx**2 &
                       + (phi_buff(ib, jb-1, kb) + phi_buff(ib, jb+1, kb) - 2*phi_buff(ib, jb, kb)) / dx**2 &
                       + (phi_buff(ib, jb, kb-1) + phi_buff(ib, jb, kb+1) - 2*phi_buff(ib, jb, kb)) / dx**2
            
      END DO
    END DO
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  res = maxval(ABS(lap_phi - 4*pi*G*rho))/maxval(4*pi*G*rho)
 


END SUBROUTINE residu_jac

!-------------------------------------------------------------------------------------------------------------
!                                              Gradient Conjugué
!-------------------------------------------------------------------------------------------------------------

  SUBROUTINE solve_gradconj(CB)
    !!! Calcule la grille de potentiel en fonction de la grille de densité
    ! CB = "per"/"mul"/"nul" => Indique les conditions aux bords à utiliser (p2riodiques/multipolaires/nulles)
    
    ! ============================= Déclarations ============================
    IMPLICIT NONE
    CHARACTER(LEN=3) :: CB
    INTEGER :: i, j, k, ib, jb, kb, n_iter, n_iter_max = 3000
    REAL(KIND=xp) :: eps = 1e-10, flag, alpha, beta, res_old
    REAL(KIND=xp), DIMENSION(:,:,:), ALLOCATABLE :: d_buff, phi_buff
    REAL(KIND=xp), DIMENSION(:,:,:), ALLOCATABLE :: res_flag, res, res_new, d_new

    ALLOCATE(  d_buff(N_grille+2, N_grille+2, N_grille+2))
    ALLOCATE(phi_buff(N_grille+2, N_grille+2, N_grille+2))
    ALLOCATE(res_flag(N_grille, N_grille, N_grille))
    ALLOCATE(     res(N_grille, N_grille, N_grille))
    ALLOCATE( res_new(N_grille, N_grille, N_grille))
    ALLOCATE(   d_new(N_grille, N_grille, N_grille))
    ! ============================= Initialisation ============================
    phi_buff = 0.
    IF (CB == 'mul') THEN 
      PRINT *, 'Calcul du potentiel via Gradient Conjugué (CB: Multipolaire)'
      CALL buffer_multipol(phi, phi_buff) !phi: Définition de la zone buffer

    ELSEIF (CB == 'per') THEN
      PRINT *, 'Calcul du potentiel via Gradient Conjugué (CB: Périodique)'
      CALL buffer_per(phi, phi_buff) !phi: Définition de la zone buffer

    ELSEIF (CB == 'nul') THEN
      PRINT *, 'Calcul du potentiel via Gradient Conjugué (CB: Nulle)'
      CALL buffer_nul(phi, phi_buff) !phi: Définition de la zone buffer

    ELSE
      PRINT *, 'Erreur: Conditions au bord invalides: ', CB
      PRINT *, "Entrez 'mul' pour Multipolaire, 'per' pour Périodique ou 'nul' pour Nulles"
      PRINT *, 'Arrêt de la subroutine'
      STOP
      
    END IF

    
    CALL residu(phi_buff, res_flag) ! Calcul du résidu initial
    res = res_flag
    CALL buffer_nul(res, d_buff) ! d: Définition de la zone buffer
    
    alpha = 0.
    beta = 0.
    n_iter = 0.

    flag = 1.- (MAXVAL(ABS(4*pi*G*rho)) - MAXVAL(ABS(res_flag)))/MAXVAL(ABS(4*pi*G*rho))

    ! ============================= Itération ============================
    PRINT *, 'Iterating ...'
    DO WHILE  ( (n_iter < n_iter_max) .AND. (flag > eps) ) 
      n_iter = n_iter+1

      !OPERATION 1 :: alpha = res_T*res / ( d_T*A*d ) 
      CALL coef_alpha(res, d_buff, alpha)

      !OPERATION 2 :: phi_new = phi + alpha*d 
      DO i = 1, N_grille
        ib = i+1
        DO j = 1, N_grille
          jb = j+1
          DO k = 1, N_grille
            kb = k+1
            !OPERATION 2 :: phi_new = phi + alpha*d 
            phi(i, j, k) = phi(i, j, k) + alpha * d_buff(ib, jb, kb)

            !OPERATION 3 :: res_new = res - alpha*A*d 
            res_new(i, j, k) = res(i, j, k) - alpha * &
                              ( (d_buff(ib-1, jb, kb) + d_buff(ib+1, jb, kb) - 2*d_buff(ib, jb, kb)) / dx**2 &
                              + (d_buff(ib, jb-1, kb) + d_buff(ib, jb+1, kb) - 2*d_buff(ib, jb, kb)) / dx**2 &
                              + (d_buff(ib, jb, kb-1) + d_buff(ib, jb, kb+1) - 2*d_buff(ib, jb, kb)) / dx**2 ) 
          END DO
        END DO
      END DO

      !OPERATION 4 :: beta = res_new_T*res_new / ( res_T*res ) 
      CALL coef_beta(res, res_new, beta)

      !OPERATION 5 :: d = r_new + beta*d
      DO i = 1, N_grille
        ib = i+1
        DO j = 1, N_grille
          jb = j+1
          DO k = 1, N_grille
            kb = k+1
            !OPERATION 5 :: d = r_new + beta*d
            d_new(i, j, k) = res_new(i, j, k) + beta * d_buff(ib, jb, kb)
          END DO
        END DO
      END DO
  
      ! ===================== Actualisation de d et de res=====================
      CALL buffer_nul(d_new, d_buff) 
      res = res_new
        
      ! ============================= Résidu et flag ============================
      IF (CB == 'per') THEN
        CALL buffer_per(phi, phi_buff) !phi: Mise à jour de la zone utile et buffer
      ELSE
        phi_buff(2:N_grille+1, 2:N_grille+1, 2:N_grille+1) = phi !phi: Mise à jour de la zone utile
      END IF
      
      res_old = MAXVAL( ABS(res_flag) )
      CALL residu(phi_buff, res_flag)
      flag =  ABS( MAXVAL(ABS(res_flag)) - res_old ) / res_old

    END DO

    ! ============================= Données finales ============================
    PRINT *, 'Nb iterations:', n_iter, 'sur', n_iter_max, 'maximum'
    PRINT *, 'Résidu final:', MAXVAL(ABS(res_flag))
    PRINT *, 'Fin de la subroutine solve_gradconj('//CB//')'
    Print *, ' '

  END SUBROUTINE solve_gradconj

!=============================================================================================================
!=============================================================================================================
!                                        Subroutines secondaires
!=============================================================================================================
!=============================================================================================================

!-------------------------------------------------------------------------------------------------------------
!                                               Résidu
!-------------------------------------------------------------------------------------------------------------

  SUBROUTINE residu(phi_buff, res) 
    ! Calcule le résidu à partir du potentiel

    IMPLICIT NONE
    INTEGER :: i, j, k, ib, jb, kb
    REAL(KIND=xp), DIMENSION(:,:,:), INTENT(IN) :: phi_buff
    REAL(KIND=xp), DIMENSION(:,:,:), INTENT(INOUT) :: res
    REAL(KIND=xp), DIMENSION(:,:,:), ALLOCATABLE :: lap_phi

    ALLOCATE( lap_phi(N_grille, N_grille, N_grille))

    DO i = 1, N_grille
      ib = i+1
      DO j = 1, N_grille
        jb = j+1
        DO k = 1, N_grille
          kb = k+1

          lap_phi(i,j,k) = (phi_buff(ib-1, jb, kb) + phi_buff(ib+1, jb, kb) - 2*phi_buff(ib, jb, kb)) / dx**2 &
                         + (phi_buff(ib, jb-1, kb) + phi_buff(ib, jb+1, kb) - 2*phi_buff(ib, jb, kb)) / dx**2 &
                         + (phi_buff(ib, jb, kb-1) + phi_buff(ib, jb, kb+1) - 2*phi_buff(ib, jb, kb)) / dx**2

        END DO
      END DO
    END DO

    res = 4*pi*G*rho - lap_phi ! b - Ax  => b=4piGrho; Ax=lap_phi; x=phi 

  END SUBROUTINE residu


!-------------------------------------------------------------------------------------------------------------
!                                       Coefficient alpha (Gradient conjugué)
!-------------------------------------------------------------------------------------------------------------

  SUBROUTINE coef_alpha(res, d_buff, alpha) 
    ! Calcule le le coefficient alpha à partir du residu et de d

    IMPLICIT NONE
    INTEGER :: i, j, k, ib, jb, kb
    REAL(KIND=xp) :: num, den
    REAL(KIND=xp), INTENT(INOUT) :: alpha
    REAL(KIND=xp), DIMENSION(:,:,:), INTENT(IN) :: res
    REAL(KIND=xp), DIMENSION(:,:,:), INTENT(IN) :: d_buff

    !OPERATION :: alpha = res_T*res / ( d_T*A*d )

    num = 0 ! numérateur
    den = 0 ! dénominateur
    
    DO i = 1, N_grille
      ib = i+1
      DO j = 1, N_grille
        jb = j+1
        DO k = 1, N_grille
          kb = k+1
          
          num = num + res(i, j, k)**2

          den = den + d_buff(ib, jb, kb) &
                    * ( (d_buff(ib-1, jb, kb) + d_buff(ib+1, jb, kb) - 2*d_buff(ib, jb, kb)) / dx**2 &
                      + (d_buff(ib, jb-1, kb) + d_buff(ib, jb+1, kb) - 2*d_buff(ib, jb, kb)) / dx**2 &
                      + (d_buff(ib, jb, kb-1) + d_buff(ib, jb, kb+1) - 2*d_buff(ib, jb, kb)) / dx**2 ) 
          
          
        END DO
      END DO
    END DO 
    
    alpha =  num / den

  END SUBROUTINE coef_alpha

!-------------------------------------------------------------------------------------------------------------
!                                       Coefficient Beta (Gradient conjugué)
!-------------------------------------------------------------------------------------------------------------

  SUBROUTINE coef_beta(res, res_new, beta) 
    ! Calcule le le coefficient beta à partir de l'ancien et du nouveau residu

    IMPLICIT NONE
    INTEGER :: i, j, k
    REAL(KIND=xp) :: num, den
    REAL(KIND=xp), INTENT(INOUT) :: beta
    REAL(KIND=xp), DIMENSION(:,:,:), INTENT(IN) :: res, res_new

    !OPERATION 4 :: beta = res_new_T*res_new / ( res_T*res ) 

    num = 0
    den = 0

    DO i = 1, N_grille
      DO j = 1, N_grille
        DO k = 1, N_grille

          num = num + res_new(i, j, k)**2
          den = den + res(i, j, k)**2

        END DO
      END DO
    END DO 
    
    beta =  num / den

  END SUBROUTINE coef_beta

!-------------------------------------------------------------------------------------------------------------
!                                       Buffer periodique
!-------------------------------------------------------------------------------------------------------------

  SUBROUTINE buffer_per(mat, mat_buff)
    ! Crée une zone buffer periodidique de mat dans mat_buff

    IMPLICIT NONE
    REAL(KIND=xp), DIMENSION(:,:,:), INTENT(IN) :: mat
    REAL(KIND=xp), DIMENSION(:,:,:), INTENT(INOUT) :: mat_buff

    ! ===================== Définition des zones buffer =====================
    ! Initialisation (arêtes du cubes non utilisées => =0)
    mat_buff = 0
    ! Zone utile:
    mat_buff(2:N_grille+1, 2:N_grille+1, 2:N_grille+1) = mat
    ! Zone fantome:
    ! Définition des surfaces exterieures du cube
    mat_buff(1, 2:N_grille+1, 2:N_grille+1) = mat(N_grille,:,:)
    mat_buff(N_grille+2, 2:N_grille+1, 2:N_grille+1) = mat(1,:,:)
    mat_buff(2:N_grille+1, 1, 2:N_grille+1) = mat(:,N_grille,:)
    mat_buff(2:N_grille+1, N_grille+2, 2:N_grille+1) = mat(:,1,:)
    mat_buff(2:N_grille+1, 2:N_grille+1, 1) = mat(:,:,N_grille)
    mat_buff(2:N_grille+1, 2:N_grille+1, N_grille+2) = mat(:,:,1)

  END SUBROUTINE buffer_per

!-------------------------------------------------------------------------------------------------------------
!                                       Buffer multipolaire
!-------------------------------------------------------------------------------------------------------------

SUBROUTINE buffer_multipol(mat, mat_buff)
    ! Crée une zone buffer periodidique de mat dans mat_buff

    ! =================================== ============ ==============================================
    ! =================================== Déclarations ==============================================
    ! =================================== ============ ==============================================

    IMPLICIT NONE
    INTEGER :: i, j, k, count_in, count_out

    ! Matrices
    REAL(KIND=xp), DIMENSION(:,:,:), INTENT(IN)    :: mat      !matrice initiale
    REAL(KIND=xp), DIMENSION(:,:,:), INTENT(INOUT) :: mat_buff !matrice buffed

    ! Paramètres du système
    REAL(KIND=xp)               :: mtot ! masse totale
    REAL(KIND=xp), DIMENSION(3) :: cdg  !distance des particules, centre de gravité et vecteur de coeffs pour dipole

    ! Coordonnées
    REAL(KIND=xp), DIMENSION(:), ALLOCATABLE         :: rp_norm, rp_all_norm 
    REAL(KIND=xp), DIMENSION(:,:), ALLOCATABLE       :: rp, rp_all, x_in, x_out   
    REAL(KIND=xp), DIMENSION(:,:,:), ALLOCATABLE     :: R_norm  
    REAL(KIND=xp), DIMENSION(:,:,:,:), ALLOCATABLE   :: R, R_out_norm     
    REAL(KIND=xp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: R_out      !coords des cellules du buffer (ligne, colonne, face du cube, x/y/z)

    ! Coefficients multipolaires
    REAL(KIND=xp), DIMENSION(3)    :: P !vecteur de coeffs pour dipole
    REAL(KIND=xp), DIMENSION(3, 3) :: Q !Matrice 3x3 de coeffs pour quadripole
    
    ! Potentiels
    REAL(KIND=xp), DIMENSION(:,:,:), ALLOCATABLE :: phi0, phi1, phi2, phi_out  

    ! =================================== ========== ==============================================
    ! =================================== Algorithme ==============================================
    ! =================================== ========== ==============================================

    cdg = 0.
    DO i = 1, Np
      cdg = cdg + M(i)*x(i,:)
    END DO
    cdg = cdg / SUM(M)

    ! =================================== Coordonnées de toutes les particules ===================================
    ! Allocations
    ALLOCATE(rp_all_norm(Np))    !distance de toutes les particules par rapport au centre de gravité global
    ALLOCATE(     rp_all(Np, 3)) !coordonnées de toutes les particules par rapport au centre de gravité global

    ! position de chacune des particules par rapport au centre de gravité global
    rp_all = 0.
    DO i = 1, Np
      rp_all(i,:) = x(i,:) - cdg
    END DO
    rp_all_norm = SQRT( rp_all(:,1)**2 + rp_all(:,2)**2 + rp_all(:,3)**2 )

    ! =================================== Tri: particules multipolaire vs évaluées individuellement ===================================
    count_out = 0 ! compte des particules évaluées individuellement              
    count_in = 0  ! compte des particules évaluées par developpement multipolaire
    DO i = 1, Np  
      IF ( rp_all_norm(i) > lim ) THEN
        count_out = count_out + 1
      ELSE
        count_in = count_in + 1
      END IF
    END DO
  
    PRINT *, 'Nombre de particules hors developpement multipolaire:', count_out, '/', Np

    ! =================================== Coordonnées des cellules du buffer ===================================
    ! Allocations
    ALLOCATE(   x_in( count_in, 3)) !coords des particules du développement multipolaire
    ALLOCATE(  x_out(count_out, 3)) !coords des particules évaluées individuellement
    ALLOCATE(     rp( count_in, 3)) !coords des particules du développement multipolaire par rapport au centre de gravité multipolaire
    ALLOCATE(rp_norm( count_in))    !distance des particules du développement multipolaire par rapport au centre de gravité multipolaire
    
    count_out = 0     
    count_in = 0 
    mtot = 0. ! Masse totale multipolaire
    cdg = 0. ! centre de gravité multipolaire
    DO i = 1, Np  
      IF ( rp_all_norm(i) > lim ) THEN
        count_out = count_out + 1
        x_out(count_out,:) = x(i,:)
      ELSE
        count_in = count_in + 1
        x_in(count_in,:) = x(i,:)
        mtot = mtot + M(i)
        cdg = cdg + M(i)*x(i,:)
      END IF
    END DO
    cdg = cdg / mtot

    ! Coordonnées des particules multipolaires par rapport au centre de gravité multipolaire
    rp = 0.
    DO i = 1, count_in
      rp(i,:) = x_in(i,:) - cdg
    END DO
    rp_norm = SQRT( rp(:,1)**2 + rp(:,2)**2 + rp(:,3)**2 )

    ! =================================== Coordonnées des cellules du buffer ===================================
    ! Allocations
    ALLOCATE(    R_norm(N_grille, N_grille, 6)) !distance des cellules du buffer au centre de gravité multipolaire (ligne, colonne, face)
    ALLOCATE(         R(N_grille, N_grille, 6, 3)) !coords des cellules du buffer au centre de gravité multipolaire (ligne, colonne, face, x/y/z)
    ALLOCATE(R_out_norm(N_grille, N_grille, count_out, 6)) !distance des cellules du buffer aux particules évaluées individuelement (ligne, colonne, particule, face)
    ALLOCATE(     R_out(N_grille, N_grille, count_out, 6, 3)) !coords des cellules du buffer aux particules évaluées individuelement (ligne, colonne, particule, face, x/y/z)

    ! Coordonnées du centre des cellules du buffer (d'abord les faces impaires, puis paire à partir des précédentes)
    DO i = 1, N_grille
      DO j = 1, N_grille
          
          ! par rapport au dévelopement multipolaire
          R(i,j,1,1) = (REAL(i)-0.5)*dx - cdg(1)
          R(i,j,1,2) = (REAL(j)-0.5)*dx - cdg(2)
          R(i,j,1,3) = (REAL(N_grille+1)-0.5)*dx - cdg(3)

          R(i,j,3,1) = (REAL(i)-0.5)*dx - cdg(1)
          R(i,j,3,2) = (REAL(N_grille+1)-0.5)*dx - cdg(2)
          R(i,j,3,3) = (REAL(j)-0.5)*dx - cdg(3)

          R(i,j,5,1) = (REAL(N_grille+1)-0.5)*dx - cdg(1)
          R(i,j,5,2) = (REAL(i)-0.5)*dx - cdg(2)
          R(i,j,5,3) = (REAL(j)-0.5)*dx - cdg(3)

          !par rapport aux particules évaluées individuellement
          DO k= 1, count_out ! ici k est utilisé comme numéro de particule
            R_out(i,j,k,1,1) = (REAL(i)-0.5)*dx -          x_out(k,1)
            R_out(i,j,k,1,2) = (REAL(j)-0.5)*dx -          x_out(k,2)
            R_out(i,j,k,1,3) = (REAL(N_grille+1)-0.5)*dx - x_out(k,3)

            R_out(i,j,k,3,1) = (REAL(i)-0.5)*dx - x_out(k,1)
            R_out(i,j,k,3,2) = (REAL(N_grille+1)-0.5)*dx - x_out(k,2)
            R_out(i,j,k,3,3) = (REAL(j)-0.5)*dx - x_out(k,3)

            R_out(i,j,k,5,1) = (REAL(N_grille+1)-0.5)*dx - x_out(k,1)
            R_out(i,j,k,5,2) = (REAL(i)-0.5)*dx - x_out(k,2)
            R_out(i,j,k,5,3) = (REAL(j)-0.5)*dx - x_out(k,3)

            R_out(i,j,k,2,3) = -0.5*dx - x_out(k,3)
            R_out(i,j,k,4,2) = -0.5*dx - x_out(k,2)
            R_out(i,j,k,6,1) = -0.5*dx - x_out(k,1)
          END DO
      END DO          
    END DO

    R(:,:,2,1) = R(:,:,1,1)
    R(:,:,2,2) = R(:,:,1,2)
    R(:,:,2,3) = -0.5*dx - cdg(3)
    R(:,:,4,1) = R(:,:,3,1)
    R(:,:,4,2) = -0.5*dx - cdg(2)
    R(:,:,4,3) = R(:,:,3,3)
    R(:,:,6,1) = -0.5*dx - cdg(1)
    R(:,:,6,2) = R(:,:,5,2)
    R(:,:,6,3) = R(:,:,5,3)
    R_norm = SQRT( R(:, :, :, 1)**2 + R(:, :, :, 2)**2 + R(:, :, :, 3)**2 ) 

    R_out(:,:,:,2,1) = R_out(:,:,:,1,1)
    R_out(:,:,:,2,2) = R_out(:,:,:,1,2)
    R_out(:,:,:,4,1) = R_out(:,:,:,3,1)
    R_out(:,:,:,4,3) = R_out(:,:,:,3,3)
    R_out(:,:,:,6,2) = R_out(:,:,:,5,2)
    R_out(:,:,:,6,3) = R_out(:,:,:,5,3)
    R_out_norm(:,:,:,:) = SQRT(R_out(:, :, :, :, 1)**2 + R_out(:, :, :, :, 2)**2 + R_out(:, :, :, :, 3)**2 ) 
    ! =================================== Coefficients multipolaires ===================================
    ! Coefficient dipolaire
    P = 0.
    DO i = 1, count_in
        P = P + M(i)*rp(i,:)
    END DO

    ! Matrice de coefficient quadripolaire
    Q = 0.
    DO i = 1, 3
      DO j = 1, 3
        DO k = 1, count_in

          IF (i == j) THEN

            Q(i,j) = Q(i,j) + M(k)* (1.5*rp(k,i)*rp(k,j) - 0.5*rp_norm(k)**2)
          ELSE

            Q(i,j) = Q(i,j) + M(k)* 1.5*rp(k,i)*rp(k,j)
          END IF
        END DO
      END DO
    END DO


    ! =================================== Calcul des potentiels ===================================
    ! Allocations
    ALLOCATE(   phi0(N_grille, N_grille, 6)) !monopole (ligne, colonne, face du cube)
    ALLOCATE(   phi1(N_grille, N_grille, 6)) !dipole (ligne, colonne, face du cube)
    ALLOCATE(   phi2(N_grille, N_grille, 6)) !quadripole (ligne, colonne, face du cube)
    ALLOCATE(phi_out(N_grille, N_grille, 6)) !potentiel des particules évaluées individuelement (ligne, colonne, face du cube)

    ! Calcul du monopole
    phi0 = -G*mtot/R_norm

    ! Calcul du dipole
    phi1 =  -G*mtot/R_norm**3 * ( R(:,:,:,1)*P(1) + R(:,:,:,2)*P(2) + R(:,:,:,3)*P(3))   

    ! Calcul du quadrupole
    phi2 = 0
    DO i = 1, 3
      DO j = 1, 3
        phi2 = phi2 - G*mtot/R_norm**4 * R(:,:,:,i)*R(:,:,:,j) * Q(i,j)    
      END DO
    END DO

    ! Calcul du potentiel des particules évaluées individuellement
    phi_out = 0.
    DO i = 1, count_out
      DO k = 1, 6
      phi_out(:, :, k) = phi_out(:, :, k) - G *M(i)/ R_out_norm(:, :, i, k)
      END DO
    END DO

    ! =================================== Définition des zones buffer ===================================
    ! Zone utile:
    mat_buff(2:N_grille+1, 2:N_grille+1, 2:N_grille+1) = mat

    ! Zone fantome:
    ! Définition des surfaces exterieures du cube
    mat_buff(           1,   2:N_grille+1, 2:N_grille+1) = phi0(:,:,6) + phi_out(:,:,6) + phi1(:,:,6) + phi2(:,:,6)
    mat_buff(  N_grille+2,   2:N_grille+1, 2:N_grille+1) = phi0(:,:,5) + phi_out(:,:,5) + phi1(:,:,5) + phi2(:,:,5)
    mat_buff(2:N_grille+1,              1, 2:N_grille+1) = phi0(:,:,4) + phi_out(:,:,4) + phi1(:,:,4) + phi2(:,:,4)
    mat_buff(2:N_grille+1,     N_grille+2, 2:N_grille+1) = phi0(:,:,3) + phi_out(:,:,3) + phi1(:,:,3) + phi2(:,:,3)
    mat_buff(2:N_grille+1,   2:N_grille+1,            1) = phi0(:,:,2) + phi_out(:,:,2) + phi1(:,:,2) + phi2(:,:,2)
    mat_buff(2:N_grille+1,   2:N_grille+1,   N_grille+2) = phi0(:,:,1) + phi_out(:,:,1) + phi1(:,:,1) + phi2(:,:,1)

  END SUBROUTINE buffer_multipol

!-------------------------------------------------------------------------------------------------------------
!                                       Buffer nul
!-------------------------------------------------------------------------------------------------------------

  SUBROUTINE buffer_nul(mat, mat_buff)
    ! Crée une zone buffer nulle de mat dans mat_buff
    
    IMPLICIT NONE
    REAL(KIND=xp), DIMENSION(:,:,:), INTENT(IN) :: mat
    REAL(KIND=xp), DIMENSION(:,:,:), INTENT(INOUT) :: mat_buff

    ! ===================== Définition des zones buffer =====================
    mat_buff = 0.
    ! Zone utile:
    mat_buff(2:N_grille+1, 2:N_grille+1, 2:N_grille+1) = mat

  END SUBROUTINE buffer_nul

END MODULE SOLV