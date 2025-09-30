MODULE CONST

    IMPLICIT NONE
    
    ! Définition des constantes en majuscules
    
    INTEGER, PARAMETER, PUBLIC :: XP = SELECTED_REAL_KIND(15)  ! Précision des réels
    
    REAL(KIND = XP), PARAMETER :: G = 6.67E-11                 ! Constante gravitationnelle
    
    INTEGER(KIND = XP), PARAMETER :: NP = 100000               ! Nombre de particules
    
    INTEGER(KIND = XP), PARAMETER :: N_GRILLE = 128           ! Taille de la grille
    
    REAL(KIND = XP), PARAMETER :: DX = 30000*3.08E16/N_GRILLE            ! Taille d'une cellule
    
    REAL(KIND = XP) :: MP = 2E30*1E5                           ! Masse d'une particule
    
    REAL(KIND = XP) :: PI = 4*ATAN(1.D0)                       ! Pi
    
    REAL(KIND = XP) :: EPS = 1E-10                             ! Pour le calcul des résidus, volontairement très bas pour forcer un maximum d'itération au premier pas de temps

    REAL(KIND = xp) :: LIM = 0.15*N_GRILLE*DX                   ! Pour le développement multipolaire de gradient conjugué

END MODULE CONST

