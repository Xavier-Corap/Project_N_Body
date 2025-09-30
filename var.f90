module var
    ! Utilisation du module contenant les constantes
    use const

    implicit none

    ! Déclaration des variables pour les positions (x), les vitesses (v), et les accélérations (a)
    real(kind=xp), dimension(:,:), allocatable :: x
    real(kind=xp), dimension(:,:), allocatable :: v
    real(kind=xp), dimension(:,:), allocatable :: a

    ! Déclaration de la masse de chaque particule
    real(kind=xp) :: M(Np)

    ! Déclaration du pas de temps
    real(kind=xp) :: dt

    ! Déclaration du résidu
    real(kind=xp) :: res = 1

    ! Déclaration des tableaux pour la densité (rho), le potentiel gravitationnel (phi),
    ! et le champ de gravité (g_field)
    real(kind=xp), dimension(:,:,:), allocatable :: rho
    real(kind=xp), dimension(:,:,:), allocatable :: phi
    real(kind=xp), dimension(:,:,:,:), allocatable :: g_field
end module var

