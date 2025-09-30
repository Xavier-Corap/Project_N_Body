program main

    use var
    use const
    use pos
    use init
    use solv
    use interp 
    use CIC
    use iso_c_binding


    implicit none

    real(kind = xp) :: Vmax

    real(kind = xp) :: t = 0, xx, yy

    integer(kind = xp) :: pas = 0, i, it_num, j

    real(kind = xp) :: C_factor = 0.4, i_real
    
    character(len=8) :: fmt
    character(len=5) :: x1
    character(len=20) :: filename
    
    fmt = '(I5.5)'

    ! Allocation des tableaux
    allocate( a(Np, 3) )
    allocate( x(Np, 3) )
    allocate( v(Np, 3) )
    allocate( phi(N_grille, N_grille, N_grille) )
    allocate( rho(N_grille, N_grille, N_grille) )
    allocate( g_field(N_grille, N_grille, N_grille, 3) )

    ! Initialisation des masses
    M = mp

    ! Initialisation des positions
    x = 0.5

    ! Génération aléatoire de positions
    DO i = 2, Np
        CALL RANDOM_NUMBER( xx )
        CALL RANDOM_NUMBER( yy ) 
        x(i,1) = xx
        x(i,2) = yy
    END DO

    ! Initialisation des variables
    call Plummer()

    ! Décalage des positions
    x(:,1) = x(:,1) -  N_grille/8*dx
    x(:,2) = x(:,2) -  N_grille/8*dx

    ! Calcul des densités des cellules
    call calculate_cell_density()

    ! Ouverture des fichiers de sortie
    open(1, file="pos_0.dat")
    open(2, file="vit_0.dat")

    ! Calcul de la vitesse maximale initiale
    Vmax = maxval(sqrt(v(:,1)**2 + v(:,2)**2 + v(:,3)**2))
    v(:, 2) = v(:, 2) + 150000
    v(1,1) = 0
    v(1,2) = - 150000
    v(1,3) = 0
    
    ! Affichage de quelques informations
    print *, maxval(rho),  sum(abs(sqrt(v(:,1)**2 + v(:,2)**2 + v(:,3)**2)))/Np/1000

    ! Calcul du pas de temps initial
    dt = min(C_factor*sqrt(3*pi/(32*G*maxval(rho))), C_factor*dx/(Vmax))
    
    ! Écriture des données initiales dans les fichiers
    DO i= 1,Np
        write(1, *) x(i,:)
        write(2, *) v(i,:)
    END DO

    ! Boucle principale
    do while(pas < 3000)
        print *, "iter", pas
        pas = pas+1
        t = t + dt
        print *, 't, dt :', t/(60*60*24*365)/1e6,"M.années" ,dt/(60*60*24*365)/1e6, "M.années"
        
        ! Calcul des accélérations
        call pos_euler()

        ! Calcul de la nouvelle vitesse maximale
        Vmax = maxval(sqrt(v(:,1)**2 + v(:,2)**2 + v(:,3)**2))

        ! Calcul du nouveau pas de temps
        dt = min(sqrt(C_factor*3*pi/(32*G*maxval(rho))),C_factor*dx/(Vmax))
   
        ! Mise à jour des vitesses et des positions
        v = v + a*dt/2 
        x = x + dt*v 

        ! Écriture des positions dans un fichier à intervalles réguliers
        if (MOD(pas, 25) == 0 .or. pas == 1) then
            it_num = pas
            write(x1, fmt) it_num
            filename = 'pos_' // trim(x1) // '.dat'
            open(6, file="images/"//trim(filename))
            do i=1,Np
                write(6, *) x(i,:)
            end do
            close(6)

            ! Écriture des densités dans un fichier
            OPEN(UNIT=7, FILE='images_density/'//trim(filename))
            DO i = 1, N_grille
                DO j = 1, N_grille
                    WRITE(7, '(ES15.6)', ADVANCE='NO') SUM(rho(i, j, :))
                END DO
                WRITE(7, *)
            END DO
            
            ! Écriture des vitesses dans un fichier
            filename = 'vit_' // trim(x1) // '.dat'
            open(8, file="vitesses/"//trim(filename))
            do i=1,Np
                write(8, *) v(i,:)
            end do
            close(8)

            CLOSE(10)
        end if


    end do

    ! Affichage du temps écoulé
    print *, 't:', t/(60*60*24*365)/1e6,"millions d'années" 
 
    close(1)
    close(2)

    
    contains
    
    ! Fonction pour convertir un entier en chaîne de caractères
    character(len=10) function int2str(i)
        integer, intent(in) :: i
        write(int2str, '(I0)') i
    end function int2str

end program main

