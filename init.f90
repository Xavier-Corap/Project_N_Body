module init

    use var
    use const


    implicit none

    contains


    ! Fonction qui servait initialement à tester notre projet à plus petite échelle.

    subroutine two_bodies()

        !call RANDOM_NUMBER(x) 
        integer(kind = xp) :: i
        real(kind = xp) :: r

            
            ! Conditions initiales simples : vx = GM/R , vy = vz = 0, ax = 0, ay= -GM/R², az = 0.


            r = sqrt((x(2,1) - x(1,1))**2 + (x(2,2) - x(1,2))**2 + (x(2,3) - x(1,3))**2)

            v(1,:) = 0
            
            a(1,:) = 0

            v(2,:) = 0

            v(2,2) = (G*M(1)/r)**0.5


    end

    

    subroutine Plummer()

        real(kind =xp), dimension(Np) :: y, r, theta, phi_
        real(kind =xp) :: xal, Vn, u, Mtot
        integer(kind =xp) :: i
        real(kind = xp) :: a=0.05*15000*3.08e16

        ! On tire selon une loi uniforme 

        call random_number(y) 
        call random_number(theta)
        call random_number(phi_)

        ! On remplace les expression de manière à bien avoir une distribution uniforme sur une surface sphérique.

        theta = theta*2.0_xp*pi
        phi_ = acos(2.0_xp*phi_-1.0_xp)

        ! On tire r aléatoirement à l'aide de y 
        r = a*(1.0_xp/(y**(-2.0_xp/3.0_xp)-1.0_xp))**0.5


        ! On remplit les positions tout en recentrant dans notre cube.
        x(:,1) = r*sin(phi_)*cos(theta) + N_grille/2*dx
        x(:,2) = r*sin(phi_)*sin(theta) + N_grille/2*dx
        x(:,3) = r*cos(phi_) + N_grille/2*dx

        call random_number(theta)
        call random_number(phi_)

        theta = theta*2*pi
        phi_ = acos(2*phi_-1)

        ! ajoute d'une masse centrale si nécessaire

        !x(1,1) = N_grille/2*dx
        !x(1,2) = N_grille/2*dx
        !x(1,3) = N_grille/2*dx
!
        !M(1) = 10**4 * Mp

        Mtot = sum(M)

        do i=1, Np
            call random_number(xal)
            call random_number(u)
            ! On applique la méthode de rejet avec la loi normalisée du poly fourni et en divisant par un majorant 2.15 calculé numériquement en python.
            do while (u >= (1.0_xp-xal**2.0_xp)**(7.0_xp/2.0_xp)*xal**2.0_xp/0.04295_xp/2.15_xp) 

                call random_number(xal)
                call random_number(u)
            
            end do


            ! On multiplie par la vitesse de libération pour obtenir la vitesse non uniformisée.
            xal = (r(i)**2 + a**2)**(-1.0_xp/4.0_xp)*xal*sqrt(2*G*Mtot)
            
            v(i,1) = xal*cos(theta(i))*sin(phi_(i))
            v(i,2) = xal*sin(theta(i))*sin(phi_(i))
            v(i,3) = xal*cos(phi_(i))

        end do

            !Uniquement avec masse centrale
            !v(1,:) = 0

    end subroutine

    !Nécesité de cette fonction pour avoir une distribution gaussienne utile plus tard pour le calcul des dispersions

    function gauss(mean, variance) result(x)
        real(kind=xp), intent(in) :: mean, variance
        real(kind=xp) :: x, rv, u

        call random_number(u)
        call random_number(rv) 
        
        x = mean + variance * sqrt(-2 * Log(u)) * Cos(2 * pi * rv) 



    end function 


    ! On définit la fonction de densité de Miyamoto-Nagai dont nous aurons besoin plus tard pour la méthode de rejet.
    function density_miy(r, z, a, b, M) result(f)
        real(kind=xp), intent(in) :: r,z,a,b,M
        real(kind=xp) :: f

        f = b**2._xp*M/(4._xp*pi)*(a*r**2._xp+(a+3._xp*sqrt(z**2._xp + b**2._xp)*(a + sqrt(z**2._xp + b**2._xp))**2._xp)) &
        /((r**2._xp + (a + sqrt(z**2._xp + b**2._xp))**2._xp)**(5._xp/2._xp)*(z**2._xp + b**2._xp)**(3._xp/2._xp))



    end function 

    subroutine Miyamoto()

        real(kind =xp), dimension(Np) ::  theta
        real(kind =xp) :: Mtot, z, r, Vc, K, Omega, s_density, v_r, v_theta, v_z, sig_theta, xal
        integer(kind =xp) :: i
        real(kind = xp) :: a = 0.05*15000*3.08e16, z_0, b, Q = 1.2_xp, s_0

        ! On définit les constantes caractéristiques

        s_0 = 10_xp**8._xp*2._xp*10._xp**30*(3.08e19)**(-2._xp)
        z_0 = 0.2*a
        b = 0.1*a

        Mtot = sum(M)

        !On tire aléatoirement theta selon une loi uniforme 
        
        call random_number(theta)

        theta = theta*2.0_xp*pi
        
        do i = 1, Np

            !on applique ici la méthode de réjet comme explicité dans notre rapport

            call random_number(z)
            call random_number(r)
            call random_number(xal)

            z = -z*20._xp*b + (1._xp - z)*20._xp*b
            r = r*100._xp*a

            ! On définit bien un majorant à la loi de densité rho(r,z) comme étant rho(0,0)
            do while (xal >= density_miy(r,z,a,b,Mtot)/(density_miy(0._xp,0._xp,a,b,Mtot))) 
                call random_number(z)
                call random_number(r)
                call random_number(xal)


                z = -z*20._xp*b + (1._xp - z)*20._xp*b
                r = r*100._xp*a
            end do


            !On remplit les positions en recentrant dans notre cube.
            x(i,1) = r*cos(theta(i)) + N_grille/2*dx
            x(i,2) = r*sin(theta(i)) + N_grille/2*dx
            x(i,3) = z + N_grille/2*dx


            ! On calcule la vitesse circulaire 
            Vc = sqrt(G*mtot*r**2._xp*(r**2._xp+(a+ sqrt(b**2._xp + z**2._xp))**2._xp)**(-3.0_xp/2.0_xp))
            ! On définit la densité surfacique
            s_density = s_0*exp(-r/a)
            ! On calcule la vitesse angulaire et la fréquence épicyclique
            Omega = sqrt(G*Mtot*(r**2._xp + (a+b)**2._xp)**(-3.0_xp/2.0_xp)) 
            K = sqrt(-3._xp*G*mtot*r*(r**2._xp + (a+ b)**2._xp)**(-5._xp/2.0_xp) + 4._xp*omega)


            ! On calcule la dispersion v_r tout en vérifient que cela ne mènera pas à une vitesse selon theta imaginaire.
            v_r = Q*3.36_xp * G*s_density/K

            if (v_r**2._xp*(1._xp-K/(4._xp*Omega)-2._xp*r/a) + Vc**2._xp < 0) then

                ! On remplace donc par une autre expression corrigée de la densité surfacique et on recalcule v_r
                s_density = s_0*exp(-sqrt(r**2._xp + 2*(a/4._xp)**2._xp)/(a))
                v_r = Q*3.36_xp * G*s_density/K
                
                if (v_r**2._xp*(1._xp-K/(4._xp*Omega)-2._xp*r/a) + Vc**2._xp < 0) then

                    print *, "Erreur, v_theta non réel."
    
                end if

            end if
            
            ! On calcule la vitesse v_theta et les dispersions correspondantes suivant une gaussienne de moyenne 0 et d'écart type les dispersions calculées.
            v_theta = sqrt(v_r**2._xp*(1-K/(4._xp*Omega)-2._xp*r/a) + Vc**2._xp)
            v_z = sqrt(pi*s_density*G*z_0)
            sig_theta = v_r**2._xp*sqrt(K/(4._xp*Omega))
            sig_theta = gauss(0._xp,sig_theta)
            v_r = gauss(0._xp, v_r)
            v_z = gauss(0._xp, v_z)

            ! On rajoute à chaque composante la valeur de vitesse selon theta ainsi que toutes les dispersions.
            v(i,1) = - (v_theta + sig_theta)*sin(theta(i)) + v_r*cos(theta(i))
            v(i,2) = (v_theta + sig_theta)*cos(theta(i)) + v_r*sin(theta(i))
            v(i,3) = v_z

            
        end do


    end subroutine


        



end module init

