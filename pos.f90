MODULE POS

   USE var, only : x, v, a, dt
   USE solv
   USE grad
   USE interp
   USE CIC

   IMPLICIT NONE


CONTAINS
!=============================================================================================================
!                                            EULER SEMI-IMPLICIT
!=============================================================================================================
   SUBROUTINE pos_euler()
     !!! Calcule x(t+1) et v(t+1)
     ! Euler Semi-implicite: 
     ! --> v(t+1) à partir de x(t) et v(t)
     ! --> x(t+1) à partir de x(t) et v(t+1)
      
     IMPLICIT NONE
     


     call calculate_cell_density()
     CALL solve_jac()  ! modifie phi      à partir de  rho, à modifier si on veut utiliser le gradient conjugué
     CALL grad_func()      ! modifie g_field  à partir de  phi  
     CALL interpol()    ! modifie a        à partir de  g_field
     
     v = v + dt/2*a


     END SUBROUTINE pos_euler

END MODULE POS
