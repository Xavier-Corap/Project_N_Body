MODULE grad
    USE const
    USE var
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE grad_func() 
        INTEGER(KIND = xp) :: i, j, k
        INTEGER(KIND = xp) :: xi1, yj1, zk1, xi2, yj2, zk2

        ! Calcul du gradient de la fonction de potentiel gravitationnel
        !$OMP PARALLEL PRIVATE(i, j, k, xi1, yj1, zk1, xi2, yj2, zk2)
        !$OMP DO SCHEDULE(DYNAMIC)

        DO i = 1, N_grille
            DO j = 1, N_grille
                DO k = 1, N_grille

                    xi1 = i
                    yj1 = j
                    zk1 = k
                    xi2 = i
                    yj2 = j
                    zk2 = k

                    ! Conditions périodiques aux limites
                    IF (i == 1) THEN
                        xi1 = xi1 + N_grille
                    END IF
                    IF (k == 1) THEN
                        zk1 = zk1 + N_grille
                    END IF
                    IF (j == 1) THEN 
                        yj1 = yj1 + N_grille
                    END IF
                    IF (i == N_grille) THEN
                        xi2 = xi2 - (N_grille)
                    END IF
                    IF (j == N_grille) THEN
                        yj2 = yj2 - (N_grille)
                    END IF
                    IF (k == N_grille) THEN
                        zk2 = zk2 - (N_grille)
                    END IF

                    ! Calcul des composantes du gradient
                    g_field(i, j, k, 1) = -(phi(xi2 + 1, j, k) - phi(xi1 - 1, j, k)) / (2 * dx)
                    g_field(i, j, k, 2) = -(phi(i, yj2 + 1, k) - phi(i, yj1 - 1, k)) / (2 * dx)
                    g_field(i, j, k, 3) = -(phi(i, j, zk2 + 1) - phi(i, j, zk1 - 1)) / (2 * dx)
                
                END DO
            END DO
        END DO

        !$OMP END DO
        !$OMP END PARALLEL
    END SUBROUTINE grad_func
END MODULE grad

