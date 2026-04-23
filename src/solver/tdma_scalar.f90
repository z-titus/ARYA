MODULE TDMA_SCALAR
    ! This module contains the subroutine the tridiagonal matrix algorithm (TDMA) solver
    ! This is built for 1D geometries and can be integrated in iterative solvers for 2D and 3D geometries

    ! This solver uses the nomenclature and procedure outline in Versteeg and Malaskera (VM) - Section 7.3

    ! This is a modified version of the original tdma.f90 code, in that the nodes don't intersect the boundaries
    ! which is common for scalars. This is important for staggered grid implementation.
    IMPLICIT NONE

    CONTAINS

        SUBROUTINE TDMA_SCALAR_MAIN(A, x, b)
            ! TDMA algorithm must accept 2 boundary conditions (BC1 is starting point and BC2 is end)
            ! A = coefficient matrix
            ! b = RHS of Ax=b (constants)
            ! x = solved for
            INTEGER, PARAMETER  :: dp = KIND(1.0D0)
            REAL(dp), DIMENSION(:,:), INTENT(IN)    :: A ! coefficient matrix
            REAL(dp), DIMENSION(:),   INTENT(OUT)   :: x ! output solved for
            REAL(dp), DIMENSION(:),   INTENT(IN)    :: b ! constant Ax = b

            REAL(dp), DIMENSION(SIZE(b))    :: alpha ! Refer to VM nomenclature
            REAL(dp), DIMENSION(SIZE(b))    :: beta
            REAL(dp), DIMENSION(SIZE(b))    :: D
            REAL(dp), DIMENSION(SIZE(b))    :: C

            REAL(dp), DIMENSION(SIZE(alpha))    :: A_frwd ! Denoted A in VM (bad naming convention)
            REAL(dp), DIMENSION(SIZE(alpha))    :: C_frwd ! Denoted C' in VM nomenclature

            REAL(dp), DIMENSION(SIZE(x))     ::  phi ! Denoted C' in VM nomenclature

            ! get alpha, beta, D, and C column vectors (VM nomenclature)
            CALL GET_COEFF_VECTORS(A, b, alpha, beta, D, C)
            CALL GET_FRWD_ELIM_VECTORS(alpha, beta, D, C, A_frwd, C_frwd)
            CALL BACK_SUB(A_frwd, C_frwd, phi)

            x = phi

        END SUBROUTINE TDMA_SCALAR_MAIN

        SUBROUTINE GET_COEFF_VECTORS(A, b, alpha, beta, D, C)

            INTEGER, PARAMETER  :: dp = KIND(1.0D0)
            REAL(dp), DIMENSION(:,:), INTENT(IN)    :: A ! coefficient matrix
            REAL(dp), DIMENSION(:),   INTENT(IN)    :: b ! constant Ax = b

            REAL(dp), DIMENSION(:),   INTENT(OUT)    :: alpha ! Refer to VM nomenclature
            REAL(dp), DIMENSION(:),   INTENT(OUT)    :: beta
            REAL(dp), DIMENSION(:),   INTENT(OUT)    :: D
            REAL(dp), DIMENSION(:),   INTENT(OUT)    :: C

            INTEGER :: idx

            ! initialize
            beta(1)  = 0
            alpha(SIZE(alpha))   = 0

            ! populate column vectors corresponding to coefficients in A
            DO idx = 1, (SIZE(alpha))

                IF (idx /= 1) THEN
                    beta(idx)  = -A(idx, idx - 1)
                END IF

                IF (idx /= SIZE(alpha)) THEN
                    alpha(idx) = -A(idx, idx + 1)
                END IF

                D(idx)     = A(idx, idx)
                C(idx)     = b(idx)
             
            END DO
    
        END SUBROUTINE GET_COEFF_VECTORS

        SUBROUTINE GET_FRWD_ELIM_VECTORS(alpha, beta, D, C, A_frwd, C_frwd)
            
            ! get vector of coefficients arising from the forward elimination step
            INTEGER, PARAMETER  :: dp = KIND(1.0D0)
            REAL(dp), DIMENSION(:),   INTENT(IN)    :: alpha ! Refer to VM nomenclature
            REAL(dp), DIMENSION(:),   INTENT(IN)    :: beta
            REAL(dp), DIMENSION(:),   INTENT(IN)    :: D
            REAL(dp), DIMENSION(:),   INTENT(IN)    :: C

            REAL(dp), DIMENSION(:),   INTENT(OUT)    :: A_frwd ! Denoted A in VM (bad naming convention)
            REAL(dp), DIMENSION(:),   INTENT(OUT)    :: C_frwd ! Denoted C' in VM nomenclature

            INTEGER :: idx

            ! populate forward elimination vectors
            DO idx = 1, SIZE(A_frwd)
                
                IF (idx == 1) THEN
                    A_frwd(idx)  = alpha(idx)/(D(idx))
                    C_frwd(idx)  = (C(idx))/(D(idx))
                ELSE
                    A_frwd(idx)  = alpha(idx)/(D(idx)-(beta(idx)*A_frwd(idx-1)))
                    C_frwd(idx)  = (beta(idx)*C_frwd(idx-1) + C(idx))/(D(idx) - (beta(idx)*A_frwd(idx-1)))
                END IF

            END DO

        END SUBROUTINE GET_FRWD_ELIM_VECTORS

        SUBROUTINE BACK_SUB(A_frwd, C_frwd, phi)
            
            INTEGER, PARAMETER    :: dp=KIND(1.0D0)
            REAL(dp), DIMENSION(:),   INTENT(IN)    :: A_frwd ! Denoted A in VM (bad naming convention)
            REAL(dp), DIMENSION(:),   INTENT(IN)    :: C_frwd ! Denoted C' in VM nomenclature

            REAL(dp), DIMENSION(:),   INTENT(OUT)   :: phi    ! solved for (basically x)

            INTEGER :: idx

            ! initialize
            phi(SIZE(phi))   = C_frwd(SIZE(C_frwd))   ! the last phi is equal to the end boundary condition

            ! back substitution - loop backwards from the last value of phi, which is known
            DO idx = (SIZE(phi)-1), 1, -1
                
                phi(idx)    = A_frwd(idx)*phi(idx+1) + C_frwd(idx)
                
            END DO


        END SUBROUTINE BACK_SUB

END MODULE TDMA_SCALAR