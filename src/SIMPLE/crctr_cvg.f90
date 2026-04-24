MODULE CRCTR_CVG

    ! This module contains a subroutine to apply the corrections to pressure and velocity based 
    ! on results from the main solver and a subroutine to check the convergence on b prime (mass conservation)
    IMPLICIT NONE

    CONTAINS
        ! subroutine to apply corrections to velocity and pressure
        SUBROUTINE CORRECTER(Nx, Ny, Nz, p, u, v, w, p_star, p_prime, u_star, v_star, w_star, d_u, d_v, d_w, &
                             relax_vel, relax_p)

            INTEGER, PARAMETER :: dp = KIND(1.0D0)

            INTEGER :: Nx, Ny, Nz

            REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE, INTENT(IN)   :: u_star, v_star, w_star, p_star ! guesses
            REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE, INTENT(IN)   :: p_prime ! pressure correction

            REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE, INTENT(INOUT)  :: u,v,w,p ! updated values

            REAL(dp), DIMENSION(:,:,:), INTENT(IN)     :: d_u, d_v, d_w

            REAL(dp), INTENT(IN)    :: relax_vel, relax_p ! relaxation factors on corrections

            INTEGER  :: i,j,k

            DO i = 1,Ny
                DO j = 1,Nx
                    DO k = 1,Nz
                        p(i,j,k) = p_star(i,j,k) + p_prime(i,j,k)*relax_p
                    END DO
                END DO
            END DO

            DO i = 1,Ny
                DO j = 2,Nx
                    DO k = 1,Nz
                        u(i,j,k) = u_star(i,j,k) + d_u(i,j,k)*(p_prime(i,j-1,k) - p_prime(i,j,k))*relax_vel
                    END DO
                END DO
            END DO

            DO i = 2,Ny
                DO j = 1,Nx
                    DO k = 1,Nz
                        v(i,j,k) = v_star(i,j,k) + d_v(i,j,k)*(p_prime(i-1,j,k) - p_prime(i,j,k))*relax_vel
                    END DO
                END DO
            END DO

            DO i = 1,Ny
                DO j = 1,Nx
                    DO k = 2,Nz
                        w(i,j,k) = w_star(i,j,k) + d_w(i,j,k)*(p_prime(i,j,k-1) - p_prime(i,j,k))*relax_vel
                    END DO
                END DO
            END DO
        
        END SUBROUTINE CORRECTER

        ! subroutine to check the convergence on b' from the pressure correction equation
        SUBROUTINE CONTINUITY_CVG(Nx, Ny, Nz, b3D_prime, main_it, cnty_tol, cnty_ref, cnvrged)

            INTEGER, PARAMETER :: dp = KIND(1.0D0)

            INTEGER, INTENT(IN) :: Nx, Ny, Nz

            LOGICAL,  INTENT(OUT)    :: cnvrged
            REAL(dp), INTENT(IN)     :: cnty_tol ! cnty refers to continuity
            REAL(dp), INTENT(INOUT)    :: cnty_ref ! store this value after a certain amount of main iterations
            
            INTEGER, INTENT(IN) :: main_it ! iteration number of the main loop. Used to designate cnty_ref

            REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE, INTENT(IN)   :: b3D_prime

            REAL(dp)    :: cnty_check

            INTEGER  :: i,j,k

            cnty_check = 0

            ! use sum of b' terms as the check
            DO i = 1,Ny
                DO j = 1,Nx
                    DO k = 1,Nz

                        cnty_check = cnty_check + b3D_prime(i,j,k)

                    END DO
                END DO
            END DO

            ! store reference continiuty check after 3 main iterations
            IF (main_it == 3) THEN
                cnty_ref = cnty_check
            END IF

            ! if the b prime is within a designated tolerance, consider the simulation converged
            IF (((cnty_check/cnty_ref) < cnty_tol) .AND. (main_it > 3)) THEN
                cnvrged = .TRUE.
            ELSE
                cnvrged = .FALSE.
            END IF
        
        END SUBROUTINE

END MODULE CRCTR_CVG