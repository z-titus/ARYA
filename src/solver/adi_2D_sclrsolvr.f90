MODULE ADI_2D_SCLRSOLVR

    ! This module implements the TDMA solver in alternating directions (ADI) to solve for scalars on a 2D grid.
    ! Sweeps west-east, east-west
    
    ! INPUTS to main subroutine -
    ! phi_old, res = the solution (phi) of the previous iteration passed from MAIN, residual from current iteration
    ! a_{xx}, phi2D = updated coefficient matrices and solution matrix of the current iteration in MAIN
    ! b2D = constant matrix for TDMA solver (Ax = b)
    ! main_i = iteration number of MAIN
    USE TDMA_SCALAR

    IMPLICIT NONE

    CONTAINS

        SUBROUTINE ADI_2D_SCLRSOLVR_MAIN(res_tol, max_it, aP, aN, aE, aS, aW, phi2D, b2D)
            ! SUBROUTINE DESCRIPTION
            ! each node has the discretized eqn a_p*phi_p = sum(a_nb*phi_nb) + b [EQN 1]
            ! this solver iteratively solves for phi on a 2D grid
            ! input aP, a_nb, and b based on problem physics in the form of EQN 1
            ! considered converged when residual of first iteration is reduced by a factor of res_tol by the mth iteration

            ! INPUTS
            ! aP, aN, aE, aS, aW (:,:) == 2D arrays for the point, north, east, south, west coefficients of [1]
            ! phi2D (:,:) == initialized 2D array of the spatial scalar solution
            ! b2D (:,:)   == initialized source term of eqn [1] for each node on the 2D grid

            ! NOTE that the row indexing starts from the bottom (the south most points) row 1 ^^ row N
            ! The column indexing starts from the west most point column 1 --> column N

            INTEGER, PARAMETER  :: dp = KIND(1.0D0)

            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: aP   
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: aN   
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: aE   
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: aS   
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: aW   
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: b2D

            REAL(dp), INTENT(IN)    :: res_tol
            INTEGER, INTENT(IN)     :: max_it

            REAL(dp), DIMENSION(:,:)   , INTENT(OUT)   :: phi2D ! output solved for

            INTEGER  :: i  ! index of the main loop

            ! FOR residuals
            REAL(dp), DIMENSION(max_it+1) :: res      ! residual for each iteration
            ! store previous iterations phi solution for residual calc
            REAL(dp), ALLOCATABLE, DIMENSION(:,:)  :: phi_old 
            
            ! imlicit variables
            INTEGER     :: dimy, dimx   ! grid dimensions    

            ! initialization
            dimy = SIZE(b2D, DIM = 1)
            dimx = SIZE(b2D, DIM = 2)
            ALLOCATE(phi_old(dimy, dimx))
            phi_old = 0.0

            ! ADI algorithm
            DO i = 1,max_it     

                ! alternating sweeps. each sweep updates phi2D
                CALL SWP_WEST_EAST(aP, aN, aE, aS, aW, b2D, dimx, dimy, phi2D)
                CALL SWP_EAST_WEST(aP, aN, aE, aS, aW, b2D, dimx, dimy, phi2D)
                        
                ! calculate residual
                CALL CALC_RMS_2D(phi_old, phi2D, res(i))
                
                ! termination
                IF (i > 1 .AND. (res(i)/res(1)) < res_tol) THEN 
                    PRINT *, 'ADI Converged on iteration', i
                    EXIT
                ELSEIF (i == max_it) THEN
                    PRINT *, 'ADI Max number of iterations reached'
                END IF

                phi_old = phi2D

            END DO

            ! update the previous solution values and pass to MAIN
            phi_old = phi2D

        END SUBROUTINE ADI_2D_SCLRSOLVR_MAIN


        SUBROUTINE SWP_WEST_EAST(aP, aN, aE, aS, aW, b2D, dimx, dimy, phi2D)

            INTEGER, PARAMETER  :: dp = KIND(1.0D0)
            ! sweep east to west with lines that span north to south
            
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: aP   
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: aN   
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: aE   
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: aS   
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: aW   
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: b2D

            INTEGER,  INTENT(IN)    :: dimx, dimy

            REAL(dp), DIMENSION(:,:)   , INTENT(OUT)   :: phi2D ! output solved for

            REAL(dp), DIMENSION(dimy, dimy)    :: A ! coefficient matrix for 1D TDMA

            REAL(dp), DIMENSION(dimy, dimx)    :: b2D_temp   ! implicit temporary b constant vector that gets updated with aE and aW terms

            INTEGER     :: i, j, m, n

            ! initialize
            b2D_temp = b2D

            DO i = 1,(dimx) ! sweep ->

                DO m = 1,(dimy)
                    DO n = 1,dimy
                        IF (n == m) THEN
                            A(m,n) = aP(m,i)
                        
                        ELSEIF (n == m-1) THEN
                            A(m,n) = -1*aS(m,i)
                           
                        ELSEIF (n == m+1) THEN
                            A(m,n) = -1*aN(m,i)
                        
                        END IF
                    END DO
                END DO

                DO j = 1,(dimy)

                    IF (i == 1) THEN
                        b2D_temp(j,i) = b2D(j,i) + aE(j,i)*phi2D(j, i+1)  
                    ELSEIF (i == dimx) THEN
                        b2D_temp(j,i) = b2D(j,i) + aW(j,i)*phi2D(j, i-1)
                    ELSE 
                        b2D_temp(j,i) = b2D(j,i) + aE(j,i)*phi2D(j, i+1)  + aW(j,i)*phi2D(j, i-1)
                    END IF

                END DO


                CALL TDMA_SCALAR_MAIN(A, phi2D(:,i), b2D_temp(:,i))

            END DO
            
        END SUBROUTINE SWP_WEST_EAST


        SUBROUTINE SWP_EAST_WEST(aP, aN, aE, aS, aW, b2D, dimx, dimy, phi2D)

            INTEGER, PARAMETER  :: dp = KIND(1.0D0)
            ! sweep east to west with lines that span north to south
            
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: aP   
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: aN   
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: aE   
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: aS   
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: aW   
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: b2D
 
            INTEGER,  INTENT(IN)    :: dimx, dimy

            REAL(dp), DIMENSION(:,:)   , INTENT(OUT)   :: phi2D ! output solved for

            ! form the A matrix and x and b vectors that will get passed to the 1D TDMA algorithm
            ! TDMA(A, x, b, BC1, BC2)
            ! BC1 = BC_N, BC2 = BC_S

            REAL(dp), DIMENSION(dimy, dimy)    :: A ! coefficient matrix for 1D TDMA

            REAL(dp), DIMENSION(dimy, dimx)    :: b2D_temp   ! implicit temporary b constant vector that gets updated with aE and aW terms

            INTEGER     :: i, j, m, n

            ! initialize
            b2D_temp = b2D

            DO i = (dimx), 1, -1 ! sweep <-
                

                DO m = 1,(dimy)
                    DO n = 1,dimy
                        IF (n == m) THEN
                            A(m,n) = aP(m,i)
                        
                        ELSEIF (n == m-1) THEN
                            A(m,n) = -1*aS(m,i)
                           
                        ELSEIF (n == m+1) THEN
                            A(m,n) = -1*aN(m,i)
                        
                        END IF
                    END DO
                END DO

                DO j = 1,dimy
                    IF (i == 1) THEN
                        b2D_temp(j,i) = b2D(j,i) + aE(j,i)*phi2D(j, i+1)  
                    ELSEIF (i == dimx) THEN
                        b2D_temp(j,i) = b2D(j,i) + aW(j,i)*phi2D(j, i-1)
                    ELSE 
                        b2D_temp(j,i) = b2D(j,i) + aE(j,i)*phi2D(j, i+1)  + aW(j,i)*phi2D(j, i-1)
                    END IF
                END DO

                CALL TDMA_SCALAR_MAIN(A, phi2D(:,i), b2D_temp(:,i))

            END DO
            
        END SUBROUTINE SWP_EAST_WEST

        SUBROUTINE CALC_RMS_2D(r1, r2, rms_res)
            ! Input 2 2D arrays to calculate the element wise RMS between them

            INTEGER,  PARAMETER  :: dp = KIND(1.0D0)

            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: r1   
            REAL(dp), DIMENSION(:,:)   , INTENT(IN)   :: r2
            
            REAL(dp), INTENT(OUT)      :: rms_res

            REAL(dp)    :: sum_sqr, temp, avg

            INTEGER     :: i,j

            sum_sqr = 0.0
            temp    = 0.0
            
            DO i = 1, SIZE(r1, dim = 1)
                DO j = 1, SIZE(r1, dim = 2)
                    
                    temp = (r2(i,j) - r1(i,j))**2
                    sum_sqr = sum_sqr + temp
                   
                END DO
            END DO
           
            
            avg = sum_sqr/(SIZE(r1, dim = 1)*SIZE(r1, dim = 2))

            rms_res = SQRT(avg)
            
        END SUBROUTINE CALC_RMS_2D

END MODULE ADI_2D_SCLRSOLVR