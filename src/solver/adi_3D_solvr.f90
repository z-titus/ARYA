! make sure you do this for your future self

MODULE ADI_3D_SOLVR

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

        SUBROUTINE ADI_3D_SOLVR_MAIN(res_tol, max_it, aP, aN, aE, aS, aW, aB, aT, phi3D, b3D)
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

            REAL(dp), DIMENSION(:,:,:)   , INTENT(IN)   :: aP, aN, aE, aS, aW, aB, aT 
            REAL(dp), DIMENSION(:,:,:)   , INTENT(IN)   :: b3D

            REAL(dp), INTENT(IN)    :: res_tol
            INTEGER, INTENT(IN)     :: max_it

            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)   :: phi3D ! output solved for

            INTEGER  :: i, k, z_xsec  ! indexs in main loop

            ! FOR residuals
            REAL(dp), DIMENSION(max_it+1) :: res      ! residual for each iteration
            ! store previous iterations phi solution for residual calc
            REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)  :: phi_old 
            
            ! imlicit variables
            INTEGER     :: dimy, dimx, dimz   ! grid dimensions    

            ! initialization
            dimy = SIZE(b3D, DIM = 1)
            dimx = SIZE(b3D, DIM = 2)
            dimz = SIZE(b3D, DIM = 3)
            ALLOCATE(phi_old(dimy, dimx, dimz))
            phi_old = 0.0

            ! ADI algorithm
            DO i = 1,max_it     

                ! ADI for every x - y plane in z direction
                DO k = 1,dimz

                    z_xsec = k
                    ! alternating sweeps. each sweep updates phi2D
                    CALL SWP_WEST_EAST(aP, aN, aE, aS, aW, aB, aT, b3D, dimx, dimy, dimz, z_xsec, phi3D)
                    CALL SWP_EAST_WEST(aP, aN, aE, aS, aW, aB, aT, b3D, dimx, dimy, dimz, z_xsec, phi3D)

                    !phi_old = phi2D
                
                END DO

                CALL CALC_RMS_3D(aP, aN, aE, aS, aW, aB, aT, b3D, dimx, dimy, dimz, phi3D, res, i)
                    
                ! termination
                IF (i > 1 .AND. (res(i)/res(1)) < res_tol) THEN 
                    PRINT *, 'ADI Converged on iteration', i
                    EXIT
                ELSEIF (i == max_it) THEN
                    PRINT *, 'ADI Max number of iterations reached', 'r1=',res(1) 
                    PRINT *, 'Final normalized residual: ', res(i)/res(1)
                END IF

            END DO

        END SUBROUTINE ADI_3D_SOLVR_MAIN


        SUBROUTINE SWP_WEST_EAST(aP, aN, aE, aS, aW, aB, aT, b3D, dimx, dimy, dimz, z_xsec, phi3D)

            INTEGER, PARAMETER  :: dp = KIND(1.0D0)
            ! sweep east to west with lines that span north to south
            
            REAL(dp), DIMENSION(:,:,:)   , INTENT(IN)   :: aP, aN, aE, aS, aW, aB, aT   
            REAL(dp), DIMENSION(:,:,:)   , INTENT(IN)   :: b3D

            INTEGER,  INTENT(IN)    :: dimx, dimy, dimz, z_xsec

            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)   :: phi3D ! output solved for

            REAL(dp), DIMENSION(dimy, dimy)    :: A ! coefficient matrix for 1D TDMA

            REAL(dp), DIMENSION(dimy, dimx)    :: b2D_temp, phi2D_temp ! implicit temporary b constant vector that gets updated with aE and aW terms

            INTEGER     :: i, j, m, n, k

            ! isolate the plane of interest in z
            b2D_temp = b3D(:,:,z_xsec)
            phi2D_temp = phi3D(:,:,z_xsec)

            A = 0.0

            DO i = 1,(dimx) ! sweep ->

                DO m = 1,(dimy)
                    DO n = 1,dimy
                        IF (n == m) THEN
                            A(m,n) =    aP(m,i,z_xsec)
                        
                        ELSEIF (n == m-1) THEN
                            A(m,n) = -1*aS(m,i,z_xsec)
                           
                        ELSEIF (n == m+1) THEN
                            A(m,n) = -1*aN(m,i,z_xsec)
                        
                        END IF
                    END DO
                END DO

                DO j = 1,(dimy)

                    IF (i == 1 .AND. (z_xsec == 1)) THEN
                        b2D_temp(j,i) = b3D(j,i,z_xsec) + aE(j,i,z_xsec)*phi2D_temp(j, i+1)  + aT(j,i,z_xsec)*phi3D(j,i,z_xsec+1)

                    ELSEIF (i == 1 .AND. (z_xsec == dimz)) THEN
                        b2D_temp(j,i) = b3D(j,i,z_xsec) + aE(j,i,z_xsec)*phi2D_temp(j, i+1)  + aB(j,i,z_xsec)*phi3D(j,i,z_xsec-1)

                    ELSEIF (i == 1 .AND. (z_xsec .NE. 1) .AND. (z_xsec .NE. dimz)) THEN
                        b2D_temp(j,i) = b3D(j,i,z_xsec) + aE(j,i,z_xsec)*phi2D_temp(j, i+1)  + aT(j,i,z_xsec)*phi3D(j,i,z_xsec+1) &
                                        + aB(j,i,z_xsec)*phi3D(j,i,z_xsec-1)

                    ELSEIF (i == dimx .AND. (z_xsec == 1)) THEN
                        b2D_temp(j,i) = b3D(j,i,z_xsec) + aW(j,i,z_xsec)*phi2D_temp(j, i-1) + aT(j,i,z_xsec)*phi3D(j,i,z_xsec+1)
                    
                    ELSEIF (i == dimx .AND. (z_xsec == dimz)) THEN
                        b2D_temp(j,i) = b3D(j,i,z_xsec) + aW(j,i,z_xsec)*phi2D_temp(j, i-1) + aB(j,i,z_xsec)*phi3D(j,i,z_xsec-1)
                    
                    ELSEIF (i == dimx .AND. (z_xsec .NE. 1) .AND. (z_xsec .NE. dimz)) THEN
                        b2D_temp(j,i) = b3D(j,i,z_xsec) + aW(j,i,z_xsec)*phi2D_temp(j, i-1)  + aT(j,i,z_xsec)*phi3D(j,i,z_xsec+1) &
                                        + aB(j,i,z_xsec)*phi3D(j,i,z_xsec-1)

                    ELSE 
                        b2D_temp(j,i) = b3D(j,i,z_xsec) + aE(j,i,z_xsec)*phi2D_temp(j, i+1)  + aW(j,i,z_xsec)*phi2D_temp(j, i-1) &
                                        + aT(j,i,z_xsec)*phi3D(j,i,z_xsec+1) + aB(j,i,z_xsec)*phi3D(j,i,z_xsec-1)
                    END IF

                END DO


                CALL TDMA_SCALAR_MAIN(A, phi2D_temp(:,i), b2D_temp(:,i))

            END DO

            phi3D(:,:,z_xsec) = phi2D_temp(:,:)
            
        END SUBROUTINE SWP_WEST_EAST


        SUBROUTINE SWP_EAST_WEST(aP, aN, aE, aS, aW, aB, aT, b3D, dimx, dimy, dimz, z_xsec, phi3D)

            INTEGER, PARAMETER  :: dp = KIND(1.0D0)
            ! sweep east to west with lines that span north to south
            
            REAL(dp), DIMENSION(:,:,:)   , INTENT(IN)   :: aP, aN, aE, aS, aW, aB, aT   
            REAL(dp), DIMENSION(:,:,:)   , INTENT(IN)   :: b3D

            INTEGER,  INTENT(IN)    :: dimx, dimy, dimz, z_xsec

            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)   :: phi3D ! output solved for

            REAL(dp), DIMENSION(dimy, dimy)    :: A ! coefficient matrix for 1D TDMA

            REAL(dp), DIMENSION(dimy, dimx)    :: b2D_temp, phi2D_temp ! implicit temporary b constant vector that gets updated with aE and aW terms

            INTEGER     :: i, j, m, n, k

            ! form the A matrix and x and b vectors that will get passed to the 1D TDMA algorithm
            ! TDMA(A, x, b, BC1, BC2)
            ! BC1 = BC_N, BC2 = BC_S

            ! initialize
            ! isolate the plane of interest in z
            b2D_temp = b3D(:,:,z_xsec)
            phi2D_temp = phi3D(:,:,z_xsec)

            A = 0.0

            DO i = (dimx), 1, -1 ! sweep <-

                DO m = 1,(dimy)
                    DO n = 1,dimy
                        IF (n == m) THEN
                            A(m,n) =    aP(m,i,z_xsec)
                        
                        ELSEIF (n == m-1) THEN
                            A(m,n) = -1*aS(m,i,z_xsec)
                           
                        ELSEIF (n == m+1) THEN
                            A(m,n) = -1*aN(m,i,z_xsec)
                        
                        END IF
                    END DO
                END DO

                DO j = 1,(dimy)

                    IF (i == 1 .AND. (z_xsec == 1)) THEN
                        b2D_temp(j,i) = b3D(j,i,z_xsec) + aE(j,i,z_xsec)*phi2D_temp(j, i+1)  + aT(j,i,z_xsec)*phi3D(j,i,z_xsec+1)

                    ELSEIF (i == 1 .AND. (z_xsec == dimz)) THEN
                        b2D_temp(j,i) = b3D(j,i,z_xsec) + aE(j,i,z_xsec)*phi2D_temp(j, i+1)  + aB(j,i,z_xsec)*phi3D(j,i,z_xsec-1)

                    ELSEIF (i == 1 .AND. (z_xsec .NE. 1) .AND. (z_xsec .NE. dimz)) THEN
                        b2D_temp(j,i) = b3D(j,i,z_xsec) + aE(j,i,z_xsec)*phi2D_temp(j, i+1)  + aT(j,i,z_xsec)*phi3D(j,i,z_xsec+1) &
                                        + aB(j,i,z_xsec)*phi3D(j,i,z_xsec-1)

                    ELSEIF (i == dimx .AND. (z_xsec == 1)) THEN
                        b2D_temp(j,i) = b3D(j,i,z_xsec) + aW(j,i,z_xsec)*phi2D_temp(j, i-1) + aT(j,i,z_xsec)*phi3D(j,i,z_xsec+1)
                    
                    ELSEIF (i == dimx .AND. (z_xsec == dimz)) THEN
                        b2D_temp(j,i) = b3D(j,i,z_xsec) + aW(j,i,z_xsec)*phi2D_temp(j, i-1) + aB(j,i,z_xsec)*phi3D(j,i,z_xsec-1)
                    
                    ELSEIF (i == dimx .AND. (z_xsec .NE. 1) .AND. (z_xsec .NE. dimz)) THEN
                        b2D_temp(j,i) = b3D(j,i,z_xsec) + aW(j,i,z_xsec)*phi2D_temp(j, i-1)  + aT(j,i,z_xsec)*phi3D(j,i,z_xsec+1) &
                                        + aB(j,i,z_xsec)*phi3D(j,i,z_xsec-1)

                    ELSE 
                        b2D_temp(j,i) = b3D(j,i,z_xsec) + aE(j,i,z_xsec)*phi2D_temp(j, i+1)  + aW(j,i,z_xsec)*phi2D_temp(j, i-1) &
                                        + aT(j,i,z_xsec)*phi3D(j,i,z_xsec+1) + aB(j,i,z_xsec)*phi3D(j,i,z_xsec-1)
                    END IF

                END DO


                CALL TDMA_SCALAR_MAIN(A, phi2D_temp(:,i), b2D_temp(:,i))

            END DO

            phi3D(:,:,z_xsec) = phi2D_temp(:,:)
            
        END SUBROUTINE SWP_EAST_WEST

        SUBROUTINE CALC_RMS_3D(aP, aN, aE, aS, aW, aB, aT, b3D, dimx, dimy, dimz, phi3D, res, main_i)
            ! residuals are calculated as r = aP*phiP - (sum.anb*phinb + b)
            ! calculate rms value of r over entire grid to get residual measure

            INTEGER,  PARAMETER  :: dp = KIND(1.0D0)

            REAL(dp), DIMENSION(:,:,:)   , INTENT(IN)   :: phi3D, aP, aN, aE, aS, aW, aB, aT, b3D

            REAL(dp), DIMENSION(:), INTENT(OUT) :: res
            
            INTEGER, INTENT(IN)  :: main_i, dimx, dimy, dimz

            REAL(dp)    :: r_sqr, sum_r_sqr, avg, rms_res, r, temp

            INTEGER     :: i,j,k

            r_sqr        = 0.0
            sum_r_sqr    = 0.0
            
            DO i = 1, dimy
                DO j = 1, dimx
                    DO k = 1, dimz

                        temp = 0.0

                        IF ((i-1) .GE. 1) THEN
                            temp = temp + aW(i,j,k)*phi3D(i-1,j,k)
                        END IF

                        IF ((i+1) .LE. dimx) THEN
                            temp = temp + aE(i,j,k)*phi3D(i+1,j,k)
                        END IF

                        IF ((j-1) .GE. 1) THEN
                            temp = temp + aS(i,j,k)*phi3D(i,j-1,k)
                        END IF

                        IF ((j+1) .LE. dimy) THEN
                            temp = temp + aN(i,j,k)*phi3D(i,j+1,k)
                        END IF

                        IF ((k-1) .GE. 1) THEN
                            temp = temp + aB(i,j,k)*phi3D(i,j,k-1)
                        END IF

                        IF ((k+1) .LE. dimz) THEN
                            temp = temp + aT(i,j,k)*phi3D(i,j,k+1)
                        END IF

                        r = (temp + b3D(i,j,k)) - aP(i,j,k)*phi3D(i,j,k)

                        r_sqr = r**2
                        sum_r_sqr = sum_r_sqr + r_sqr

                    END DO
                END DO
            END DO
           
            
            avg = sum_r_sqr/(dimx*dimy*dimz)

            rms_res = SQRT(avg)

            res(main_i) = rms_res
            
        END SUBROUTINE CALC_RMS_3D

END MODULE ADI_3D_SOLVR