MODULE PRESS_CRCTR_SLVR
    ! This module adapts the code from the unsteady convection-diffusion energy solver to build the coefficient matrices
    ! for the pressure correction in the SIMPLE algorithm

    ! solve pressure correction equations
    USE ADI_3D_SOLVR

    IMPLICIT NONE

    CONTAINS
        SUBROUTINE PRESS_CRCTR_SLVR_MAIN(Nx, Ny, Nz, params, del_x, del_y, del_z, &
                                        p_prime, u_star, v_star, w_star, & 
                                        aP_ndl, aN_ndl, aE_ndl, aS_ndl, aW_ndl, aB_ndl, aT_ndl, & 
                                        d_u, d_v, d_w, & 
                                        b3D_prime, &
                                        res_tol_p, max_it_ADI_p)

            ! build coefficients
            CALL PRESS_CRCTR_COEFFS_BLDR(Nx, Ny, Nz, params, del_x, del_y, del_z, &
                                        p_prime, u_star, v_star, w_star, & 
                                        aP_ndl, aN_ndl, aE_ndl, aS_ndl, aW_ndl, aB_ndl, aT_ndl, & 
                                        d_u, d_v, d_w, & 
                                        b3D_prime)
            
            ! solve pressure correction equations                           
            CALL ADI_3D_SOLVR_MAIN(res_tol_p, max_it_ADI_p, &
                                   aP_ndl, aN_ndl, aE_ndl, aS_ndl, aW_ndl, aB_ndl, aT_ndl, p_prime, b3D_prime)
        
        END SUBROUTINE PRESS_CRCTR_SLVR_MAIN

        ! coeffiient matrix builder
        SUBROUTINE PRESS_CRCTR_COEFFS_BLDR(Nx, Ny, Nz, params, del_x, del_y, del_z, &
                                        p_prime, u_star, v_star, w_star, & 
                                        aP_ndl, aN_ndl, aE_ndl, aS_ndl, aW_ndl, aB_ndl, aT_ndl, & 
                                        d_u, d_v, d_w, & 
                                        b3D_prime)

            IMPLICIT NONE

            INTEGER, PARAMETER  :: dp = KIND(1.0D0)

            ! inputs
            INTEGER,  INTENT(IN)     :: Nx, Ny, Nz

            REAL(dp), DIMENSION(:,:,:), INTENT(IN)     :: u_star, v_star, w_star ! velocity correction values
            !REAL(dp), DIMENSION(:,:,:), INTENT(INOUT)  :: p_prime ! pressure correction values
            !REAL(dp), INTENT(IN)    :: t
            REAL(dp), INTENT(IN)    :: delx, del_y, del_z

            REAL(dp), DIMENSION(:,:,:), INTENT(IN)     :: aN_stg, aE_stg, aS_stg, aW_stg, aB_stg, aT_stg
            !REAL(dp)                :: dt

            !REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(IN)   :: conv_grid_x_stg, conv_grid_y_stg, conv_grid_z_stg
            REAL(dp), DIMENSION(:),     ALLOCATABLE, INTENT(IN)   :: params

            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)     :: aP_ndl, aN_ndl, aE_ndl, aS_ndl, aW_ndl, aB_ndl, aT_ndl
            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)     :: b3D_prime
            !REAL(dp), DIMENSION(:,:,:)   , INTENT(INOUT)   :: phi3D_vel_crctr

            ! IMPLICIT
            INTEGER     :: i,j,k

            !REAL(dp)    :: diff_coeff, gamma ! alpha an_ndld specific heat_ndl

            !REAL(dp)    :: del_x, del_y, del_z, del_V, BC_E, BC_N, BC_S, BC_W, BC_B, BC_T, Lx, Ly, Lz, area, alpha, cp, rho  ! params for building coeff mat_ndlrices

            REAL(dp)    :: area_xz, area_xy, area_yz


            ! Lx          = params
            ! Ly          = params
            ! Lz          = params
        
            area_xy = del_x*del_y
            area_xz = del_x*del_z
            area_yz = del_y*del_z

            ! calculate a prime coefficients in pressure correction eqn
            DO i = 1,Ny
                DO j = 1,Nx
                    DO k = 1,Nz

                        aS_ndl(i,j,k) = rho*area_xz*d_v(i,j,k) 
                        aN_ndl(i,j,k) = rho*area_xz*d_v(i+1,j,k)
                        aW_ndl(i,j,k) = rho*area_yz*d_u(i,j,k)
                        aE_ndl(i,j,k) = rho*area_yz*d_u(i,j+1,k)
                        aB_ndl(i,j,k) = rho*area_xy*d_w(i,j,k)
                        aT_ndl(i,j,k) = rho*area_xy*d_w(i,j,k+1)
                        
                    END DO
                END DO
            END DO

            ! Build point P coeff matrix
            DO i = 1,Ny
                DO j = 1,Nx
                    DO k = 1,Nz
                
                        aP_ndl(i,j,k) = aN_ndl(i,j,k) + aE_ndl(i,j,k) + aS_ndl(i,j,k) + aW_ndl(i,j,k)+ aB_ndl(i,j,k) + aT_ndl(i,j,k)
                            
                    END DO
                END DO
            END DO
                

            ! Build constant matrix b (b_prime for pressure correction)
            DO i = 1,Ny
                DO j = 1,Nx
                    DO k = 1,Nz

                        ! continuity at cell faces in the staggered grid
                        b3D_prime(i,j,k) = ((rho*u_star(i,j,k)) - (rho*u_star(i,j+1,k))) + &
                                        ((rho*v_star(i,j,k)) - (rho*v_star(i+1,j,k))) + &
                                        ((rho*w_star(i,j,k)) - (rho*w_star(i,j,k+1)))
                    
                    END DO
                END DO
            END DO  


        END SUBROUTINE PRESS_CRCTR_COEFFS_BLDR

END MODULE PRESS_CRCT_SLVR