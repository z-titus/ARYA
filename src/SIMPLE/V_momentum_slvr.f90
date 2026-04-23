MODULE V_MOMENTUM_SLVR
    ! This module adapts the code from the unsteady convection-diffusion energy solver to build the coefficient matrices
    ! for the momentum equations solution in the SIMPLE algorithm

    ! get coefficients to solve the u momentum eqn (u star)
    USE ADI_3D_SOLVR

    IMPLICIT NONE

    CONTAINS

        SUBROUTINE V_MOMENTUM_SLVR_MAIN(Nx, Ny, Nz, params, del_x, del_y, del_z, & 
                                        aP_stg_v, aN_stg_v, aE_stg_v, aS_stg_v, aW_stg_v, aB_stg_v, aT_stg_v, b3D_v, &
                                        u_star, v_star, w_star, p_star, &
                                        BC_N_v, BC_E_v, BC_S_v, BC_W_v, BC_T_v, BC_B_v, &
                                        d_v, &
                                        res_tol_vel, max_it_ADI_vel)

            INTEGER, PARAMETER  :: dp = KIND(1.0D0)  

            ! Stuff to keep around after
            REAL(dp), DIMENSION(:,:,:)   , INTENT(INOUT)     :: v_star
            REAL(dp), DIMENSION(:,:,:)   , INTENT(INOUT)     :: d_v

            ! domain stuff
            INTEGER,  INTENT(IN)     :: Nx, Ny, Nz
            REAL(dp), DIMENSION(:),     ALLOCATABLE, INTENT(IN)   :: params
            REAL(dp), INTENT(IN)    :: del_x, del_y, del_z

            ! boundary conditions
            REAL(dp), INTENT(IN)    :: BC_N_v, BC_E_v, BC_S_v, BC_W_v, BC_T_v, BC_B_v

            ! velocities and pressures not updated
            REAL(dp), DIMENSION(:,:,:)   , INTENT(IN)     :: u_star, w_star, p_star ! corrrection values

            ! coeff matrices
            REAL(dp), DIMENSION(:,:,:)   , INTENT(INOUT)     :: aP_stg_v, aN_stg_v, aE_stg_v, aS_stg_v, & 
                                                                aW_stg_v, aB_stg_v, aT_stg_v
            REAL(dp), DIMENSION(:,:,:)   , INTENT(INOUT)     :: b3D_v

            ! residual stuff for ADI
            REAL(dp), INTENT(IN)    :: res_tol_vel
            INTEGER, INTENT(IN)     :: max_it_ADI_vel


            CALL V_COEFF_BLDR(Nx, Ny, Nz, params, del_x, del_y, del_z, & 
                              aP_stg_v, aN_stg_v, aE_stg_v, aS_stg_v, aW_stg_v, aB_stg_v, aT_stg_v, b3D_v, &
                              u_star, v_star, w_star, p_star, &
                              BC_N_v, BC_E_v, BC_S_v, BC_W_v, BC_T_v, BC_B_v, &
                              d_v)
            
            CALL ADI_3D_SOLVR_MAIN(res_tol_vel, max_it_ADI_vel, aP_stg_v(2:Ny,:,:), aN_stg_v(2:Ny,:,:), &
                                   aE_stg_v(2:Ny,:,:), aS_stg_v(2:Ny,:,:), aW_stg_v(2:Ny,:,:), aB_stg_v(2:Ny,:,:), &
                                   aT_stg_v(2:Ny,:,:), v_star(2:Ny,:,:), b3D_v(2:Ny,:,:))

        END SUBROUTINE V_MOMENTUM_SLVR_MAIN

        ! coefficient matrices builder
        SUBROUTINE  V_COEFF_BLDR(Nx, Ny, Nz, params, del_x, del_y, del_z, & 
                                        aP_stg_v, aN_stg_v, aE_stg_v, aS_stg_v, aW_stg_v, aB_stg_v, aT_stg_v, b3D_v, &
                                        u_star, v_star, w_star, p_star, &
                                        BC_N_v, BC_E_v, BC_S_v, BC_W_v, BC_T_v, BC_B_v, &
                                        d_v)

            IMPLICIT NONE

            INTEGER, PARAMETER  :: dp = KIND(1.0D0)

            ! domain stuff
            INTEGER,  INTENT(IN)     :: Nx, Ny, Nz
            REAL(dp), DIMENSION(:),     ALLOCATABLE, INTENT(IN)   :: params
            REAL(dp), INTENT(IN)    :: del_x, del_y, del_z

            ! coeff matrices
            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)     :: aP_stg_v, aN_stg_v, aE_stg_v, aS_stg_v, & 
                                                            aW_stg_v, aB_stg_v, aT_stg_v
            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)     :: b3D_v
            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)     :: d_v ! for pressure correction

            ! boundary conditions
            REAL(dp), INTENT(IN)    :: BC_N_v, BC_E_v, BC_S_v, BC_W_v, BC_T_v, BC_B_v

            ! flow quantities
            REAL(dp), DIMENSION(:,:,:), INTENT(IN)        :: u_star, v_star, w_star ! used for convection term
            REAL(dp), DIMENSION(:,:,:), INTENT(IN)        :: p_star ! pressure 

            ! implicit stuff
            REAL(dp)                :: dt

            !REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(IN)   :: conv_grid_x_stg, conv_grid_y_stg, conv_grid_z_stg

            INTEGER     :: i,j,k
            INTEGER     :: Ny_stg !add 1 to Nx input to form the staggered grid loops

            REAL(dp)    :: diff_coeff, rho, visc
            REAL(dp)    :: del_V, Lx, Ly, Lz ! volume of element and length of domain
            REAL(dp)    :: area_xz, area_xy, area_yz

            REAL(dp)    :: D_e, D_w, D_s, D_n, D_t, D_b, aP_stg_0   ! cst reusab_stgle neighbour coeffs 

            REAL(dp), DIMENSION(Ny+1,Nx,Nz)   :: Su_D_N, Su_D_E, Su_D_S, Su_D_W, Su_D_T, Su_D_B, &   ! boundary contributions from diffusion
                                                Sp_D_N, Sp_D_E, Sp_D_S, Sp_D_W, Sp_D_T, Sp_D_B

            REAL(dp), DIMENSION(Ny+1,Nx,Nz)   :: Su_F_N, Su_F_E, Su_F_S, Su_F_W, Su_F_B, Su_F_T, &   
                                                Sp_F_N, Sp_F_E, Sp_F_S, Sp_F_W, Sp_F_B, Sp_F_T ! boundary contributions from adv

            ! 3D convection array
            REAL(dp),  DIMENSION(Ny+1,Nx,Nz)  :: F_e, F_n, F_w, F_s, F_t, F_b


            Ny_stg = Ny+1

            Lx          = params(4)
            Ly          = params(5)
            Lz          = params(6)

            dt          = params(7)

            rho         = params(10)
            visc        = params(12)

            diff_coeff  = visc

            ! diffusion contributions -
            D_e = diff_coeff/del_x
            D_w = diff_coeff/del_x

            D_n = diff_coeff/del_y
            D_s = diff_coeff/del_y

            D_b = diff_coeff/del_z
            D_t = diff_coeff/del_z

            area_xy = del_x*del_y
            area_xz = del_x*del_z
            area_yz = del_y*del_z

            del_V = del_x*del_y*del_z

            ! implicit unsteady piece contribution
            aP_stg_0 = (rho*del_V)/dt 

            ! 3D convection terms
            DO i = 1,Ny_stg
                DO j = 1,Nx ! start from the second column on the staggered grid for u
                    DO k = 1,Nz

                        IF (i==1 .OR. i==Ny_stg) THEN
                            ! set everything to zero - start from the second column on the staggered grid for u
                            F_w(i,j,k) = 0
                            F_e(i,j,k) = 0
                            F_s(i,j,k) = 0
                            F_n(i,j,k) = 0
                            F_b(i,j,k) = 0
                            F_t(i,j,k) = 0
                        
                        ELSE
                            F_w(i,j,k) = rho*0.5*(u_star(i-1,j,k)     + u_star(i,j,k))
                            F_e(i,j,k) = rho*0.5*(u_star(i-1,j+1,k)   + u_star(i,j+1,k))

                            F_s(i,j,k) = rho*0.5*(v_star(i-1,j,k) + v_star(i,j,k))
                            F_n(i,j,k) = rho*0.5*(v_star(i,j,k)   + v_star(i+1,j,k))

                            F_b(i,j,k) = rho*0.5*(w_star(i-1,j,k)   + w_star(i,j,k))
                            F_t(i,j,k) = rho*0.5*(w_star(i-1,j,k+1) + w_star(i,j,k+1))
                        END IF

                    END DO
                END DO
            END DO

            ! initialize boundary matrices for diffusion with zeros
            DO i = 1,Ny_stg
                DO j = 1,Nx
                    DO k = 1,Nz
                        Su_D_N(i,j,k) = 0
                        Sp_D_N(i,j,k) = 0
                        Su_D_E(i,j,k) = 0
                        Sp_D_E(i,j,k) = 0
                        Su_D_S(i,j,k) = 0
                        Sp_D_S(i,j,k) = 0
                        Su_D_W(i,j,k) = 0
                        Sp_D_W(i,j,k) = 0
                        Su_D_B(i,j,k) = 0
                        Sp_D_B(i,j,k) = 0
                        Su_D_T(i,j,k) = 0
                        Sp_D_T(i,j,k) = 0
                    END DO
                END DO
            END DO

            ! diffusion contributions near boundaries
            DO i = 1,Ny_stg
                DO j = 1,Nx
                    DO k = 1,Nz

                        ! boundary faces
                        IF (i == 2) THEN
                            Su_D_S(i,j,k) = 0
                            Sp_D_S(i,j,k) = -1*visc*area_xz/del_y
                        END IF

                        IF (j == 1) THEN
                            Su_D_W(i,j,k) = 0
                            Sp_D_W(i,j,k) = -2*visc*area_yz/del_x
                        END IF

                        IF (i == Ny) THEN
                            Su_D_N(i,j,k) =  0
                            Sp_D_N(i,j,k) = -1*visc*area_xz/del_y
                        END IF

                        IF (j == Nx) THEN
                            Su_D_E(i,j,k) = 0
                            Sp_D_E(i,j,k) = -2*visc*area_yz/del_x
                        END IF

                        IF (k == 1) THEN
                            Su_D_B(i,j,k) = 0
                            Sp_D_B(i,j,k) = -2*visc*area_xy/del_z
                        END IF

                        IF (k == Nz) THEN
                            Su_D_T(i,j,k) = 0
                            Sp_D_T(i,j,k) = -2*visc*area_xy/del_z
                        END IF

                    END DO
                END DO
            END DO

            ! initialize boundary matrices for advection with zeros
            DO i = 1,Ny_stg
                DO j = 1,Nx
                    DO k = 1,Nz
                        Su_F_N(i,j,k) = 0
                        Sp_F_N(i,j,k) = 0
                        Su_F_E(i,j,k) = 0
                        Sp_F_E(i,j,k) = 0
                        Su_F_S(i,j,k) = 0
                        Sp_F_S(i,j,k) = 0
                        Su_F_W(i,j,k) = 0
                        Sp_F_W(i,j,k) = 0
                        Su_F_B(i,j,k) = 0
                        Sp_F_B(i,j,k) = 0
                        Su_F_T(i,j,k) = 0
                        Sp_F_T(i,j,k) = 0
                    END DO
                END DO
            END DO
        
            ! Build neighbour coefficient mat_stgrices
            DO i = 1,Ny_stg
                DO j = 1,Nx
                    DO k = 1,Nz
                        aE_stg_v(i,j,k) = D_e + F_e(i,j,k)
                        aW_stg_v(i,j,k) = D_w + F_w(i,j,k)
                        aS_stg_v(i,j,k) = D_s + F_s(i,j,k)
                        aN_stg_v(i,j,k) = D_n + F_n(i,j,k)
                        aT_stg_v(i,j,k) = D_t + F_t(i,j,k)
                        aB_stg_v(i,j,k) = D_b + F_b(i,j,k)
                    END DO 
                END DO
            END DO

            ! replace appropriate neighbour coefficients with zeroes at boundaries 
            DO i = 1,Ny_stg
                DO j = 1,Nx
                    DO k = 1,Nz

                        IF     (i == 1 .OR. i == 2) THEN
                            aS_stg_v(i,j,k) = 0
                        END IF

                        IF (i == Ny_stg .OR. i == Ny_stg-1) THEN
                            aN_stg_v(i,j,k) = 0
                        END IF

                        IF (j == 1) THEN
                            aW_stg_v(i,j,k) = 0
                        END IF

                        IF (j == Nx) THEN
                            aE_stg_v(i,j,k) = 0
                        END IF

                        IF (k == 1) THEN
                            aB_stg_v(i,j,k) = 0
                        END IF

                        IF (k == Nz) THEN
                            aT_stg_v(i,j,k) = 0
                        END IF

                    END DO
                END DO
            END DO

            ! Build point P coeff matrix
            DO i = 1,Ny_stg
                DO j = 1,Nx
                    DO k = 1,Nz
                
                        aP_stg_v(i,j,k) =    aN_stg_v(i,j,k) + aE_stg_v(i,j,k) + aS_stg_v(i,j,k) &
                                            + aW_stg_v(i,j,k)+ aB_stg_v(i,j,k) + aT_stg_v(i,j,k) & 
                                            - Sp_D_N(i,j,k) - Sp_D_E(i,j,k) - Sp_D_S(i,j,k) - Sp_D_W(i,j,k) &
                                            - Sp_D_B(i,j,k) - Sp_D_T(i,j,k)  &
                                            - Sp_F_N(i,j,k) - Sp_F_E(i,j,k) - Sp_F_S(i,j,k) &
                                            - Sp_F_W(i,j,k) - Sp_F_B(i,j,k) - Sp_F_T(i,j,k) & ! linearized boundary
                                            + aP_stg_0
                    END DO
                END DO
            END DO
                

            ! Build constant matrix b
            DO i = 1,Ny_stg
                DO j = 1,Nx
                    DO k = 1,Nz

                        IF (i == 1) THEN
                        b3D_v(i,j,k) =     Su_D_N(i,j,k) + Su_D_E(i,j,k) + Su_D_S(i,j,k) &
                                        + Su_D_W(i,j,k) + Su_D_B(i,j,k) + Su_D_T(i,j,k) &
                                        + Su_F_N(i,j,k) + Su_F_E(i,j,k) + Su_F_S(i,j,k) &
                                        + Su_F_W(i,j,k) + Su_F_B(i,j,k) + Su_F_T(i,j,k) &
                                        + aP_stg_0*v_star(i,j,k) &   ! unsteady piece
                                        + (-p_star(i,j,k))*area_xz ! pressure piece
                        
                        ELSEIF (i == Ny_stg) THEN
                        b3D_v(i,j,k) =     Su_D_N(i,j,k) + Su_D_E(i,j,k) + Su_D_S(i,j,k) &
                                        + Su_D_W(i,j,k) + Su_D_B(i,j,k) + Su_D_T(i,j,k) &
                                        + Su_F_N(i,j,k) + Su_F_E(i,j,k) + Su_F_S(i,j,k) &
                                        + Su_F_W(i,j,k) + Su_F_B(i,j,k) + Su_F_T(i,j,k) &
                                        + aP_stg_0*v_star(i,j,k) &   ! unsteady piece
                                        + (p_star(i-1,j,k))*area_xz ! pressure piece
                        
                        ELSE
                        b3D_v(i,j,k) =     Su_D_N(i,j,k) + Su_D_E(i,j,k) + Su_D_S(i,j,k) &
                                        + Su_D_W(i,j,k) + Su_D_B(i,j,k) + Su_D_T(i,j,k) &
                                        + Su_F_N(i,j,k) + Su_F_E(i,j,k) + Su_F_S(i,j,k) &
                                        + Su_F_W(i,j,k) + Su_F_B(i,j,k) + Su_F_T(i,j,k) &
                                        + aP_stg_0*v_star(i,j,k) & ! unsteady piece
                                        + (p_star(i-1,j,k)-p_star(i,j,k))*area_xz ! pressure piece
                        END IF
                    
                    END DO
                END DO
            END DO  

            ! build d coefficients with SIMPLEC formulation to bring to pressure correcter
            DO i = 1,Ny_stg
                DO j = 1,Nx
                    DO k = 1,Nz

                        d_v(i,j,k)      =  area_xz/(aP_stg_v(i,j,k)  - (aN_stg_v(i,j,k) + aE_stg_v(i,j,k) + aS_stg_v(i,j,k) & 
                                                                        + aW_stg_v(i,j,k)+ aB_stg_v(i,j,k) + aT_stg_v(i,j,k)))
                
                    END DO
                END DO
            END DO
            

        END SUBROUTINE V_COEFF_BLDR

END MODULE V_MOMENTUM_SLVR