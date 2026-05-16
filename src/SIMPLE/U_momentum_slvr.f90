MODULE U_MOMENTUM_SLVR

    ! This module adapts the code from the unsteady convection-diffusion energy solver to build the coefficient matrices
    ! for the momentum equations solution in the SIMPLE algorithm

    USE ADI_3D_SOLVR

    IMPLICIT NONE

    CONTAINS

        ! SUBROUTINE U_MOMENTUM_SLVR_MAIN(Nx, Ny, Nz, params, del_x, del_y, del_z, & 
        !                       aP_stg_u, aN_stg_u, aE_stg_u, aS_stg_u, aW_stg_u, aB_stg_u, aT_stg_u, b3D_u, &
        !                       u_star, v_star, w_star, p_star, &
        !                       BC_N_u, BC_E_u, BC_S_u, BC_W_u, BC_T_u, BC_B_u, &
        !                       d_u, &
        !                       res_tol_vel, max_it_ADI_vel)
            
        !     INTEGER, PARAMETER  :: dp = KIND(1.0D0)                  

        !     ! Stuff to keep around after
        !     REAL(dp), DIMENSION(:,:,:)   , INTENT(INOUT)     :: u_star
        !     REAL(dp), DIMENSION(:,:,:)   , INTENT(INOUT)     :: d_u

        !     ! domain stuff
        !     INTEGER,  INTENT(IN)     :: Nx, Ny, Nz
        !     REAL(dp), DIMENSION(:),     ALLOCATABLE, INTENT(IN)   :: params
        !     REAL(dp), INTENT(IN)    :: del_x, del_y, del_z

        !     ! boundary conditions
        !     REAL(dp), INTENT(IN)    :: BC_N_u, BC_E_u, BC_S_u, BC_W_u, BC_T_u, BC_B_u

        !     ! velocities and pressures not updated
        !     REAL(dp), DIMENSION(:,:,:)   , INTENT(IN)     :: v_star, w_star, p_star ! corrrection values
        !     REAL(dp), DIMENSION(:,:,:)   , INTENT(IN)        :: u_0

        !     ! coeff matrices
        !     REAL(dp), DIMENSION(:,:,:)   , INTENT(INOUT)     :: aP_stg_u, aN_stg_u, aE_stg_u, aS_stg_u, & 
        !                                                         aW_stg_u, aB_stg_u, aT_stg_u
        !     REAL(dp), DIMENSION(:,:,:)   , INTENT(INOUT)     :: b3D_u

        !     ! residual stuff for ADI
        !     REAL(dp), INTENT(IN)    :: res_tol_vel
        !     INTEGER, INTENT(IN)     :: max_it_ADI_vel

        !     ! Build coefficient matrices
        !     CALL U_COEFF_BLDR(Nx, Ny, Nz, params, visc, del_x, del_y, del_z, & 
        !                       aP_stg_u, aN_stg_u, aE_stg_u, aS_stg_u, aW_stg_u, aB_stg_u, aT_stg_u, b3D_u, &
        !                       u_star, v_star, w_star, p_star, &
        !                       BC_N_u, BC_E_u, BC_S_u, BC_W_u, BC_T_u, BC_B_u, &
        !                       d_u, u_0)

        !     ! solve for u_star with ADI routine - don't solve at the walls
        !     CALL ADI_3D_SOLVR_MAIN(res_tol_vel, max_it_ADI_vel, aP_stg_u(:,2:Nx,:), aN_stg_u(:,2:Nx,:), &
        !                            aE_stg_u(:,2:Nx,:), aS_stg_u(:,2:Nx,:), aW_stg_u(:,2:Nx,:), aB_stg_u(:,2:Nx,:), &
        !                            aT_stg_u(:,2:Nx,:), u_star(:,2:Nx,:), b3D_u(:,2:Nx,:))

        ! END SUBROUTINE U_MOMENTUM_SLVR_MAIN


        ! get coefficients to solve the u momentum eqn (u star)
        SUBROUTINE U_COEFF_BLDR(Nx, Ny, Nz, params, visc, del_x, del_y, del_z, & 
                                aP_stg_u, aN_stg_u, aE_stg_u, aS_stg_u, aW_stg_u, aB_stg_u, aT_stg_u, b3D_u, &
                                u_star, v_star, w_star, p_star, &
                                BC_N_u, BC_E_u, BC_S_u, BC_W_u, BC_T_u, BC_B_u, &
                                d_u, u_0)

            IMPLICIT NONE

            INTEGER, PARAMETER  :: dp = KIND(1.0D0)

            ! domain stuff
            INTEGER,  INTENT(IN)     :: Nx, Ny, Nz
            REAL(dp), DIMENSION(:),     ALLOCATABLE, INTENT(IN)   :: params
            REAL(dp), INTENT(IN)    :: del_x, del_y, del_z, visc

            ! coeff matrices
            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)     :: aP_stg_u, aN_stg_u, aE_stg_u, aS_stg_u, & 
                                                            aW_stg_u, aB_stg_u, aT_stg_u
            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)     :: b3D_u
            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)     :: d_u ! for pressure correction

            ! boundary conditions
            REAL(dp), INTENT(IN)    :: BC_N_u, BC_E_u, BC_S_u, BC_W_u, BC_T_u, BC_B_u

            ! flow quantities
            REAL(dp), DIMENSION(:,:,:), INTENT(IN)        :: u_star, v_star, w_star ! velocity correction values
            REAL(dp), DIMENSION(:,:,:), INTENT(IN)        :: p_star ! pressure 
            REAL(dp), DIMENSION(:,:,:), INTENT(IN)        :: u_0

            ! implicit stuff
            REAL(dp)                :: dt

            !REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(IN)   :: conv_grid_x_stg, conv_grid_y_stg, conv_grid_z_stg

            INTEGER     :: i,j,k
            INTEGER     :: Nx_stg !add 1 to Nx input to form the staggered grid loops

            REAL(dp)    :: diff_coeff, rho
            REAL(dp)    :: del_V, Lx, Ly, Lz ! volume of element and length of domain
            REAL(dp)    :: area_xz, area_xy, area_yz
            REAL(dp)    ::  F_e_pos,  F_w_pos,  F_s_pos,  F_n_pos

            REAL(dp)    :: D_e, D_w, D_s, D_n, D_t, D_b, aP_stg_0   ! cst reusab_stgle neighbour coeffs 

            REAL(dp), DIMENSION(Ny,Nx+1,Nz)   :: Su_D_N, Su_D_E, Su_D_S, Su_D_W, Su_D_T, Su_D_B, &   ! boundary contributions from diffusion
                                                Sp_D_N, Sp_D_E, Sp_D_S, Sp_D_W, Sp_D_T, Sp_D_B

            REAL(dp), DIMENSION(Ny,Nx+1,Nz)   :: Su_F_N, Su_F_E, Su_F_S, Su_F_W, Su_F_B, Su_F_T, &   
                                                Sp_F_N, Sp_F_E, Sp_F_S, Sp_F_W, Sp_F_B, Sp_F_T ! boundary contributions from adv

            ! logical array for velocity direction at_stg each node
            REAL(dp),  DIMENSION(Ny,Nx+1,Nz)  :: F_e, F_n, F_w, F_s, F_t, F_b

            Nx_stg = Nx+1

            Lx          = params(4)
            Ly          = params(5)
            Lz          = params(6)

            dt          = params(7)

            rho         = params(10)

            diff_coeff  = visc


            D_b = 0
            D_t = 0

            area_xy = del_x*del_y
            area_xz = del_x*del_z
            area_yz = del_y*del_z

            del_V = del_x*del_y*del_z

            ! diffusion contributions
            D_e = (diff_coeff/del_x)*area_yz
            D_w = (diff_coeff/del_x)*area_yz

            D_n = (diff_coeff/del_y)*area_xz
            D_s = (diff_coeff/del_y)*area_xz


            ! implicit unsteady piece contribution
            aP_stg_0 = (rho*del_V)/dt 

            ! 3D convection terms
            DO i = 1,Ny
                DO j = 1,Nx_stg ! start from the second column on the staggered grid for u
                    DO k = 1,Nz

                        IF (j==1 .OR. j==Nx_stg) THEN
                            ! set everything to zero - start from the second column on the staggered grid for u
                            F_w(i,j,k) = 0
                            F_e(i,j,k) = 0
                            F_s(i,j,k) = 0
                            F_n(i,j,k) = 0
                            F_b(i,j,k) = 0
                            F_t(i,j,k) = 0
                        
                        ELSE
                            F_w(i,j,k) = rho*0.5*(u_star(i,j-1,k) + u_star(i,j,k))*area_yz
                            F_e(i,j,k) = rho*0.5*(u_star(i,j,k)   + u_star(i,j+1,k))*area_yz

                            F_s(i,j,k) = rho*0.5*(v_star(i,j-1,k)   + v_star(i,j,k))*area_xz
                            F_n(i,j,k) = rho*0.5*(v_star(i+1,j-1,k) + v_star(i+1,j,k))*area_xz

                            ! F_b(i,j,k) = rho*0.5*(w_star(i,j-1,k)   + w_star(i,j,k))
                            ! F_t(i,j,k) = rho*0.5*(w_star(i,j-1,k+1) + w_star(i,j,k+1))
                            F_b(i,j,k) = 0
                            F_t(i,j,k) = 0
                        END IF

                    END DO
                END DO
            END DO

            ! initialize boundary matrices for diffusion with zeros
            DO i = 1,Ny
                DO j = 1,Nx_stg
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
            DO i = 1,Ny
                DO j = 1,Nx_stg
                    DO k = 1,Nz

                        ! boundary faces
                        IF (i == 1) THEN
                            Su_D_S(i,j,k) = 0
                            Sp_D_S(i,j,k) = -2*visc*area_xz/del_y
                            !Sp_D_S(i,j,k) = -2*visc/del_y
                        END IF

                        IF (j == 2) THEN
                            Su_D_W(i,j,k) = 0
                            Sp_D_W(i,j,k) = -visc*area_yz/del_x
                            !Sp_D_W(i,j,k) = -visc/del_x
                        END IF

                        IF (i == Ny) THEN
                            !Su_D_N(i,j,k) =  2*visc*area_xz*BC_N_u/del_y ! contribution from the slip wall
                            Su_D_N(i,j,k) =  2*visc*BC_N_u*area_xz/del_y ! contribution from the slip wall
                            Sp_D_N(i,j,k) = -2*visc*area_xz/del_y
                            !Sp_D_N(i,j,k) = -2*visc/del_y
                        END IF

                        IF (j == Nx) THEN
                            Su_D_E(i,j,k) = 0
                            Sp_D_E(i,j,k) = -visc*area_yz/del_x
                            !Sp_D_E(i,j,k) = -visc/del_x
                        END IF

                        ! IF (k == 1) THEN
                        !     Su_D_B(i,j,k) = 0
                        !     Sp_D_B(i,j,k) = -2*visc*area_xy/del_z
                        ! END IF

                        ! IF (k == Nz) THEN
                        !     Su_D_T(i,j,k) = 0
                        !     Sp_D_T(i,j,k) = -2*visc*area_xy/del_z
                        ! END IF

                    END DO
                END DO
            END DO

            ! initialize boundary matrices for advection with zeros
            DO i = 1,Ny
                DO j = 1,Nx_stg
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
            DO i = 1,Ny
                DO j = 1,Nx_stg
                    DO k = 1,Nz

                        ! IF (F_w(i,j,k) .GE. 0) THEN
                        !     F_e(i,j,k) = 0
                        ! ELSE 
                        !     F_w(i,j,k) = 0
                        ! END IF

                        ! IF (F_s(i,j,k) .GE. 0) THEN
                        !     F_n(i,j,k) = 0
                        ! ELSE 
                        !     F_s(i,j,k) = 0
                        ! END IF

                        F_e_pos = MAX(-F_e(i,j,k), 0.0_dp)
                        F_w_pos = MAX( F_w(i,j,k), 0.0_dp)

                        F_n_pos = MAX(-F_n(i,j,k), 0.0_dp)
                        F_s_pos = MAX( F_s(i,j,k), 0.0_dp)

                        aE_stg_u(i,j,k) = D_e + F_e_pos
                        aW_stg_u(i,j,k) = D_w + F_w_pos
                        aS_stg_u(i,j,k) = D_s + F_s_pos
                        aN_stg_u(i,j,k) = D_n + F_n_pos
                        ! aT_stg_u(i,j,k) = D_t - F_t(i,j,k)
                        ! aB_stg_u(i,j,k) = D_b + F_b(i,j,k)
                    END DO 
                END DO
            END DO

            ! replace appropriate neighbour coefficients with zeroes at boundaries 
            DO i = 1,Ny
                DO j = 1,Nx_stg
                    DO k = 1,Nz

                        IF     (i == 1) THEN
                            aS_stg_u(i,j,k) = 0
                        END IF

                        IF (i == Ny) THEN
                            aN_stg_u(i,j,k) = 0
                        END IF

                        IF (j == 1 .OR. j ==2) THEN
                            aW_stg_u(i,j,k) = 0
                        END IF

                        IF (j == Nx_stg .OR. j == Nx_stg-1) THEN
                            aE_stg_u(i,j,k) = 0
                        END IF

                        IF (k == 1) THEN
                            aB_stg_u(i,j,k) = 0
                        END IF

                        IF (k == Nz) THEN
                            aT_stg_u(i,j,k) = 0
                        END IF

                    END DO
                END DO
            END DO

            ! Build point P coeff matrix
            DO i = 1,Ny
                DO j = 1,Nx_stg
                    DO k = 1,Nz
                
                        aP_stg_u(i,j,k) =    aN_stg_u(i,j,k) + aE_stg_u(i,j,k) + aS_stg_u(i,j,k) + & 
                                             aW_stg_u(i,j,k) &
                                             !+ aB_stg_u(i,j,k) + aT_stg_u(i,j,k) & 
                                            - Sp_D_N(i,j,k) - Sp_D_E(i,j,k) - Sp_D_S(i,j,k) - Sp_D_W(i,j,k) &
                                            !- Sp_D_B(i,j,k) - Sp_D_T(i,j,k)  &
                                            - Sp_F_N(i,j,k) - Sp_F_E(i,j,k) - Sp_F_S(i,j,k) - Sp_F_W(i,j,k) &
                                            - F_w(i,j,k) - F_s(i,j,k) + F_n(i,j,k) + F_e(i,j,k) &
                                            !- Sp_F_B(i,j,k) - Sp_F_T(i,j,k) & ! linearized boundary
                                            + aP_stg_0
                    END DO
                END DO
            END DO
                

            ! Build constant matrix b
            DO i = 1,Ny
                DO j = 1,Nx_stg
                    DO k = 1,Nz
                        ! ADD PRESSURE TIMES AREA TERM TO b

                        IF (j == 1) THEN
                        b3D_u(i,j,k) =   Su_D_N(i,j,k) + Su_D_E(i,j,k) + Su_D_S(i,j,k) &
                                        + Su_D_W(i,j,k) &
                                        !+ Su_D_B(i,j,k) + Su_D_T(i,j,k) &
                                        + Su_F_N(i,j,k) + Su_F_E(i,j,k) + Su_F_S(i,j,k) &
                                        + Su_F_W(i,j,k) &
                                        !+ Su_F_B(i,j,k) + Su_F_T(i,j,k) &
                                        + aP_stg_0*u_0(i,j,k)! &  ! unsteady piece
                                        !+ (-p_star(i,j,k))*area_yz ! pressure piece
                        
                        ELSEIF (j == Nx_stg) THEN
                        b3D_u(i,j,k) =   Su_D_N(i,j,k) + Su_D_E(i,j,k) + Su_D_S(i,j,k) &
                                        + Su_D_W(i,j,k) + Su_D_B(i,j,k) + Su_D_T(i,j,k) &
                                        + Su_F_N(i,j,k) + Su_F_E(i,j,k) + Su_F_S(i,j,k) &
                                        + Su_F_W(i,j,k) + Su_F_B(i,j,k) + Su_F_T(i,j,k) &
                                        + aP_stg_0*u_0(i,j,k)  !&  ! unsteady piece
                                        !+ (p_star(i,j-1,k))*area_yz ! pressure piece
                        
                        ELSE
                        b3D_u(i,j,k) =   Su_D_N(i,j,k) + Su_D_E(i,j,k) + Su_D_S(i,j,k) &
                                        + Su_D_W(i,j,k) &
                                        !+ Su_D_B(i,j,k) + Su_D_T(i,j,k) &
                                        + Su_F_N(i,j,k) + Su_F_E(i,j,k) + Su_F_S(i,j,k) &
                                        + Su_F_W(i,j,k) &
                                        !+ Su_F_B(i,j,k) + Su_F_T(i,j,k) &
                                        + aP_stg_0*u_0(i,j,k) &! unsteady piece
                                        + (p_star(i,j-1,k)-p_star(i,j,k))*area_yz ! pressure piece
                        END IF
                    
                    END DO
                END DO
            END DO  

            ! build d coefficients with SIMPLEC formulation to bring to pressure correcter
            DO i = 1,Ny
                DO j = 1,Nx_stg
                    DO k = 1,Nz

                        d_u(i,j,k)      =  area_yz/(aP_stg_u(i,j,k) - (aN_stg_u(i,j,k) + aE_stg_u(i,j,k) + aS_stg_u(i,j,k) & 
                                                                        + aW_stg_u(i,j,k)))
                                                                       !+ aB_stg_u(i,j,k) + aT_stg_u(i,j,k)))
                
                    END DO
                END DO
            END DO

        END SUBROUTINE U_COEFF_BLDR

        ! get coefficients to solve the u momentum eqn with an internal geometry
        SUBROUTINE U_COEFF_BLDR_INT_GEOMETRY(Nx, Ny, Nz, params, visc, del_x, del_y, del_z, & 
                                aP_stg_u, aN_stg_u, aE_stg_u, aS_stg_u, aW_stg_u, aB_stg_u, aT_stg_u, b3D_u, &
                                u_star, v_star, w_star, p_star, &
                                BC_N_u, BC_E_u, BC_S_u, BC_W_u, BC_T_u, BC_B_u, &
                                d_u, u_0)

            IMPLICIT NONE

            INTEGER, PARAMETER  :: dp = KIND(1.0D0)

            ! domain stuff
            INTEGER,  INTENT(IN)     :: Nx, Ny, Nz
            REAL(dp), DIMENSION(:),     ALLOCATABLE, INTENT(IN)   :: params
            REAL(dp), INTENT(IN)    :: del_x, del_y, del_z, visc

            ! coeff matrices
            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)     :: aP_stg_u, aN_stg_u, aE_stg_u, aS_stg_u, & 
                                                            aW_stg_u, aB_stg_u, aT_stg_u
            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)     :: b3D_u
            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)     :: d_u ! for pressure correction

            ! boundary conditions
            REAL(dp), INTENT(IN)    :: BC_N_u, BC_E_u, BC_S_u, BC_W_u, BC_T_u, BC_B_u

            ! flow quantities
            REAL(dp), DIMENSION(:,:,:), INTENT(IN)        :: u_star, v_star, w_star ! velocity correction values
            REAL(dp), DIMENSION(:,:,:), INTENT(IN)        :: p_star ! pressure 
            REAL(dp), DIMENSION(:,:,:), INTENT(IN)        :: u_0

            ! implicit stuff
            REAL(dp)                :: dt

            !REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(IN)   :: conv_grid_x_stg, conv_grid_y_stg, conv_grid_z_stg

            INTEGER     :: i,j,k
            INTEGER     :: Nx_stg !add 1 to Nx input to form the staggered grid loops

            REAL(dp)    :: diff_coeff, rho
            REAL(dp)    :: del_V, Lx, Ly, Lz ! volume of element and length of domain
            REAL(dp)    :: area_xz, area_xy, area_yz
            REAL(dp)    ::  F_e_pos,  F_w_pos,  F_s_pos,  F_n_pos

            REAL(dp)    :: D_e, D_w, D_s, D_n, D_t, D_b, aP_stg_0   ! cst reusab_stgle neighbour coeffs 

            REAL(dp), DIMENSION(Ny,Nx+1,Nz)   :: Su_D_N, Su_D_E, Su_D_S, Su_D_W, Su_D_T, Su_D_B, &   ! boundary contributions from diffusion
                                                Sp_D_N, Sp_D_E, Sp_D_S, Sp_D_W, Sp_D_T, Sp_D_B

            REAL(dp), DIMENSION(Ny,Nx+1,Nz)   :: Su_F_N, Su_F_E, Su_F_S, Su_F_W, Su_F_B, Su_F_T, &   
                                                Sp_F_N, Sp_F_E, Sp_F_S, Sp_F_W, Sp_F_B, Sp_F_T ! boundary contributions from adv

            ! logical array for velocity direction at_stg each node
            REAL(dp),  DIMENSION(Ny,Nx+1,Nz)  :: F_e, F_n, F_w, F_s, F_t, F_b

            INTEGER :: Nx_int_start, Nx_int_end, Ny_int_start, Ny_int_end

            Nx_stg = Nx+1

            Lx          = params(4)
            Ly          = params(5)
            Lz          = params(6)

            dt          = params(7)

            rho         = params(10)

            diff_coeff  = visc


            D_b = 0
            D_t = 0

            area_xy = del_x*del_y
            area_xz = del_x*del_z
            area_yz = del_y*del_z

            del_V = del_x*del_y*del_z

            ! diffusion contributions
            D_e = (diff_coeff/del_x)*area_yz
            D_w = (diff_coeff/del_x)*area_yz

            D_n = (diff_coeff/del_y)*area_xz
            D_s = (diff_coeff/del_y)*area_xz


            ! implicit unsteady piece contribution
            aP_stg_0 = (rho*del_V)/dt 

            ! points where there is an internal box
            Nx_int_start=params(37)
            Nx_int_end=params(38)
            Ny_int_start=params(39)
            Ny_int_end=params(40)

            ! 3D convection terms
            DO i = 1,Ny
                DO j = 1,Nx_stg ! start from the second column on the staggered grid for u
                    DO k = 1,Nz

                        IF (j==1 .OR. j==Nx_stg) THEN
                            ! set everything to zero - start from the second column on the staggered grid for u
                            F_w(i,j,k) = 0
                            F_e(i,j,k) = 0
                            F_s(i,j,k) = 0
                            F_n(i,j,k) = 0
                            F_b(i,j,k) = 0
                            F_t(i,j,k) = 0
                        
                        ELSE
                            F_w(i,j,k) = rho*0.5*(u_star(i,j-1,k) + u_star(i,j,k))*area_yz
                            F_e(i,j,k) = rho*0.5*(u_star(i,j,k)   + u_star(i,j+1,k))*area_yz

                            F_s(i,j,k) = rho*0.5*(v_star(i,j-1,k)   + v_star(i,j,k))*area_xz
                            F_n(i,j,k) = rho*0.5*(v_star(i+1,j-1,k) + v_star(i+1,j,k))*area_xz

                            ! F_b(i,j,k) = rho*0.5*(w_star(i,j-1,k)   + w_star(i,j,k))
                            ! F_t(i,j,k) = rho*0.5*(w_star(i,j-1,k+1) + w_star(i,j,k+1))
                            F_b(i,j,k) = 0
                            F_t(i,j,k) = 0
                        END IF

                    END DO
                END DO
            END DO

            ! initialize boundary matrices for diffusion with zeros
            DO i = 1,Ny
                DO j = 1,Nx_stg
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
            DO i = 1,Ny
                DO j = 1,Nx_stg
                    DO k = 1,Nz

                        ! boundary faces
                        IF (i == 1) THEN
                            Su_D_S(i,j,k) = 0
                            Sp_D_S(i,j,k) = -2*visc*area_xz/del_y
                            !Sp_D_S(i,j,k) = -2*visc/del_y
                        END IF

                        IF (j == 2) THEN
                            Su_D_W(i,j,k) = 0
                            Sp_D_W(i,j,k) = -visc*area_yz/del_x
                            !Sp_D_W(i,j,k) = -visc/del_x
                        END IF

                        IF (i == Ny) THEN
                            !Su_D_N(i,j,k) =  2*visc*area_xz*BC_N_u/del_y ! contribution from the slip wall
                            Su_D_N(i,j,k) =  2*visc*BC_N_u*area_xz/del_y ! contribution from the slip wall
                            Sp_D_N(i,j,k) = -2*visc*area_xz/del_y
                            !Sp_D_N(i,j,k) = -2*visc/del_y
                        END IF

                        IF (j == Nx) THEN
                            Su_D_E(i,j,k) = 0
                            Sp_D_E(i,j,k) = -visc*area_yz/del_x
                            !Sp_D_E(i,j,k) = -visc/del_x
                        END IF
                        
                        ! internal boundary stuff
                        IF ((i == Ny_int_start) .AND. (Nx_int_start .LE. j) .AND. (j .LE. Nx_int_end)) THEN
                            Sp_D_N(i-1,j,k) = -visc*area_yz/del_x
                        END IF

                        IF (j == Nx_int_start .AND. (Ny_int_start .LE. i) .AND. (i .LE. Ny_int_end)) THEN
                            Sp_D_E(i,j-1,k) = -visc*area_yz/del_x
                        END IF

                        IF (i == Ny_int_end  .AND. (Nx_int_start .LE. j) .AND. (j .LE. Nx_int_end)) THEN
                            Sp_D_S(i+1,j,k) = -visc*area_yz/del_y
                        END IF

                        IF (j == Nx_int_end .AND. (Ny_int_start .LE. i) .AND. (i .LE. Ny_int_end)) THEN
                            Sp_D_W(i,j+1,k) = -visc*area_yz/del_y
                        END IF

                        ! corners
                        IF (j == Nx_int_start .AND. i == Ny_int_start) THEN
                            Sp_D_E(i-1,j-1,k) = -visc*area_yz/del_x
                            Sp_D_N(i-1,j-1,k) = -visc*area_yz/del_y
                        END IF

                        IF (j == Nx_int_start .AND. i == Ny_int_end) THEN
                            Sp_D_E(i+1,j-1,k) = -visc*area_yz/del_x
                            Sp_D_S(i+1,j-1,k) = -visc*area_yz/del_y
                        END IF

                        IF (j == Nx_int_end .AND. i == Ny_int_start) THEN
                            Sp_D_W(i-1,j+1,k) = -visc*area_yz/del_x
                            Sp_D_N(i-1,j+1,k) = -visc*area_yz/del_y
                        END IF

                        IF (j == Nx_int_end .AND. i == Ny_int_end) THEN
                            Sp_D_W(i+1,j+1,k) = -visc*area_yz/del_x
                            Sp_D_S(i+1,j+1,k) = -visc*area_yz/del_y
                        END IF

                    END DO
                END DO
            END DO

            ! initialize boundary matrices for advection with zeros
            DO i = 1,Ny
                DO j = 1,Nx_stg
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
            DO i = 1,Ny
                DO j = 1,Nx_stg
                    DO k = 1,Nz

                        ! IF (F_w(i,j,k) .GE. 0) THEN
                        !     F_e(i,j,k) = 0
                        ! ELSE 
                        !     F_w(i,j,k) = 0
                        ! END IF

                        ! IF (F_s(i,j,k) .GE. 0) THEN
                        !     F_n(i,j,k) = 0
                        ! ELSE 
                        !     F_s(i,j,k) = 0
                        ! END IF

                        F_e_pos = MAX(-F_e(i,j,k), 0.0_dp)
                        F_w_pos = MAX( F_w(i,j,k), 0.0_dp)

                        F_n_pos = MAX(-F_n(i,j,k), 0.0_dp)
                        F_s_pos = MAX( F_s(i,j,k), 0.0_dp)

                        aE_stg_u(i,j,k) = D_e + F_e_pos
                        aW_stg_u(i,j,k) = D_w + F_w_pos
                        aS_stg_u(i,j,k) = D_s + F_s_pos
                        aN_stg_u(i,j,k) = D_n + F_n_pos
                        ! aT_stg_u(i,j,k) = D_t - F_t(i,j,k)
                        ! aB_stg_u(i,j,k) = D_b + F_b(i,j,k)
                    END DO 
                END DO
            END DO

            ! replace appropriate neighbour coefficients with zeroes at boundaries 
            DO i = 1,Ny
                DO j = 1,Nx_stg
                    DO k = 1,Nz

                        IF     (i == 1) THEN
                            aS_stg_u(i,j,k) = 0
                        END IF

                        IF (i == Ny) THEN
                            aN_stg_u(i,j,k) = 0
                        END IF

                        IF (j == 1 .OR. j ==2) THEN
                            aW_stg_u(i,j,k) = 0
                        END IF

                        IF (j == Nx_stg .OR. j == Nx_stg-1) THEN
                            aE_stg_u(i,j,k) = 0
                        END IF

                        IF (k == 1) THEN
                            aB_stg_u(i,j,k) = 0
                        END IF

                        IF (k == Nz) THEN
                            aT_stg_u(i,j,k) = 0
                        END IF

                    END DO
                END DO
            END DO

            ! Build point P coeff matrix
            DO i = 1,Ny
                DO j = 1,Nx_stg
                    DO k = 1,Nz
                
                        aP_stg_u(i,j,k) =    aN_stg_u(i,j,k) + aE_stg_u(i,j,k) + aS_stg_u(i,j,k) + & 
                                             aW_stg_u(i,j,k) &
                                             !+ aB_stg_u(i,j,k) + aT_stg_u(i,j,k) & 
                                            - Sp_D_N(i,j,k) - Sp_D_E(i,j,k) - Sp_D_S(i,j,k) - Sp_D_W(i,j,k) &
                                            !- Sp_D_B(i,j,k) - Sp_D_T(i,j,k)  &
                                            - Sp_F_N(i,j,k) - Sp_F_E(i,j,k) - Sp_F_S(i,j,k) - Sp_F_W(i,j,k) &
                                            - F_w(i,j,k) - F_s(i,j,k) + F_n(i,j,k) + F_e(i,j,k) &
                                            !- Sp_F_B(i,j,k) - Sp_F_T(i,j,k) & ! linearized boundary
                                            + aP_stg_0
                    END DO
                END DO
            END DO
                

            ! Build constant matrix b
            DO i = 1,Ny
                DO j = 1,Nx_stg
                    DO k = 1,Nz
                        ! ADD PRESSURE TIMES AREA TERM TO b

                        IF (j == 1) THEN
                        b3D_u(i,j,k) =   Su_D_N(i,j,k) + Su_D_E(i,j,k) + Su_D_S(i,j,k) &
                                        + Su_D_W(i,j,k) &
                                        !+ Su_D_B(i,j,k) + Su_D_T(i,j,k) &
                                        + Su_F_N(i,j,k) + Su_F_E(i,j,k) + Su_F_S(i,j,k) &
                                        + Su_F_W(i,j,k) &
                                        !+ Su_F_B(i,j,k) + Su_F_T(i,j,k) &
                                        + aP_stg_0*u_0(i,j,k)! &  ! unsteady piece
                                        !+ (-p_star(i,j,k))*area_yz ! pressure piece
                        
                        ELSEIF (j == Nx_stg) THEN
                        b3D_u(i,j,k) =   Su_D_N(i,j,k) + Su_D_E(i,j,k) + Su_D_S(i,j,k) &
                                        + Su_D_W(i,j,k) + Su_D_B(i,j,k) + Su_D_T(i,j,k) &
                                        + Su_F_N(i,j,k) + Su_F_E(i,j,k) + Su_F_S(i,j,k) &
                                        + Su_F_W(i,j,k) + Su_F_B(i,j,k) + Su_F_T(i,j,k) &
                                        + aP_stg_0*u_0(i,j,k)  !&  ! unsteady piece
                                        !+ (p_star(i,j-1,k))*area_yz ! pressure piece
                        
                        ELSE
                        b3D_u(i,j,k) =   Su_D_N(i,j,k) + Su_D_E(i,j,k) + Su_D_S(i,j,k) &
                                        + Su_D_W(i,j,k) &
                                        !+ Su_D_B(i,j,k) + Su_D_T(i,j,k) &
                                        + Su_F_N(i,j,k) + Su_F_E(i,j,k) + Su_F_S(i,j,k) &
                                        + Su_F_W(i,j,k) &
                                        !+ Su_F_B(i,j,k) + Su_F_T(i,j,k) &
                                        + aP_stg_0*u_0(i,j,k) &! unsteady piece
                                        + (p_star(i,j-1,k)-p_star(i,j,k))*area_yz ! pressure piece
                        END IF
                    
                    END DO
                END DO
            END DO  

            ! build d coefficients with SIMPLEC formulation to bring to pressure correcter
            DO i = 1,Ny
                DO j = 1,Nx_stg
                    DO k = 1,Nz

                        d_u(i,j,k)      =  area_yz/(aP_stg_u(i,j,k) - (aN_stg_u(i,j,k) + aE_stg_u(i,j,k) + aS_stg_u(i,j,k) & 
                                                                        + aW_stg_u(i,j,k)))
                                                                       !+ aB_stg_u(i,j,k) + aT_stg_u(i,j,k)))
                
                    END DO
                END DO
            END DO

        END SUBROUTINE U_COEFF_BLDR_INT_GEOMETRY

END MODULE U_MOMENTUM_SLVR