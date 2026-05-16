MODULE ENERGY_SLVR

    IMPLICIT NONE

    CONTAINS 
    ! advection-diffusion coefficient matrix builder
    SUBROUTINE ENERGY_COEFFBLDR(Nx, Ny, Nz, params, u, v, w, & 
                                    aP_ndl, aN_ndl, aE_ndl, aS_ndl, aW_ndl, &
                                    aB_ndl, aT_ndl, b3D_energy, phi3D_energy)

        IMPLICIT NONE

        INTEGER, PARAMETER  :: dp = KIND(1.0D0)

        !REAL(dp), PARAMETER  :: coeff_test = 0.5

        INTEGER,  INTENT(IN)     :: Nx, Ny, Nz
        REAL(dp), DIMENSION (:,:,:), INTENT(IN)     :: u, v, w

        REAL(dp), DIMENSION(:),     ALLOCATABLE, INTENT(IN)   :: params

        REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)     :: aP_ndl, aN_ndl, aE_ndl, aS_ndl, aW_ndl, aB_ndl, aT_ndl
        REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)     :: b3D_energy
        REAL(dp), DIMENSION(:,:,:)   , INTENT(INOUT)   :: phi3D_energy

        ! IMPLICIT
        INTEGER     :: i,j,k

        REAL(dp)                :: dt

        REAL(dp)    :: diff_coeff, gamma ! alpha and specific heat

        REAL(dp)    :: del_x, del_y, del_z, del_V, Lx, Ly, Lz, alpha, cp, rho  ! params for building coeff matrices

        REAL(dp)    :: area_xz, area_xy, area_yz

        REAL(dp)    :: bN, bE, bS, bW, bB, bT   ! modified boundary conditons from BC input

        REAL(dp)    :: BC_N, BC_E, BC_S, BC_W, BC_B, BC_T   ! modified boundary conditons from BC input

        REAL(dp)    :: flux_N, flux_E, flux_S, flux_W, flux_B, flux_T   ! modified boundary conditons from BC input

        REAL(dp)    :: aX_D, aY_D, aZ_D, aP_0   ! cst reusable neighbour coeffs 

        REAL(dp), DIMENSION(Ny, Nx, Nz)    :: F_e, F_s, F_n, F_w, F_b, F_t

        REAL(dp), DIMENSION(Ny,Nx,Nz)   :: Su_D_N, Su_D_E, Su_D_S, Su_D_W, Su_D_T, Su_D_B, &   ! boundary contributions from diffusion
                                        Sp_D_N, Sp_D_E, Sp_D_S, Sp_D_W, Sp_D_T, Sp_D_B

        ! logical array for velocity direction at each node
        LOGICAL,  DIMENSION(Ny,Nx,Nz)   :: uwind_N, uwind_E, uwind_S, uwind_W, uwind_B, uwind_T

        REAL(dp), DIMENSION(Ny,Nx,Nz)   :: Su_F_N, Su_F_E, Su_F_S, Su_F_W, Su_F_B, Su_F_T, &   
                                        Sp_F_N, Sp_F_E, Sp_F_S, Sp_F_W, Sp_F_B, Sp_F_T ! boundary contributions from adv

        Lx          = params(4)
        Ly          = params(5)
        Lz          = params(6)

        dt          = params(7)

        alpha       = params(9)
        rho         = params(10)
        cp          = params(11)

        BC_N        = params(13)
        BC_E        = params(14)
        BC_S        = params(15)
        BC_W        = params(16)
        BC_T        = params(17)
        BC_B        = params(18)

        flux_N      = params(19)
        flux_E      = params(20)
        flux_S      = params(21)
        flux_W      = params(22)
        flux_T      = params(23)
        flux_B      = params(24)

        diff_coeff  = alpha*rho ! k /cp
        gamma = 1.4

        ! convert boundary conditions (temp deg C) to energy values (J/kg)
        CALL BC_CALC(BC_N, BC_E, BC_S, BC_W, BC_T, BC_B, cp, bN, bE, bS, bW, bT, bB)

        ! neighbours in x have cst coefficient diffusion -
        del_x = Lx/Nx
        aX_D = (diff_coeff/del_x)

        del_y = Ly/Ny
        del_z = Lz/Nz

        area_xy = del_x*del_y
        area_xz = del_x*del_z
        area_yz = del_y*del_z

        ! reusable constant  diffusion contributions
        aX_D = (diff_coeff/del_x)*area_yz
        aY_D = (diff_coeff/del_y)*area_xz
        aZ_D = (diff_coeff/del_z)*area_xy

        del_V = del_x*del_y*del_z

        ! implicit unsteady piece contribution
        aP_0 = ((rho/gamma)*del_V)/dt

        ! convection terms
        DO i = 1,Ny
            DO j = 1,Nx
                DO k = 1,Nz

                    F_w(i,j,k) = rho*u(i,j,k)*area_yz
                    F_e(i,j,k) = rho*u(i,j+1,k)*area_yz

                    F_s(i,j,k) = rho*v(i,j,k)*area_xz
                    F_n(i,j,k) = rho*v(i+1,j,k)*area_xz

                    ! F_b(i,j,k) = rho*0.5*(w_star(i-1,j,k)   + w_star(i,j,k))
                    ! F_t(i,j,k) = rho*0.5*(w_star(i-1,j,k+1) + w_star(i,j,k+1))

                    F_b(i,j,k) = 0
                    F_t(i,j,k) = 0
                END DO
            END DO
        END DO

        ! '' boolean array for upwind differencing
        DO i = 1,Ny
            DO j = 1,Nx
                DO k = 1,Nz

                    IF (u(i,j,k) .GE. 0) THEN
                        uwind_W(i,j,k) = .TRUE.
                    ELSE
                        uwind_W(i,j,k) = .FALSE.
                    END IF

                    IF (u(i,j+1,k) .LE. 0) THEN
                        uwind_E(i,j,k) = .TRUE.
                    ELSE
                        uwind_E(i,j,k) = .FALSE.
                    END IF

                    IF (v(i,j,k) .GE. 0) THEN
                        uwind_S(i,j,k) = .TRUE.
                    ELSE
                        uwind_S(i,j,k) = .FALSE.
                    END IF

                    IF (v(i+1,j,k) .LE. 0) THEN
                        uwind_N(i,j,k) = .TRUE.
                    ELSE
                        uwind_N(i,j,k) = .FALSE.
                    END IF

                    IF (w(i,j,k) .GE. 0) THEN
                        uwind_B(i,j,k) = .TRUE.
                        uwind_T(i,j,k) = .FALSE.
                    ELSE
                        uwind_B(i,j,k) = .FALSE.
                        uwind_T(i,j,k) = .TRUE.
                    END IF
                END DO
            END DO
        END DO

        ! initialize boundary matrices for diffusion with zeros
        DO i = 1,Ny
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

        ! update boundary contributions from diffusion
        DO i = 1,Ny
            DO j = 1,Nx
                DO k = 1,Nz

                    ! boundary faces
                    IF (i == 1) THEN
                        Su_D_S(i,j,k) = flux_S*area_xz
                        Sp_D_S(i,j,k) = 0
                    END IF

                    IF (j == 1) THEN
                        Su_D_W(i,j,k) = flux_W*area_yz
                        Sp_D_W(i,j,k) = 0
                    END IF

                    IF (i == Ny) THEN
                        Su_D_N(i,j,k) = 2*aY_D*bN*area_xz
                        Sp_D_N(i,j,k) = -2*aY_D*area_xz
                    END IF

                    IF (j == Nx) THEN
                        Su_D_E(i,j,k) = flux_E*area_yz
                        Sp_D_E(i,j,k) = 0
                    END IF

                    ! IF (k == 1) THEN
                    !     Su_D_B(i,j,k) = 2*aZ_D*bB
                    !     Sp_D_B(i,j,k) = -2*aZ_D
                    ! END IF

                    ! IF (k == Nz) THEN
                    !     Su_D_T(i,j,k) = 2*aZ_D*bT
                    !     Sp_D_T(i,j,k) = -2*aZ_D
                    ! END IF

                END DO
            END DO
        END DO

        ! initialize boundary matrices for advection with zeros
        DO i = 1,Ny
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

        ! update boundary contributions from advection based on upwind scheme
        DO i = 1,Ny
            DO j = 1,Nx
                DO k = 1,Nz

                    ! boundary faces
                    IF     (i == 1) THEN
                        ! advection contribution from boundary if upwind cell is adjacent
                        IF      (uwind_S(i,j,k)) THEN
                            Su_F_S(i,j,k) = F_s(i,j,k)*bS*area_xz
                            Sp_F_S(i,j,k) = -F_s(i,j,k)*area_xz
                        END IF
                    END IF
                
                    IF (j == 1) THEN
                        IF     (uwind_W(i,j,k)) THEN
                            Su_F_W(i,j,k) = -F_w(i,j,k)*bW*area_yz
                            Sp_F_W(i,j,k) = F_w(i,j,k)*area_yz
                        END IF
                    END IF
                
                    IF (i == Ny) THEN
                        IF     (uwind_N(i,j,k)) THEN
                            Su_F_N(i,j,k) = -F_n(i,j,k)*bN*area_xz
                            Sp_F_N(i,j,k) = F_n(i,j,k)*area_xz
                        END IF
                    END IF

                    IF (j == Nx) THEN
                        IF     (uwind_E(i,j,k)) THEN
                            Su_F_E(i,j,k) = F_e(i,j,k)*bE*area_yz
                            Sp_F_E(i,j,k) = -F_e(i,j,k)*area_yz
                        END IF
                    END IF
                    
                    ! ELSEIF (k == 1) THEN
                    !     IF     (uwind_B(i,j,k)) THEN
                    !         Su_F_B(i,j,k) = aZ_F*bB
                    !         Sp_F_B(i,j,k) = -aZ_F
                    !     END IF

                    ! ELSEIF (k == Nz) THEN
                    !     IF     (uwind_T(i,j,k)) THEN
                    !         Su_F_T(i,j,k) = aZ_F*bT
                    !         Sp_F_T(i,j,k) = -aZ_F
                    !     END IF

                    ! END IF
                END DO
            END DO
        END DO

        ! Update source terms near outflow or adiabatic boundaries
        ! CALL   OUTFLW_AD(Su_F_S, Sp_F_S, Su_D_S, Sp_D_S, Su_F_W, Sp_F_W, Su_D_W, Sp_D_W, Su_F_N, Sp_F_N, Su_D_N, Sp_D_N, &
        !                 Su_F_E, Sp_F_E, Su_D_E, Sp_D_E, Su_F_B, Sp_F_B, Su_D_B, Sp_D_B, Su_F_T, Sp_F_T, Su_D_T, Sp_D_T, &
        !                 otflw_ad_S, otflw_ad_W, otflw_ad_N, otflw_ad_E, otflw_ad_B, otflw_ad_T, Nx, Ny, Nz)

        ! Build neighbour coefficient matrices
        DO i = 1,Ny
            DO j = 1,Nx
                DO k = 1,Nz

                    IF  (uwind_E(i,j,k)) THEN
                        aE_ndl(i,j,k) = aX_D - F_e(i,j,k)
                    ELSE
                        aE_ndl(i,j,k) = aX_D
                    END IF

                    IF  (uwind_W(i,j,k)) THEN
                        aW_ndl(i,j,k) = aX_D + F_w(i,j,k)
                    ELSE
                        aW_ndl(i,j,k) = aX_D
                    END IF
                    
                    IF  (uwind_N(i,j,k)) THEN
                        aN_ndl(i,j,k) = aY_D - F_n(i,j,k)
                    ELSE 
                        aN_ndl(i,j,k) = aY_D
                    END IF

                    IF  (uwind_S(i,j,k)) THEN
                        aS_ndl(i,j,k) = aY_D + F_s(i,j,k)
                    ELSE
                        aS_ndl(i,j,k) = aY_D
                    END IF

                    ! IF      (uwind_T(i,j,k)) THEN
                    !     aT(i,j,k) = aZ_D - aZ_F
                    !     aB(i,j,k) = aZ_D
                    ! ELSEIF  (uwind_B(i,j,k)) THEN
                    !     aT(i,j,k) = aZ_D
                    !     aB(i,j,k) = aZ_D + aZ_F
                    ! END IF
                END DO 
            END DO
        END DO

        ! replace appropriate neighbour coefficients with zeroes at boundaries 
        DO i = 1,Ny
            DO j = 1,Nx
                DO k = 1,Nz

                    IF     (i == 1) THEN
                        aS_ndl(i,j,k) = 0
                    END IF

                    IF (i == Ny) THEN
                        aN_ndl(i,j,k) = 0
                    END IF

                    IF (j == 1) THEN
                        aW_ndl(i,j,k) = 0
                    END IF

                    IF (j == Nx) THEN
                        aE_ndl(i,j,k) = 0
                    END IF

                    ! IF (k == 1) THEN
                    !     aB(i,j,k) = 0
                    ! END IF

                    ! IF (k == Nz) THEN
                    !     aT(i,j,k) = 0
                    ! END IF


                END DO
            END DO
        END DO

        ! Build point P coeff matrix
        DO i = 1,Ny
            DO j = 1,Nx
                DO k = 1,Nz
            
                    aP_ndl(i,j,k) = aN_ndl(i,j,k) + aE_ndl(i,j,k) + aS_ndl(i,j,k) + aW_ndl(i,j,k) &
                                !+ aB(i,j,k) + aT(i,j,k) & 
                                - Sp_D_N(i,j,k) - Sp_D_E(i,j,k) - Sp_D_S(i,j,k) - Sp_D_W(i,j,k) &
                                !- Sp_D_B(i,j,k) - Sp_D_T(i,j,k)  &
                                - Sp_F_N(i,j,k) - Sp_F_E(i,j,k) - Sp_F_S(i,j,k) - Sp_F_W(i,j,k) &
                                !- Sp_F_B(i,j,k) - Sp_F_T(i,j,k) & ! linearized boundary
                                !+ (aX_F+aY_F+aZ_F) - (aX_F+aY_F+aZ_F) &  
                                + F_e(i,j,k) - F_w(i,j,k) + F_n(i,j,k) - F_s(i,j,k) &! change in advection quantity F
                                + aP_0
                END DO
            END DO
        END DO
            

        ! Build constant matrix b
        DO i = 1,Ny
            DO j = 1,Nx
                DO k = 1,Nz

                    b3D_energy(i,j,k) = Su_D_N(i,j,k) + Su_D_E(i,j,k) + Su_D_S(i,j,k) + Su_D_W(i,j,k) &
                                !+ Su_D_B(i,j,k) + Su_D_T(i,j,k) &
                                + Su_F_N(i,j,k) + Su_F_E(i,j,k) + Su_F_S(i,j,k) + Su_F_W(i,j,k) &
                                !+ Su_F_B(i,j,k) + Su_F_T(i,j,k) &
                                + aP_0*phi3D_energy(i,j,k) ! unsteady piece
                
                END DO
            END DO
        END DO  


    END SUBROUTINE ENERGY_COEFFBLDR

    SUBROUTINE ENERGY_COEFFBLDR_INT_GEOMETRY(Nx, Ny, Nz, params, u, v, w, & 
                                    aP_ndl, aN_ndl, aE_ndl, aS_ndl, aW_ndl, &
                                    aB_ndl, aT_ndl, b3D_energy, phi3D_energy)

        IMPLICIT NONE

        INTEGER, PARAMETER  :: dp = KIND(1.0D0)

        !REAL(dp), PARAMETER  :: coeff_test = 0.5

        INTEGER,  INTENT(IN)     :: Nx, Ny, Nz
        REAL(dp), DIMENSION (:,:,:), INTENT(IN)     :: u, v, w

        REAL(dp), DIMENSION(:),     ALLOCATABLE, INTENT(IN)   :: params

        REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)     :: aP_ndl, aN_ndl, aE_ndl, aS_ndl, aW_ndl, aB_ndl, aT_ndl
        REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)     :: b3D_energy
        REAL(dp), DIMENSION(:,:,:)   , INTENT(INOUT)   :: phi3D_energy

        ! IMPLICIT
        INTEGER     :: i,j,k

        REAL(dp)                :: dt

        REAL(dp)    :: diff_coeff, gamma ! alpha and specific heat

        REAL(dp)    :: del_x, del_y, del_z, del_V, Lx, Ly, Lz, alpha, cp, rho  ! params for building coeff matrices

        REAL(dp)    :: area_xz, area_xy, area_yz

        REAL(dp)    :: bN, bE, bS, bW, bB, bT   ! modified boundary conditons from BC input

        REAL(dp)    :: BC_N, BC_E, BC_S, BC_W, BC_B, BC_T   ! modified boundary conditons from BC input

        REAL(dp)    :: flux_N, flux_E, flux_S, flux_W, flux_B, flux_T   ! modified boundary conditons from BC input

        REAL(dp)    :: aX_D, aY_D, aZ_D, aP_0   ! cst reusable neighbour coeffs 

        REAL(dp), DIMENSION(Ny, Nx, Nz)    :: F_e, F_s, F_n, F_w, F_b, F_t

        REAL(dp), DIMENSION(Ny,Nx,Nz)   :: Su_D_N, Su_D_E, Su_D_S, Su_D_W, Su_D_T, Su_D_B, &   ! boundary contributions from diffusion
                                        Sp_D_N, Sp_D_E, Sp_D_S, Sp_D_W, Sp_D_T, Sp_D_B

        ! logical array for velocity direction at each node
        LOGICAL,  DIMENSION(Ny,Nx,Nz)   :: uwind_N, uwind_E, uwind_S, uwind_W, uwind_B, uwind_T

        REAL(dp), DIMENSION(Ny,Nx,Nz)   :: Su_F_N, Su_F_E, Su_F_S, Su_F_W, Su_F_B, Su_F_T, &   
                                        Sp_F_N, Sp_F_E, Sp_F_S, Sp_F_W, Sp_F_B, Sp_F_T ! boundary contributions from adv

        INTEGER :: Nx_int_start, Nx_int_end, Ny_int_start, Ny_int_end
        REAL(dp) :: block_temp, block_phi

        Lx          = params(4)
        Ly          = params(5)
        Lz          = params(6)

        dt          = params(7)

        alpha       = params(9)
        rho         = params(10)
        cp          = params(11)

        BC_N        = params(13)
        BC_E        = params(14)
        BC_S        = params(15)
        BC_W        = params(16)
        BC_T        = params(17)
        BC_B        = params(18)

        flux_N      = params(19)
        flux_E      = params(20)
        flux_S      = params(21)
        flux_W      = params(22)
        flux_T      = params(23)
        flux_B      = params(24)

        Nx_int_start=params(37)
        Nx_int_end=params(38)
        Ny_int_start=params(39)
        Ny_int_end=params(40)

        block_temp = params(41)
        block_phi  = (block_temp+273)*cp

        diff_coeff  = alpha*rho ! k /cp
        gamma = 1.4

        ! convert boundary conditions (temp deg C) to energy values (J/kg)
        CALL BC_CALC(BC_N, BC_E, BC_S, BC_W, BC_T, BC_B, cp, bN, bE, bS, bW, bT, bB)

        ! neighbours in x have cst coefficient diffusion -
        del_x = Lx/Nx
        aX_D = (diff_coeff/del_x)

        del_y = Ly/Ny
        del_z = Lz/Nz

        area_xy = del_x*del_y
        area_xz = del_x*del_z
        area_yz = del_y*del_z

        ! reusable constant  diffusion contributions
        aX_D = (diff_coeff/del_x)*area_yz
        aY_D = (diff_coeff/del_y)*area_xz
        aZ_D = (diff_coeff/del_z)*area_xy

        del_V = del_x*del_y*del_z

        ! implicit unsteady piece contribution
        aP_0 = ((rho/gamma)*del_V)/dt

        ! convection terms
        DO i = 1,Ny
            DO j = 1,Nx
                DO k = 1,Nz

                    F_w(i,j,k) = rho*u(i,j,k)*area_yz
                    F_e(i,j,k) = rho*u(i,j+1,k)*area_yz

                    F_s(i,j,k) = rho*v(i,j,k)*area_xz
                    F_n(i,j,k) = rho*v(i+1,j,k)*area_xz

                    ! F_b(i,j,k) = rho*0.5*(w_star(i-1,j,k)   + w_star(i,j,k))
                    ! F_t(i,j,k) = rho*0.5*(w_star(i-1,j,k+1) + w_star(i,j,k+1))

                    F_b(i,j,k) = 0
                    F_t(i,j,k) = 0
                END DO
            END DO
        END DO

        ! '' boolean array for upwind differencing
        DO i = 1,Ny
            DO j = 1,Nx
                DO k = 1,Nz

                    IF (u(i,j,k) .GE. 0) THEN
                        uwind_W(i,j,k) = .TRUE.
                    ELSE
                        uwind_W(i,j,k) = .FALSE.
                    END IF

                    IF (u(i,j+1,k) .LE. 0) THEN
                        uwind_E(i,j,k) = .TRUE.
                    ELSE
                        uwind_E(i,j,k) = .FALSE.
                    END IF

                    IF (v(i,j,k) .GE. 0) THEN
                        uwind_S(i,j,k) = .TRUE.
                    ELSE
                        uwind_S(i,j,k) = .FALSE.
                    END IF

                    IF (v(i+1,j,k) .LE. 0) THEN
                        uwind_N(i,j,k) = .TRUE.
                    ELSE
                        uwind_N(i,j,k) = .FALSE.
                    END IF

                    IF (w(i,j,k) .GE. 0) THEN
                        uwind_B(i,j,k) = .TRUE.
                        uwind_T(i,j,k) = .FALSE.
                    ELSE
                        uwind_B(i,j,k) = .FALSE.
                        uwind_T(i,j,k) = .TRUE.
                    END IF
                END DO
            END DO
        END DO

        ! initialize boundary matrices for diffusion with zeros
        DO i = 1,Ny
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

        ! update boundary contributions from diffusion
        DO i = 1,Ny
            DO j = 1,Nx
                DO k = 1,Nz

                    ! boundary faces
                    IF (i == 1) THEN
                        Su_D_S(i,j,k) = flux_S*area_xz
                        Sp_D_S(i,j,k) = 0
                    END IF

                    IF (j == 1) THEN
                        Su_D_W(i,j,k) = flux_W*area_yz
                        Sp_D_W(i,j,k) = 0
                    END IF

                    IF (i == Ny) THEN
                        Su_D_N(i,j,k) = flux_N*area_xz
                        Sp_D_N(i,j,k) = 0
                    END IF

                    IF (j == Nx) THEN
                        Su_D_E(i,j,k) = 2*aX_D*bN*area_yz + flux_E*area_yz
                        Sp_D_E(i,j,k) = -2*aX_D*area_yz
                    END IF

                    ! Force the source term to be extremely large inside the internal geometry
                    IF ((Ny_int_start .LE. i) .AND. (i .LE. Ny_int_end) .AND. &
                        (Nx_int_start .LE. j) .AND. (j .LE. Nx_int_end)) THEN
                        Sp_D_N(i,j,k)   = -1.0E30
                        Su_D_N(i,j,k)   = 1.0E30*block_phi
                    END IF

                END DO
            END DO
        END DO

        !PRINT *, Sp

        ! initialize boundary matrices for advection with zeros
        DO i = 1,Ny
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

        ! update boundary contributions from advection based on upwind scheme
        ! DO i = 1,Ny
        !     DO j = 1,Nx
        !         DO k = 1,Nz

        !             ! boundary faces
        !             IF     (i == 1) THEN
        !                 ! advection contribution from boundary if upwind cell is adjacent
        !                 IF      (uwind_S(i,j,k)) THEN
        !                     Su_F_S(i,j,k) = F_s(i,j,k)*bS*area_xz
        !                     Sp_F_S(i,j,k) = -F_s(i,j,k)*area_xz
        !                 END IF
        !             END IF
                
        !             IF (j == 1) THEN
        !                 IF     (uwind_W(i,j,k)) THEN
        !                     Su_F_W(i,j,k) = -F_w(i,j,k)*bW*area_yz
        !                     Sp_F_W(i,j,k) = F_w(i,j,k)*area_yz
        !                 END IF
        !             END IF
                
        !             IF (i == Ny) THEN
        !                 IF     (uwind_N(i,j,k)) THEN
        !                     Su_F_N(i,j,k) = -F_n(i,j,k)*bN*area_xz
        !                     Sp_F_N(i,j,k) = F_n(i,j,k)*area_xz
        !                 END IF
        !             END IF

        !             IF (j == Nx) THEN
        !                 IF     (uwind_E(i,j,k)) THEN
        !                     Su_F_E(i,j,k) = F_e(i,j,k)*bE*area_yz
        !                     Sp_F_E(i,j,k) = -F_e(i,j,k)*area_yz
        !                 END IF
        !             END IF
                    
        !             ! ELSEIF (k == 1) THEN
        !             !     IF     (uwind_B(i,j,k)) THEN
        !             !         Su_F_B(i,j,k) = aZ_F*bB
        !             !         Sp_F_B(i,j,k) = -aZ_F
        !             !     END IF

        !             ! ELSEIF (k == Nz) THEN
        !             !     IF     (uwind_T(i,j,k)) THEN
        !             !         Su_F_T(i,j,k) = aZ_F*bT
        !             !         Sp_F_T(i,j,k) = -aZ_F
        !             !     END IF

        !             ! END IF
        !         END DO
        !     END DO
        ! END DO

        ! Update source terms near outflow or adiabatic boundaries
        ! CALL   OUTFLW_AD(Su_F_S, Sp_F_S, Su_D_S, Sp_D_S, Su_F_W, Sp_F_W, Su_D_W, Sp_D_W, Su_F_N, Sp_F_N, Su_D_N, Sp_D_N, &
        !                 Su_F_E, Sp_F_E, Su_D_E, Sp_D_E, Su_F_B, Sp_F_B, Su_D_B, Sp_D_B, Su_F_T, Sp_F_T, Su_D_T, Sp_D_T, &
        !                 otflw_ad_S, otflw_ad_W, otflw_ad_N, otflw_ad_E, otflw_ad_B, otflw_ad_T, Nx, Ny, Nz)

        ! Build neighbour coefficient matrices
        DO i = 1,Ny
            DO j = 1,Nx
                DO k = 1,Nz

                    IF  (uwind_E(i,j,k)) THEN
                        aE_ndl(i,j,k) = aX_D - F_e(i,j,k)
                    ELSE
                        aE_ndl(i,j,k) = aX_D
                    END IF

                    IF  (uwind_W(i,j,k)) THEN
                        aW_ndl(i,j,k) = aX_D + F_w(i,j,k)
                    ELSE
                        aW_ndl(i,j,k) = aX_D
                    END IF
                    
                    IF  (uwind_N(i,j,k)) THEN
                        aN_ndl(i,j,k) = aY_D - F_n(i,j,k)
                    ELSE 
                        aN_ndl(i,j,k) = aY_D
                    END IF

                    IF  (uwind_S(i,j,k)) THEN
                        aS_ndl(i,j,k) = aY_D + F_s(i,j,k)
                    ELSE
                        aS_ndl(i,j,k) = aY_D
                    END IF

                    ! IF      (uwind_T(i,j,k)) THEN
                    !     aT(i,j,k) = aZ_D - aZ_F
                    !     aB(i,j,k) = aZ_D
                    ! ELSEIF  (uwind_B(i,j,k)) THEN
                    !     aT(i,j,k) = aZ_D
                    !     aB(i,j,k) = aZ_D + aZ_F
                    ! END IF
                END DO 
            END DO
        END DO

        ! replace appropriate neighbour coefficients with zeroes at boundaries 
        DO i = 1,Ny
            DO j = 1,Nx
                DO k = 1,Nz

                    IF     (i == 1) THEN
                        aS_ndl(i,j,k) = 0
                    END IF

                    IF (i == Ny) THEN
                        aN_ndl(i,j,k) = 0
                    END IF

                    IF (j == 1) THEN
                        aW_ndl(i,j,k) = 0
                    END IF

                    IF (j == Nx) THEN
                        aE_ndl(i,j,k) = 0
                    END IF

                    ! IF (k == 1) THEN
                    !     aB(i,j,k) = 0
                    ! END IF

                    ! IF (k == Nz) THEN
                    !     aT(i,j,k) = 0
                    ! END IF


                END DO
            END DO
        END DO

        ! Build point P coeff matrix
        DO i = 1,Ny
            DO j = 1,Nx
                DO k = 1,Nz
            
                    aP_ndl(i,j,k) = aN_ndl(i,j,k) + aE_ndl(i,j,k) + aS_ndl(i,j,k) + aW_ndl(i,j,k) &
                                !+ aB(i,j,k) + aT(i,j,k) & 
                                - Sp_D_N(i,j,k) - Sp_D_E(i,j,k) - Sp_D_S(i,j,k) - Sp_D_W(i,j,k) &
                                !- Sp_D_B(i,j,k) - Sp_D_T(i,j,k)  &
                                - Sp_F_N(i,j,k) - Sp_F_E(i,j,k) - Sp_F_S(i,j,k) - Sp_F_W(i,j,k) &
                                !- Sp_F_B(i,j,k) - Sp_F_T(i,j,k) & ! linearized boundary
                                !+ (aX_F+aY_F+aZ_F) - (aX_F+aY_F+aZ_F) &  
                                !+ F_e(i,j,k) - F_w(i,j,k) + F_n(i,j,k) - F_s(i,j,k) &! change in advection quantity F
                                + aP_0
                END DO
            END DO
        END DO
            

        ! Build constant matrix b
        DO i = 1,Ny
            DO j = 1,Nx
                DO k = 1,Nz

                    b3D_energy(i,j,k) = Su_D_N(i,j,k) + Su_D_E(i,j,k) + Su_D_S(i,j,k) + Su_D_W(i,j,k) &
                                !+ Su_D_B(i,j,k) + Su_D_T(i,j,k) &
                                + Su_F_N(i,j,k) + Su_F_E(i,j,k) + Su_F_S(i,j,k) + Su_F_W(i,j,k) &
                                !+ Su_F_B(i,j,k) + Su_F_T(i,j,k) &
                                + aP_0*phi3D_energy(i,j,k) ! unsteady piece
                
                END DO
            END DO
        END DO  


    END SUBROUTINE ENERGY_COEFFBLDR_INT_GEOMETRY

    ! SUBROUTINE for calculating boundary conditions - especially useful for unsteady BCs
    SUBROUTINE BC_CALC(BC_N, BC_E, BC_S, BC_W, BC_T, BC_B, cp, bN, bE, bS, bW, bT, bB)

        IMPLICIT NONE

        INTEGER, PARAMETER  :: dp = KIND(1.0D0)

        REAL(dp), INTENT(IN)    :: BC_N, BC_E, BC_S, BC_W, BC_T, BC_B, cp
        REAL(dp), INTENT(OUT)   :: bN, bE, bS, bW, bT, bB

            bN      = (BC_N+273)*cp
            bE      = (BC_E+273)*cp
            bS      = (BC_S+273)*cp
            bW      = (BC_W+273)*cp
            bT      = (BC_T+273)*cp
            bB      = (BC_B+273)*cp

    END SUBROUTINE BC_CALC

END MODULE 