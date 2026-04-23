PROGRAM MAIN

    USE ADI_3D_SOLVR
    USE READ_INPUT
    USE WRITE_OUTPUT

    IMPLICIT NONE

    INTEGER, PARAMETER  :: dp = KIND(1.0D0)

    REAL(dp), DIMENSION(:), ALLOCATABLE   :: params   ! (Nx, Ny, Nz, Lx, Ly, Lz Area, alpha, rho, cp,
                                                      !  T_N, T_E, T_S, T_W, T_T, T_B, mdptx, mdpty, mdotz
                                                      ! max_iterations, residual tolerance)

    ! pointers for necessary parameters in main routine
    REAL(dp)    ::  cp, BC_N, BC_E, BC_S, BC_W, BC_T, BC_B, bN, bE, bS, bW, bT, bB

    CHARACTER(256)             :: in_file
    CHARACTER(20) :: str           ! explicit buffer for integer-to-string conversion
    CHARACTER(:), ALLOCATABLE  :: out_file
    CHARACTER(:), ALLOCATABLE  :: out_file_tmp

    REAL(dp)        :: res_tol
    INTEGER         :: Nx, Ny, Nz, max_it_ADI

    ! unsteady stuff
    INTEGER     :: n_iter
    REAL(dp)    :: t, dt ! current time and time step
    LOGICAL     :: unst_N, unst_E, unst_S, unst_W, unst_T, unst_B

    ! initialize
    REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE   :: aP, aN, aE, aS, aW, aB, aT   ! coefficient matrices
    REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE   :: phi3D, b3D ! output solved for

    INTEGER     :: i, j, k ! loop index

    ! residuals
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: res      ! residual for each iteration
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)  :: phi_old 

    ! ADVECTION 
    REAL(dp)    :: mdotx, mdoty, mdotz
    ! at each cell face assign a mass flow (constant for a now)
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)   :: adv_clfc_x, adv_clfc_y, adv_clfc_z


    ! Interface block for external subroutine with assumed-shape arguments
    INTERFACE
        SUBROUTINE UNST_CV_DIFF_COEFFBLDR(Nx, Ny, Nz, params, mdotx, mdoty, mdotz, adv_clfc_x, adv_clfc_y, adv_clfc_z, & 
                                          aP, aN, aE, aS, aW, aB, aT, b3D, phi3D, t)
            IMPLICIT NONE

            INTEGER, PARAMETER  :: dp = KIND(1.0D0)

            INTEGER,  INTENT(IN)     :: Nx, Ny, Nz
            REAL(dp), INTENT(IN)     :: mdotx, mdoty, mdotz

            REAL(dp), INTENT(IN)    :: t

            REAL(dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(IN)   :: adv_clfc_x, adv_clfc_y, adv_clfc_z
            REAL(dp), DIMENSION(:),       ALLOCATABLE, INTENT(IN)   :: params

            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)   :: aP, aN, aE, aS, aW, aB, aT
            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)   :: b3D
            REAL(dp), DIMENSION(:,:,:)   , INTENT(INOUT)   :: phi3D

        END SUBROUTINE UNST_CV_DIFF_COEFFBLDR
    END INTERFACE
    
    CALL GET_COMMAND_ARGUMENT(1, in_file)   ! Recieve an input file from command line
    CALL READ_INPUT_MAIN(params, in_file, out_file)   ! Store paramaters from input file into params array

    ! PRELIMINARIES
    !unpack parameters needed in main
    Nx =    params(1)
    Ny =    params(2)   
    Nz =    params(3)

    max_it_ADI = params(SIZE(params)-1)
    res_tol    = params(SIZE(params))

    n_iter      = params(7)
    dt          = params(8)

    ! boundary conditions to write to output file
    cp          = params(11)
    BC_N        = params(12)
    BC_E        = params(13)
    BC_S        = params(14)
    BC_W        = params(15)
    BC_T        = params(16)
    BC_B        = params(17)

    ! allocate array sizes
    ALLOCATE(aP(Ny,Nx,Nz))
    ALLOCATE(aN(Ny,Nx,Nz))
    ALLOCATE(aE(Ny,Nx,Nz))
    ALLOCATE(aS(Ny,Nx,Nz))
    ALLOCATE(aW(Ny,Nx,Nz))
    ALLOCATE(aB(Ny,Nx,Nz))
    ALLOCATE(aT(Ny,Nx,Nz))
    ALLOCATE(phi3D(Ny,Nx,Nz))
    ALLOCATE(b3D(Ny,Nx,Nz))
    ALLOCATE(res(max_it_ADI+1))
    ALLOCATE(adv_clfc_x(Ny+1,Nx+1,Nz+1))
    ALLOCATE(adv_clfc_y(Ny+1,Nx+1,Nz+1))
    ALLOCATE(adv_clfc_z(Ny+1,Nx+1,Nz+1))

    ! UPWIND ADVECTION 
    mdotx       = params(24)
    mdoty       = params(25)
    mdotz       = params(26)
    ! construct cell face advection matrix to be updated in main (all constant for a now)
    DO i = 1,Ny+1
        DO j = 1,Nx+1
            DO k = 1, Nz+1

            adv_clfc_x(i,j,k) = mdotx
            adv_clfc_y(i,j,k) = mdoty
            adv_clfc_z(i,j,k) = mdotz
            
            END DO
        END DO
    END DO

    ! initialize
    phi3D = 0.0
    t = 0
    ! Main routine
    DO i = 1,n_iter
        ! update cofficient matrices
        CALL UNST_CV_DIFF_COEFFBLDR(Nx, Ny, Nz, params, mdotx, mdoty, mdotz, adv_clfc_x, adv_clfc_y, adv_clfc_z, & 
                                    aP, aN, aE, aS, aW, aB, aT, b3D, phi3D, t)


        ! Solve current set of coefficients using ADI
        CALL ADI_3D_SOLVR_MAIN(res_tol, max_it_ADI, aP, aN, aE, aS, aW, aB, aT, phi3D, b3D)

        ! update BCs to write out
        CALL BC_CALC(BC_N, BC_E, BC_S, BC_W, BC_T, BC_B, cp, bN, bE, bS, bW, bT, bB)

        ! format iteration number into a string
        WRITE(str, '(I0)') i
        out_file_tmp = TRIM(out_file) // TRIM(str)
        ! write output for current timestep
        CALL WRITE3D_OUTPUT_MAIN(bN, bE, bS, bW, bB, bT, phi3D, out_file_tmp)

        t = t + dt
    END DO

    DEALLOCATE(aP, aN, aE, aS, aW, aB, aT, b3D, res, adv_clfc_x, adv_clfc_y, adv_clfc_z)

END PROGRAM MAIN

! advection-diffusion coefficient matrix builder
SUBROUTINE UNST_CV_DIFF_COEFFBLDR(Nx, Ny, Nz, params, mdotx, mdoty, mdotz, adv_clfc_x, adv_clfc_y, adv_clfc_z, & 
                                  aP, aN, aE, aS, aW, aB, aT, b3D, phi3D, t)

    IMPLICIT NONE

    INTEGER, PARAMETER  :: dp = KIND(1.0D0)

    !REAL(dp), PARAMETER  :: coeff_test = 0.5

    INTEGER,  INTENT(IN)     :: Nx, Ny, Nz
    REAL(dp), INTENT(IN)     :: mdotx, mdoty, mdotz

    REAL(dp), INTENT(IN)    :: t
    REAL(dp)                :: dt

    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(IN)   :: adv_clfc_x, adv_clfc_y, adv_clfc_z
    REAL(dp), DIMENSION(:),     ALLOCATABLE, INTENT(IN)   :: params

    REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)   :: aP, aN, aE, aS, aW, aB, aT
    REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)   :: b3D
    REAL(dp), DIMENSION(:,:,:)   , INTENT(INOUT)   :: phi3D

    ! IMPLICIT
    INTEGER     :: i,j,k

    REAL(dp)    :: diff_coeff, gamma ! alpha and specific heat

    REAL(dp)    :: del_x, del_y, del_z, del_V, BC_E, BC_N, BC_S, BC_W, BC_B, BC_T, Lx, Ly, Lz, area, alpha, cp, rho  ! params for building coeff matrices

    REAL(dp)    :: area_xz, area_xy, area_yz

    REAL(dp)    :: bN, bE, bS, bW, bB, bT   ! modified boundary conditons from BC input

    REAL(dp)    :: aX_D, aY_D, aZ_D, aX_F, aY_F, aZ_F, aP_0   ! cst reusable neighbour coeffs 

    REAL(dp), DIMENSION(Ny,Nx,Nz)   :: Su_D_N, Su_D_E, Su_D_S, Su_D_W, Su_D_T, Su_D_B, &   ! boundary contributions from diffusion
                                       Sp_D_N, Sp_D_E, Sp_D_S, Sp_D_W, Sp_D_T, Sp_D_B

    ! logical array for velocity direction at each node
    LOGICAL,  DIMENSION(Ny,Nx,Nz)   :: uwind_N, uwind_E, uwind_S, uwind_W, uwind_B, uwind_T

    REAL(dp), DIMENSION(Ny,Nx,Nz)   :: Su_F_N, Su_F_E, Su_F_S, Su_F_W, Su_F_B, Su_F_T, &   
                                       Sp_F_N, Sp_F_E, Sp_F_S, Sp_F_W, Sp_F_B, Sp_F_T ! boundary contributions from adv

    ! logical values for boundaries which are adiabatic or outflow                                    
    LOGICAL     :: otflw_ad_S, otflw_ad_W, otflw_ad_N, otflw_ad_E, otflw_ad_B, otflw_ad_T

    Lx          = params(4)
    Ly          = params(5)
    Lz          = params(6)

    dt          = params(8)

    alpha       = params(9)
    rho         = params(10)
    cp          = params(11)

    BC_N        = params(12)
    BC_E        = params(13)
    BC_S        = params(14)
    BC_W        = params(15)
    BC_T        = params(16)
    BC_B        = params(17)
    otflw_ad_N      = (params(18) /= 0.0)
    otflw_ad_E      = (params(19) /= 0.0)
    otflw_ad_S      = (params(20) /= 0.0)
    otflw_ad_W      = (params(21) /= 0.0)
    otflw_ad_T      = (params(22) /= 0.0)
    otflw_ad_B      = (params(23) /= 0.0)

    diff_coeff  = alpha 
    gamma = 1.4

    ! convert boundary conditions (temp deg C) to energy values (J/kg)
    CALL BC_CALC(BC_N, BC_E, BC_S, BC_W, BC_T, BC_B, cp, bN, bE, bS, bW, bT, bB)

    ! neighbours in x have cst coefficient diffusion -
    del_x = Lx/Nx
    aX_D = diff_coeff/del_x

    ! neighbouts in y have cst coefficient diffusion-
    del_y = Ly/Ny
    aY_D = diff_coeff/del_y

    del_z = Lz/Nz
    aZ_D = diff_coeff/del_z

    area_xy = del_x*del_y
    area_xz = del_x*del_z
    area_yz = del_y*del_z

    del_V = del_x*del_y*del_z

    ! '' convection -
    aX_F = (mdotx/(rho*area_yz))
    aY_F = (mdoty/(rho*area_xz))
    aZ_F = (mdotz/(rho*area_xy))

    ! implicit unsteady piece contribution
    aP_0 = ((rho/gamma)*del_V)/dt

    ! '' boolean array for upwind differencing
    DO i = 1,Ny
        DO j = 1,Nx
            DO k = 1,Nz

                IF (adv_clfc_x(i,j,k) .GE. 0) THEN
                    uwind_W(i,j,k) = .TRUE.
                    uwind_E(i,j,k) = .FALSE.
                ELSE
                    uwind_W(i,j,k) = .FALSE.
                    uwind_E(i,j,k) = .TRUE.
                END IF

                IF (adv_clfc_y(i,j,k) .GE. 0) THEN
                    uwind_S(i,j,k) = .TRUE.
                    uwind_N(i,j,k) = .FALSE.
                ELSE
                    uwind_S(i,j,k) = .FALSE.
                    uwind_N(i,j,k) = .TRUE.
                END IF

                IF (adv_clfc_z(i,j,k) .GE. 0) THEN
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
                    Su_D_S(i,j,k) = 2*aY_D*bS
                    Sp_D_S(i,j,k) = -2*aY_D
                END IF

                IF (j == 1) THEN
                    Su_D_W(i,j,k) = 2*aX_D*bW
                    Sp_D_W(i,j,k) = -2*aX_D
                END IF

                IF (i == Ny) THEN
                    Su_D_N(i,j,k) = 2*aY_D*bN
                    Sp_D_N(i,j,k) = -2*aY_D
                END IF

                IF (j == Nx) THEN
                    Su_D_E(i,j,k) = 2*aX_D*bE
                    Sp_D_E(i,j,k) = -2*aX_D
                END IF

                IF (k == 1) THEN
                    Su_D_B(i,j,k) = 2*aZ_D*bB
                    Sp_D_B(i,j,k) = -2*aZ_D
                END IF

                IF (k == Nz) THEN
                    Su_D_T(i,j,k) = 2*aZ_D*bT
                    Sp_D_T(i,j,k) = -2*aZ_D
                END IF

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
                        Su_F_S(i,j,k) = aY_F*bS
                        Sp_F_S(i,j,k) = -aY_F
                    END IF
            
                ELSEIF (j == 1) THEN

                    IF     (uwind_W(i,j,k)) THEN
                        Su_F_W(i,j,k) = aX_F*bW
                        Sp_F_W(i,j,k) = -aX_F
                    END IF
            

                ELSEIF (i == Ny) THEN
                    IF     (uwind_N(i,j,k)) THEN
                        Su_F_N(i,j,k) = aX_F*bN
                        Sp_F_N(i,j,k) = -aX_F
                    END IF

                ELSEIF (j == Nx) THEN
                    IF     (uwind_E(i,j,k)) THEN
                        Su_F_E(i,j,k) = aX_F*bE
                        Sp_F_E(i,j,k) = -aX_F
                    END IF
                
                ELSEIF (k == 1) THEN
                    IF     (uwind_B(i,j,k)) THEN
                        Su_F_B(i,j,k) = aZ_F*bB
                        Sp_F_B(i,j,k) = -aZ_F
                    END IF

                ELSEIF (k == Nz) THEN
                    IF     (uwind_T(i,j,k)) THEN
                        Su_F_T(i,j,k) = aZ_F*bT
                        Sp_F_T(i,j,k) = -aZ_F
                    END IF

                END IF
            END DO
        END DO
    END DO

    ! Update source terms near outflow or adiabatic boundaries
    CALL   OUTFLW_AD(Su_F_S, Sp_F_S, Su_D_S, Sp_D_S, Su_F_W, Sp_F_W, Su_D_W, Sp_D_W, Su_F_N, Sp_F_N, Su_D_N, Sp_D_N, &
                     Su_F_E, Sp_F_E, Su_D_E, Sp_D_E, Su_F_B, Sp_F_B, Su_D_B, Sp_D_B, Su_F_T, Sp_F_T, Su_D_T, Sp_D_T, &
                     otflw_ad_S, otflw_ad_W, otflw_ad_N, otflw_ad_E, otflw_ad_B, otflw_ad_T, Nx, Ny, Nz)

    ! Build neighbour coefficient matrices
    DO i = 1,Ny
        DO j = 1,Nx
            DO k = 1,Nz

                IF      (uwind_E(i,j,k)) THEN
                    aE(i,j,k) = aX_D - aX_F
                    aW(i,j,k) = aX_D
                ELSEIF  (uwind_W(i,j,k)) THEN
                    aE(i,j,k) = aX_D
                    aW(i,j,k) = aX_D + aX_F
                END IF
                
                IF      (uwind_N(i,j,k)) THEN
                    aN(i,j,k) = aY_D - aY_F
                    aS(i,j,k) = aY_D
                ELSEIF  (uwind_S(i,j,k)) THEN
                    aN(i,j,k) = aY_D
                    aS(i,j,k) = aY_D + aY_F
                END IF

                IF      (uwind_T(i,j,k)) THEN
                    aT(i,j,k) = aZ_D - aZ_F
                    aB(i,j,k) = aZ_D
                ELSEIF  (uwind_B(i,j,k)) THEN
                    aT(i,j,k) = aZ_D
                    aB(i,j,k) = aZ_D + aZ_F
                END IF
            END DO 
        END DO
    END DO

    ! replace appropriate neighbour coefficients with zeroes at boundaries 
    DO i = 1,Ny
        DO j = 1,Nx
            DO k = 1,Nz

                IF     (i == 1) THEN
                    aS(i,j,k) = 0
                END IF

                IF (i == Ny) THEN
                    aN(i,j,k) = 0
                END IF

                IF (j == 1) THEN
                    aW(i,j,k) = 0
                END IF

                IF (j == Nx) THEN
                    aE(i,j,k) = 0
                END IF

                IF (k == 1) THEN
                    aB(i,j,k) = 0
                END IF

                IF (k == Nz) THEN
                    aT(i,j,k) = 0
                END IF


            END DO
        END DO
    END DO

    ! Build point P coeff matrix
    DO i = 1,Ny
        DO j = 1,Nx
            DO k = 1,Nz
          
                aP(i,j,k) = aN(i,j,k) + aE(i,j,k) + aS(i,j,k) + aW(i,j,k)+ aB(i,j,k) + aT(i,j,k) & 
                            - Sp_D_N(i,j,k) - Sp_D_E(i,j,k) - Sp_D_S(i,j,k) - Sp_D_W(i,j,k) - Sp_D_B(i,j,k) - Sp_D_T(i,j,k)  &
                            - Sp_F_N(i,j,k) - Sp_F_E(i,j,k) - Sp_F_S(i,j,k) - Sp_F_W(i,j,k) - Sp_F_B(i,j,k) - Sp_F_T(i,j,k) & ! linearized boundary
                            + (aX_F+aY_F+aZ_F) - (aX_F+aY_F+aZ_F) &  ! change in advection quantity F
                            + aP_0
            END DO
        END DO
    END DO
         

    ! Build constant matrix b
    DO i = 1,Ny
        DO j = 1,Nx
            DO k = 1,Nz

                b3D(i,j,k) = Su_D_N(i,j,k) + Su_D_E(i,j,k) + Su_D_S(i,j,k) + Su_D_W(i,j,k) + Su_D_B(i,j,k) + Su_D_T(i,j,k) &
                             + Su_F_N(i,j,k) + Su_F_E(i,j,k) + Su_F_S(i,j,k) + Su_F_W(i,j,k) + Su_F_B(i,j,k) + Su_F_T(i,j,k) &
                             + aP_0*phi3D(i,j,k) ! unsteady piece
            
            END DO
        END DO
    END DO  


END SUBROUTINE UNST_CV_DIFF_COEFFBLDR

SUBROUTINE OUTFLW_AD(Su_F_S, Sp_F_S, Su_D_S, Sp_D_S, Su_F_W, Sp_F_W, Su_D_W, Sp_D_W, Su_F_N, Sp_F_N, Su_D_N, Sp_D_N, &
                     Su_F_E, Sp_F_E, Su_D_E, Sp_D_E, Su_F_B, Sp_F_B, Su_D_B, Sp_D_B, Su_F_T, Sp_F_T, Su_D_T, Sp_D_T, &
                     otflw_ad_S, otflw_ad_W, otflw_ad_N, otflw_ad_E, otflw_ad_B, otflw_ad_T, Nx, Ny, Nz)

    INTEGER, PARAMETER :: dp = KIND(1.0D0)

    REAL(dp), DIMENSION(Ny,Nx,Nz), INTENT(INOUT)   :: Su_F_S, Sp_F_S, Su_D_S, Sp_D_S, &
                                                      Su_F_W, Sp_F_W, Su_D_W, Sp_D_W, &
                                                      Su_F_N, Sp_F_N, Su_D_N, Sp_D_N, &
                                                      Su_F_E, Sp_F_E, Su_D_E, Sp_D_E, &
                                                      Su_F_B, Sp_F_B, Su_D_B, Sp_D_B, &
                                                      Su_F_T, Sp_F_T, Su_D_T, Sp_D_T
    
    LOGICAL, INTENT(IN) :: otflw_ad_S, otflw_ad_W, otflw_ad_N, otflw_ad_E, otflw_ad_B, otflw_ad_T

    INTEGER, INTENT(IN) :: Nx, Ny, Nz

    INTEGER :: i,j,k

    ! Cells adjacent boundaries which are either outflows or adiabatic need to be have a source term of zero
    DO i = 1,Ny
        DO j = 1,Nx
            DO k = 1,Nz

                ! boundary faces
                IF     (i == 1 .AND. otflw_ad_S) THEN
                   
                    Su_F_S(i,j,k) = 0
                    Sp_F_S(i,j,k) = 0
                    Su_D_S(i,j,k) = 0
                    Sp_D_S(i,j,k) = 0

                END IF
                    
                IF (j == 1 .AND. otflw_ad_W) THEN

                    Su_F_W(i,j,k) = 0
                    Sp_F_W(i,j,k) = 0
                    Su_D_W(i,j,k) = 0
                    Sp_D_W(i,j,k) = 0

                END IF
            

                IF (i == Ny .AND. otflw_ad_N) THEN
                    
                    Su_F_N(i,j,k) = 0
                    Sp_F_N(i,j,k) = 0
                    Su_D_N(i,j,k) = 0
                    Sp_D_N(i,j,k) = 0
                   
                END IF

                IF (j == Nx .AND. otflw_ad_E) THEN

                    Su_F_E(i,j,k) = 0
                    Sp_F_E(i,j,k) = 0
                    Su_D_E(i,j,k) = 0
                    Sp_D_E(i,j,k) = 0

                END IF
                
                IF (k == 1 .AND. otflw_ad_B) THEN

                    Su_F_B(i,j,k) = 0
                    Sp_F_B(i,j,k) = 0
                    Su_D_B(i,j,k) = 0
                    Sp_D_B(i,j,k) = 0

                END IF

                IF (k == Nz .AND. otflw_ad_T) THEN

                    Su_F_T(i,j,k) = 0
                    Sp_F_T(i,j,k) = 0
                    Su_D_T(i,j,k) = 0
                    Sp_D_T(i,j,k) = 0

                END IF

            END DO
        END DO
    END DO

END SUBROUTINE OUTFLW_AD

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

! function used to generate unsteady T profile
REAL FUNCTION T_PWRLAW(t, c, T_unst)

    INTEGER, PARAMETER  :: dp = KIND(1.0D0)
    REAL(dp)            :: t, c, T_unst
    
    T_unst = -((0.1+t)**(-1.3)) + c
    RETURN

END FUNCTION