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
    CHARACTER(:), ALLOCATABLE  :: out_file

    REAL(dp)        :: res_tol
    INTEGER         :: Nx, Ny, Nz, max_it

    ! initialize
    REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE   :: aP, aN, aE, aS, aW, aB, aT   ! coefficient matrices
    REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE   :: phi3D, b3D ! output solved for

    INTEGER     :: i, j, k ! loop index

    ! residuals
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: res      ! residual for each iteration
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)  :: phi_old 

    ! ADVECTION 
    REAL(dp)    :: mdotx, mdoty, mdotz
    ! at each cell face assign a mass flow (constant for a now)
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)   :: adv_clfc_x, adv_clfc_y, adv_clfc_z


    ! Interface block for external subroutine with assumed-shape arguments
    INTERFACE
        SUBROUTINE ADV_DF_COEFFBLDR(Nx, Ny, Nz, params, mdotx, mdoty, mdotz, adv_clfc_x, adv_clfc_y, adv_clfc_z, & 
                                    aP, aN, aE, aS, aW, aB, aT, b3D, phi3D)
            IMPLICIT NONE

            INTEGER, PARAMETER  :: dp = KIND(1.0D0)

            INTEGER,  INTENT(IN)     :: Nx, Ny, Nz
            REAL(dp), INTENT(IN)     :: mdotx, mdoty, mdotz

            REAL(dp), DIMENSION(:, :, :), ALLOCATABLE, INTENT(IN)   :: adv_clfc_x, adv_clfc_y, adv_clfc_z
            REAL(dp), DIMENSION(:),       ALLOCATABLE, INTENT(IN)   :: params

            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)   :: aP, aN, aE, aS, aW, aB, aT
            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)   :: b3D
            REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)   :: phi3D

        END SUBROUTINE ADV_DF_COEFFBLDR
    END INTERFACE
    
    CALL GET_COMMAND_ARGUMENT(1, in_file)   ! Recieve an input file from command line
    CALL READ_INPUT_MAIN(params, in_file, out_file)   ! Store paramaters from input file into params array

    ! PRELIMINARIES
    !unpack parameters needed in main
    Nx =    params(1)
    Ny =    params(2)   
    Nz =    params(3)
    max_it = params(SIZE(params)-1)
    res_tol = params(SIZE(params))

    ! boundary conditions to write to output file
    cp          = params(10)
    BC_N        = params(11)
    BC_E        = params(12)
    BC_S        = params(13)
    BC_W        = params(14)
    BC_T        = params(15)
    BC_B        = params(16)
    ! convert boundary conditions (temp deg C) to energy values (J/kg)
    bN        = (BC_N+273)*cp
    bE        = (BC_E+273)*cp
    bS        = (BC_S+273)*cp
    bW        = (BC_W+273)*cp
    bB        = (BC_B+273)*cp
    bT        = (BC_T+273)*cp

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
    ALLOCATE(res(max_it+1))
    ALLOCATE(adv_clfc_x(Ny+1,Nx+1,Nz+1))
    ALLOCATE(adv_clfc_y(Ny+1,Nx+1,Nz+1))
    ALLOCATE(adv_clfc_z(Ny+1,Nx+1,Nz+1))

    ! UPWIND ADVECTION 
    mdotx       = params(17)
    mdoty       = params(18)
    mdotz       = params(19)
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

    ! Main routine

    ! update cofficient matrices
    CALL ADV_DF_COEFFBLDR(Nx, Ny, Nz, params, mdotx, mdoty, mdotz, adv_clfc_x, adv_clfc_y, adv_clfc_z, & 
                            aP, aN, aE, aS, aW, aB, aT, b3D, phi3D)

    ! Solve current set of coefficients using ADI
    CALL ADI_3D_SOLVR_MAIN(res_tol, max_it, aP, aN, aE, aS, aW, aB, aT, phi3D, b3D)

    DEALLOCATE(aP, aN, aE, aS, aW, aB, aT, b3D, res, adv_clfc_x, adv_clfc_y, adv_clfc_z)

    CALL WRITE3D_OUTPUT_MAIN(bN, bE, bS, bW, bB, bT, phi3D, out_file)

END PROGRAM MAIN

! advection-diffusion coefficient matrix builder
SUBROUTINE ADV_DF_COEFFBLDR(Nx, Ny, Nz, params, mdotx, mdoty, mdotz, adv_clfc_x, adv_clfc_y, adv_clfc_z, & 
                            aP, aN, aE, aS, aW, aB, aT, b3D, phi3D)

    IMPLICIT NONE

    INTEGER, PARAMETER  :: dp = KIND(1.0D0)

    !REAL(dp), PARAMETER  :: coeff_test = 0.5

    INTEGER,  INTENT(IN)     :: Nx, Ny, Nz
    REAL(dp), INTENT(IN)     :: mdotx, mdoty, mdotz

    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(IN)   :: adv_clfc_x, adv_clfc_y, adv_clfc_z
    REAL(dp), DIMENSION(:),     ALLOCATABLE, INTENT(IN)   :: params

    REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)   :: aP, aN, aE, aS, aW, aB, aT
    REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)   :: b3D
    REAL(dp), DIMENSION(:,:,:)   , INTENT(OUT)   :: phi3D

    ! IMPLICIT
    INTEGER     :: i,j,k

    REAL(dp)    :: diff_coeff

    REAL(dp)    :: del_x, del_y, del_z, BC_E, BC_N, BC_S, BC_W, BC_B, BC_T, Lx, Ly, Lz, area, alpha, cp, rho  ! params for building coeff matrices

    REAL(dp)    :: area_xz, area_xy, area_yz

    REAL(dp)    :: bN, bE, bS, bW, bB, bT   ! modified boundary conditons from BC input

    REAL(dp)    :: aX_D, aY_D, aZ_D, aX_F, aY_F, aZ_F   ! cst reusable neighbour coeffs 

    REAL(dp), DIMENSION(Ny,Nx,Nz)   :: Su_D_N, Su_D_E, Su_D_S, Su_D_W, Su_D_T, Su_D_B, &   ! boundary contributions from diffusion
                                       Sp_D_N, Sp_D_E, Sp_D_S, Sp_D_W, Sp_D_T, Sp_D_B

    ! logical array for velocity direction at each node
    LOGICAL,  DIMENSION(Ny,Nx,Nz)   :: uwind_N, uwind_E, uwind_S, uwind_W, uwind_B, uwind_T

    REAL(dp), DIMENSION(Ny,Nx,Nz)   :: Su_F_N, Su_F_E, Su_F_S, Su_F_W, Su_F_B, Su_F_T, &   
                                       Sp_F_N, Sp_F_E, Sp_F_S, Sp_F_W, Sp_F_B, Sp_F_T ! boundary contributions from adv

    Lx          = params(4)
    Ly          = params(5)
    Lz          = params(6)
    area        = params(7)
    alpha       = params(8)
    rho         = params(9)
    cp          = params(10)
    BC_N        = params(11)
    BC_E        = params(12)
    BC_S        = params(13)
    BC_W        = params(14)
    BC_T        = params(15)
    BC_B        = params(16)

    diff_coeff  = alpha

    ! convert boundary conditions (temp deg C) to energy values (J/kg)
    bN        = (BC_N+273)*cp
    bE        = (BC_E+273)*cp
    bS        = (BC_S+273)*cp
    bW        = (BC_W+273)*cp
    bB        = (BC_B+273)*cp
    bT        = (BC_T+273)*cp

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

    ! '' convection -
    aX_F = (mdotx/(rho*area_yz))
    aY_F = (mdoty/(rho*area_xz))
    aZ_F = (mdotz/(rho*area_xy))

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
                IF     (i == 1) THEN
                    Su_D_S(i,j,k) = 2*aY_D*bS
                    Sp_D_S(i,j,k) = -2*aY_D
                    
                ELSEIF (j == 1) THEN
                    Su_D_W(i,j,k) = 2*aX_D*bW
                    Sp_D_W(i,j,k) = -2*aX_D
                    
                ELSEIF (i == Ny) THEN
                    Su_D_N(i,j,k) = 2*aY_D*bN
                    Sp_D_N(i,j,k) = -2*aY_D           

                ELSEIF (j == Nx) THEN
                    Su_D_E(i,j,k) = 2*aX_D*bE
                    Sp_D_E(i,j,k) = -2*aX_D
                
                ELSEIF (k == 1) THEN
                    Su_D_B(i,j,k) = 2*aZ_D*bB
                    Sp_D_B(i,j,k) = -2*aZ_D 

                ELSEIF (k == Nz) THEN
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
                            + (aX_F+aY_F+aZ_F) - (aX_F+aY_F+aZ_F) ! change in advection quantity F
            END DO
        END DO
    END DO
         

    ! Build matrix for contribution of BCs (total Su source term contribution)
    DO i = 1,Ny
        DO j = 1,Nx
            DO k = 1,Nz

                b3D(i,j,k) = Su_D_N(i,j,k) + Su_D_E(i,j,k) + Su_D_S(i,j,k) + Su_D_W(i,j,k) + Su_D_B(i,j,k) + Su_D_T(i,j,k) &
                             + Su_F_N(i,j,k) + Su_F_E(i,j,k) + Su_F_S(i,j,k) + Su_F_W(i,j,k) + Su_F_B(i,j,k) + Su_F_T(i,j,k)
            
            END DO
        END DO
    END DO  


END SUBROUTINE ADV_DF_COEFFBLDR