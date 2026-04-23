PROGRAM MAIN

    USE ADI_2D_SCLRSOLVR
    USE READ_INPUT
    USE WRITE2D_OUTPUT

    IMPLICIT NONE

    INTEGER, PARAMETER  :: dp = KIND(1.0D0)

    REAL(dp), DIMENSION(:), ALLOCATABLE   :: params   ! (Nx, Ny, Lx, Ly, Area, alpha, rho, cp,
                                                      ! mdot_N, T_N, mdot_E, T_E,  mdot_S, T_S, mdot_W, T_W,
                                                      ! max_iterations, residual tolerance)

    CHARACTER(256)             :: in_file
    CHARACTER(:), ALLOCATABLE  :: out_file

    REAL(dp)        :: res_tol
    INTEGER         :: Nx, Ny, max_it

    ! initialize
    REAL(dp), DIMENSION(:,:)   , AllOCATABLE   :: aP, aN, aE, aS, aW   ! coefficient matrices
    REAL(dp), DIMENSION(:,:)   , AllOCATABLE   :: phi2D, b2D ! output solved for

    INTEGER     :: i, j ! loop index

    ! residuals
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: res      ! residual for each iteration
    REAL(dp), ALLOCATABLE, DIMENSION(:,:)  :: phi_old 

    ! ADVECTION 
    REAL(dp)    :: mdotx, mdoty
    ! at each cell face assign a mass flow (constant for a now)
    REAL(dp), ALLOCATABLE, DIMENSION(:,:)   :: adv_clfc_x, adv_clfc_y


    ! Interface block for external subroutine with assumed-shape arguments
    INTERFACE
        SUBROUTINE ADV_DF_COEFFBLDR(Nx, Ny, params, mdotx, mdoty, adv_clfc_x, adv_clfc_y, & 
                                    aP, aN, aE, aS, aW, b2D, phi2D)
            IMPLICIT NONE

            INTEGER, PARAMETER  :: dp = KIND(1.0D0)

            INTEGER,  INTENT(IN)     :: Nx, Ny
            REAL(dp), INTENT(IN)     :: mdotx, mdoty

            REAL(dp), DIMENSION(:, :), ALLOCATABLE, INTENT(IN)   :: adv_clfc_x, adv_clfc_y
            REAL(dp), DIMENSION(:),    ALLOCATABLE, INTENT(IN)   :: params

            REAL(dp), DIMENSION(:,:)   , INTENT(OUT)   :: aP, aN, aE, aS, aW
            REAL(dp), DIMENSION(:,:)   , INTENT(OUT)   :: b2D
            REAL(dp), DIMENSION(:,:)   , INTENT(OUT)   :: phi2D

        END SUBROUTINE ADV_DF_COEFFBLDR
    END INTERFACE
    
    CALL GET_COMMAND_ARGUMENT(1, in_file)   ! Recieve an input file from command line
    CALL READ_INPUT_MAIN(params, in_file, out_file)   ! Store paramaters from input file into params array

    ! PRELIMINARIES
    !unpack parameters needed in main
    Nx =    params(1)
    Ny =    params(2)   
    max_it = params(SIZE(params)-1)
    res_tol = params(SIZE(params))
    ! allocate array sizes
    ALLOCATE(aP(Ny,Nx))
    ALLOCATE(aN(Ny,Nx))
    ALLOCATE(aE(Ny,Nx))
    ALLOCATE(aS(Ny,Nx))
    ALLOCATE(aW(Ny,Nx))
    ALLOCATE(phi2D(Ny,Nx))
    ALLOCATE(b2D(Ny,Nx))
    ALLOCATE(res(max_it+1))
    ALLOCATE(adv_clfc_x(Ny+1,Nx+1))
    ALLOCATE(adv_clfc_y(Ny+1,Nx+1))

    ! UPWIND ADVECTION 
    mdotx       = params(13)
    mdoty       = params(14)
    ! construct cell face advection matrix to be updated in main (all constant for a now)
    DO i = 1,Ny+1
        DO j = 1,Nx+1

            adv_clfc_x(i,j) = mdotx
            adv_clfc_y(i,j) = mdoty

        END DO
    END DO

    ! initialize
    phi2D = 0.0

    ! Main routine

    ! update cofficient matrices
    CALL ADV_DF_COEFFBLDR(Nx, Ny, params, mdotx, mdoty, adv_clfc_x, adv_clfc_y, & 
                            aP, aN, aE, aS, aW, b2D, phi2D)

    ! Solve current set of coefficients using ADI
    CALL ADI_2D_SCLRSOLVR_MAIN(res_tol, max_it, aP, aN, aE, aS, aW, phi2D, b2D)

    DEALLOCATE(aP, aN, aE, aS, aW, b2D, res, adv_clfc_x, adv_clfc_y)

    CALL WRITE2D_OUTPUT_MAIN(TRANSPOSE(phi2D), out_file)

END PROGRAM MAIN

! advection-diffusion coefficient matrix builder
SUBROUTINE ADV_DF_COEFFBLDR(Nx, Ny, params, mdotx, mdoty, adv_clfc_x, adv_clfc_y, & 
                            aP, aN, aE, aS, aW, b2D, phi2D)

    IMPLICIT NONE

    INTEGER, PARAMETER  :: dp = KIND(1.0D0)

    !REAL(dp), PARAMETER  :: coeff_test = 0.5

    INTEGER,  INTENT(IN)     :: Nx, Ny
    REAL(dp), INTENT(IN)     :: mdotx, mdoty

    REAL(dp), DIMENSION(:, :), ALLOCATABLE, INTENT(IN)   :: adv_clfc_x, adv_clfc_y
    REAL(dp), DIMENSION(:),    ALLOCATABLE, INTENT(IN)   :: params

    REAL(dp), DIMENSION(:,:)   , INTENT(OUT)   :: aP, aN, aE, aS, aW
    REAL(dp), DIMENSION(:,:)   , INTENT(OUT)   :: b2D
    REAL(dp), DIMENSION(:,:)   , INTENT(OUT)   :: phi2D

    ! IMPLICIT
    INTEGER     :: i,j

    REAL(dp)    :: diff_coeff

    REAL(dp)    :: del_x, del_y, BC_E, BC_N, BC_S, BC_W, Lx, Ly, area, alpha, cp, rho  ! params for building coeff matrices

    REAL(dp)    :: bN, bE, bS, bW   ! modified boundary conditons from BC input

    REAL(dp)    :: aX_D, aY_D, aX_F, ay_F   ! cst reusable neighbour coeffs 

    REAL(dp), DIMENSION(Ny  ,  Nx)   :: Su_D_N, Su_D_E, Su_D_S, Su_D_W, &   ! boundary contributions from diffusion
                   Sp_D_N, Sp_D_E, Sp_D_S, Sp_D_W

    ! logical array for velocity direction at each node
    LOGICAL,  DIMENSION(Ny  ,  Nx)   :: uwind_N, uwind_E, uwind_S, uwind_W 

    REAL(dp), DIMENSION(Ny  ,  Nx)   :: Su_F_N, Su_F_E, Su_F_S, Su_F_W, &   
                                        Sp_F_N, Sp_F_E, Sp_F_S, Sp_F_W ! boundary contributions from adv

    Lx          = params(3)
    Ly          = params(4)
    area        = params(5)
    alpha       = params(6)
    rho         = params(7)
    cp          = params(8)
    BC_N        = params(9)
    BC_E        = params(10)
    BC_S        = params(11)
    BC_W        = params(12)

    diff_coeff  = alpha

    ! convert boundary conditions (temp deg C) to energy values (J/kg)
    bN        = (BC_N+273)*cp
    bE        = (BC_E+273)*cp
    bS        = (BC_S+273)*cp
    bW        = (BC_W+273)*cp

    ! neighbours in x have cst coefficient diffusion -
    del_x = Lx/Nx
    aX_D = diff_coeff/del_x

    ! neighbouts in y have cst coefficient diffusion-
    del_y = Ly/Ny
    aY_D = diff_coeff/del_y

    ! '' convection -
    aX_F = (mdotx/(rho*area))
    aY_F = (mdoty/(rho*area))

    ! '' boolean array for upwind differencing
    DO i = 1,Ny
        DO j = 1,Nx

            IF (adv_clfc_x(i,j) .GE. 0) THEN
                uwind_W(i,j) = .TRUE.
                uwind_E(i,j) = .FALSE.
            ELSE
                uwind_W(i,j) = .FALSE.
                uwind_E(i,j) = .TRUE.
            END IF

            IF (adv_clfc_y(i,j) .GE. 0) THEN
                uwind_S(i,j) = .TRUE.
                uwind_N(i,j) = .FALSE.
            ELSE
                uwind_S(i,j) = .FALSE.
                uwind_N(i,j) = .TRUE.
            END IF

        END DO
    END DO

    ! initialize boundary matrices for diffusion with zeros
    DO i = 1,Ny
        DO j = 1,Nx
            Su_D_N(i,j) = 0
            Sp_D_N(i,j) = 0
            Su_D_E(i,j) = 0
            Sp_D_E(i,j) = 0
            Su_D_S(i,j) = 0
            Sp_D_S(i,j) = 0
            Su_D_W(i,j) = 0
            Sp_D_W(i,j) = 0
        END DO
    END DO

    ! update boundary contributions from diffusion
    DO i = 1,Ny
        DO j = 1,Nx

            ! boundary faces
            IF     (i == 1) THEN
                Su_D_S(i,j) = 2*aY_D*bS
                Sp_D_S(i,j) = -2*aY_D
                
            ELSEIF (j == 1) THEN
                Su_D_W(i,j) = 2*aX_D*bW
                Sp_D_W(i,j) = -2*aX_D
                
            ELSEIF (i == Ny) THEN
                Su_D_N(i,j) = 2*aY_D*bN
                Sp_D_N(i,j) = -2*aY_D
                

            ELSEIF (j == Nx) THEN
                Su_D_E(i,j) = 2*aX_D*bE
                Sp_D_E(i,j) = -2*aX_D    
            
            END IF
        END DO
    END DO

    ! initialize boundary matrices for advection with zeros
    DO i = 1,Ny
        DO j = 1,Nx
            Su_F_N(i,j) = 0
            Sp_F_N(i,j) = 0
            Su_F_E(i,j) = 0
            Sp_F_E(i,j) = 0
            Su_F_S(i,j) = 0
            Sp_F_S(i,j) = 0
            Su_F_W(i,j) = 0
            Sp_F_W(i,j) = 0
        END DO
    END DO

    ! update boundary contributions from advection based on upwind scheme
    DO i = 1,Ny
        DO j = 1,Nx

            ! boundary faces
            IF     (i == 1) THEN
                ! advection contribution from boundary if upwind cell is adjacent
                IF      (uwind_S(i,j)) THEN
                    Su_F_S(i,j) = aY_F*bS
                    Sp_F_S(i,j) = -aY_F
                END IF
        
            ELSEIF (j == 1) THEN

                IF     (uwind_W(i,j)) THEN
                    Su_F_W(i,j) = aX_F*bW
                    Sp_F_W(i,j) = -aX_F
                END IF
        

            ELSEIF (i == Ny) THEN
                IF     (uwind_N(i,j)) THEN
                    Su_F_N(i,j) = aX_F*bN
                    Sp_F_N(i,j) = -aX_F
                END IF

            ELSEIF (j == Nx) THEN
                IF     (uwind_E(i,j)) THEN
                    Su_F_E(i,j) = aX_F*bE
                    Sp_F_E(i,j) = -aX_F
                END IF
            
            END IF
        END DO
    END DO

    ! Build neighbour coefficient matrices
    DO i = 1,Ny
        DO j = 1,Nx

            IF      (uwind_E(i,j)) THEN
                aE(i,j) = aX_D - aX_F
                aW(i,j) = aX_D
            ELSEIF  (uwind_W(i,j)) THEN
                aE(i,j) = aX_D
                aW(i,j) = aX_D + aX_F
            END IF
            
            IF      (uwind_N(i,j)) THEN
                aN(i,j) = aY_D - aY_F
                aS(i,j) = aY_D
            ELSEIF  (uwind_S(i,j)) THEN
                aN(i,j) = aY_D
                aS(i,j) = aY_D + aY_F
            END IF
           
        END DO
    END DO

    ! replace appropriate neighbour coefficients with zeroes at boundaries 
    DO i = 1,Ny
        DO j = 1,Nx

            IF     (i == 1) THEN
                aS(i,j) = 0
            END IF

            IF (i == Ny) THEN
                aN(i,j) = 0
            END IF

            IF (j == 1) THEN
                aW(i,j) = 0
            END IF

            IF (j == Nx) THEN
                aE(i,j) = 0
            END IF
        END DO
    END DO

    ! Build point P coeff matrix
    DO i = 1,Ny
        DO j = 1,Nx
          
            aP(i,j) = aN(i,j) + aE(i,j) + aS(i,j) + aW(i,j) & 
                        - Sp_D_N(i,j) - Sp_D_E(i,j) - Sp_D_S(i,j) - Sp_D_W(i,j) &
                        - Sp_F_N(i,j) - Sp_F_E(i,j) - Sp_F_S(i,j) - Sp_F_W(i,j) & ! linearized boundary
                        + (aX_F+aY_F) - (aX_F+aY_F) ! change in advection quantity F
                 
        END DO
    END DO
         

    ! Build matrix for contribution of BCs (total Su source term contribution)
    DO i = 1,Ny
        DO j = 1,Nx

            b2D(i,j) = Su_D_N(i,j) + Su_D_E(i,j) + Su_D_S(i,j) + Su_D_W(i,j) + &
                       Su_F_N(i,j) + Su_F_E(i,j) + Su_F_S(i,j) + Su_F_W(i,j)

        END DO
    END DO  


END SUBROUTINE ADV_DF_COEFFBLDR