PROGRAM MAIN

    USE U_MOMENTUM_SLVR
    USE V_MOMENTUM_SLVR
    USE W_MOMENTUM_SLVR
    USE PRESS_CRCT_SLVR
    USE CRCTR_CVG
    USE READ_INPUT
    USE WRITE_OUTPUT

    IMPLICIT NONE

    INTEGER, PARAMETER  :: dp = KIND(1.0D0)

    ! declarations for necessary parameters in main routine
    REAL(dp), DIMENSION(:), ALLOCATABLE   :: params ! read from input

    ! BCs - velocity stuff
    REAL(dp)    ::  BC_N_u, BC_E_u, BC_S_u, BC_W_u, BC_T_u, BC_B_u, &
                    BC_N_v, BC_E_v, BC_S_v, BC_W_v, BC_T_v, BC_B_v, &
                    BC_N_w, BC_E_w, BC_S_w, BC_W_w, BC_T_w, BC_B_w

    ! BCs - energy stuff

    ! file handling
    CHARACTER(256)             :: in_file
    CHARACTER(20)              :: str           ! buffer for integer-to-string conversion
    CHARACTER(:), ALLOCATABLE  :: out_file, out_file_tmp

    ! residuals stuff
    REAL(dp)        :: res_tol_p, res_tol_vel ! residual tolerance momentum and pressure solvers
    INTEGER         :: max_it_ADI_p, max_it_ADI_vel, max_it_main
    REAL(dp)        :: res_main, res_tol_main ! residual and residual tolerance on main continuity measure

    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: res_temp_p, res_temp_vel  ! reusable array to store residuals ine each ADI call

    ! unsteady stuff
    INTEGER     :: n_iter
    REAL(dp)    :: t, dt ! current time and time step

    ! grid stuff
    INTEGER         :: Nx, Ny, Nz

    ! Declare coefficient matrices for each grid
    REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE   :: aP_ndl, aN_ndl, aE_ndl, & 
                                                    aS_ndl, aW_ndl, aB_ndl, aT_ndl   ! coefficient matrices for scalar

    REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE   :: aP_stg_u, aN_stg_u, aE_stg_u, & 
                                                    aS_stg_u, aW_stg_u, aB_stg_u, aT_stg_u   ! coefficient matrices for u momentum

    REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE   :: aP_stg_v, aN_stg_v, aE_stg_v, & 
                                                    aS_stg_v, aW_stg_v, aB_stg_v, aT_stg_v   ! coefficient matrices for v momentum

    REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE   :: aP_stg_w, aN_stg_w, aE_stg_w, & 
                                                    aS_stg_w, aW_stg_w, aB_stg_w, aT_stg_w   ! coefficient matrices for w momentum

    ! velocities and pressure solved for
    REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE   :: u_star, v_star, w_star, p_star ! guesses
    REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE   :: p_prime ! pressure correction

    REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE   :: u,v,w,p ! updated values

    ! lumped d coefficients to store and transfer between solvers
    REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE   :: d_u, d_v, d_w

    ! constant b matrices
    REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE   :: b3D_u, b3D_v, b3D_w, b3D_prime

    ! implicit stuff
    INTEGER     :: i, j, k ! loop index
    INTEGER     :: n_iter_simple ! how many times to loop through the SIMPLE algorithm

    ! residuals
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: res      ! residual for each iteration
    !REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)  :: phi_old 

    
    CALL GET_COMMAND_ARGUMENT(1, in_file)   ! Recieve an input file from command line
    CALL READ_INPUT_MAIN(params, in_file, out_file)   ! Store paramaters from input file into params array

    ! PRELIMINARIES
    !unpack parameters needed in main
    Nx =    params
    Ny =    params   
    Nz =    params
    Lx =    params
    Ly =    params   
    Lz =    params

    max_it_ADI_vel = params
    max_it_ADI_p   = params
    max_it_main    = params

    res_tol_vel  = params
    res_tol_p    = params
    res_tol_main = params

    n_iter      = params
    dt          = params

    n_iter_simple = 5

    relax_vel = 0.8 ! relaxation factors
    relax_p = 0.4

    ! boundary conditions 
    BC_N_u      = 1
    BC_E_u      = 0
    BC_S_u      = 0
    BC_W_u      = 0
    BC_T_u      = 0
    BC_B_u      = 0

    BC_N_v      = 0
    BC_E_v      = 0
    BC_S_v      = 0
    BC_W_v      = 0
    BC_T_v      = 0
    BC_B_v      = 0 

    BC_N_w      = 0
    BC_E_w      = 0
    BC_S_w      = 0
    BC_W_w      = 0
    BC_T_w      = 0
    BC_B_w      = 0

    ! allocate array sizes - coefficients
    ALLOCATE(aP_ndl(Ny,Nx,Nz))
    ALLOCATE(aN_ndl(Ny,Nx,Nz))
    ALLOCATE(aE_ndl(Ny,Nx,Nz))
    ALLOCATE(aS_ndl(Ny,Nx,Nz))
    ALLOCATE(aW_ndl(Ny,Nx,Nz))
    ALLOCATE(aB_ndl(Ny,Nx,Nz))
    ALLOCATE(aT_ndl(Ny,Nx,Nz))

    ALLOCATE(aP_stg_u(Ny,Nx+1,Nz))
    ALLOCATE(aN_stg_u(Ny,Nx+1,Nz))
    ALLOCATE(aE_stg_u(Ny,Nx+1,Nz))
    ALLOCATE(aS_stg_u(Ny,Nx+1,Nz))
    ALLOCATE(aW_stg_u(Ny,Nx+1,Nz))
    ALLOCATE(aB_stg_u(Ny,Nx+1,Nz))
    ALLOCATE(aT_stg_u(Ny,Nx+1,Nz))

    ALLOCATE(aP_stg_v(Ny+1,Nx,Nz))
    ALLOCATE(aN_stg_v(Ny+1,Nx,Nz))
    ALLOCATE(aE_stg_v(Ny+1,Nx,Nz))
    ALLOCATE(aS_stg_v(Ny+1,Nx,Nz))
    ALLOCATE(aW_stg_v(Ny+1,Nx,Nz))
    ALLOCATE(aB_stg_v(Ny+1,Nx,Nz))
    ALLOCATE(aT_stg_v(Ny+1,Nx,Nz))

    ALLOCATE(aP_stg_w(Ny,Nx,Nz+1))
    ALLOCATE(aN_stg_w(Ny,Nx,Nz+1))
    ALLOCATE(aE_stg_w(Ny,Nx,Nz+1))
    ALLOCATE(aS_stg_w(Ny,Nx,Nz+1))
    ALLOCATE(aW_stg_w(Ny,Nx,Nz+1))
    ALLOCATE(aB_stg_w(Ny,Nx,Nz+1))
    ALLOCATE(aT_stg_w(Ny,Nx,Nz+1))
    
    ! guessed values and values iterated on
    ALLOCATE(p_star(Ny,Nx,Nz))
    ALLOCATE(u_star(Ny,Nx+1,Nz))
    ALLOCATE(v_star(Ny+1,Nx,Nz))
    ALLOCATE(w_star(Ny,Nx,Nz+1))

    ! Actual values
    ALLOCATE(p(Ny,Nx,Nz))
    ALLOCATE(u(Ny,Nx+1,Nz))
    ALLOCATE(v(Ny+1,Nx,Nz))
    ALLOCATE(w(Ny,Nx,Nz+1))

    ! d coefficients to transfer between solvers
    ALLOCATE(d_u(Ny,Nx+1,Nz))
    ALLOCATE(d_v(Ny+1,Nx,Nz))
    ALLOCATE(d_w(Ny,Nx,Nz+1))

    ALLOCATE(res_temp_p(max_it_ADI_p+1))
    ALLOCATE(res_temp_vel(max_it_ADI_vel+1))

    ! Apply appropriate boundary conditons to u_star to initialize
    DO i = 1,Ny
        DO j = 1,Nx_stg ! start from the second column on the staggered grid for u
            DO k = 1,Nz

                IF (i==Ny) THEN
                    u(i,j,k) = BC_N_u
                ELSE
                    u(i,j,k) = 0
                END IF

            END DO
        END DO
    END DO

    ! initialize everything else to zero
    v = 0.0
    w = 0.0
    p = 0.0
    t = 0

    ! Main routine
    DO i = 1,n_iter_time

        DO j = 1, n_iter_simple
            ! change to updated values
            u_star = u
            v_star = v
            w_star = w
            p_star = p

            !solve momentum equations
            CALL U_MOMENTUM_SLVR_MAIN(Nx, Ny, Nz, params, del_x, del_y, del_z, & 
                                    aP_stg_u, aN_stg_u, aE_stg_u, aS_stg_u, aW_stg_u, aB_stg_u, aT_stg_u, b3D_u, &
                                    u_star, v_star, w_star, p_star, &
                                    BC_N_u, BC_E_u, BC_S_u, BC_W_u, BC_T_u, BC_B_u, &
                                    d_u, &
                                    res_tol_vel, max_it_ADI_vel)

            CALL V_MOMENTUM_SLVR_MAIN(Nx, Ny, Nz, params, del_x, del_y, del_z, & 
                                    aP_stg_v, aN_stg_v, aE_stg_v, aS_stg_v, aW_stg_v, aB_stg_v, aT_stg_v, b3D_v, &
                                    u_star, v_star, w_star, p_star, &
                                    BC_N_v, BC_E_v, BC_S_v, BC_W_v, BC_T_v, BC_B_v, &
                                    d_v, &
                                    res_tol_vel, max_it_ADI_vel)

            CALL W_MOMENTUM_SLVR_MAIN(Nx, Ny, Nz, params, del_x, del_y, del_z, & 
                                    aP_stg_w, aN_stg_w, aE_stg_w, aS_stg_w, aW_stg_w, aB_stg_w, aT_stg_w, b3D_w, &
                                    u_star, v_star, w_star, p_star, &
                                    BC_N_w, BC_E_w, BC_S_w, BC_W_w, BC_T_w, BC_B_w, &
                                    d_w, &
                                    res_tol_vel, max_it_ADI_vel)

            !solve pressure correction equations
            CALL PRESS_CRCTR_SLVR_MAIN(Nx, Ny, Nz, params, del_x, del_y, del_z, &
                                        p_prime, u_star, v_star, w_star, & 
                                        aP_ndl, aN_ndl, aE_ndl, aS_ndl, aW_ndl, aB_ndl, aT_ndl, & 
                                        d_u, d_v, d_w, & 
                                        b3D_prime, &
                                        res_tol_p, max_it_ADI_p)


            !correct pressure and velocity
            CALL CORRECTER(p, u, v, w, p_star, p_prime, u_star, v_star, w_star, d_u, d_v, d_w, &
                           relax_vel, relax_p)
        END DO

        !check if b term is convered
        CALL CONTINUITY_CVG(Nx, Ny, Nz, b3D_prime, main_it, cnty_tol, cnty_ref, cnvrged)

        IF (cnvrged == .TRUE.) THEN
            PRINT *, 'SIMPLE algorithm converged on main iteration ', i
            EXIT
        END IF

        t = t + dt

    END DO

    ! write out each velocity component
    IF (cnvrged == .TRUE.) THEN
        CALL WRITE3D_OUTPUT_MAIN(BC_N_u, BC_E_u, BC_S_u, BC_W_u, BC_B_u, BC_T_u, u(:,2:Nx,:), out_file_u)
        CALL WRITE3D_OUTPUT_MAIN(BC_N_v, BC_E_v, BC_S_v, BC_W_v, BC_B_v, BC_T_v, v(2:Ny,:,:), out_file_v)
        CALL WRITE3D_OUTPUT_MAIN(BC_N_w, BC_E_w, BC_S_w, BC_W_w, BC_B_w, BC_T_w, w(:,:,2:Nz), out_file_w)
    ELSE
        PRINT *, 'Solver did not converge'
    END IF

END PROGRAM MAIN