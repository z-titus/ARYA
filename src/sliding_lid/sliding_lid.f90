PROGRAM MAIN

    USE U_MOMENTUM_SLVR
    USE V_MOMENTUM_SLVR
    USE W_MOMENTUM_SLVR
    USE PRESS_CRCTR_SLVR
    USE CRCTR_CVG
    USE READ_INPUT
    USE WRITE_OUTPUT
    USE ADI_2D_SCLRSOLVR

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
    INTEGER     :: n_iter_time
    REAL(dp)    :: t, dt ! current time and time step

    ! grid stuff
    INTEGER         :: Nx, Ny, Nz
    REAL(dp)        :: Lx, Ly, Lz
    REAL(dp)        :: del_x, del_y, del_z

    ! flow properties
    REAL(dp)    :: Re, rho, visc

    ! relaxation factors
    REAL(dp)    :: relax_p, relax_vel

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
    REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE   :: u_0, v_0, w_0  ! old values used for unsteady part

    ! lumped d coefficients to store and transfer between solvers
    REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE   :: d_u, d_v, d_w

    ! constant b matrices
    REAL(dp), DIMENSION(:,:,:)   , AllOCATABLE   :: b3D_u, b3D_v, b3D_w, b3D_prime

    ! implicit stuff
    INTEGER     :: i, j, k, l, m, n ! loop index
    INTEGER     :: n_iter_simple ! how many times to loop through the SIMPLE algorithm

    ! residuals
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: res      ! residual for each iteration

    REAL(dp)    :: cnty_ref, cnty_tol
    LOGICAL     :: cnvrged

    !output
    CHARACTER(:), ALLOCATABLE  :: out_file_u, out_file_v, out_file_w, out_file_p

    
    CALL GET_COMMAND_ARGUMENT(1, in_file)   ! Recieve an input file from command line
    CALL READ_INPUT_MAIN(params, in_file, out_file)   ! Store paramaters from input file into params array

    ! PRELIMINARIES
    !unpack parameters needed in main
    Nx =    params(1)
    Ny =    params(2)  
    Nz =    params(3)
    Lx =    params(4)
    Ly =    params(5) 
    Lz =    params(6)

    del_x = Lx/Nx
    del_y = Ly/Ny
    del_z = Lz/Nz

    dt = params(7)

    Re = params(8)

    rho         = params(10)

    max_it_ADI_vel = params(25)
    max_it_ADI_p   = params(26)

    res_tol_vel    = params(27)
    res_tol_p      = params(28)

    n_iter_time    = params(29)
    n_iter_simple  = params(30)
    cnty_tol       = params(31)

    relax_vel      = params(32)
    relax_p        = params(33)

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

    ! viscosity is the 'controlled' variable - to keep time scales constant with changing Re
    visc = (rho*BC_N_u*Ly)/Re

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

    ! pressure correction
    ALLOCATE(p_prime(Ny,Nx,Nz))

    ! Actual values
    ALLOCATE(p(Ny,Nx,Nz))
    ALLOCATE(u(Ny,Nx+1,Nz))
    ALLOCATE(v(Ny+1,Nx,Nz))
    ALLOCATE(w(Ny,Nx,Nz+1))

    ! Old values
    ALLOCATE(u_0(Ny,Nx+1,Nz))
    ALLOCATE(v_0(Ny+1,Nx,Nz))
    ALLOCATE(w_0(Ny,Nx,Nz+1))

    ! d coefficients to transfer between solvers
    ALLOCATE(d_u(Ny,Nx+1,Nz))
    ALLOCATE(d_v(Ny+1,Nx,Nz))
    ALLOCATE(d_w(Ny,Nx,Nz+1))

    ALLOCATE(b3D_prime(Ny,Nx,Nz))
    ALLOCATE(b3D_u(Ny,Nx+1,Nz))
    ALLOCATE(b3D_v(Ny+1,Nx,Nz))
    ALLOCATE(b3D_w(Ny,Nx,Nz+1))

    ALLOCATE(res_temp_p(max_it_ADI_p+1))
    ALLOCATE(res_temp_vel(max_it_ADI_vel+1))

    ! initialize everything else to zero
    u_0 = 0.0
    v_0 = 0.0
    w_0 = 0.0

    u = 0.0
    v = 0.0
    w = 0.0
    p = 0.0
    t = 0

    ! Main routine
    DO i = 1,n_iter_time

        PRINT *, 'in time loop', i

        DO j = 1, n_iter_simple
            ! anchor the pressure
            ! DO l = 1,Ny
            !     DO m = 1,Nx
            !         DO n = 1,Nz

            !             p(l,m,n) = p(l,m,n) - p(1,1,1)
                
            !         END DO
            !     END DO
            ! END DO
            p(1,1,1) = 0

            ! change to updated values
            u_star = u
            v_star = v
            w_star = 0.
            p_star = p

            !build coeffs for the momentum eqns
            CALL U_COEFF_BLDR(Nx, Ny, Nz, params, visc, del_x, del_y, del_z, & 
                                aP_stg_u, aN_stg_u, aE_stg_u, aS_stg_u, aW_stg_u, aB_stg_u, aT_stg_u, b3D_u, &
                                u_star, v_star, w_star, p_star, &
                                BC_N_u, BC_E_u, BC_S_u, BC_W_u, BC_T_u, BC_B_u, &
                                d_u, u_0)

            CALL V_COEFF_BLDR(Nx, Ny, Nz, params, visc, del_x, del_y, del_z, & 
                                aP_stg_v, aN_stg_v, aE_stg_v, aS_stg_v, aW_stg_v, aB_stg_v, aT_stg_v, b3D_v, &
                                u_star, v_star, w_star, p_star, &
                                BC_N_v, BC_E_v, BC_S_v, BC_W_v, BC_T_v, BC_B_v, &
                                d_v, v_0)

            CALL W_COEFF_BLDR(Nx, Ny, Nz, params, visc, del_x, del_y, del_z, & 
                                aP_stg_w, aN_stg_w, aE_stg_w, aS_stg_w, aW_stg_w, aB_stg_w, aT_stg_w, b3D_w, &
                                u_star, v_star, w_star, p_star, &
                                BC_N_w, BC_E_w, BC_S_w, BC_W_w, BC_T_w, BC_B_w, &
                                d_w, w_0)

            ! solve momentum eequations with coeffs
                             
            ! CALL ADI_3D_SOLVR_MAIN(res_tol_vel, max_it_ADI_vel, aP_stg_u(:,2:Nx,:), aN_stg_u(:,2:Nx,:), &
            !                        aE_stg_u(:,2:Nx,:), aS_stg_u(:,2:Nx,:), aW_stg_u(:,2:Nx,:), aB_stg_u(:,2:Nx,:), &
            !                        aT_stg_u(:,2:Nx,:), u_star(:,2:Nx,:), b3D_u(:,2:Nx,:))
            PRINT *, 'Solving u momentum...'
            CALL ADI_2D_SCLRSOLVR_MAIN(res_tol_vel, max_it_ADI_vel, aP_stg_u(:,2:Nx,1), aN_stg_u(:,2:Nx,1), &
                                   aE_stg_u(:,2:Nx,1), aS_stg_u(:,2:Nx,1), aW_stg_u(:,2:Nx,1), u_star(:,2:Nx,1), b3D_u(:,2:Nx,1))
            
            ! CALL ADI_3D_SOLVR_MAIN(res_tol_vel, max_it_ADI_vel, aP_stg_v(2:Ny,:,:), aN_stg_v(2:Ny,:,:), &
            !                        aE_stg_v(2:Ny,:,:), aS_stg_v(2:Ny,:,:), aW_stg_v(2:Ny,:,:), aB_stg_v(2:Ny,:,:), &
            !                        aT_stg_v(2:Ny,:,:), v_star(2:Ny,:,:), b3D_v(2:Ny,:,:))

            PRINT *, 'Solving v momentum...'
            CALL ADI_2D_SCLRSOLVR_MAIN(res_tol_vel, max_it_ADI_vel, aP_stg_v(2:Ny,:,1), aN_stg_v(2:Ny,:,1), &
                                   aE_stg_v(2:Ny,:,1), aS_stg_v(2:Ny,:,1), aW_stg_v(2:Ny,:,1), v_star(2:Ny,:,1), b3D_v(2:Ny,:,1))
            
            ! CALL ADI_3D_SOLVR_MAIN(res_tol_vel, max_it_ADI_vel, aP_stg_w(:,:,2:Nz), aN_stg_w(:,:,2:Nz), &
            !                         aE_stg_w(:,:,2:Nz), aS_stg_w(:,:,2:Nz), aW_stg_w(:,:,2:Nz), aB_stg_w(:,:,2:Nz), &
            !                         aT_stg_w(:,:,2:Nz), w_star(:,:,2:Nz), b3D_w(:,:,2:Nz))
                   
            !solve pressure correction equations
            CALL PRESS_CRCTR_COEFFS_BLDR(Nx, Ny, Nz, params, del_x, del_y, del_z, &
                                        u_star, v_star, w_star, & 
                                        aP_ndl, aN_ndl, aE_ndl, aS_ndl, aW_ndl, aB_ndl, aT_ndl, & 
                                        d_u, d_v, d_w, & 
                                        b3D_prime)
            PRINT *, 'Solving pressure correction...'
            CALL ADI_2D_SCLRSOLVR_MAIN(res_tol_p, max_it_ADI_p, aP_ndl(:,:,1), aN_ndl(:,:,1), aE_ndl(:,:,1), aS_ndl(:,:,1), &
                                      aW_ndl(:,:,1), p_prime(:,:,1), b3D_prime(:,:,1))


            ! CALL ADI_3D_SOLVR_MAIN(res_tol_p, max_it_ADI_p, &
            !                        aP_ndl, aN_ndl, aE_ndl, aS_ndl, aW_ndl, aB_ndl, aT_ndl, p_prime, b3D_prime)
            
            
            out_file_u = TRIM(out_file) // 'u'
            out_file_v = TRIM(out_file) // 'v'
            out_file_p = TRIM(out_file) // 'p'
            
            
            CALL WRITE2D_OUTPUT_MAIN(BC_N_u, BC_E_u, BC_S_u, BC_W_u, u(:,2:Nx,1), out_file_u)
            CALL WRITE2D_OUTPUT_MAIN(BC_N_v, BC_E_v, BC_S_v, BC_W_v, v(2:Ny,:,1), out_file_v)
            CALL WRITE2D_OUTPUT_MAIN(BC_N_v, BC_E_v, BC_S_v, BC_W_v, p(:,:,1), out_file_p)
            
            !correct pressure and velocity
            CALL CORRECTER(Nx, Ny, Nz, p, u, v, w, p_star, p_prime, u_star, v_star, w_star, d_u, d_v, d_w, &
                             relax_vel, relax_p)
                    
            
        END DO

        !check if b term is convered
        CALL CONTINUITY_CVG(Nx, Ny, Nz, b3D_prime, i, cnty_tol, cnty_ref, cnvrged, out_file)
       

        IF (cnvrged) THEN
            PRINT *, 'SIMPLE algorithm converged on main iteration ', i
            EXIT
        ELSE
            PRINT *, 'Not yet converged'
        END IF

        ! update previous timesteps values
        u_0 = u
        v_0 = v
        w_0 = w

        t = t + dt

    END DO

    ! write out each velocity component
    IF (cnvrged) THEN
        PRINT *, 'Solver converged'
    ELSE
        PRINT *, 'Solver did not converge'
    END IF

    ! solve for energy based on the velocity field

    ! Write out velocity fields
    out_file_u = TRIM(out_file) // 'u'
    out_file_v = TRIM(out_file) // 'v'
    !out_file_w = TRIM(out_file) // 'w'
    
    CALL WRITE2D_OUTPUT_MAIN(BC_N_u, BC_E_u, BC_S_u, BC_W_u, u(:,2:Nx,1), out_file_u)
    CALL WRITE2D_OUTPUT_MAIN(BC_N_v, BC_E_v, BC_S_v, BC_W_v, v(2:Ny,:,1), out_file_v)
    !CALL WRITE3D_OUTPUT_MAIN(BC_N_w, BC_E_w, BC_S_w, BC_W_w, BC_B_w, BC_T_w, w(:,:,2:Nz), out_file_w)


END PROGRAM MAIN