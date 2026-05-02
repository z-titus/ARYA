MODULE WRITE_OUTPUT
    ! Module that takes in a 2D or 3D output from a simulation
    ! and writes appropriate data files
    IMPLICIT NONE

    CONTAINS

        SUBROUTINE WRITE2D_OUTPUT_MAIN(bN, bE, bS, bW, out_data, out_file)

            ! PRECISION
            INTEGER, PARAMETER :: dp = KIND(1.0D0) ! double point

            ! INPUTS
            CHARACTER(LEN=:), AllOCATABLE, INTENT(IN) ::   out_file
        
            REAL(dp), DIMENSION(:, :)  ::   out_data
            LOGICAL :: exists
            INTEGER :: io, stat
            INTEGER :: i, j, k, nrows, ncols, ydim, xdim   ! num of rows and columns in 2D data array

            REAL(dp), INTENT(IN)    :: bN, bE, bS, bW

            nrows = SIZE(out_data, 1)+2
            ncols = SIZE(out_data, 2)+2

            ydim = SIZE(out_data, 1)
            xdim = SIZE(out_data, 2)

            ! Delete file if it already exists
            INQUIRE(file=out_file, exist=exists)

            IF (exists) THEN
                OPEN(file=out_file, newunit=io, iostat=stat)
                IF (stat == 0) CLOSE(io, status="delete", iostat=stat)
            END IF
            
            ! Write to file now
            OPEN(newunit=io, file=out_file, status="new", action="write", iostat=stat)

            PRINT *, stat

            ! write a row of south
            DO i = 1, ncols
                WRITE(io,*) bS
            END DO

            DO k = 1, ydim

                ! write west
                WRITE(io,*) bW

               
                    DO j = 1, xdim
                        
                        WRITE(io,*) out_data(k,j)
                        
                    END DO
                

                ! write east
                WRITE(io,*) bE  
            END DO
            
            ! write a row of north
            DO i = 1, ncols
                WRITE(io,*) bN
            END DO
            
            CLOSE(io)

        END SUBROUTINE WRITE2D_OUTPUT_MAIN

        SUBROUTINE WRITE3D_OUTPUT_MAIN(bN, bE, bS, bW, bB, bT, out_data, out_file) ! ADD BC arguments

            ! PRECISION
            INTEGER, PARAMETER :: dp = KIND(1.0D0) ! double point

            ! INPUTS
            CHARACTER(LEN=:), AllOCATABLE, INTENT(IN) ::   out_file

            REAL(dp), INTENT(IN)    :: bN, bE, bS, bW, bB, bT
        
            REAL(dp), DIMENSION(:, :, :)  ::   out_data
            LOGICAL :: exists
            INTEGER :: io, stat
            INTEGER :: i, j, k, l, nrows, ncols, nz, xdim, ydim, zdim  ! num of rows and columns in 3D data array

            nrows = SIZE(out_data, 1)+2  ! add two for boundaries
            ncols = SIZE(out_data, 2)+2
            nz    = SIZE(out_data, 3)+2

            ydim = SIZE(out_data, 1)
            xdim = SIZE(out_data, 2)
            zdim = SIZE(out_data, 3)

            ! Delete file if it already exists
            INQUIRE(file=out_file, exist=exists)

            IF (exists) THEN
                OPEN(file=out_file, newunit=io, iostat=stat)
                IF (stat == 0) CLOSE(io, status="delete", iostat=stat)
            END IF
            
            ! Write to file now
            OPEN(newunit=io, file=out_file, status="new", action="write", iostat=stat)

            ! write a plane of bottom
            DO i = 1, nrows
                DO j = 1, ncols
                    WRITE(io,*) bB
                END DO  
            END DO

            ! write internal to BC phi solution
            DO k = 1,zdim
                ! write a row of south
                DO l = 1, ncols
                    WRITE(io,*) bS
                END DO

                DO i = 1, ydim
                    ! write west
                    WRITE(io,*) bW

                    DO j = 1, xdim
                        
                        WRITE(io,*) out_data(i,j,k)
                        
                    END DO
                    ! write east
                    WRITE(io,*) bE
                END DO

                ! write a row of north
                DO l = 1, ncols
                    WRITE(io,*) bN
                END DO
                
            END DO

            ! write a plane of top
            DO i = 1, nrows
                DO j = 1, ncols
                    WRITE(io,*) bT
                END DO
            END DO

            CLOSE(io)

        END SUBROUTINE WRITE3D_OUTPUT_MAIN

END MODULE WRITE_OUTPUT