MODULE READ_INPUT
    ! Module that reads input values from file
    ! From the input file, the format should be
    ! '! <comments>' preface comments by a !
    ! '<var>=<number>' number = int or float
    ! directory to the output file should be labeled output=<dir> and should be last
    IMPLICIT NONE

    CONTAINS

        SUBROUTINE READ_INPUT_MAIN(params, in_file, out_file)
            ! INPUTS
            CHARACTER(256)  :: in_file
            ! PRECISION
            INTEGER, PARAMETER :: dp = KIND(1.0D0) ! double point
            ! DUMMYS
            INTEGER :: ios
            CHARACTER(256) :: string ! name of the parameter
            REAL(dp) :: num     ! numerical value assigned to parameter
            INTEGER :: dim  ! Number of parameters (dim of params)
            INTEGER :: idx
            ! OUTPUTS
            REAL(dp) , ALLOCATABLE, DIMENSION(:)    :: params
            CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT)  :: out_file ! read directory to output
           
            OPEN(unit=9, file=in_file, status="old", action="read", IOSTAT = ios)
            
            dim = 0
            DO ! Count the number of parameters
                READ(9, '(A)', iostat = ios) string
                IF (string(1:1) == '!') THEN 
                    CYCLE
                ELSEIF (string(1:6) == 'Output') THEN ! store output data directory
                    CALL DIR_PARSE(string, out_file)
                    EXIT
                ELSEIF (ios /= 0) THEN 
                    EXIT
                ELSE
                    dim = dim + 1
                END IF
               
            END DO
            

            CLOSE(unit=9)

            ALLOCATE(params(dim))

            OPEN(unit=9, file=in_file, status="old", action="read", IOSTAT = ios)

            idx = 1
            DO  ! read input line by line and parse parameter string for its numerical value
                READ(9, *, iostat = ios) string
                IF (string(1:1) == '!') THEN 
                    CYCLE
                ELSEIF (ios /= 0 .OR. string(1:6) == 'Output') THEN 
                    EXIT
                ELSE
                    CALL STRING_PARSE(string, num)
                   
                END IF

                params(idx) = num   ! store the parameter value in array
                idx = idx + 1
            END DO
           
            CLOSE(unit=9)

            RETURN

        END SUBROUTINE READ_INPUT_MAIN

        SUBROUTINE STRING_PARSE(string, num)
            ! PRECISION
            INTEGER, PARAMETER :: dp = KIND(1.0D0) ! double point
            ! INPUTS
            CHARACTER(256) :: string
            ! DUMMYS
            CHARACTER(1) :: delim = '='        
            INTEGER :: idx
            CHARACTER(16) :: num_str
            INTEGER :: io_status
            ! OUTPUTS
            REAL(dp) :: num

            idx = SCAN(string, delim)           ! Find numerical value after deliminator
            
            num_str = string(idx+1:LEN(string))

            READ(num_str, *, iostat=io_status) num ! Convert string to dp float

            RETURN
        
        END SUBROUTINE STRING_PARSE

        SUBROUTINE DIR_PARSE(string, out_file)
            ! scan the output=<> declaration for the directory to output data file

            ! PRECISION
            INTEGER, PARAMETER :: dp = KIND(1.0D0) ! double point
            ! INPUTS
            CHARACTER(256) :: string
          
            CHARACTER(1) :: delim = '='        
            INTEGER :: idx

            ! OUTPUTS
            CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT)  :: out_file

            idx = SCAN(string, delim)           ! Find numerical value after deliminator
            
            ALLOCATE(CHARACTER(LEN=LEN(string) - idx) :: out_file)
            out_file = string(idx+1:LEN(string))
    

            !WRITE(string(idx+1:LEN(string)), *) out_file ! Convert string to dp float

            ! RETURN
        
        END SUBROUTINE DIR_PARSE

END MODULE READ_INPUT