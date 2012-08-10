MODULE variables
  !---------------------------------------------------------------------------!
  ! This module deals with setting up and initializing variables as well as
  ! reading and writing variables to output files. It also reads and writes
  ! Cholesky decompositions. 
  !---------------------------------------------------------------------------!
  USE grid
  IMPLICIT NONE

  ! In what follows, the last index refers to the grid level, first 1 (or 2) 
  ! indices to point in space
  
  REAL(KIND(0.D0)), DIMENSION(:,:,:), ALLOCATABLE :: omega    ! vorticity 
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: fb           ! forces
  REAL(KIND(0.D0)), DIMENSION(:,:,:), ALLOCATABLE :: rhs_old  ! nonlinear terms at previous time-step (for AB2 integration)
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: q, q0      ! fluxes, q0 is the underlying potential flow flux
  
  ! NOTE: Additional potential flow to be added
  ! (currrently uniform flow if unif_flow=T)

  ! Variables for Cholesky (stationary body)
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: cholmat
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: cholvec

CONTAINS
  
!=============================================================================!
  SUBROUTINE setup_variables
    !-------------------------------------------------------------------------!
    ! This subroutine first allocates space for all the pertinent grid 
    ! quantities. If istart = 0, then it will initialize a zero flow 
    ! condition (or uniform flow if appropriate). Otherwise, it will search
    ! for a restart and read the appropriate variables and setup the 
    ! regularization operator again. If the body is stationary, then a
    ! Cholesky factorization is used to speed up the solution to the modified
    ! Poisson's equation for the boundary forces.
    !-------------------------------------------------------------------------!
    USE parameters
    USE myfft
    INTEGER :: i,j,k
    
    ! SETUP THE FFT
    CALL setup_fft

    ! ALLOCATE ENOUGH SPACE FOR ALL QUANTITIES
    ALLOCATE( omega(2:m,2:n,mgridlev), fb(Nf), rhs_old(2:m,2:n,mgridlev) )

    ALLOCATE( q(Nq,mgridlev), q0(Nq,mgridlev) )

    IF (stationary) THEN
       ALLOCATE( cholmat(Nf,Nf), cholvec(Nf) )
    END IF
    
    ! IF ISTART = 0, THEN INITIAL CONDITION IS ZERO FLOW
    IF (istart==0) THEN
       q = 0.d0
       fb = 0.d0
       omega = 0.d0
       rhs_old = 0.d0
       IF (unif_flow) THEN ! OR UNIFORM FLOW
          DO k=1,mgridlev
             q0(1:(m+1)*n,k) = delta * 2.d0**(k-1)
             q0((m+1)*n+1:Nq,k) = 0.d0
          END DO
       ELSE
          q0 = 0.d0
       END IF
    ! ELSE READ THE APPROPRIATE VARIABLE   
    ELSE
       CALL read_variables(istart)
       CALL setup_reg
       IF (stationary) THEN
          CALL read_cholesky
       END IF
    END IF

  END SUBROUTINE setup_variables

!=============================================================================!
  SUBROUTINE destroy_variables
    !-------------------------------------------------------------------------!
    ! Deallocates space!!
    !-------------------------------------------------------------------------!
    USE parameters
    USE myfft

    CALL destroy_fft
    DEALLOCATE( omega, fb, rhs_old, q, q0)
    IF (stationary) THEN
       DEALLOCATE( cholmat, cholvec )
    END IF

  END SUBROUTINE destroy_variables

!=============================================================================!
  SUBROUTINE destroy_grid
    !-------------------------------------------------------------------------!
    ! Deallocate the variables assocaited with the grid including index 
    ! pointers and smear variables.
    !-------------------------------------------------------------------------!
    USE parameters
    INTEGER :: i

    DEALLOCATE( f, u, v )
    DEALLOCATE( xb, yb, bdyno, smear, ismear, xb_og, yb_og )

  END SUBROUTINE destroy_grid

!=============================================================================!
  SUBROUTINE write_cholesky
    !-------------------------------------------------------------------------!
    ! Writes the computed Cholesky decomposition to input/ibfs.chd.
    !-------------------------------------------------------------------------!
    USE parameters

    OPEN(unit=100,file="input/ibfs.chd",form="unformatted",status="unknown")
    WRITE(100) cholmat, cholvec
    CLOSE(100)
    WRITE(111,*) "Saving cholesky decomposition..."

  END SUBROUTINE write_cholesky
  
!=============================================================================!
  SUBROUTINE read_cholesky
    !-------------------------------------------------------------------------!
    ! Reads a premade Cholesky factorization given in ibfs.chd.
    !-------------------------------------------------------------------------!
    USE parameters
    USE grid

    OPEN(unit=100,file="input/ibfs.chd",form="unformatted",status="unknown")
    READ(100) cholmat, cholvec
    CLOSE(100)
    WRITE(111,*) "Reading Cholesky decompopsition..."
    WRITE(*,*) "Reading Cholesky premade decomposition..."

  END SUBROUTINE read_cholesky
  

!=============================================================================!
  SUBROUTINE write_variables(it)
    !-------------------------------------------------------------------------!
    ! Writes restart variables to a file. These files can be post processed to
    ! give tecplot amenable ASCII data.
    !-------------------------------------------------------------------------!
    USE parameters
    CHARACTER(7) :: charit
    INTEGER :: it

    WRITE(charit,"(I7.7)") it
    OPEN(unit=100,file="output/ibfs"//charit//".var",form="unformatted",status="unknown")
    WRITE(100) omega,fb,rhs_old, q, q0
    CLOSE(100)
    WRITE(111,*) "Writing restart file at itime=",it

  END SUBROUTINE write_variables
  
!=============================================================================!
  SUBROUTINE read_variables(it)
    !-------------------------------------------------------------------------!
    ! This subroutine reads variables from time step it. It is mainly used to
    ! restart a simulation.
    !-------------------------------------------------------------------------!
    USE parameters
    CHARACTER(7) :: charit
    INTEGER :: it

    WRITE(charit,"(I7.7)") it
    OPEN(unit=100,file="output/ibfs"//charit//".var",form="unformatted",status="unknown")
    READ(100) omega, fb, rhs_old, q, q0
    CLOSE(100)
    WRITE(111,*) "Reading restart file at itime=",it
    WRITE(*,*) "Reading restart file at itime =",it

  END SUBROUTINE read_variables

END MODULE variables
