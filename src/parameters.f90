MODULE parameters
  !---------------------------------------------------------------------------!
  ! MODULE parameters
  !
  ! This module contains the default run case as well as a subroutine to read
  ! in the input file given at ./input/ibfs.inp.
  !---------------------------------------------------------------------------!
  
  IMPLICIT NONE
  
  ! PARAMETERS (AND DEFAULT VALUES)
  INTEGER :: istart = 0									! initial time index
  INTEGER :: istop = 400                ! last time index
  INTEGER :: isave = 25                 ! save a restart every isave steps
  INTEGER :: m  = 100                   ! cells in x
  INTEGER :: n  = 100                   ! cells in y
  REAL(KIND(0.0D0)) :: dt = 0.005d0     ! time step
  REAL(KIND(0.0D0)) :: Re = 300.d0      ! Reynolds number
  REAL(KIND(0.0D0)) :: cgtol = 1.D-7    ! tol. for cg convergence (poission eq)
  INTEGER :: cg_max_iter=3000           ! max. iterations for any cg iteration
  INTEGER :: n_act = 0                  ! number of actuators 
  INTEGER :: n_act_pt = 0               ! total number of points for actuators
                                        ! (will be counted in grid.f90) 
  REAL(KIND(0.D0)) :: len = 2.d0				! length scale for grid
  REAL(KIND(0.D0)) :: offsetx = 1.d0    ! offset for grid in x
  REAL(KIND(0.D0)) :: offsety = 1.d0    ! offset for grid in y
  INTEGER :: mgridlev = 4								! number of MGRID levels

  LOGICAL :: unif_flow = .TRUE.         ! add a uniform flow
  LOGICAL :: make_vortex = .FALSE.      ! initialize a vortex for ic
  LOGICAL :: stationary = .FALSE.       ! if stationary body, TRUE
  
  ! NON-INPUT FILE VALUES
  REAL(KIND(0.D0)) :: cputime           ! computer run time
  REAL(KIND(0.D0)) :: pi = 3.1415926535897932384d0

CONTAINS
  
  SUBROUTINE input
    !-------------------------------------------------------------------------!
    ! This subroutine reads the input file given at ./input/ibfs.inp and if 
    ! the input file does not exist, then the default case given in the 
    ! parameters above are substituted and an input file is created.
    !-------------------------------------------------------------------------!
    
    LOGICAL :: readinput 
    
    NAMELIST /read_parameters/ istart,istop,isave,m,n,dt,Re,cgtol, &
                                cg_max_iter,unif_flow,len, offsetx, &
                                offsety,mgridlev,make_vortex,stationary

    !! READ INPUT
    INQUIRE(file='input/ibfs.inp',exist=readinput)

    IF (readinput) THEN ! If exist, then read input file
      OPEN(unit=3,file='input/ibfs.inp',form='formatted',status='old')
      READ(unit=3,nml=read_parameters)
      CLOSE(3)
      WRITE(111,*) "... Read parameters from ibfs.inp"
      WRITE(111,*) "... Input parameters follow"
      WRITE(unit=111,nml=read_parameters)
    ELSE ! Else, use the default case
      OPEN(unit=3,file='input/ibfs.inp',form='formatted',status='new')
      WRITE(unit=3,nml=read_parameters)
      CLOSE(3)
      WRITE(111,*) "... loading default run"
      WRITE(unit=111,nml=read_parameters)
    END IF

  END SUBROUTINE input
  
END MODULE parameters


