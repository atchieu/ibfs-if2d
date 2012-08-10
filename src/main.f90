PROGRAM main
  !---------------------------------------------------------------------------!
  ! (c) Kunihiko (Sam) Taira and Tim Colonius
  !
  ! Immersed Boundary Fractional Step Method w/ MGRID
  !    - 2D, Incompressible Navier-Stokes
  !    - 2nd-order-accurate except for near IB
  !    - 1st-order-accurate IB
  !    - Requires uniform grid using multi-grid method
  !---------------------------------------------------------------------------!

  USE parameters
  USE grid
  USE variables
  USE operators

  IMPLICIT NONE
  INTEGER :: it
  CHARACTER(10):: date, time

  ! START INFO FILE
  OPEN(unit=111,file='output/ibfs.inf',status='unknown',position='append', &
           form='formatted')
  WRITE(111,*) "==================================="
  WRITE(111,*) "|   Welcome to IBFS MGRID v1.0    |"
  WRITE(111,*) "==================================="

  ! TIME AND DATE
  CALL DATE_AND_TIME(date,time) 
  WRITE(111,*) "Executed at time:",time," on date:",date
  
  ! WRITE TO SCREEN
  WRITE(*,*) " "
  WRITE(*,*) "======================================="
  WRITE(*,*) "|     Welcome to IBFS MGRID v1.0      |"
  WRITE(*,*) "| (c) Kunihiko Taira and Tim Colonius |"
  WRITE(*,*) "======================================="
  WRITE(*,*) " "
  WRITE(*,*) "Executed at time: ",time," on date: ",date

  ! SETUP SIMULATION
  CALL input
  CALL setup
  CALL setup_variables
  
  ! BEGIN TIME MARCH
  it = istart
  DO WHILE (it < istop)
     CALL advance(it)
  END DO
  
  ! CLEAN
  CALL destroy_variables
  CALL destroy_grid
  
  ! TOTAL CPU RUNTIME
  CALL CPU_TIME(cputime)
  WRITE(*,*) " "
  WRITE(*,*) "Total run time: ",cputime," s"
  WRITE(111,*) "Total run time: ",cputime," s"
  WRITE(*,*) " " 
  
END PROGRAM main
