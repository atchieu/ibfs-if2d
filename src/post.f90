PROGRAM post
  !---------------------------------------------------------------------------!
  ! Post processing driver program. Writes files to Tecplot readable ASCII
  ! files for all grids.
  !---------------------------------------------------------------------------!
  USE parameters
  USE grid
  USE variables
  USE operators

  IMPLICIT NONE

  INTEGER :: it,k
  CHARACTER(10):: date, time

  ! WRITE GENERAL RUN INFO TO FILE
  OPEN(unit=111,file='output/ibfs.inf',status='unknown',position='append', &
           form='formatted')
  WRITE(111,*) "================================="
  WRITE(111,*) "|   Welcome to IBFSpost v1.0    |"
  WRITE(111,*) "================================="

  CALL date_and_time(date,time)
  WRITE(111,*) "Executed at time:",time," on date:",date

  CALL input
  CALL setup

  it = istart
 
  DO WHILE (it <= istop )
    IF ((MOD(it,isave).eq.0).or.(it==istop)) THEN
       istart = it
       WRITE(*,*) "Post-processing: itime = ",it 
       CALL setup_variables
       DO k=1,mgridlev
          CALL write_tecplot(it,k)
       END DO
       CALL destroy_variables
    END IF
    it=it+1
 END DO
 CALL destroy_grid

END PROGRAM post
