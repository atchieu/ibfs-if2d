MODULE grid
  !---------------------------------------------------------------------------!
  ! This module sets up the grid, sets up the body, allows motion of the body,
  ! and also contains some information about actuation.
  !---------------------------------------------------------------------------!
  
  USE parameters
  IMPLICIT NONE

  ! VARIABLES FOR IB GEOMETRY AND ITS MOTION (IF ANY)
  REAL(KIND(0.0D0)) :: support = 3.0D0  ! support for smearing delta functions
  REAL(KIND(0.D0)) :: delta             ! near field grid spacing
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: xb,yb ! body coordinates, fixed frame
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: xb_og, yb_og ! original body coordinates, this will NEVER change once read
  INTEGER, DIMENSION(:), ALLOCATABLE :: bdyno ! body number used for assigning motion
  
  ! ARRAYS FOR SMEARING COEFFICIENTS
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE :: smear
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ismear
  INTEGER :: nsmear
  
  ! NUMBER OF CELLS, EDGES, ETC.
  INTEGER :: nq, nb, nf, ns

  ! INTEGER POINTER FOR FIELD VARIABLES 
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: f ! f(i,k) gives kth-component of force at ith point on body
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: u ! u(i,j) gives point in the q array where u at face of cell i,j lives
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: v ! v(i,j) gives point in the q array where v at face of cell i,j lives

  ! FOR BCs
  INTEGER :: top,bottom,left,right

  ! A SPECIAL TYPE THAT IS USED IN ASSEMBLING GEOMETRIES, USED BY SETUP_GEOMETRY
  INTEGER, PARAMETER :: maxpts = 10000 ! a large integer
  TYPE bodylist
     INTEGER :: ntot=0
     REAL(KIND(0.D0)), DIMENSION(maxpts) :: xbody,ybody
     INTEGER, DIMENSION(maxpts) :: bodyno
  END TYPE bodylist

CONTAINS

!=============================================================================!
  SUBROUTINE setup
    !-------------------------------------------------------------------------!
    ! Initial setup of the simulation. This subroutine does:
    ! 
    ! SETUP GRID - Checks that the grid is divisible by 4, determines the 
    !   indicies of the boundaries, determines max number of grid levels
    !   possible
    !
    ! SETUP BODY - Calls setup_geometry to get body and put it in "geom," 
    !   allocate and assign position of body to xb, yb, and bdyno
    !
    ! SETUP INDEXING - Determine the indecies corresponding to the
    ! streamfuction, forces, and velocity
    !
    ! SETUP FORCE REGULARIZATION OPERATORS - See subroutine
    !-------------------------------------------------------------------------!
    INTEGER :: i,j,next,nlevmax
    TYPE( bodylist ) :: geom

    !!! SETUP GRID !!!
    delta = len/REAL(m) 

    ! FOR MGRID, CHECK THAT GRID IS DIVISIBLE BY 4
    IF ((MOD(m,4).ne.0).or.(MOD(n,4).ne.0)) THEN
       STOP 'HEY! Grid nodes must be divisible by 4!!'
    END IF
    
    ! FOR MGRID BCs
    ! Get the indecies for the boundaries
    left = 0
    right = n+1
    bottom = 2*(n+1)
    top = 2*(n+1) + m+1
    
    ! NUMBER OF GRID LEVELS, MAX AND USED
    nlevmax = FLOOR(  LOG( REAL( MIN(m,n) /2)  ) / LOG( 2.d0) ) - 1
    WRITE(*,*) 'Maximum:',nlevmax,'levels in MGRID'
    mgridlev = MIN(mgridlev,nlevmax)
    WRITE(*,*) 'Using:',mgridlev,'levels in MGRID'

    !!! SETUP BODY !!!
    CALL setup_geometry( delta, geom ) ! CALL and Put everything into "geom" and get "delta"
    
    nb = geom%ntot
    ALLOCATE( xb(nb),  yb(nb), bdyno(nb)  ) ! Allocate space for the Lagrangian points for each body  
    ALLOCATE( xb_og(nb), yb_og(nb) ) ! Allocate space for this original points
    DO i=1,nb ! for each body
       xb(i) = geom%xbody(i)
       yb(i) = geom%ybody(i)
       bdyno(i) = geom%bodyno(i)
    END DO
    
    xb_og = xb
    yb_og = yb

    !!! SETUP INDEXING !!!

    ! INDEXING FOR STREAMFUNCTION
    ns = (m-1)*(n-1) 

    ! INDEXING FOR FORCES
    nf  = 2*nb ! number of forces
    ALLOCATE( f(1:nb,1:2)  )      
    next = 0
    DO i=1,nb
       next = next + 1
       f(i,1) = next
    END DO
    DO i = 1,nb
       next = next + 1
       f(i,2) = next
    END DO
    IF (next.ne.nf) STOP "Error in setup_parms - f (forces)..."

    ! INDEXING FOR VELOCITIES (FLUX COMPONENTS)
    nq = (m+1)*n + (n+1)*m
    ALLOCATE( u(1:m+1,1:n),  v(1:m,1:n+1) )
    next = 0
    DO j=1,n ! I think this is poor mans way of keeping indecies so you dont have to deal with them later...
       DO i=1,m+1
          next = next+1
          u(i,j) = next
       END DO
    END DO
    DO j=1,n+1
       DO i=1,m
          next = next+1
          v(i,j) = next
       END DO
    END DO
    IF (next.ne.nq) STOP "Error in setup_parms - u (velocities)..."

    !!! SETUP FORCE REGULARIZATION OPERATORS !!!
    CALL setup_reg
    WRITE(111,*) "Set up force regularization operators..."
  END SUBROUTINE setup

!=============================================================================!
  SUBROUTINE setup_reg
    !-------------------------------------------------------------------------!
    ! This subroutine sets up the regularization operator.
    !-------------------------------------------------------------------------!
    INTEGER :: ijsmear
    REAL(KIND(0.D0)) :: x,y,d2,sup
    INTEGER :: i,j,k

    ! CHECK TO SEE IF ISMEAR AND SMEAR ARE ALREADY ALLOCATED
    IF (ALLOCATED(ismear)) THEN
       DEALLOCATE( ismear, smear )
    END IF

    d2 = 0.5d0*delta
    sup = support*d2

    ! FIRST, COUNT UP HOW MANY BOUNDARY POINTS SMEAR OVER A GIVEN INTERIOR POINT
    Nsmear = 0
    DO i=1,m+1
       DO j=1,n
          x = delta*REAL(i-1)-offsetx
          y = delta*REAL(j-1)-offsety+d2
          ijsmear = 0
          DO k=1,Nb
             ! x-directional vector
             IF ((ABS(x-xb(k)) < sup).AND.(ABS(y-yb(k)) < sup)) THEN
                ijsmear = ijsmear + 1
             END IF
          END DO
          Nsmear = MAX(Nsmear,ijsmear)
       END DO
    END DO
    DO i=1,m
       DO j=1,n+1
          x = delta*REAL(i-1)-offsetx+d2
          y = delta*REAL(j-1)-offsety
          ijsmear = 0
          DO k=1,Nb
             ! y-directional vector
             IF ((ABS(x-xb(k)) < sup).AND.(ABS(y-yb(k)) < sup)) THEN
                ijsmear = ijsmear + 1
              END IF
          END DO
          Nsmear = MAX(Nsmear,ijsmear)
       END DO
    END DO

    ALLOCATE( ismear(Nq,Nsmear), smear(Nq,Nsmear) )

    ! NOW GET THE POINTS OVER WHICH THE FORCE IS SMEARED
    ! Default is to set the coefficient to 0 and the influence point to 1.  
    ! This avoids having to have any logic in the routines reg and regT.
    ismear = 1
    smear = 0.d0

    DO i=1,m+1
       DO j=1,n
          x = delta*REAL(i-1)-offsetx
          y = delta*REAL(j-1)-offsety+d2
          ijsmear = 0
          DO k=1,Nb
             ! x-directional vector
             IF ((ABS(x-xb(k)) < sup).AND.(ABS(y-yb(k)) < sup)) THEN
                ijsmear = ijsmear + 1
                ismear(u(i,j),ijsmear) = f(k,1)
                smear(u(i,j),ijsmear) =  (delta*delta)* deltafnc(x, xb(k), delta) & 
                                                      * deltafnc(y, yb(k), delta) 
             END IF
          END DO
       END DO
    END DO
    DO i=1,m
       DO j=1,n+1
          x = delta*REAL(i-1)-offsetx + d2
          y = delta*REAL(j-1)-offsety
          ijsmear = 0
          DO k=1,Nb
             ! y-directional vector
             IF ((ABS(x-xb(k))< sup).AND.(ABS(y-yb(k))< sup)) THEN
                ijsmear = ijsmear+1
                ismear(v(i,j),ijsmear) = f(k,2)
                smear(v(i,j),ijsmear) = (delta*delta)* deltafnc(x, xb(k), delta) &
                                                     * deltafnc(y, yb(k), delta) 
             END IF
          END DO
       END DO
    END DO

  END SUBROUTINE setup_reg

!=============================================================================!
  FUNCTION deltafnc( r,r0, dr ) 
    !-------------------------------------------------------------------------!
    ! This function outputs the value of the delta function given the grid
    ! spacing, dr. r and r0 are just some distances.
    !-------------------------------------------------------------------------! 
    REAL(KIND(0.D0)) :: r,r0,dr,deltafnc

    ! FROM Roma, Peskin, & Berger (JCP 1999)
    IF (ABS(r-r0)<=(0.5D0*dr)) THEN
       deltafnc = (1.D0+SQRT(-3.D0*((r-r0)/dr)**2+1.D0))/(3.0D0*dr)
    ELSEIF (ABS(r-r0)<=(1.5D0*dr)) THEN
       deltafnc = (5.D0-3.D0*ABS((r-r0)/dr)-SQRT(-3.D0*(1.D0-ABS((r-r0)/dr))**2+1.D0))/(6.D0*dr)
    ELSE 
       deltafnc = 0.D0
    END IF

  END FUNCTION deltafnc

!=============================================================================!
  SUBROUTINE setup_geometry( ds, geom )
    !-------------------------------------------------------------------------!
    ! This subroutine sets up the geometry used for the simulation. It will
    ! either read a list of Lagrangian points in the file geom.inp OR create a
    ! rudimentary geometry using some primative shapes (functions given below)
    !-------------------------------------------------------------------------!
    LOGICAL :: readinput, readact
    TYPE(bodylist) :: geom
    REAL(KIND(0.D0)) :: ds
    INTEGER :: i, i_act, n_pt
    CHARACTER(3) :: act_num
    REAL(KIND(0.D0)) ::  x_j,y_j,h_j,u_j,u_j_osc,theta_j,omg_j
    NAMELIST /read_act_parameters/ x_j,y_j,h_j,u_j,u_j_osc,theta_j,omg_j

    ! CHECK IF FILE EXIST
    INQUIRE(file='input/geom.inp',exist=readinput)
    
    ! IF FILE EXIST, READ INPUT
    IF (readinput) THEN
       OPEN(unit=8,file='input/geom.inp',form='formatted',status='old')
       READ(8,*) geom%ntot
       DO i=1,geom%ntot
          READ(8,*) geom%xbody(i), geom%ybody(i), geom%bodyno(i)
       END DO
       CLOSE(8)
       WRITE(111,*) "Reading geometry from geom.inp"
       WRITE(*,*) "Reading geometry from geom.inp"
       
       !----------------------------------------------------------------------!
       ! NOTE: There should really be some check of the supplied geom.inp
       ! for example to insure that the point spacing is about equal to ds.  
       ! For now, we just assume it is okay.
       !----------------------------------------------------------------------!
        
       ! IF BODYNO < 0, THEN ACTUATORS ARE ALREADY DEFINED
       IF (ANY(geom%bodyno<0)) THEN 

          DO i_act = 1,100 ! max no of actuator = 100.
             WRITE(act_num,"(I3.3)") i_act
             INQUIRE(file="input/actuator."//act_num//".inp",exist=readact)

             IF (readact) THEN
                OPEN(unit=9,file="input/slot."//act_num//".inp",form='formatted',status='old')
                READ(9,*) n_pt 
                CLOSE(9)
                n_act    = n_act + 1  ! counting number of actuators
                n_act_pt = n_act_pt + n_pt 
             END IF

          END DO
          
       ! ELSE, CHECK IF WE NEED TO APPEND ACTUATORS AND DO IT   
       ELSE 
          DO i_act = 1,100 ! max no of actuator = 100.

             WRITE(act_num,"(I3.3)") i_act
             INQUIRE(file="input/actuator."//act_num//".inp",exist=readact)
             
             IF (readact) THEN
                n_act = n_act + 1  ! counting number of actuators
                OPEN(unit=7,file="input/actuator."//act_num//".inp",form='formatted',status='old')
                READ(unit=7,nml=read_act_parameters)
                CLOSE(7)
                geom = actuator_slot( x_j, y_j, h_j, ds, theta_j, geom, -i_act )
             END IF
             
          END DO
          OPEN(unit=8,file='input/geom.inp',form='formatted',status='replace')
          WRITE(8,*) geom%ntot
          DO i=1,geom%ntot
             WRITE(8,*) geom%xbody(i), geom%ybody(i), geom%bodyno(i)
          END DO
          CLOSE(8)
          WRITE(111,*) "...appended actuator info to geom.inp" 
       END IF
       
    ! ELSE CREATE A GEOMETRY FROM THE PRIMATIVES GIVEN
    ! All routines append the geometry onto the list geom. The final argument
    ! in each is the body number for a given piece of geometry. These
    ! correspond to different equations of motion for different bodies if the
    ! body will be stationary, use body number 0.
    
    ELSE
!      ! EXAMPLE 1: EQUILATERAL TRIANGLE WITH BASE LENGTH 1
!      geom = vertex(-0.5d0, 0.d0, geom, 0 )
!      geom = segment( ds, -0.5d0, 0.d0, 0.5d0, 0.577350269d0, geom, 0 )
!      geom = vertex( 0.5d0,  0.577350269d0, geom, 0 )
!      geom = segment( ds, 0.5d0, 0.577350269d0, 0.5d0, -0.577350269d0, geom, 0 )
!      geom = vertex( 0.5d0,  -0.577350269d0, geom, 0 )
!      geom = segment( ds, 0.5d0,  -0.577350269d0, -0.5d0, 0.d0, geom, 0 )
 
!      geom = vertex( -0.5d0, 0.5d0, geom, 0) 
!      geom = segment( ds, -0.5d0, 0.5d0, 0.5d0, 0.5d0, geom, 0) 
!      geom = vertex( 0.5d0, 0.5d0, geom, 0) 
!      geom = segment( ds, 0.5d0, 0.5d0, 0.5d0, -0.5d0, geom, 0) 
!      geom = vertex( 0.5d0, -0.5d0, geom, 0) 
!      geom = segment( ds, 0.5d0, -0.5d0, -0.5d0, -0.5d0, geom, 0) 
!      geom = vertex( -0.5d0, -0.5d0, geom, 0) 
!      geom = segment( ds, -0.5d0, -0.5d0, -0.5d0, 0.5d0, geom, 0) 

      ! EXAMPLE 2: A CIRCLE OF DIAMETER 1 AT THE ORIGIN
      geom = circle( ds, 1.d0, 0.0d0, 0.d0, geom, 0)

      ! EXAMPLE 3: A LINE SEGEMENT (INCLUDING END POINTS): LENGTH 1, 30 DEGREE AOA
!      geom = line( ds, 1.d0, 10.d0, 0.5d0, 0.d0, geom, 0 )
       
!      ! EXAMPLE 4:A LINE SEGEMENT (BASED ON SPECIFYING ITS END POINTS)
!      geom = segment( ds, 0.0d0, -0.0d0, 1.0d0, 0.0d0, geom, 0 )

!      ! EXAMPLE 5: A SQUARE WITH SIDE-LENGTH 1
!      geom = vertex( -0.5d0, 0.5d0, geom, 0) 
!      geom = segment( ds, -0.5d0, 0.5d0, 0.5d0, 0.5d0, geom, 0) 
!      geom = vertex( 0.5d0, 0.5d0, geom, 0) 
!      geom = segment( ds, 0.5d0, 0.5d0, 0.5d0, -0.5d0, geom, 0) 
!      geom = vertex( 0.5d0, -0.5d0, geom, 0) 
!      geom = segment( ds, 0.5d0, -0.5d0, -0.5d0, -0.5d0, geom, 0) 
!      geom = vertex( -0.5d0, -0.5d0, geom, 0) 
!      geom = segment( ds, -0.5d0, -0.5d0, -0.5d0, 0.5d0, geom, 0) 
       
!       ! EXAMPLE 6: A CIRCLE "FILLED IN"
!       geom = circle( ds, 1.00D0, 0.D0, 0.D0, geom, 0)
!       geom = circle( ds, 0.90D0, 0.D0, 0.D0, geom, 0)
!       geom = circle( ds, 0.80D0, 0.D0, 0.D0, geom, 0)
!       geom = circle( ds, 0.60D0, 0.D0, 0.D0, geom, 0)
!       geom = circle( ds, 0.50D0, 0.D0, 0.D0, geom, 0)
!       geom = circle( ds, 0.40D0, 0.D0, 0.D0, geom, 0)
!       geom = circle( ds, 0.30D0, 0.D0, 0.D0, geom, 0)
!       geom = circle( ds, 0.20D0, 0.D0, 0.D0, geom, 0)
!       geom = circle( ds, 0.10D0, 0.D0, 0.D0, geom, 0)
                     
       ! CHECK IF WE NEED ACTUATION (I really have no clue about actuation)
       DO i_act = 1,100 ! max no of actuator = 100.
          WRITE(act_num,"(I3.3)") i_act
          INQUIRE(file="input/actuator."//act_num//".inp",exist=readact)

          IF (readact) THEN
             n_act = n_act + 1  ! counting number of actuators
             OPEN(unit=7,file="input/actuator."//act_num//".inp",form='formatted',status='old')
             READ(unit=7,nml=read_act_parameters)
             CLOSE(7)
             geom = actuator_slot( x_j, y_j, h_j, ds, theta_j, geom, -i_act )
          END IF
       END DO
       
       ! WRITE BODY TO GEOM.INP
       OPEN(unit=8,file='input/geom.inp',form='formatted',status='new')
       WRITE(8,*) geom%ntot
       DO i=1,geom%ntot
          WRITE(8,*) geom%xbody(i), geom%ybody(i), geom%bodyno(i)
       END DO
       CLOSE(8)
       WRITE(111,*) "Generated geometry from primatives."
       WRITE(*,*) "Generated geometry from primatives."
    END IF
  END SUBROUTINE setup_geometry
 
!=============================================================================!
  FUNCTION vertex( x, y, geom, num )
    !-------------------------------------------------------------------------!
    ! Create a vertex for the geometry with body number, num.
    !-------------------------------------------------------------------------!
    TYPE(bodylist) :: geom, vertex
    REAL(KIND(0.D0)) :: x,y
    INTEGER :: num,i

    DO i=1,geom%ntot
       vertex%xbody(i) = geom%xbody(i)
       vertex%ybody(i) = geom%ybody(i)
       vertex%bodyno(i) = geom%bodyno(i)
    END DO

    vertex%xbody(geom%ntot+1) = x
    vertex%ybody(geom%ntot+1) = y
    vertex%bodyno(geom%ntot+1) = num
    vertex%ntot = geom%ntot+1

  END FUNCTION vertex

!=============================================================================!
  FUNCTION line( ds, len, aoa, x,y , geom, num)
    !-------------------------------------------------------------------------!
    ! Create a line including end points for the geometry with body number, 
    ! num.
    !-------------------------------------------------------------------------!
    TYPE(bodylist) :: geom, line
    REAL(KIND(0.D0)) :: ds,x,y, len, aoa, angle, seg
    REAL(KIND(0.D0)) :: x_LE, y_LE
    INTEGER :: i
    INTEGER :: num, npt

    angle = -aoa*pi/180.D0
    x_LE  = -0.5*len*COS(angle)
    y_LE  = -0.5*len*SIN(angle)

    DO i=1,geom%ntot
       line%xbody(i) = geom%xbody(i)
       line%ybody(i) = geom%ybody(i)
       line%bodyno(i) = geom%bodyno(i)
    END DO

    ! Achtung... angle input in degrees
    npt = INT( len/ds )
    DO i=1,npt+1
       seg = REAL(i-1)*len/REAL(npt)
       line%xbody(geom%ntot+i) = x_LE + seg*COS(angle)
       line%ybody(geom%ntot+i) = y_LE + seg*SIN(angle)
       line%bodyno(geom%ntot+i) = num
    END DO
    line%ntot = geom%ntot+npt+1

  END FUNCTION line


!=============================================================================!
  FUNCTION segment( ds, x1, y1, x2, y2, geom, num )
    !-------------------------------------------------------------------------!
    ! Create a segment based on end points for the geometry with body number, 
    ! num.
    !-------------------------------------------------------------------------!
    TYPE(bodylist) :: geom, segment
    REAL(KIND(0.D0)) :: ds,x1,x2,y1,y2, len, angle, seg
    INTEGER :: i
    INTEGER :: num, npt

    DO i=1,geom%ntot
       segment%xbody(i) = geom%xbody(i)
       segment%ybody(i) = geom%ybody(i)
       segment%bodyno(i) = geom%bodyno(i)
    END DO

    len = SQRT( (x2-x1)**2 + (y2-y1)**2 )
    angle = ATAN2( (y2-y1), (x2-x1) )
    npt = INT( len/ds )
    DO i=1,npt-1
       seg = REAL(i)*len/REAL(npt)
       segment%xbody(geom%ntot+i) = x1 + seg * COS(angle)
       segment%ybody(geom%ntot+i) = y1 + seg * SIN(angle)
       segment%bodyno(geom%ntot+i) = num
    END DO
    segment%ntot = geom%ntot+npt-1

  END FUNCTION segment

!=============================================================================!
  FUNCTION circle( ds, dia, x, y, geom, num )
    !-------------------------------------------------------------------------!
    ! Create a cicle with diameter, dia, at location (x,y) with body number,
    ! num.
    !-------------------------------------------------------------------------!
    TYPE(bodylist) :: geom, circle
    REAL(KIND(0.D0)) :: ds, dia, x, y !, pi
    INTEGER :: i
    INTEGER :: num, npt

    DO i=1,geom%ntot
       circle%xbody(i) = geom%xbody(i)
       circle%ybody(i) = geom%ybody(i)
       circle%bodyno(i) = geom%bodyno(i)
    END DO

    npt = INT( pi*dia/ds )
    DO i=1,npt
       circle%bodyno(geom%ntot+i) = num
       circle%xbody(geom%ntot+i) = x + 0.5D0*dia*COS(2.*pi*REAL(i-1)/REAL(npt))
       circle%ybody(geom%ntot+i) = y + 0.5D0*dia*SIN(2.*pi*REAL(i-1)/REAL(npt))
    END DO
    circle%ntot = geom%ntot+npt

  END FUNCTION circle

!=============================================================================!
  FUNCTION move_body(itime) RESULT (motion)
    !-------------------------------------------------------------------------!
    ! This function returns the motion of the body for each Lagrangian point 
    ! that is not an actuator. One must also update the body location as well
    ! to be consistent.
    !-------------------------------------------------------------------------!
    REAL(KIND(0.D0)), DIMENSION(Nf) ::  motion
    REAL(KIND(0.D0)) :: angular = 0.5D0
    INTEGER :: itime, i
    motion = 0.d0

!     ! EXAMPLE 1: Impulsively start the cylinder moving to the left
!     DO i=1,(nb-n_act)
!       motion(f(i,1)) = - 1.D0 
!       motion(f(i,2)) =  0.D0 
!       xb(i) = xb(i) - 1.D0* (dt*REAL(itime) - dt*REAL(itime-1));
!       yb(i) = yb(i)
!     END DO
!
!    ! EXAMPLE 2a: Gradually rotate a cylinder from rest
!    DO i=1,(nb-n_act)
!       motion(f(i,1)) = - 1.0*exp( - (dt*REAL(itime)-4.d0*0.125)**2 / 2.d0/0.125**2 ) &
!                               *sin(2.d0*pi*REAL(i-1)/REAL(nb-n_act))
!       motion(f(i,2)) =   1.0*exp( - (dt*REAL(itime)-4.d0*0.125)**2 / 2.d0/0.125**2 ) &
!                               *cos(2.d0*pi*REAL(i-1)/REAL(nb-n_act)) 
!    END DO

!    ! EXAMPLE 2b: Impulsively rotate a cylinder
!    DO i=1,(nb-n_act)
!       motion(f(i,1)) = - 1.0*angular*yb(i)
!       motion(f(i,2)) =   1.0*angular*xb(i)
!       xb(i) = xb_og(i)*COS(angular*dt*itime) - yb_og(i)*SIN(angular*dt*itime)
!       yb(i) = xb_og(i)*SIN(angular*dt*itime) + yb_og(i)*COS(angular*dt*itime)
!    END DO

!    ! EXAMPLE 3: Sinusoidally oscillating cylinder in the y-direction
!    !
!    ! NOTE: Turn the freestream on! The cylinder will oscillate with
!    !
!    ! y(t)  = sin(t)
!    ! y'(t) = cos(t)
!    
!    DO i=1,(nb-n_act)
!      motion(f(i,1)) = 0.D0
!      motion(f(i,2)) = 0.5D0*COS(dt*DBLE(itime))
!      xb(i) = xb(i)
!      yb(i) = yb(i) + 0.5D0*(SIN(dt*DBLE(itime)) - SIN(dt*DBLE(itime-1)))
!    END DO

    ! Scaling required
    motion(f(1:(nb-n_act),1)) = motion(f(1:(nb-n_act),1)) * delta
    motion(f(1:(nb-n_act),2)) = motion(f(1:(nb-n_act),2)) * delta

  END FUNCTION move_body

!=============================================================================!  
  FUNCTION actuation(itime,old_motion) RESULT(motion)
    !-------------------------------------------------------------------------!
    ! I am not going to deal with actuation at this time!
    !-------------------------------------------------------------------------!
    REAL(KIND(0.D0)), DIMENSION(Nf) :: motion,old_motion
    INTEGER :: itime, npt, count, i_act, i_pt
    CHARACTER(3) :: act_num

    motion = old_motion
    count = 0

    DO i_act = 1,n_act

       WRITE(act_num,"(I3.3)") i_act
       OPEN(unit=210,file="input/slot."//act_num//".inp",form='formatted',status='old')
       READ(210,*) npt
       CLOSE(210)

       DO i_pt = 1,npt
          count = count + 1
          motion(f(nb-n_act_pt+count,1)) = u_jet(itime,i_act,1) * delta
          motion(f(nb-n_act_pt+count,2)) = u_jet(itime,i_act,2) * delta
       END DO

    END DO

  END FUNCTION actuation

!=============================================================================!
  FUNCTION u_jet( it, i_act, dir ) 
  !---------------------------------------------------------------------------!
  ! This function returnes the time-varying component of the jet velocity
  !---------------------------------------------------------------------------!
    IMPLICIT NONE

    INTEGER :: it, i_act, dir
    REAL(KIND(0.D0)) :: u_jet
    CHARACTER(3) :: act_num
    REAL(KIND(0.D0)) :: x_j,y_j,h_j,u_j,u_j_osc,theta_j,omg_j
    NAMELIST /read_act_parameters/ x_j,y_j,h_j,u_j,u_j_osc,theta_j,omg_j

    WRITE(act_num,"(I3.3)") i_act
    OPEN(unit=7,file="input/actuator."//act_num//".inp",form="formatted",status="old")
    READ(unit=7,nml=read_act_parameters)
    CLOSE(7)

    IF (dir==1) THEN
       u_jet = ( u_j + u_j_osc * SIN( omg_j*dt*DBLE(it) ) ) * COS( pi/180.d0*theta_j )
    ELSEIF (dir==2) THEN
       u_jet = ( u_j + u_j_osc * SIN( omg_j*dt*DBLE(it) ) ) * SIN( pi/180.d0*theta_j )
    ELSE
       STOP ' undefined direction in u_jet'
    END IF

  END FUNCTION u_jet

!=============================================================================!
  FUNCTION actuator_slot( x, y, hs, ds, angle, geom, num )
  !---------------------------------------------------------------------------!
  ! Actuation slot?!? What the hell is this?
  !
  ! This is probably some kind of model...
  !---------------------------------------------------------------------------!

    TYPE(bodylist) :: geom, actuator_slot
    REAL(KIND(0.D0)) :: x,y,seg,angle,hs,ds
    INTEGER :: num,i,npt
    LOGICAL :: readslot,readslot2
    CHARACTER(3) :: act_num

    ! Closest number of points to approximate the actuator slot
    ! with point forces
    npt = INT(hs/ds)

    ! INITIALIZE
    DO i=1,geom%ntot
       actuator_slot%xbody(i)  = geom%xbody(i)
       actuator_slot%ybody(i)  = geom%ybody(i)
       actuator_slot%bodyno(i) = geom%bodyno(i)
    END DO

    ! Define slot model location
    IF (npt==0) THEN
       ! For this case, the slot is approximated by a single point
       actuator_slot%xbody(geom%ntot+1)  = x
       actuator_slot%ybody(geom%ntot+1)  = y
       actuator_slot%bodyno(geom%ntot+1) = num
       actuator_slot%ntot                = geom%ntot+1
    ELSE
       ! For this case, we place finite number of point forces
       ! to represent the slot
       DO i = 1,npt+1
          seg = REAL(i-1)*hs/REAL(npt)
          actuator_slot%xbody(geom%ntot+i) = x - hs/2.0d0*COS(pi*(angle+90.)/180.) + seg*COS(pi*(angle+90.)/180.)
          actuator_slot%ybody(geom%ntot+i) = y - hs/2.0d0*SIN(pi*(angle+90.)/180.) + seg*SIN(pi*(angle+90.)/180.)
          actuator_slot%bodyno(geom%ntot+i)= num
       END DO
       actuator_slot%ntot = geom%ntot+npt+1
    END IF
    
    ! Writing slot setup information to file
    WRITE(act_num,"(I3.3)") -num

    IF (num==-1) THEN
       OPEN(unit=116,file="input/slot.dat",form="formatted",status="replace")
       OPEN(unit=118,file="input/slot."//act_num//".inp",form="formatted",status="replace")
    ELSE
       INQUIRE(file="input/slot.dat",exist=readslot)
       INQUIRE(file="input/slot."//act_num//".inp",exist=readslot2)
       IF (readslot) THEN 
          OPEN(unit=116,file="input/slot.dat",form="formatted",status="old",position="append")
       ELSE
          OPEN(unit=116,file="input/slot.dat",form="formatted",status="new")
       END IF
       IF (readslot2) THEN
          OPEN(unit=118,file="input/slot."//act_num//".inp",form="formatted",status="old",position="append")
       ELSE
          OPEN(unit=118,file="input/slot."//act_num//".inp",form="formatted",status="new")
       END IF
    END IF

    ! The following is written to input/slot.dat for users
    WRITE(116,*) " ----------- ACTUATOR No. ",-num
    WRITE(116,*) "   o Desired slot width               = ",hs
    IF (npt==0)  THEN
    WRITE(116,*) "   o Actual slot width in simulation  = ",ds
    ELSE
    WRITE(116,*) "   o Actual slot width in simulation  = ",hs
    END IF
    WRITE(116,*) "   o Number of pts used to model slot = ",(npt+1)
    CLOSE(116)
    
    ! The following is written for the code
    WRITE(118,*) (npt+1)
    CLOSE(118)

    ! Counting total number of points used for actuators
    n_act_pt = n_act_pt + (npt + 1)

  END FUNCTION actuator_slot
  
END MODULE grid

