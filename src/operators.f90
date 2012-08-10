MODULE operators

  !******************************************************************!
  !*                               IBFS 2.0                         *!
  !*                 c. Tim Colonius & Kunihiko (Sam) Taira         *!
  !*                 California Institute of Technology             *!
  !*                             March 1, 2007                      *!
  !******************************************************************!
  !*                                                                *!
  !*   New features                                                 *!
  !*     - completely revamped solution method                      *!
  !*     - discrete streamfunction/vorticity formulation            *!
  !*     - super grid solution of Lapalce eqn. for better           *!
  !*       far field boundary conditions                            *!
  !*     - Fast sine transform to invert Laplacians                 *!
  !*     - Cholesky factorization for solution of stationary        *!
  !*       bodies                                                   *!
  !*     - 2D for now                                               *!
  !*                                                                *!
  !******************************************************************!

  USE parameters
  USE grid
  IMPLICIT NONE
  
CONTAINS

!=============================================================================!
  SUBROUTINE preprocess
    !-------------------------------------------------------------------------!
    ! Performs the Cholesky factorization of EC(C^t C)^-1 C^t E^t. It is 
    ! performed once only and only for stationary bodies. 
    !-------------------------------------------------------------------------!
    USE variables

    REAL(KIND(0.D0)), DIMENSION(Nf) :: z   ! temporary force vector
    INTEGER :: i
    LOGICAL :: readchol

    ! HAS IT ALREADY BEEN DONE? IF NOT DO IT!
    INQUIRE(file='input/ibfs.chd',exist=readchol)
    IF (readchol) THEN
       CALL read_cholesky
    ELSE
       ! BUILD MATRIX ONE COLUMN AT A TIME
       WRITE(*,*) 'Precomputing body matrix... i ='
       WRITE(111,*) 'Precomputing body matrix...'         
       DO i=1,nf
          WRITE(*,"(A,I3)",ADVANCE="NO") '...',i
          z = 0.d0
          z(i) = 1.d0
          cholmat(1:Nf,i) = a_times( z )
       END DO
       
       WRITE(*,*) 'Performing Cholesky decomposition...'
       CALL choldc       
       CALL write_cholesky
    END IF
  END SUBROUTINE preprocess

!=============================================================================!
  SUBROUTINE advance(itime)
    !-------------------------------------------------------------------------!
    ! Advances simulation to the next time level.
    !-------------------------------------------------------------------------!
    USE myfft
    USE variables
    INTEGER :: itime

    ! INTERMEDIATE RESULTS
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: rhs, rhsbc
    REAL(KIND(0.D0)), DIMENSION(Nf) :: motion, rhsf
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1)) :: omega_bc
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1),mgridlev) :: lastbc
    INTEGER :: i,j,k

    itime = itime + 1

    !! STEP 0: PREPROCESSING AND SETUP ----------------------------------------
    
    ! WRITE DIAGNOSTICS TO SCREEN
    CALL CPU_TIME(cputime)
    WRITE(*,*) " "
    WRITE(*,*) "-------------------------------------------"
    WRITE(*,*) "* For itime =",itime
    WRITE(*,"(A,EN12.3,A)") " * cputime =",cputime," s"
    WRITE(111,*) "...Advancing to itime =",itime

    ! PREPROCESS IF STATIONARY BODY TO GET CHOLESKY DECOMPOSITION
    IF (itime==1.and.stationary) THEN
       CALL preprocess
    END IF

    !! STEP 1: SOLVE INTERMEDIATE CURL(MOMENTUM) EQ. ON EACH GRID LEVEL -------
    
    ! SAVE OMEGA BC's FOR EACH LEVEL BEFORE TELESCOPING IN
    DO k=1,mgridlev-1
       CALL get_bc(omega(:,:,k+1), lastbc(:,k), 0.25d0)
    END DO
    lastbc(:,mgridlev) = 0.d0 ! Largest grid level BC's set to zero
    
    ! LOOP OVER ALL GRID LEVELS GOING BACKWARD (EQ. 35)
    DO k = mgridlev,1,-1
       IF (k.lt.mgridlev) THEN
          CALL get_bc( omega(:,:,k+1), omega_bc, 0.25d0 )
       ELSE
          omega_bc = 0.d0 ! Largest grid level BC's set to zero
       END IF

       ! COMPUTE THE NONLINEAR TERM
       ! (This is not yet multiplied by the AB2 3/2 and 1/2 "stuff" yet)
       rhs = rot( nonlinear(omega(:,:,k), q(:,k), q0(:,k), lastbc(:,k) ))
    
       ! USE EXPLICIT EULER IF FIRST TIME STEP
       IF (itime==1) THEN
          rhs_old(:,:,k) = rhs(:,:)
       END IF
      
       ! GET BC'S ON OMEGA AND ASSIGN TO A SINGLE VECTOR, RHSBC
       ! (Note: This averages "lastbc" and "omega_bc")
       rhsbc = 0.d0
       CALL apply_bc( rhsbc, lastbc(:,k), vfac(k) )
       CALL apply_bc( rhsbc, omega_bc, vfac(k) )

       ! SOLVE FOR INTERMEDIATE OMEGA (using sine transforms)
       omega(:,:,k) = dst( lam1i(:,:,k) * &
                    ( dst( con1(k)*rhs(:,:)       + &
                           con2(k)*rhs_old(:,:,k) + &
                                   rhsbc(:,:)   ) + & 
                      lam1(:,:,k)*dst(omega(:,:,k)) ) )

      ! SAVE FOR NEXT TIME STEP
       rhs_old(:,:,k) = rhs(:,:)
    END DO

    ! We now have the intermediate vorticity evolution in all domains. 
    ! For the smallest domain we add the vorticity created by the body.

    !! STEP 2: SOLVE THE MODIFIED POISSON'S EQ. 36 ----------------------------
    
    ! BOUNDARY MOTION? UPDATE REG KERNEL IF NECESSARY
    IF (ANY(bdyno>0)) THEN
       motion = move_body(itime)
       CALL setup_reg
    ELSE
       motion = 0.d0
    END IF
    
    ! ACTUATION?
    IF (n_act>0) THEN
       motion = actuation(itime,motion)
    END IF
    
    ! ADD MOTION AND UNDERLYING POTENTIAL FLOW
    motion = motion - regT(q0(:,1))

    ! SOLVE MODIFIED POISSON'S EQ. 36
    
    ! Calculate the RHS of EQ. 36,
    ! q is re-used as temp variable here, the rhs (omega) is coarsified and
    ! Laplacian inverted we only need resulting velocity on grid 1 because
    ! q is not needed on other grid levels.
    CALL vort2flux( q(:,1), omega, 1 )
    rhsf = regT( q(:,1) ) - motion
    
    ! Solve the linear system 
    IF (stationary) THEN ! With Cholesky factorization (precomputed A)
       fb = cholsl( rhsf )
    ELSE ! Or by conjugate gradient method (uses a_times to get A(t))
       CALL cjgr( fb, rhsf ) 
    END IF
  
    !! STEP 3: PROJECTION STEP - COMPLETE OMEGA ON GRID 1 ---------------------
    omega(:,:,1) = omega(:,:,1) - ainv( rot( reg( fb ) ) )
    
    !! STEP 4: COARSIFY FINAL OMEGA AND GET NEW VELOCITIES ON ALL GRIDS -------
    CALL vort2flux( q, omega, mgridlev )

    ! SAVE A RESTART AND WRITE OTHER STUFF TO FILES
    IF ((MOD(itime,isave).eq.0).or.(itime==istop)) THEN
       CALL write_variables(itime)
    END IF

    CALL write_force( itime, fb/dt )
    CALL write_circ( itime )
    CALL write_slip( itime, q, motion )

  END SUBROUTINE advance

!=============================================================================!
  SUBROUTINE vort2flux( vel, omega, nlev )
  !---------------------------------------------------------------------------!
  ! Multiscale method to solve C^T C s = omega and return the velocity, C s.
  ! Results are returned in vel on each of the first nlev grids.
  !
  ! Warning: The vorticity field on all but the finest mesh is modified by the
  ! the routine in the following way: the value in the center of the domain is
  ! interpolated from the next finer mesh (the value near the edge is not 
  ! changed).
  !---------------------------------------------------------------------------!
    USE myfft
    INTEGER :: nlev
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n,mgridlev) :: omega
    REAL(KIND(0.D0)), DIMENSION(Nq,nlev) :: vel

    REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: s, vort ! streamfunction & vort
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1)) :: sbc  ! streamfunction bcs

    INTEGER :: k
    
    ! FIND COARSIFIED OMEGAS, EQ. 30
    DO k=2,mgridlev
       omega(:,:,k) = coarsify(omega(:,:,k), omega(:,:,k-1))
    END DO

    ! INVERT LAPLACIAN ON LARGEST GRID
    ! Zero bc's on s assumed for largest grid, EQ. 31
    sbc = 0.d0
    
    ! COMPUTE S, EQ. 32
    vort = omega(:,:,mgridlev)
    s = ctci( vort )
    
    ! COMPUTE VEL IF DESIRED BY TAKING CURL
    IF (nlev.ge.mgridlev) THEN
       vel(:,mgridlev) = curl(  s, sbc )
    END IF

    ! REPEAT EQ. 32 GOING FROM LARGEST GRID LEVEL IN
    DO k=mgridlev-1,1,-1
       CALL get_bc( s, sbc, 1.d0)       ! Get BC's from pervious s
       vort = omega(:,:,k)              ! Temporarily put in vort vector
       CALL apply_bc( vort, sbc, 1.d0)  ! Apply BC's
       s = ctci( vort )                 ! Compute new s
       IF (nlev.ge.k) THEN              ! Compute vel if desired
          vel(:,k) = curl( s, sbc )
       END IF
    END DO

  END SUBROUTINE vort2flux

!=============================================================================!
  FUNCTION curl( x, xbc )
    !-------------------------------------------------------------------------!
    ! Returns curl(x) given x and the values of the streamfunction on the 
    ! body.
    !-------------------------------------------------------------------------!
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: x
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1)) :: xbc
    REAL(KIND(0.D0)), DIMENSION(Nq) :: curl
    INTEGER :: i,j

    DO i=2,m
       DO j=2,n-1
          curl(u(i,j)) = x(i,j+1) - x(i,j)
       END DO
       j=1
       curl(u(i,j)) = x(i,j+1) - xbc( bottom + i )
       j=n
       curl(u(i,j)) = xbc( top + i ) - x(i,j)
    END DO
    i=1
    DO j=1,n
       curl(u(i,j)) = xbc( left + j + 1 ) -xbc( left + j )
    END DO
    i=m+1
    DO j=1,n
       curl(u(i,j)) = xbc( right + j + 1 ) -xbc( right + j )
    END DO

    DO j=2,n
       DO i=2,m-1
          curl(v(i,j)) = x(i,j) - x(i+1,j)
       END DO
       i=1
       curl(v(i,j)) = - x(i+1,j) + xbc( left + j )
       i=m
       curl(v(i,j)) = - xbc( right + j) +  x(i,j)
    END DO
    j=1
    DO i=1,m
       curl(v(i,j)) = xbc( bottom + i ) -xbc( bottom + i + 1 )
    END DO
    j=n+1
    DO i=1,m
       curl(v(i,j)) = xbc( top + i ) -xbc( top + i + 1 )
    END DO
  END  FUNCTION curl

!=============================================================================!
  FUNCTION rot( x )
    !-------------------------------------------------------------------------!
    ! Returns the transpose of the curl for x.
    !-------------------------------------------------------------------------!
    REAL(KIND(0.D0)), DIMENSION(Nq) :: x
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: rot
    INTEGER :: i,j

    DO i=2,m
       DO j=2,n
          rot(i,j) = x(v(i,j))-x(v(i-1,j)) - x(u(i,j)) + x(u(i,j-1))
       END DO
    END DO
  END FUNCTION rot

!=============================================================================!
  FUNCTION coarsify( crhs, rhs ) RESULT( arhs )
    !-------------------------------------------------------------------------!
    ! Given vorticity on a smaller, fine mesh, (rhs) this function 
    ! interpolates values to the center region of a larger, coarser 
    ! mesh (crhs). The values outside the center region are not
    ! modified. Result is placed in arhs.
    !-------------------------------------------------------------------------!
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: crhs, rhs
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: arhs
    INTEGER :: i,j,ii,jj
    
    arhs = crhs ! So you don't modify crhs inadvertently

    DO i=m/4+2,3*m/4
       ii = 2*i-m/2-1

      DO j=n/4+2,3*n/4 
        jj = 2*j-n/2-1

        arhs(i,j) = rhs(ii,jj) + &
                    0.50d0*( rhs(ii-1,jj  ) + rhs(ii+1,jj  )   + &
                             rhs(ii  ,jj-1) + rhs(ii  ,jj+1) ) + &
                    0.25d0*( rhs(ii-1,jj-1) + rhs(ii-1,jj+1)   + &
                             rhs(ii+1,jj-1) + rhs(ii+1,jj+1) )
      END DO
    END DO          
  END FUNCTION coarsify

!=============================================================================!
  SUBROUTINE get_bc( r, rbc, fac)
    !-------------------------------------------------------------------------!
    ! Given vorticity on a larger, coarser mesh, this subroutine
    ! interpolates its values to the edge of a smaller, finer mesh.
    !-------------------------------------------------------------------------!

    REAL(KIND(0.D0)), DIMENSION(:,:) :: r
    REAL(KIND(0.D0)), DIMENSION(:) :: rbc
    REAL(KIND(0.D0)) :: fac
    INTEGER :: i,j
    
    ! get interpolated boundary conditions on finer grid
    DO i=0,m,2
       rbc(bottom + i+1) = r(m/4+i/2,n/4)
       rbc(top + i+1) = r(m/4+i/2,3*n/4)
    END DO
    DO i=1,m-1,2
       rbc(bottom +i+1)  = 0.5d0*( r(m/4+(i+1)/2,n/4) + r(m/4-1+(i+1)/2,n/4) )
       rbc(top + i+1) = 0.5d0*( r(m/4+(i+1)/2,3*n/4) + r(m/4-1+(i+1)/2,3*n/4) )
    END DO

    DO j=0,n,2
       rbc(left + j+1) = r(m/4, n/4+j/2)
       rbc(right + j+1) = r(3*m/4, n/4+j/2)
    END DO
    DO j=1,n-1,2
       rbc(left + j+1) = 0.5d0*( r(m/4, n/4+(j+1)/2) + r(m/4, n/4-1+(j+1)/2) )
       rbc(right + j+1) = 0.5d0*( r(3*m/4, n/4+(j+1)/2) + r(3*m/4, n/4-1+(j+1)/2) )
    END DO

    rbc = rbc*fac
  END SUBROUTINE get_bc

!=============================================================================!
  SUBROUTINE apply_bc( r, rbc, fac)
    !-------------------------------------------------------------------------!
    ! Given vorticity at the edges of the domain, rbc (from larger, coarser
    ! mesh), this subroutine adds the values to correct the Laplacian of 
    ! vorticity on the smaller, finer domain, r.
    !-------------------------------------------------------------------------!
    REAL(KIND(0.D0)), DIMENSION(:,:) :: r
    REAL(KIND(0.D0)), DIMENSION(:) :: rbc
    REAL(KIND(0.D0)) :: fac
    INTEGER :: i,j
    
    ! Add bc's from coarser grid
    DO i=1,m-1
       r(i,1) = r(i,1) + fac* rbc( bottom + i+1 )
       r(i,n-1) = r(i,n-1) + fac* rbc( top + i+1 )
    END DO
    DO j=1,n-1
       r(1,j) = r(1,j) + fac* rbc( left + j+1 )
       r(m-1,j) = r(m-1,j) + fac* rbc( right + j+1 )
    END DO
  END SUBROUTINE apply_bc

!=============================================================================!
  SUBROUTINE cjgr( x, b ) 
    !-------------------------------------------------------------------------!
    ! Conjugate gradient solver for A x = b.
    ! Uses the routine a_times(x) to represent the "A matrix." x is the
    ! the returned solution. No preconditioner (yet).
    !-------------------------------------------------------------------------!
    IMPLICIT NONE

    REAL(KIND(0.d0)), DIMENSION(:), INTENT(IN)  :: b
    REAL(KIND(0.d0)), DIMENSION(:), INTENT(OUT) :: x
    REAL(KIND(0.d0)), DIMENSION(SIZE(b)) :: r,d,q,s
    REAL(KIND(0.d0)) :: delta_new,delta_old,eps, alpha, beta
    INTEGER :: iter

    iter = 0

    r = b - a_times(x) 

    d = r ! m_inverse_times(r)
    delta_new = DOT_PRODUCT(r,d) ! same as above for check
    eps = cgtol*cgtol

    DO WHILE ((iter.lt.cg_max_iter).and.(delta_new.gt.eps))
       iter=iter+1
       q = a_times( d )
       alpha = delta_new / DOT_PRODUCT( d, q )
       x = x + alpha * d

       IF (MOD(iter,50).eq.0) THEN
          r = b - a_times( x ) 
       ELSE
          r = r - alpha * q
       END IF

       s = r ! m_inverse_times(r)
       delta_old = delta_new
       delta_new = DOT_PRODUCT(r,s)

       beta = delta_new/delta_old
       d = s + beta * d
    END DO

    IF (iter.eq.cg_max_iter) THEN
       WRITE(111,*)  "WARNING!! CJGR used max iterations!"
       WRITE(*,*)    "=> itmax = ",iter,", res = ",delta_new
    END IF
    WRITE(111,*) "CJGR used ",iter," iterations." ! Just to check later
    WRITE(*,*) "* CJGR iterations = ",iter
  END SUBROUTINE cjgr

!=============================================================================!
  FUNCTION a_times(x) 
    !-------------------------------------------------------------------------!
    ! Gives the LHS of the equation that must be solved using either the 
    ! conjugate gradient method or Cholesky decomposition.
    !-------------------------------------------------------------------------!
    USE myfft
    REAL(KIND(0.D0)), DIMENSION(Nf) :: x, a_times
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n,mgridlev) :: vort
    REAL(KIND(0.D0)), DIMENSION(Nq) :: vel

    ! Zero all levels of temporary vorticity
    vort = 0.d0
    ! Regularize forces for level 1 vorticity
    vort(:,:,1) = ainv( rot( reg(x) ) )
    ! Coarsify vorticity and find curl(laplace^-1)
    CALL vort2flux( vel, vort, 1 )
    ! Regularize
    a_times = regT( vel )
  END FUNCTION a_times

!=============================================================================!
  FUNCTION nonlinear( omega, q, q0, bc ) RESULT(fq)
    !-------------------------------------------------------------------------!
    ! Result gives the nonlinear terms in rotational form, fq. Inputs are
    ! vorticity, omega, fluxes, q and q0, and boundary conditions, bc.
    !-------------------------------------------------------------------------!
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n) :: omega
    REAL(KIND(0.D0)), DIMENSION(Nq) :: q, q0, fq
    REAL(KIND(0.D0)), DIMENSION(2*(m+1)+2*(n+1)) :: bc
    REAL(KIND(0.D0)), DIMENSION(1:m+1,2:n) :: uavg
    REAL(KIND(0.D0)), DIMENSION(2:m,1:n+1) :: vavg
    INTEGER :: i,j

    DO j=2,n
       DO i=1,m+1
          uavg(i,j) = 0.5D0*( q(u(i,j))+q(u(i,j-1)) + q0(u(i,j))+q0(u(i,j-1)))
       END DO
    END DO
    DO j=1,n+1
       DO i=2,m
          vavg(i,j) = 0.5D0*( q(v(i,j))+q(v(i-1,j)) + q0(v(i,j))+q0(v(i-1,j)))
       END DO
    END DO

    DO i=2,m
       DO j=2,n-1
           fq(u(i,j)) = 0.5D0*( vavg(i,j+1)*omega(i,j+1) + vavg(i,j)*omega(i,j))
       END DO
       j = 1
           fq(u(i,j)) = 0.5D0*( vavg(i,j+1)*omega(i,j+1) + vavg(i,j)*bc(bottom+i))
       j = n
           fq(u(i,j)) = 0.5D0*( vavg(i,j+1)*bc(top+i) + vavg(i,j)*omega(i,j) )
    END DO
    ! NOTE...we don't need result for i=1 or 1=m+1 since it isn't needed by rot

    DO j=2,n
       DO i=2,m-1
           fq(v(i,j)) = -0.5D0*( uavg(i+1,j)*omega(i+1,j) + uavg(i,j)*omega(i,j))
       END DO
       i = 1
           fq(v(i,j)) = -0.5D0*( uavg(i+1,j)*omega(i+1,j) + uavg(i,j)*bc(left+j))
       i = m
           fq(v(i,j)) = -0.5D0*( uavg(i+1,j)*bc(right+j) + uavg(i,j)*omega(i,j) )
    END DO
    ! NOTE...we don't need result for j=1 or j=n+1 since it isn't needed by rot

  END FUNCTION nonlinear

!=============================================================================!
  FUNCTION reg( h0 ) RESULT( h )
    !-------------------------------------------------------------------------!
    ! Result gives the regularization of the body force, h, on the grid given 
    ! the body force, h0.
    !-------------------------------------------------------------------------!
    REAL(KIND(0.D0)), DIMENSION(Nf), INTENT(IN) :: h0
    REAL(KIND(0.D0)), DIMENSION(Nq)             :: h
    INTEGER :: i,k 

    h = 0.D0
    DO i=1,Nq
       DO k=1,Nsmear
          h(i) = h(i) + smear(i,k)*h0(ismear(i,k))
       END DO
    END DO
  END FUNCTION reg

!=============================================================================!
  FUNCTION regT( h ) RESULT( h0 ) 
    !-------------------------------------------------------------------------!
    ! Result interpolates h to a body point. This is the transpose of reg.
    !-------------------------------------------------------------------------!
    REAL(KIND(0.D0)), DIMENSION(Nq), INTENT(IN) :: h
    REAL(KIND(0.D0)), DIMENSION(Nf)             :: h0
    INTEGER :: i,k 

    h0 = 0.D0

    DO i=1,Nq
       DO k=1,Nsmear
          h0(ismear(i,k)) = h0(ismear(i,k)) + smear(i,k)*h(i)
       END DO
    END DO
  END FUNCTION regT

!=============================================================================!
  SUBROUTINE write_slip( it, x, xb )
    !-------------------------------------------------------------------------!
    ! This subroutine writes the slip velocity at time it, at the boundary,
    ! xb. This is to monitor the error, which is usually largest near the 
    ! body. It is written to output/slip.dat in ASCII format. The input x
    ! represents the fluxes.
    !-------------------------------------------------------------------------!
    USE parameters
    INTEGER          :: it
    REAL(KIND(0.D0)), DIMENSION(Nq) :: x
    REAL(KIND(0.D0)), DIMENSION(Nf) :: xb
    REAL(KIND(0.D0)) :: slip
    LOGICAL :: readslip

    ! Computing slip 
    slip = MAXVAL( ABS (regT(x) - xb) ) 

    INQUIRE(file='output/slip.dat',exist=readslip)

    IF (it==1) THEN
       OPEN(unit=106,file="output/slip.dat",form="formatted",status="replace")
    ELSE
       IF (readslip) THEN 
          OPEN(unit=106,file="output/slip.dat",form="formatted",status="old",position="append")
       ELSE
          OPEN(unit=106,file="output/slip.dat",form="formatted",status="new")
       END IF
    END IF
    WRITE(106,*) it, slip
    CLOSE(106)
    WRITE(*,*) "SLIP =", slip
  END SUBROUTINE write_slip

!=============================================================================!
  SUBROUTINE write_force( it, xg )
    !-------------------------------------------------------------------------!
    ! This subroutine writes forces at time it to a file in ASCII format.
    !-------------------------------------------------------------------------!
    USE parameters
    INTEGER          :: it, i, j, k

    REAL(KIND(0.D0)) :: forcex(minval(bdyno):maxval(bdyno)) 
    REAL(KIND(0.D0)) :: forcey(minval(bdyno):maxval(bdyno)) 
    REAL(KIND(0.D0)), DIMENSION(Nf) :: xg
    REAL(KIND(0.D0)), DIMENSION(Nf) :: gtemp ! temporary storage 
    REAL(KIND(0.D0)), DIMENSION(Nq) :: qtemp
    LOGICAL :: readforce

    forcex = 0.D0; forcey = 0.D0

    ! Loop over the number of bodies
    DO k = minval(bdyno),maxval(bdyno)

       gtemp = 0.D0
       DO i = 1,nb
          IF (bdyno(i)==k) THEN
             gtemp(i)    = xg(i)    ! x-directional force
             gtemp(i+nb) = xg(i+nb) ! y-directional force
          END IF
       END DO
       qtemp = reg( gtemp )
    
       ! Assuming that the immersed object is correctly placed 
       ! within the uniform mesh region.

       ! The forces are integrated with a simple 2-D rectangular
       ! quadrature rule
       ! ...x-force
       DO i = 1,m+1
          DO j = 1,n
             forcex(k) = forcex(k) + qtemp(u(i,j))
          END DO
       END DO
       ! ...y-force
       DO j = 1,n+1
          DO i = 1,m
             forcey(k) = forcey(k) + qtemp(v(i,j))
          END DO
       END DO

       ! The force is normalized by (NEEDS TO BE FIXED)
       !
       !   f --> f / (1/2 * rho * U^2 * L)
       !
       ! where we take rho = U = L = 1.
       forcex(k) = forcex(k) * 2.d0 * delta
       forcey(k) = forcey(k) * 2.d0 * delta
   
    END DO
    
    ! Force over the bodies are saved with the following format
    !
    ! (time index) (fx:body1) (fx:body2) ... (fy:body1) (fy:body2) ...
    !
    ! If a simulation is restarted, force calculations are appended to 
    ! the pre-existing file.

    INQUIRE(file='output/force.dat',exist=readforce)
    IF (it==1) THEN
       OPEN(unit=100,file="output/force.dat",form="formatted",status="replace")
    ELSE
       IF (readforce) THEN
          OPEN(unit=100,file="output/force.dat",form="formatted",status="old",position="append")
       ELSE
          OPEN(unit=100,file="output/force.dat",form="formatted",status="new")
       END IF
    END IF
    WRITE(100,*) it,(forcex(k), k=(minval(bdyno)),(maxval(bdyno))),(forcey(k), k=(minval(bdyno)),(maxval(bdyno)))
    CLOSE(100)
    WRITE(111,*) "Computed x- & y-forces at itime=",it
  END SUBROUTINE write_force

!=============================================================================!
  SUBROUTINE write_circ( it )
    !-------------------------------------------------------------------------!
    ! This subroutine writes the circulation of the bodies to a file in 
    ! ASCII format.
    !-------------------------------------------------------------------------!
    USE parameters
    USE variables
    INTEGER          :: it, i, j, k

    REAL(KIND(0.D0)) :: circ(1:mgridlev) 

    LOGICAL :: readcirc

    circ = 0.D0

    ! Loop over the number of bodies
    DO k = 1,mgridlev
       DO i = 2,m
          DO j = 2,n
             circ(k) = circ(k) + omega(i,j,k)
          END DO
       END DO
    END DO

    INQUIRE(file='output/circ.dat',exist=readcirc)
    IF (it==1) THEN
       OPEN(unit=100,file="output/circ.dat",form="formatted",status="replace")
    ELSE
       IF (readcirc) THEN
          OPEN(unit=100,file="output/circ.dat",form="formatted",status="old",position="append")
       ELSE
          OPEN(unit=100,file="output/circ.dat",form="formatted",status="new")
       END IF
    END IF
    WRITE(100,*) it,(circ(k), k=1,mgridlev) 
    CLOSE(100)

  END SUBROUTINE write_circ

!=============================================================================!
  SUBROUTINE write_tecplot(it, igrid)
    !-------------------------------------------------------------------------!
    ! Writes the variables of the igrid level of MGRID to a Tecplot readable 
    ! ASCII file. This is considered a post processing routine.
    !-------------------------------------------------------------------------!
    USE myfft
    USE variables

    CHARACTER(7) :: charit
    CHARACTER(1) :: chargrid
    INTEGER :: it,igrid
    REAL(KIND(0.D0)), DIMENSION(2:m,2:n,5) :: fout   ! output at vertices
    REAL(KIND(0.D0)), DIMENSION(Nq) :: etf
    INTEGER :: i,j,k
    REAL(KIND(0.D0)) :: fac, del, offx,offy

    fac = REAL( 2**(igrid-1) )
    offx = (fac-1.d0)/2.d0 * len
    offy = (fac-1.d0)/2.d0 * delta*REAL(n)
    del = fac*delta

    WRITE(charit,"(I7.7)") it
    WRITE(chargrid,"(I1.1)") igrid

    etf = reg(fb)

    DO i=2,m
       DO j=2,n
          fout(i,j,5) = omega(i,j,igrid)/del/del
          fout(i,j,1) = 0.5d0*(q(u(i,j),igrid)+q(u(i,j-1),igrid)+ &
                               q0(u(i,j),igrid)+q0(u(i,j-1),igrid))/del
          fout(i,j,2) = 0.5d0*(q(v(i,j),igrid)+q(v(i-1,j),igrid)+ &
                               q0(v(i,j),igrid)+q0(v(i-1,j),igrid))/del
          fout(i,j,3) = delta*(etf(u(i,j))+etf(u(i,j-1)))/dt
          fout(i,j,4) = delta*(etf(v(i,j))+etf(v(i-1,j)))/dt
       END DO
    END DO

    WRITE(charit,"(I7.7)") it
    OPEN(unit=100,file="post/ibfs_g"//chargrid//"t"//charit//"_vel.dat",form="formatted",status="unknown")

    WRITE(100,*) 'TITLE     = "IBFS V1.0 AT ITIME=',it,'"'
    WRITE(100,*) 'VARIABLES = "x"'
    WRITE(100,*) '"y"'
    WRITE(100,*) '"u"'
    WRITE(100,*) '"v"'
    WRITE(100,*) '"etfx"'
    WRITE(100,*) '"etfy"'
    WRITE(100,*) '"vort"'
    WRITE(100,*) 'ZONE T="Rectangular zone"'
    WRITE(100,*) 'I=',m-1,', J=',n-1,', K=1, ZONETYPE=Ordered'
    WRITE(100,*) 'DATAPACKING=POINT'
    WRITE(100,*) 'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)'
    DO j=2,n
       DO i=2,m
          write(100,"(8E14.5)") fac*delta*REAL(i-1) - offx - offsetx , &
                                fac*delta*REAL(j-1) - offy - offsety, &
               ( fout(i,j,k), k=1,5)
       END DO
    END DO
    CLOSE(100)
  END SUBROUTINE write_tecplot

!=============================================================================!
  SUBROUTINE choldc
    !-------------------------------------------------------------------------!
    ! Performs the Cholesky factorization for the body to quickly find the
    ! unknown boundary forces.
    !-------------------------------------------------------------------------!
    USE variables
    REAL(KIND(0.D0)) :: sum

    INTEGER :: i,j,k

    DO i=1,Nf
      WRITE(*,"(A,I4)",ADVANCE="NO") '....',i
       DO j=i,Nf
          sum=cholmat(i,j)
          DO k=i-1,1,-1
             sum=sum-cholmat(i,k)*cholmat(j,k)
          END DO
          IF(i.EQ.j)THEN
            IF(sum.LE.0.) STOP 'choldc failed'
            cholvec(i)=SQRT(sum)
          ELSE
            cholmat(j,i)=sum/cholvec(i)
          ENDIF
       END DO
    END DO
  END SUBROUTINE choldc

!=============================================================================!
  FUNCTION cholsl(b) RESULT(x)
    !-------------------------------------------------------------------------!
    ! Solves A x = b given its Cholesky decomposition.
    !-------------------------------------------------------------------------!
    USE variables
    REAL(KIND(0.D0)), DIMENSION(:) :: b
    REAL(KIND(0.D0)), DIMENSION(SIZE(b)) :: x

    INTEGER :: i,k
    REAL(KIND(0.D0)) ::  sum

    DO i=1,Nf
       sum=b(i)
       DO  k=i-1,1,-1
          sum=sum-cholmat(i,k)*x(k)
       END DO
       x(i)=sum/cholvec(i)
    END DO
    DO i=Nf,1,-1
       sum=x(i)
       DO k=i+1,Nf
          sum=sum-cholmat(k,i)*x(k)
       END DO
       x(i)=sum/cholvec(i)
    END DO
  END FUNCTION cholsl
!=============================================================================!
  
END MODULE operators

