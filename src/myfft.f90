MODULE myfft
  !---------------------------------------------------------------------------!
  ! A simple module for performing FFTs and IFFTs. This also includes various
  ! variables related to the simulation using MGRID. This includes various
  ! constant factors and normalizations, eigenvalues, and inverses of 
  ! eigenvalue matrices.
  !---------------------------------------------------------------------------!
  IMPLICIT NONE

  INCLUDE 'fftw3.f' ! INCLUDE the FFTW C library  
  
  REAL(KIND(0.D0)), DIMENSION(:),     ALLOCATABLE :: vfac, con1, con2
  REAL(KIND(0.D0)), DIMENSION(:,:),   ALLOCATABLE :: in, lam, laminv
  REAL(KIND(0.D0)), DIMENSION(:,:,:), ALLOCATABLE :: lam1, lam1i
  INTEGER*8 :: forward
  REAL(KIND(0.D0)) :: normalize
  INTEGER*8 :: mm,nn

CONTAINS
!=============================================================================!
  SUBROUTINE setup_fft
    !-------------------------------------------------------------------------!
    ! Sets up the FFT for use with the current code. Also gets some useful
    ! factors necessary for the simulation in addition to eigenvalues to 
    ! populate the \Lambda matrix (and its mutations).
    !-------------------------------------------------------------------------!
    USE grid
    USE parameters
    INTEGER :: i,j,k
    REAL(KIND(0.D0)) :: del2, del22

    mm = m
    nn = n

    ALLOCATE( con1(mgridlev), con2(mgridlev) )
    ALLOCATE( vfac(mgridlev) )
    ALLOCATE( in(mm-1,nn-1), laminv(mm-1,nn-1), lam(mm-1,nn-1) )
    ALLOCATE( lam1(mm-1,nn-1,mgridlev), lam1i(mm-1,nn-1,mgridlev) )
    
    ! CREATE A FFTW PLAN
    CALL dfftw_plan_r2r_2d(forward, mm-1,nn-1,in,in, FFTW_RODFT00, FFTW_RODFT00,  FFTW_FORWARD, FFTW_EXHAUSTIVE)
    normalize = 4.d0*REAL(mm*nn) ! Normalization factor for FFT pair

    ! SOME USEFUL FACTORS
    del2 = delta*delta ! Grid spacing squared on the smallest grid level
    DO k=1,mgridlev
      del22 = del2*4.d0**(k-1) ! Grid spacing on the kth grid level
      vfac(k) =  0.5d0*dt/Re/del22/normalize ! Factor used to get vorticity BC's from coarser mesh
      con1(k) =  1.5d0*dt/   del22/normalize ! Factor for AB2 on convective terms
      con2(k) = -0.5d0*dt/   del22/normalize ! Factor for AB2 on convective terms
    END DO

    ! GET EIGENVALUES FOR C^T C (AND I + C^T C)
    DO j=1,nn-1
      DO i=1,mm-1
        ! Eigenvalues of C^T C
        lam(i,j) = -2.d0*( COS( pi*REAL(i)/REAL(mm) ) + &
                             COS( pi*REAL(j)/REAL(nn) ) - 2.d0 )
        laminv(i,j) = 1.d0/lam(i,j)/normalize
          
        DO k=1,mgridlev
          del22 = del2* 4.d0**(k-1)
          ! For (I - 0.5 b dt C^T C)
          lam1(i,j,k) =       (1.d0 - 0.5d0*dt*lam(i,j)/del22/Re)/normalize
          ! For (I + 0.5 b dt C^T C)^(-1)
          lam1i(i,j,k) = 1.d0/(1.d0 + 0.5d0*dt*lam(i,j)/del22/Re)
        END DO
      END DO
    END DO
  END SUBROUTINE setup_fft

!=============================================================================!
  SUBROUTINE destroy_fft
    !-------------------------------------------------------------------------!
    ! Deallocates the variables involved in the FFT. Also destroys the plans.
    !-------------------------------------------------------------------------!    
    CALL dfftw_destroy_plan(forward)
    DEALLOCATE(in, laminv, lam, vfac, con1, con2)
    DEALLOCATE(lam1, lam1i)
  END SUBROUTINE destroy_fft

!=============================================================================!
  FUNCTION dst( psi )
    !-------------------------------------------------------------------------!
    ! Performs the discrete sine transform on psi. CAREFUL! Two calls of dst 
    ! need to be divided by "normalize" to return the original vector. 
    !-------------------------------------------------------------------------!
    REAL(KIND(0.D0)), DIMENSION(:,:) :: psi
    REAL(KIND(0.D0)), DIMENSION(2:mm,2:nn) :: dst

    in =  psi
    CALL dfftw_execute(forward)
    dst = in
  END FUNCTION dst

!=============================================================================!
  FUNCTION ctci( omega ) 
    !-------------------------------------------------------------------------!
    ! This represents (C^T C)^{-1}. See FUNCTION dst( psi ) for more info.
    !-------------------------------------------------------------------------!    
    REAL(KIND(0.D0)), DIMENSION(:,:) :: omega
    REAL(KIND(0.D0)), DIMENSION(2:mm,2:nn) :: ctci
        
    in =  omega
    CALL dfftw_execute(forward)
    in = laminv * in
    CALL dfftw_execute(forward)
    ctci =  in
  END  FUNCTION ctci

!=============================================================================!
  FUNCTION ainv( omega ) 
    !-------------------------------------------------------------------------!
    ! This represents (I + 0.5 b dt C^T C)^{-1}.
    !-------------------------------------------------------------------------!
    REAL(KIND(0.D0)), DIMENSION(:,:) :: omega
    REAL(KIND(0.D0)), DIMENSION(2:mm,2:nn) :: ainv
        
    in =  omega
    CALL dfftw_execute(forward)
    in = lam1i(:,:,1) * in / normalize
    CALL dfftw_execute(forward)
    ainv =  in
  END  FUNCTION ainv

END MODULE myfft
