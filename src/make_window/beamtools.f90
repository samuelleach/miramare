module beamtools
  
  use healpix_types, ONLY: SP,DP,I4B,SPC
  use pix_tools
  use mpi_stop
  use Amlutils

contains

  subroutine ReadBeam(beam, filename, lmax_beam)
    real(dp), pointer :: beam(:)
    character(LEN=*), intent(in) :: filename
    integer, intent(out) :: lmax_beam
    integer l,lmax
    real(dp) inB
    character(LEN=256) :: inline
    real(dp), pointer :: temp_beam(:)

    lmax = 2500
    allocate(temp_beam(0:lmax))
    
    call OpenTxtFile(filename,1)
    do
       read (1,'(a)', err = 500, end =500) inline
       if (inline(1:1) /= '#' .and. inline/='') then
          read(inline,*) l, inB
          if (l<=lmax) temp_beam(l) = inB
       end if
    end do
    
500 close(1)
    lmax_beam = l
    allocate(beam(0:lmax_beam))
    do l=0,lmax_beam
       beam(l)=temp_beam(l)
    end do
    
    deallocate(temp_beam)
    
  end subroutine ReadBeam

  subroutine WriteBeam(beam, filename, lmax)
    real(dp), pointer            :: beam(:)
    character(LEN=*), intent(in) :: filename
    integer, intent(in)          :: lmax
    integer l,aunit

    aunit = 100

    open(unit=aunit,file=filename,form='formatted',&
         status='replace',err=500)
    
    do l = 0, lmax
       write(aunit,*) l,beam(l)
    end do

500 close(aunit)  
    
  end subroutine WriteBeam

  Subroutine variancebeam(signal_beam, lmax_beam,nside, variance_beam,lmax_variancebeam)
    integer(kind=I4B), intent(in)  :: lmax_beam
    real(kind=DP),     intent(in)  :: signal_beam(0:lmax_beam)
    integer(kind=I4B), intent(in)  :: nside
    integer(kind=I4B), intent(in)  :: lmax_variancebeam
    real(kind=DP),     intent(out) :: variance_beam(0:lmax_variancebeam)
    
    integer lhalf,nint,ll,ii
    real(kind=DP), allocatable :: rad(:),x(:),sinrad(:)    
    real(kind=DP), allocatable :: llarray(:),signal_beam2(:)
    real(kind=DP), allocatable :: conva(:)
    real(kind=DP), allocatable :: mult(:)
    real(kind=DP), allocatable :: temp(:)
    real(kind=DP), allocatable :: lgndr(:,:)    

    real(kind=DP), parameter :: threshold =4d-4

    !----------------------------------------------------------------------
    ! Make Lagrange-polynomial transform of amp to get beam window function
    !----------------------------------------------------------------------

    !----------------------------------------------------------------
    ! Choose scale to match size of beam. First find half-power point
    !----------------------------------------------------------------
    do ll = 0, lmax_beam
       if(signal_beam(ll) .gt. 0.5) lhalf = ll
    end do

    allocate(llarray(0:lmax_variancebeam))
    allocate(signal_beam2(0:lmax_variancebeam))
    do ll = 0, lmax_variancebeam
       llarray(ll) = float(ll)
       if (ll .le. lmax_beam) then
          signal_beam2(ll) = signal_beam(ll)
       else
          signal_beam2(ll) = 0.
       endif
    end do

    !-----------------------------------------------------------------------
    ! Calculate beam out to about 20 * half-power radius, roughly (note that
    ! this gives 50 points out to the half-power point).
    !-----------------------------------------------------------------------
    nint = 1000
    allocate(rad(0:nint),sinrad(0:nint),x(0:nint))
    do ii = 0, nint
       rad(ii)    = dfloat(ii)*10.d0*pi/dfloat(lhalf*nint)
       x(ii)      = cos(rad(ii))
       sinrad(ii) = sin(rad(ii))
    end do
    
    allocate(lgndr(0:nint,0:lmax_variancebeam))
    do ii = 0, nint
       call lgnd(lmax_variancebeam,x(ii),lgndr(ii,:))
    end do

    !-------------------------------------------
    ! Generate radial profile of convolving beam:
    !-------------------------------------------
    allocate(conva(0:nint))
    do ii = 0, nint
       conva(ii) = sum((llarray(:)+0.5d0)*signal_beam2(:)*lgndr(ii,:))
    end do
    do ii = 0, nint
       conva(ii) = conva(ii) / (2.0*pi)
    end do
    PRINT * , 'Peak of convolving beam is ', conva(0), MAXVAL(conva)

    !------------------------------------------------------------
    ! Define variance beam amplitude array of same size as convbl

    ! Square convolving beam and convert back to window function
    !------------------------------------------------------------
    allocate(mult(0:nint),temp(0:nint))

    do ii=0,nint
       mult(ii) = sinrad(ii)*conva(ii)**2
    end do

    do ll = 0, lmax_variancebeam
       temp(:) = mult(:)*lgndr(:,ll)
       call avint(temp,rad,nint+1,rad(0),rad(nint),variance_beam(ll))
    end do
    do ll = 0, lmax_variancebeam
       variance_beam(ll) = variance_beam(ll)*2.*pi
    end do
    
    !----------------------
    ! Multiply by Omega_pix
    !----------------------
    do ll = 0, lmax_variancebeam
       variance_beam(ll) = variance_beam(ll)*4.*pi/nside2npix(nside)
    end do

    PRINT *, ' CVBL[0] =', variance_beam(0)

  end subroutine variancebeam


  subroutine CleanBeam(beam,lmax_beam, threshold, lthreshold)
    integer lmax_beam
    real(dp), intent(inout) :: beam(0:lmax_beam)
    integer, intent(out)    :: lthreshold
    real(dp) threshold
    logical past_threshold


    integer ll

    past_threshold = .false.

    !--------------------------
    ! Clean up beam at high ell
    !--------------------------
    do ll=1,lmax_beam
       if ( beam(ll) .lt. beam(0)*threshold) then
          beam(ll)       = 0.0
          past_threshold = .true.
       else
          if (.not. past_threshold) then
             lthreshold  = ll
          end if
       end if

       if ( beam(ll) .lt. 0.1*beam(0)) then
          if (beam(ll) .gt. beam(ll-1) ) then
             beam(ll) = beam(0)*threshold
          endif
       endif
    end do

  end subroutine CleanBeam

    

!!!!!!!!!!!!!!!!!!!!!!!!!!!   Program 5.4   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
! Please Note:                                                          !
!                                                                       !
! (1) This computer program is written by Tao Pang in conjunction with  !
!     his book, "An Introduction to Computational Physics," published   !
!     by Cambridge University Press in 1997.                            !
!                                                                       !
! (2) No warranties, express or implied, are made for this program.     !
!                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE LGND (LMAX,X,P)
  !
  ! Subroutine to generate Legendre polynomials P_L(X)
  ! for L = 0,1,...,LMAX with given X.
  !
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: LMAX
  INTEGER :: L
  real(kind=dp), INTENT (IN) :: X
  REAL(kind=dp), INTENT (OUT), DIMENSION (0:LMAX) :: P
  real(kind=dp) :: ll
  !
  P(0) = 1.0d0
  P(1) = X
  DO L = 1, LMAX-1
     ll=dfloat(l)
    P(L+1) = ((2.0d0*ll+1d0)*X*P(L)-ll*P(L-1))/(ll+1d0)
  END DO
END SUBROUTINE LGND


SUBROUTINE lpoly (n,x,pl)
!
! FUNCTIONAL DESCRIPTION:
!
!srk this subroutine calculates the ordinary Legendre Polynomials of
!srk order 0 to n-1 of argument  x  and stores them in the vector
!srk pl.  They are calculated by recursion relation from the first two
!srk polynomials.
!
!srk written by A. J. Sierk   LANL  T-9  February,1984
!
!srk NOTE:  pl and x must be REAL (KIND=8) :: on 32-bit computers!
!
! DUMMY ARGUMENTS:
!
!    input: x (REAL (KIND=4) ::) - argument of Legendre Polynomial
!           n (INTEGER (kind=4) ::) - n-1 is highest order of polynomial
!    output: pl(1:20) (REAL (KIND=4) ::) - value of Legendre Polynomials
!

INTEGER, INTENT(IN) :: n
REAL (KIND=dp), INTENT(IN)  :: x
REAL (KIND=dp), DIMENSION(0:n), INTENT(OUT) :: pl

INTEGER :: i

! pl(1) = 1.0
! pl(2) = x
! DO i=3,n
!   pl(i) = (REAL(2*i-3,KIND=r8)*x*pl(i-1) &
!   - REAL(i-2,KIND=r8)*pl(i-2))/REAL(i-1,KIND=r8)
! END DO
pl(0) = 1.0
pl(1) = x
DO i=2,n
  pl(i) = (REAL(2*i-3,KIND=dp)*x*pl(i-1) &
  - REAL(i-2,KIND=dp)*pl(i-2))/REAL(i-1,KIND=dp)
END DO

RETURN
END SUBROUTINE lpoly


subroutine avint ( ftab, xtab, ntab, a, b, result )
  !
  !***********************************************************************
  !! AVINT estimates the integral of unevenly spaced data.
  !
  !
  !  Discussion:
  !
  !    The method uses overlapping parabolas and smoothing.
  !
  !  Reference:
  !
  !    Philip Davis and Philip Rabinowitz,
  !    Methods of Numerical Integration,
  !    Blaisdell Publishing, 1967.
  !
  !    P E Hennion, &
  !    Algorithm 77, &
  !    Interpolation, Differentiation and Integration,
  !    Communications of the Association for Computing Machinery,
  !    Volume 5, page 96, 1962.
  !
  !  Modified:
  !
  !    30 October 2000
  !
  !  Parameters:
  !
  !    Input, real FTAB(NTAB), the function values, &
  !    FTAB(I) = F(XTAB(I)).
  !
  !    Input, real XTAB(NTAB), the abscissas at which the
  !    function values are given.  The XTAB's must be distinct
  !    and in ascending order.
  !
  !    Input, integer NTAB, the number of entries in FTAB and
  !    XTAB.  NTAB must be at least 3.
  !
  !    Input, real A, the lower limit of integration.  A should
  !    be, but need not be, near one endpoint of the interval
  !    (X(1), X(NTAB)).
  !
  !    Input, real B, the upper limit of integration.  B should
  !    be, but need not be, near one endpoint of the interval
  !    (X(1), X(NTAB)).
  !
  !    Output, real RESULT, the approximate value of the integral.
  !
  implicit none
  !
  integer ntab
  !
  real(kind=dp) a
  real atemp
  real(kind=dp) b
  real btemp
  real ca
  real cb
  real cc
  real ctemp
  real(kind=dp) ftab(ntab)
  integer i
  integer ihi
  integer ilo
  integer ind
  real(kind=dp) result
  real sum1
  real syl
  real term1
  real term2
  real term3
  real x1
  real x2
  real x3
  real(kind=dp) xtab(ntab)
  !
  if ( ntab < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AVINT - Fatal error!'
    write ( *, '(a,i6)' ) '  NTAB is less than 3.  NTAB = ', ntab
    stop
  end if

  do i = 2, ntab

    if ( xtab(i) <= xtab(i-1) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'AVINT - Fatal error!'
      write ( *, '(a)' ) '  XTAB(I) is not greater than XTAB(I-1).'
      write ( *, '(a,i6)' ) '  Here, I = ', I
      write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ', xtab(i-1)
      write ( *, '(a,g14.6)' ) '  XTAB(I) =   ', xtab(i)
      stop
    end if

  end do

  result = 0.0E+00

  if ( a == b ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AVINT - Warning!'
    write ( *, '(a)' ) '  A = B, integral=0.'
    return
  end if
  !
  !  If A > B, temporarily switch A and B, and store sign.
  !
  if ( a > b ) then
    syl = b
    b = a
    a = syl
    ind = -1
  else
    syl = a
    ind = 1
  end if
  !
  !  Bracket A and B between XTAB(ILO) and XTAB(IHI).
  !
  ilo = 1
  ihi = ntab

  do i = 1, ntab
    if ( xtab(i) >= a ) then
      exit
    end if
    ilo = ilo + 1
  end do

  ilo = max ( 2, ilo )
  ilo = min ( ilo, ntab-1 )

  do i = 1, ntab
    if ( b >= xtab(i) ) then
      exit
    end if
    ihi = ihi - 1
  end do

  ihi = min ( ihi, ntab-1 )
  ihi = max ( ilo, ihi-1 )
  !
  !  Carry out approximate integration from XTAB(ILO) to XTAB(IHI).
  !
  sum1 = 0.0E+00

  do i = ilo, ihi

    x1 = xtab(i-1)
    x2 = xtab(i)
    x3 = xtab(i+1)

    term1 = ftab(i-1) / ((x1-x2)*(x1-x3))
    term2 = ftab(i) / ((x2-x1)*(x2-x3))
    term3 = ftab(i+1) / ((x3-x1)*(x3-x2))

    atemp = term1 + term2 + term3
    btemp = -(x2+x3)*term1-(x1+x3)*term2-(x1+x2)*term3
    ctemp = x2*x3*term1+x1*x3*term2+x1*x2*term3

    if ( i <= ilo ) then
      ca = atemp
      cb = btemp
      cc = ctemp
    else
      ca = 0.5E+00 * ( atemp + ca )
      cb = 0.5E+00 * ( btemp + cb )
      cc = 0.5E+00 * ( ctemp + cc )
    end if

    sum1 = sum1 &
&     + ca * ( x2**3 - syl**3 ) / 3.0E+00 &
&     + cb * 0.5E+00 * ( x2**2 - syl**2 ) &
&     + cc * ( x2 - syl )

    ca = atemp
    cb = btemp
    cc = ctemp

    syl = x2

  end do

  result = sum1 &
&   + ca * ( b**3 - syl**3 ) / 3.0E+00 &
&   + cb * 0.5E+00 * ( b**2 - syl**2 ) &
&   + cc * ( b - syl )
  !
  !  Restore original values of A and B, reverse sign of integral
  !  because of earlier switch.
  !
  if ( ind /= 1 ) then
    ind = 1
    syl = b
    b = a
    a = syl
    result = -result
  end if

  return
end subroutine avint



end module beamtools
