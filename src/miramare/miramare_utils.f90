module miramare_utils

  use healpix_types
  
  implicit none
  

contains
  
function r8vec_ascends_strictly( n, x )

!*****************************************************************************80
!
!! R8VEC_ASCENDS_STRICTLY determines if an R8VEC is strictly ascending.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!    Notice the effect of entry number 6 in the following results:
!
!      X = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.4, 9.8 )
!      Y = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.5, 9.8 )
!      Z = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.6, 9.8 )
!
!      R8VEC_ASCENDS_STRICTLY ( X ) = FALSE
!      R8VEC_ASCENDS_STRICTLY ( Y ) = FALSE
!      R8VEC_ASCENDS_STRICTLY ( Z ) = TRUE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, real ( kind = 8 ) X(N), the array to be examined.
!
!    Output, logical R8VEC_ASCENDS_STRICTLY, is TRUE if the
!    entries of X strictly ascend.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical r8vec_ascends_strictly
  real    ( kind = 8 ) x(n)

  do i = 1, n - 1
    if ( x(i+1) <= x(i) ) then
      r8vec_ascends_strictly = .false.
      return
    end if
  end do

  r8vec_ascends_strictly = .true.

  return
end function r8vec_ascends_strictly


subroutine interp_linear ( dim_num, data_num, t_data, p_data, interp_num, &
  t_interp, p_interp )

!*****************************************************************************80
!
!! INTERP_LINEAR applies piecewise linear interpolation to data.
!
!  Discussion:
!
!    From a space of DIM_NUM dimensions, we are given a sequence of
!    DATA_NUM points, which are presumed to be successive samples
!    from a curve of points P.
!
!    We are also given a parameterization of this data, that is,
!    an associated sequence of DATA_NUM values of a variable T.
!    The values of T are assumed to be strictly increasing.
!
!    Thus, we have a sequence of values P(T), where T is a scalar,
!    and each value of P is of dimension DIM_NUM.
!
!    We are then given INTERP_NUM values of T, for which values P
!    are to be produced, by linear interpolation of the data we are given.
!
!    Note that the user may request extrapolation.  This occurs whenever
!    a T_INTERP value is less than the minimum T_DATA or greater than the
!    maximum T_DATA.  In that case, linear extrapolation is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
!    Input, real ( kind = 8 ) T_DATA(DATA_NUM), the value of the
!    independent variable at the sample points.  The values of T_DATA
!    must be strictly increasing.
!
!    Input, real ( kind = 8 ) P_DATA(DIM_NUM,DATA_NUM), the value of the
!    dependent variables at the sample points.
!
!    Input, integer ( kind = 4 ) INTERP_NUM, the number of points
!    at which interpolation is to be done.
!
!    Input, real ( kind = 8 ) T_INTERP(INTERP_NUM), the value of the
!    independent variable at the interpolation points.
!
!    Output, real ( kind = 8 ) P_INTERP(DIM_NUM,DATA_NUM), the interpolated
!    values of the dependent variables at the interpolation points.
!
  implicit none

  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) interp_num

  integer ( kind = 4 ) interp
  integer ( kind = 4 ) left
  real    ( kind = 8 ) p_data(dim_num,data_num)
  real    ( kind = 8 ) p_interp(dim_num,interp_num)
!  logical              r8vec_ascends_strictly
  integer ( kind = 4 ) right
  real    ( kind = 8 ) t
  real    ( kind = 8 ) t_data(data_num)
  real    ( kind = 8 ) t_interp(interp_num)

  if ( .not. r8vec_ascends_strictly ( data_num, t_data ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INTERP_LINEAR - Fatal error!'
    write ( *, '(a)' ) '  Independent variable array T_DATA is not strictly increasing.'
    stop
  end if

  do interp = 1, interp_num

    t = t_interp(interp)
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
!  nearest to, TVAL.
!
    call r8vec_bracket ( data_num, t_data, t, left, right )

    p_interp(1:dim_num,interp) = &
      ( ( t_data(right) - t                ) * p_data(1:dim_num,left)   &
      + (                 t - t_data(left) ) * p_data(1:dim_num,right) ) &
      / ( t_data(right)     - t_data(left) )

  end do

  return
end subroutine interp_linear


subroutine r8vec_bracket ( n, x, xval, left, right )

!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, real ( kind = 8 ) X(N), an array sorted into ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real    ( kind = 8 ) x(n)
  real    ( kind = 8 ) xval

  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

  return
end subroutine r8vec_bracket


subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine timestamp


subroutine read_spectrum(filename, spectrum, nu_ref, numpoint)
  ! Based on read_spectrum by Hans Kristian Eriksen, Modified by SL.
  implicit none
  
  character(len=256),                         intent(in) :: filename
  real(dp),           pointer, dimension(:,:)            :: spectrum
  real(dp),                                   intent(in) :: nu_ref
  integer(i4b),        intent(out)                       :: numpoint
  
  real(dp)            :: S_ref(1),freq_ref(1)
  integer(i4b)        :: i
  character(len=128)  :: nu, val, string

  integer(i4b) :: unit
  
  unit = 123
  freq_ref(1) = nu_ref

  open(unit, file=trim(filename))
  
  ! Find the number of entries
  numpoint = 0
  do while (.true.)
     read(unit,*,end=1) string
     
     if (string(1:1)=='#') cycle
     
     backspace(unit)
     read(unit,*,end=1) nu, val
     numpoint = numpoint + 1
  end do
1 close(unit)
  
  if (numpoint == 0) then
     write(*,*) 'No valid data entries in spectrum file ', trim(filename)
     stop
  end if
  
  allocate(spectrum(numpoint,2))
  i = 0
  open(unit, file=trim(filename))
  do while (.true.)
     read(unit,*,end=2) string
     
     if (string(1:1)=='#') cycle
     
     backspace(unit)
     read(unit,*,end=1) nu, val
     i = i+1
     
     read(nu,*)  spectrum(i,1)
     read(val,*) spectrum(i,2)
     
  end do
2 close(unit)

  !----------------
  !Rescale function 
  !----------------
  call interp_linear(1,numpoint,spectrum(:,1),spectrum(:,2),1,freq_ref,s_ref)
  do i = 1, numpoint
     spectrum(i,2) = spectrum(i,2)/s_ref(1)
  end do


end subroutine read_spectrum


end module miramare_utils
