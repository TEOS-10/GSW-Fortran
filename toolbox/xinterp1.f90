!==========================================================================
function xinterp1(x,y,n,x0)
!==========================================================================

! Linearly interpolate a real array   
!
! x      : y array (Must be monotonic)               
! y      : y array     
! n      : length of X and Y arrays
! x0     : value to be interpolated
!
! xinterp1 : Linearly interpolated value

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n, k

real (r14), dimension(n) :: x, y
real (r14) :: x0, r, xinterp1

call indx(x,n,x0,k)
r = (x0-x(k))/(x(k+1)-x(k))
xinterp1 = y(k) + r*(y(k+1)-y(k))

return
end function

!--------------------------------------------------------------------------
