!==========================================================================
pure function gsw_util_xinterp1 (x, y, n, x0)
!==========================================================================
!
! Linearly interpolate a real array   
!
! x      : x array (must be monotonic)               
! y      : y array     
! n      : length of X and Y arrays
! x0     : value to be interpolated
!
! xinterp1 : Linearly interpolated value
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_util_indx

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

integer, intent(in) :: n
real (r14), intent(in) :: x0
real (r14), dimension(n), intent(in) :: x, y

real (r14) :: gsw_util_xinterp1

integer :: k

real (r14) :: r

call gsw_util_indx(x,n,x0,k)
r = (x0-x(k))/(x(k+1)-x(k))
gsw_util_xinterp1 = y(k) + r*(y(k+1)-y(k))

return
end function

!--------------------------------------------------------------------------


