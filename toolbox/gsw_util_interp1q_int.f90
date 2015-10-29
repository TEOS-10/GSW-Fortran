!==========================================================================
pure function gsw_util_interp1q_int (x, iy, x_i) result(y_i)
!==========================================================================
! Returns the value of the 1-D function iy (integer) at the points of column
! vector x_i using linear interpolation. The vector x specifies the
! coordinates of the underlying interval.
!==========================================================================

use gsw_mod_toolbox, only : gsw_util_sort_real

use gsw_mod_kinds

implicit none

integer, intent(in) :: iy(:)
real (r8), intent(in) :: x(:), x_i(:)

real (r8) :: y_i(size(x_i))

logical, allocatable :: in_rng(:)
integer, allocatable :: j(:), k(:), r(:), jrev(:), ki(:)
real (r8), allocatable :: xi(:), xxi(:), u(:)

integer :: imax_x(1), imin_x(1), i, n, m, ii, nx, nxi
real (r8) :: max_x, min_x

min_x = minval(x); imin_x = minloc(x)
max_x = maxval(x); imax_x = maxloc(x)

where (x_i .le. min_x) y_i = iy(imin_x(1))
where (x_i .ge. max_x) y_i = iy(imax_x(1))

nx = size(x)
nxi = size(x_i)

allocate (in_rng(nxi))
in_rng = (x_i .gt. min_x .and. x_i .lt. max_x)

n = count(in_rng)
if (n .eq. 0) return

allocate (xi(n), k(n), ki(n), r(n), u(n))
m = nx + n
allocate (xxi(m), j(m), jrev(m))

ii = 0
do i = 1, nxi
    if (in_rng(i)) then
        ii = ii + 1
        xi(ii) = x_i(i)
        ki(ii) = i
    end if
end do

k = gsw_util_sort_real(xi)
xxi = (/ x, xi(k) /)
j = gsw_util_sort_real(xxi)

jrev(j) = (/ (i, i=1,size(j)) /)
r(k) = jrev(nx+1:) - (/ (i, i=1,size(xi)) /)

u = (xi-x(r))/(x(r+1)-x(r))
y_i(ki) = iy(r) + (iy(r+1)-iy(r))*u

return

end function

!--------------------------------------------------------------------------
