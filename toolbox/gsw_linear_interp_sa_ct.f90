!==========================================================================
pure subroutine gsw_linear_interp_sa_ct (sa, ct, p, p_i, sa_i, ct_i)
!==========================================================================
! This function interpolates the cast with respect to the interpolating 
! variable p. This function finds the values of SA, CT at p_i on this cast.
!
! VERSION NUMBER: 3.05 (27th January 2015)
!
! This fuction was adapted from Matlab's interp1q.
!==========================================================================

use gsw_mod_toolbox, only : gsw_util_sort_real

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: sa(:), ct(:), p(:), p_i(:)
real (r8), intent(out) :: sa_i(:), ct_i(:)

logical, allocatable :: in_rng(:)
integer, allocatable :: j(:), k(:), r(:), jrev(:), ki(:)
real (r8), allocatable :: xi(:), xxi(:), u(:)

integer :: imax_p(1), imin_p(1), i, n, m, ii, np, npi
real (r8) :: max_p, min_p

np = size(p)
npi = size(p_i)

min_p = minval(p); imin_p = minloc(p)
max_p = maxval(p); imax_p = maxloc(p)

where (p_i .le. min_p)
    sa_i = sa(imin_p(1))  ! set equal to the shallowest bottle.
    ct_i = ct(imin_p(1))
end where
where (p_i .ge. max_p)
    sa_i = sa(imax_p(1))  ! set equal to the deepest bottle.
    ct_i = ct(imax_p(1))
end where

allocate (in_rng(npi))
in_rng = (p_i .gt. min_p .and. p_i .lt. max_p)

n = count(in_rng)
if (n .eq. 0) return

allocate (xi(n), k(n), ki(n), r(n), u(n))
m = np + n
allocate (xxi(m), j(m), jrev(m))

ii = 0
do i = 1, npi
    if (in_rng(i)) then
        ii = ii + 1
        xi(ii) = p_i(i)
        ki(ii) = i
    end if
end do

k = gsw_util_sort_real(xi)
xxi = (/ p, xi(k) /)
j = gsw_util_sort_real(xxi)

jrev(j) = (/ (i, i=1,size(j)) /)
r(k) = jrev(np+1:) - (/ (i, i=1,size(xi)) /)

u = (xi-p(r))/(p(r+1)-p(r))
sa_i(ki) = sa(r) + (sa(r+1)-sa(r))*u
ct_i(ki) = ct(r) + (ct(r+1)-ct(r))*u

return

end subroutine

!--------------------------------------------------------------------------
