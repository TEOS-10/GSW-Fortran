!==========================================================================
module gsw_mod_baltic_data
!==========================================================================
!
! Coordinate data for the Baltic Sea

implicit none

integer, parameter :: gbd_r14 = selected_real_kind(14,30)

real (gbd_r14), dimension(2) :: xb_right, yb_right
real (gbd_r14), dimension(3) :: xb_left, yb_left

data xb_left/12.6d0, 7.d0, 26.d0/, yb_left/50.d0, 59.d0, 69.d0/
data xb_right/45.d0, 26.d0/, yb_right/50.d0, 69.d0/

end module

!--------------------------------------------------------------------------



