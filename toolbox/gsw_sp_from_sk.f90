!==========================================================================
elemental function gsw_sp_from_sk (sk)       
!==========================================================================
!
! Calculates Practical Salinity, SP, from SK
!
!  SK    : Knudsen Salinity                        [parts per thousand, ppt]
!
! gsw_sp_from_sk  : Practical Salinity                              [unitless]
!--------------------------------------------------------------------------

use gsw_mod_teos10_constants, only : gsw_soncl

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sk       

real (r14) :: gsw_sp_from_sk

gsw_sp_from_sk = (sk - 0.03d0)*(gsw_soncl/1.805d0) 

return
end function

!--------------------------------------------------------------------------
