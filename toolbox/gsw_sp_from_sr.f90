!==========================================================================
elemental function gsw_sp_from_sr (sr)  
!==========================================================================
!
! Calculates Practical Salinity, sp, from Reference Salinity, sr. 
!
! sr     : Reference Salinity                              [g/kg]
!
! gsw_sp_from_sr  : Practical Salinity                     [unitless]
!--------------------------------------------------------------------------

use gsw_mod_teos10_constants, only : gsw_ups

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sr  

real (r14) :: gsw_sp_from_sr

gsw_sp_from_sr = sr/gsw_ups

return
end function

!--------------------------------------------------------------------------



