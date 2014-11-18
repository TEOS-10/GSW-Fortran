!==========================================================================
elemental function gsw_sr_from_sp (sp) 
!==========================================================================
!
! Calculates Reference Salinity, SR, from Practical Salinity, SP. 
!
! sp     : Practical Salinity                              [unitless]
!
! gsw_sr_from_sp : Reference Salinity                      [g/kg]
!--------------------------------------------------------------------------

use gsw_mod_teos10_constants, only : gsw_ups

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sp 

real (r14) :: gsw_sr_from_sp

gsw_sr_from_sp = sp*gsw_ups

return
end function

!--------------------------------------------------------------------------



