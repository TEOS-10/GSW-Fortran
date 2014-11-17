!==========================================================================
function gsw_sr_from_sp(sp) 
!==========================================================================

! Calculates Reference Salinity, SR, from Practical Salinity, SP. 
!
! sp     : Practical Salinity                              [unitless]
!
! gsw_sr_from_sp : Reference Salinity                      [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sp, gsw_sr_from_sp

gsw_sr_from_sp = 1.004715428571429d0*sp;

if (gsw_sr_from_sp.ge.1d10) then
    gsw_sr_from_sp = 9d15
end if

return
end function

!--------------------------------------------------------------------------

