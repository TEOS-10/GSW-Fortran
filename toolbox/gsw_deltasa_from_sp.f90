!==========================================================================
function gsw_deltasa_from_sp(sp,p,long,lat) 
!==========================================================================

! Calculates Absolute Salinity Anomaly, deltaSA, from Practical Salinity, SP. 
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_deltasa_from_sp : Absolute Salinty Anomaly           [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sp, long, lat, p, gsw_sa_from_sp, gsw_sr_from_sp
real (r14) :: gsw_deltasa_from_sp

gsw_deltasa_from_sp = gsw_sa_from_sp(sp,p,long,lat) - gsw_sr_from_sp(sp)

if (gsw_deltasa_from_sp.gt.1d10) then
    gsw_deltasa_from_sp = 9d15
end if

return
end function

!--------------------------------------------------------------------------

