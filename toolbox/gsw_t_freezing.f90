!==========================================================================
function gsw_t_freezing(sa,p,saturation_fraction) 
!==========================================================================

! Calculates the in-situ temperature at which of seawater freezes 
! from Absolute Salinity and pressure.
!
! sa     : Absolute Salinity                                 [g/kg]
! p      : sea pressure                                      [dbar]
! saturation_fraction : saturation fraction
!
! gsw_t_freezing : in-situ temperature freezing point  [deg C]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, p, saturation_fraction, ct_freezing, gsw_ct_freezing
real (r14) :: gsw_t_from_ct, t_freezing, gsw_t_freezing

ct_freezing = gsw_CT_freezing(sa,p,saturation_fraction)
t_freezing = gsw_t_from_ct(sa,ct_freezing,p)

if (ct_freezing.gt.9d10) then
 t_freezing = 9d15
end if

gsw_t_freezing = t_freezing

return
end function

!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! isobaric melting enthalpy and isobaric evaporation enthalpy
!--------------------------------------------------------------------------

