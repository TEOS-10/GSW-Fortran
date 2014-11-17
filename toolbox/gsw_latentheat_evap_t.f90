!==========================================================================
function gsw_latentheat_evap_t(sa,t)  
!==========================================================================

! Calculates latent heat, or enthalpy, of evaporation.
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! 
! gsw_latentheat_evap_t : latent heat of evaporation       [J/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, t, ct, gsw_ct_from_pt, gsw_latentheat_evap_ct
real (r14) :: gsw_latentheat_evap_t

ct = gsw_ct_from_pt(sa,t)

gsw_latentheat_evap_t = gsw_latentheat_evap_ct(sa,ct)

return
end

!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! planet Earth properties
!--------------------------------------------------------------------------

