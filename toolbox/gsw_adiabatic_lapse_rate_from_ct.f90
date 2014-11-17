!==========================================================================
function gsw_adiabatic_lapse_rate_from_ct(sa,ct,p) 
!==========================================================================

! Calculates the adiabatic lapse rate from Conservative Temperature
!
! sa     : Absolute Salinity                                 [g/kg]
! ct     : Conservative Temperature                          [deg C]
! p      : sea pressure                                      [dbar]
! 
! gsw_adiabatic_lapse_rate_from_ct : adiabatic lapse rate    [K/Pa]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n0, n1, n2 

real (r14) :: sa, ct, p, gsw_adiabatic_lapse_rate_from_ct, gsw_gibbs
real (r14) :: gsw_pt_from_ct, gsw_pt_from_t, pt0, t, pr0

n0 = 0
n1 = 1
n2 = 2

pr0 = 0d0
pt0 = gsw_pt_from_ct(sa,ct)
t = gsw_pt_from_t(sa,pt0,pr0,p)

gsw_adiabatic_lapse_rate_from_ct = -gsw_gibbs(n0,n1,n1,sa,t,p)/(gsw_gibbs(n0,n2,n0,SA,t,p))

return
end

!--------------------------------------------------------------------------


!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! density and enthalpy, based on the 48-term expression for density
!--------------------------------------------------------------------------

