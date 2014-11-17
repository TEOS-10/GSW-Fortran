!==========================================================================
function gsw_sigma1(sa,ct) 
!==========================================================================

!  Calculates potential density anomaly with reference pressure of 1000 dbar,
!  this being this particular potential density minus 1000 kg/m^3.  This
!  function has inputs of Absolute Salinity and Conservative Temperature.
!  This function uses the computationally-efficient 48-term expression for 
!  density in terms of SA, CT and p (IOC et al., 2010).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! 
! gsw_sigma1  : potential density anomaly with reference pressure of 1000
!                                                      (48 term equation)

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, ct, gsw_sigma1, gsw_rho

gsw_sigma1 = gsw_rho(sa,ct,1000d0) - 1000

return
end function

!--------------------------------------------------------------------------

