!==========================================================================
elemental function gsw_sigma4 (sa, ct) 
!==========================================================================
!
!  Calculates potential density anomaly with reference pressure of 4000 dbar,
!  this being this particular potential density minus 1000 kg/m^3.  This
!  function has inputs of Absolute Salinity and Conservative Temperature.
!  This function uses the computationally-efficient 48-term expression for 
!  density in terms of SA, CT and p (IOC et al., 2010).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! 
! gsw_sigma4  : potential density anomaly with reference pressure of 4000
!                                                      (48 term equation)
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_rho

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, ct 

real (r14) :: gsw_sigma4

gsw_sigma4 = gsw_rho(sa,ct,4d3) - 1d3

return
end function

!--------------------------------------------------------------------------
