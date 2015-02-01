!==========================================================================
elemental function gsw_specvol (sa, ct, p) 
!==========================================================================
!
!  Calculates specific volume of seawater using the computationally-
!  efficient 48-term expression for density in terms of SA, CT and p
!  (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_specvol  :  specific volume of seawater (48 term equation)
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_rho

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: sa, ct, p 

real (r8) :: gsw_specvol

gsw_specvol = 1.0_r8/gsw_rho(sa,ct,p)

return
end function

!--------------------------------------------------------------------------
