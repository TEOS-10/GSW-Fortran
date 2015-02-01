!==========================================================================
elemental function gsw_specvol_anom (sa, ct, p)  
!==========================================================================
!
!  Calculates specific volume anomaly of seawater using the computationally
!  efficient 48-term expression for density in terms of SA, CT and p
!  (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_specvol_anom  :  specific volume anomaly of seawater (48 term equation)
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_specvol, gsw_specvol_sso_0_p

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: sa, ct, p  

real (r8) :: gsw_specvol_anom

gsw_specvol_anom = gsw_specvol(sa,ct,p) - gsw_specvol_sso_0_p(p)

return
end function

!--------------------------------------------------------------------------
