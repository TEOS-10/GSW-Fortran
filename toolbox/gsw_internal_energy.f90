!==========================================================================
elemental function gsw_internal_energy (sa, ct, p)  
!==========================================================================
!
!  Calculates internal energy of seawater using the computationally
!  efficient 48-term expression for density in terms of SA, CT and p
!  (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_internal_energy  :  internal_energy of seawater (48 term equation)
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_enthalpy, gsw_specvol

use gsw_mod_teos10_constants, only : gsw_p0, db2pa

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: sa, ct, p  

real (r8) :: gsw_internal_energy

gsw_internal_energy = gsw_enthalpy(sa,ct,p) - &
				(gsw_p0 + db2pa*p)*gsw_specvol(sa,ct,p)
return
end function

!--------------------------------------------------------------------------


