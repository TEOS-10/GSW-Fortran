!==========================================================================
elemental function gsw_z_from_p (p, lat) 
!==========================================================================
!
! Calculates the height z from pressure p
!
! p      : sea pressure                                    [dbar]
! lat    : latitude                                        [deg]
! 
! gsw_z_from_p : height                                    [m]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_enthalpy_sso_0

use gsw_mod_teos10_constants, only : deg2rad, gamma

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: p, lat 

real (r8) :: gsw_z_from_p

real (r8) :: sin2, b, c, a

sin2 = sin(lat*deg2rad)**2
b = 9.780327_r8*(1.0_r8 + (5.2792e-3_r8 + (2.32e-5_r8*sin2))*sin2) 
a = -0.5_r8*gamma*b 
c = gsw_enthalpy_sso_0(p)

gsw_z_from_p = -2.0_r8*c/(b + sqrt(b*b - 4.0_r8*a*c))

return 
end function

!--------------------------------------------------------------------------
