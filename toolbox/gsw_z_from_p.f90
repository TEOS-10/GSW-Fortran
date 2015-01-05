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

use gsw_mod_toolbox, only : gsw_enthalpy_sso_0_p

use gsw_mod_teos10_constants, only : deg2rad, gamma

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: p, lat 

real (r14) :: gsw_z_from_p

real (r14) :: sin2, b, c, a

sin2 = sin(lat*deg2rad)**2
b = 9.780327d0*(1d0 + (5.2792d-3 + (2.32d-5*sin2))*sin2) 
a = -0.5d0*gamma*b 
c = gsw_enthalpy_sso_0_p(p)

gsw_z_from_p = -2d0*c/(b + sqrt(b*b - 4d0*a*c))

return 
end function

!--------------------------------------------------------------------------
