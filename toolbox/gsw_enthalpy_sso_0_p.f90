!==========================================================================
elemental function gsw_enthalpy_sso_0_p (p) 
!==========================================================================
!
!  This function calculates enthalpy at the Standard Ocean Salinity, SSO, 
!  and at a Conservative Temperature of zero degrees C, as a function of
!  pressure, p, in dbar, using a streamlined version of the 48-term CT 
!  version of the Gibbs function, that is, a streamlined version of the 
!  code "gsw_enthalpy(SA,CT,p)".
!
! p      : sea pressure                                    [dbar]
! 
! gsw_enthalpy_sso_0_p : enthalpy(sso,0,p)
!--------------------------------------------------------------------------

use gsw_mod_rho_coefficients

use gsw_mod_teos10_constants, only : gsw_sso, gsw_sqrtsso, db2pa

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: p 

real (r8) :: gsw_enthalpy_sso_0_p
   
real (r8) :: a0, a1, a2, a3
real (r8) :: b0, b1, b2, b1sq, sqrt_disc, n, m, a, b, part

a0 = v21 + gsw_sso*(v26 + v36*gsw_sso + v31*gsw_sqrtsso)
 
a1 = v37 + v41*gsw_sso

a2 = v43

a3 = v47

b0 = v01 + gsw_sso*(v05 + v08*gsw_sqrtsso)
 
b1 = 0.5_r8*(v12 + v15*gsw_sso)

b2 = v17 + v20*gsw_sso

b1sq = b1*b1
sqrt_disc = sqrt(b1sq - b0*b2)

n = a0 + (2.0_r8*a3*b0*b1/b2 - a2*b0)/b2

m = a1 + (4.0_r8*a3*b1sq/b2 - a3*b0 - 2.0_r8*a2*b1)/b2

a = b1 - sqrt_disc
b = b1 + sqrt_disc

part = (n*b2 - m*b1)/(b2*(b - a))

gsw_enthalpy_sso_0_p = db2pa*(p*(a2 - 2.0_r8*a3*b1/b2 + 0.5_r8*a3*p)/b2 + &
                       (m/(2.0_r8*b2))*log(1.0_r8 + p*(2.0_r8*b1 + b2*p)/b0) + &
                       part*log(1.0_r8 + (b2*p*(b - a))/(a*(b + b2*p))))

return
end function

!--------------------------------------------------------------------------
