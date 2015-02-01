!==========================================================================
elemental function gsw_specvol_sso_0_p (p) 
!==========================================================================
!
!  This function calculates specifc volume at the Standard Ocean Salinty,
!  SSO, and at a Conservative Temperature of zero degrees C, as a function 
!  of pressure, p, in dbar, using a streamlined version of the 48-term CT
!  version of specific volume, that is, a streamlined version of the code
!  "gsw_specvol(SA,CT,p)".
!
! p      : sea pressure                                    [dbar]
! 
! gsw_specvol_sso_0_p : specvol(sso,0,p)
!--------------------------------------------------------------------------

use gsw_mod_rho_coefficients

use gsw_mod_teos10_constants, only : gsw_sso, gsw_sqrtsso

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: p 

real (r8) :: gsw_specvol_sso_0_p

gsw_specvol_sso_0_p = (v21 + gsw_sso*(v26 + v36*gsw_sso + v31*gsw_sqrtsso)  &
                  + p*(v37 + v41*gsw_sso + p*(v43 + v47*p )))/ &
                      (v01 + gsw_sso*(v05 + v08*gsw_sqrtsso) &
                  + p*(v12 + v15*gsw_sso + p*(v17 + v20*gsw_sso)))

return
end function

!--------------------------------------------------------------------------
