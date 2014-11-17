!==========================================================================
function gsw_specvol_sso_0_p(p) 
!==========================================================================

!  This function calculates specifc volume at the Standard Ocean Salinty,
!  SSO, and at a Conservative Temperature of zero degrees C, as a function 
!  of pressure, p, in dbar, using a streamlined version of the 48-term CT
!  version of specific volume, that is, a streamlined version of the code
!  "gsw_specvol(SA,CT,p)".
!
! p      : sea pressure                                    [dbar]
! 
! gsw_specvol_sso_0_p : specvol(sso,0,p)

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: v01 =  9.998420897506056d2, v05 = -6.698001071123802d0
real (r14), parameter :: v08 = -3.988822378968490d-2, v12 = -2.233269627352527d-2
real (r14), parameter :: v15 = -1.806789763745328d-4, v17 = -3.087032500374211d-7
real (r14), parameter :: v20 =  1.550932729220080d-10, v21 =  1.0d0;
real (r14), parameter :: v26 = -7.521448093615448d-3, v31 = -3.303308871386421d-5
real (r14), parameter :: v36 =  5.419326551148740d-6, v37 = -2.742185394906099d-5
real (r14), parameter :: v41 = -1.105097577149576d-7, v43 = -1.119011592875110d-10
real (r14), parameter :: v47 = -1.200507748551599d-15

real (r14) :: sso, sqrtsso, p, gsw_specvol_sso_0_p

sso = 35.16504d0
sqrtsso = 5.930011804372737d0     ! sqrt(SSO) = 5.930011804372737

gsw_specvol_sso_0_p = (v21 + sso*(v26 + v36*sso + v31*sqrtsso)  &
             + p*(v37 + v41*sso + p*(v43 + v47*p )))/ &
             (v01 + sso*(v05 + v08*sqrtsso) &
             + p*(v12 + v15*sso + p*(v17 + v20*sso)))

return
end function

!--------------------------------------------------------------------------

