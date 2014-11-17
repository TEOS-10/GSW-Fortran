!==========================================================================
function gsw_enthalpy_SSO_0_p(p) 
!==========================================================================

!  This function calculates enthalpy at the Standard Ocean Salinity, SSO, 
!  and at a Conservative Temperature of zero degrees C, as a function of
!  pressure, p, in dbar, using a streamlined version of the 48-term CT 
!  version of the Gibbs function, that is, a streamlined version of the 
!  code "gsw_enthalpy(SA,CT,p)".
!
! p      : sea pressure                                    [dbar]
! 
! gsw_enthalpy_SSO_0_p : enthalpy(sso,0,p)

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)
   
real (r14), parameter :: v01 =  9.998420897506056d2, v05 = -6.698001071123802d0
real (r14), parameter :: v08 = -3.988822378968490d-2, v12 = -2.233269627352527d-2
real (r14), parameter :: v15 = -1.806789763745328d-4, v17 = -3.087032500374211d-7
real (r14), parameter :: v20 =  1.550932729220080d-10, v21 =  1.0d0;
real (r14), parameter :: v26 = -7.521448093615448d-3, v31 = -3.303308871386421d-5
real (r14), parameter :: v36 =  5.419326551148740d-6, v37 = -2.742185394906099d-5
real (r14), parameter :: v41 = -1.105097577149576d-7, v43 = -1.119011592875110d-10
real (r14), parameter :: v47 = -1.200507748551599d-15, db2pa = 1d4  
real (r14), parameter :: sso = 35.16504d0, sqrtsso = 5.930011804372737d0

real (r14) :: p, gsw_enthalpy_sso_0_p, a0, a1, a2, a3
real (r14) :: b0, b1, b2, b1sq, sqrt_disc, n, m, a, b, part

a0 = v21 + sso*(v26 + v36*sso + v31*sqrtsso)
 
a1 = v37 + v41*sso

a2 = v43

a3 = v47

b0 = v01 + sso*(v05 + v08*sqrtsso)
 
b1 = 0.5*(v12 + v15*sso)

b2 = v17 + v20*sso

b1sq = b1*b1
sqrt_disc = sqrt(b1sq - b0*b2)

n = a0 + (2d0*a3*b0*b1/b2 - a2*b0)/b2

m = a1 + (4d0*a3*b1sq/b2 - a3*b0 - 2*a2*b1)/b2

a = b1 - sqrt_disc
b = b1 + sqrt_disc

part = (n*b2 - m*b1)/(b2*(b - a))

gsw_enthalpy_sso_0_p = db2pa*(p*(a2 - 2d0*a3*b1/b2 + 0.5d0*a3*p)/b2 + &
          (m/(2d0*b2))*log(1d0 + p*(2d0*b1 + b2*p)/b0) + &
           part*log(1d0 + (b2*p*(b - a))/(a*(b + b2*p))))

return
end function

!--------------------------------------------------------------------------

