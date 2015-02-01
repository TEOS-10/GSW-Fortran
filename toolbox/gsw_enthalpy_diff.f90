!==========================================================================
elemental function gsw_enthalpy_diff (sa, ct, p_shallow, p_deep)
!==========================================================================
!
!  Calculates the difference of the specific enthalpy of seawater between 
!  two different pressures, p_deep (the deeper pressure) and p_shallow
!  (the shallower pressure), at the same values of SA and CT.  This 
!  function uses the computationally-efficient 48-term expression for
!  density in terms of SA, CT and p (IOC et al., 2010).  The output
!  (enthalpy_diff_CT) is the specific enthalpy evaluated at (SA,CT,p_deep)
!  minus the specific enthalpy at (SA,CT,p_shallow). 
!
!  SA         =  Absolute Salinity                                 [ g/kg ]
!  CT         =  Conservative Temperature (ITS-90)                [ deg C ]
!  p_shallow  =  upper sea pressure                                [ dbar ]
!                ( i.e. shallower absolute pressure - 10.1325 dbar ) 
!  p_deep     =  lower sea pressure                                [ dbar ]
!                ( i.e. deeper absolute pressure - 10.1325 dbar )
!
!  enthalpy_diff_CT  =  difference of specific enthalpy            [ J/kg ]
!                       (deep minus shallow)
!--------------------------------------------------------------------------

use gsw_mod_rho_coefficients

use gsw_mod_teos10_constants, only : db2pa

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: sa, ct, p_shallow, p_deep

real (r8) :: gsw_enthalpy_diff

real (r8) :: a, a0, a1, a2, a3, b, b0, b1, b1sq, b2, delta_p
real (r8) :: m, n, part1, part2, part3, p_sum, sqrt_disc, sqrtsa

sqrtsa = sqrt(sa)

a0 = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct))) &
         + sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa &
     + sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))

a1 = v37 + ct*(v38 + ct*(v39 + v40*ct)) + sa*(v41 + v42*ct)

a2 = v43 + ct*(v44 + v45*ct + v46*sa)

a3 = v47 + v48*ct

b0 = v01 + ct*(v02 + ct*(v03 + v04*ct))  &
         + sa*(v05 + ct*(v06 + v07*ct) &
     + sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))))
 
b1 = 0.5_r8*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct))

b2 = v17 + ct*(v18 + v19*ct) + v20*sa

b1sq = b1*b1 
sqrt_disc = sqrt(b1sq - b0*b2)

n = a0 + (2.0_r8*a3*b0*b1/b2 - a2*b0)/b2

m = a1 + (4.0_r8*a3*b1sq/b2 - a3*b0 - 2.0_r8*a2*b1)/b2

a = b1 - sqrt_disc
b = b1 + sqrt_disc
delta_p = p_deep - p_shallow
p_sum = p_deep + p_shallow
part1 = b0 + p_shallow*(2.0_r8*b1 + b2*p_shallow)

part2 = (b + b2*p_deep)*(a + b2*p_shallow)

part3 = (n*b2 - m*b1)/(b2*(b - a))

gsw_enthalpy_diff = db2pa* &
         (delta_p*(a2 - 2.0_r8*a3*b1/b2 + 0.5_r8*a3*p_sum)/b2 + &
         (m/(2.0_r8*b2))*log(1.0_r8 + delta_p*(2.0_r8*b1 + b2*p_sum)/part1) + & 
         part3*log(1.0_r8 + delta_p*b2*(b - a)/part2))

return
end function

!--------------------------------------------------------------------------
