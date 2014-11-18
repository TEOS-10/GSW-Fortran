!==========================================================================
elemental function gsw_dynamic_enthalpy (sa, ct, p) 
!==========================================================================
!
! Calculates dynamic enthalpy of seawater using the computationally
! efficient 48-term expression for density in terms of SA, CT and p
! (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_dynamic_enthalpy  :  dynamic enthalpy of seawater (48 term equation)
!--------------------------------------------------------------------------

use gsw_mod_teos10_constants, only : db2pa

use gsw_mod_rho_coefficients

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, ct, p 

real (r14) :: gsw_dynamic_enthalpy

real (r14) :: sqrtsa, a0, a1, a2, a3, b0, b1, b2, b1sq
real (r14) :: sqrt_disc, ca, cb, cn, cm, part

sqrtsa = sqrt(sa)

a0 = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct))) &
         + sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa  &
     + sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))
 
a1 = v37 + ct*(v38 + ct*(v39 + v40*ct)) + sa*(v41 + v42*ct)

a2 = v43 + ct*(v44 + v45*ct + v46*sa)

a3 = v47 + v48*ct

b0 = v01 + ct*(v02 + ct*(v03 + v04*ct))  &
         + sa*(v05 + ct*(v06 + v07*ct)  &
     + sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))))
 
b1 = 0.5d0*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct))

b2 = v17 + ct*(v18 + v19*ct) + v20*sa

b1sq = b1*b1 
sqrt_disc = sqrt(b1sq - b0*b2)

cn = a0 + (2d0*a3*b0*b1/b2 - a2*b0)/b2

cm = a1 + (4d0*a3*b1sq/b2 - a3*b0 - 2d0*a2*b1)/b2

ca = b1 - sqrt_disc
cb = b1 + sqrt_disc

part = (cn*b2 - cm*b1)/(b2*(cb - ca))

gsw_dynamic_enthalpy = db2pa*(p*(a2 - 2d0*a3*b1/b2 + 0.5d0*a3*p)/b2  &
                       + (cm/(2d0*b2))*log(1d0 + p*(2d0*b1 + b2*p)/b0)  &
                       + part*log(1d0 + (b2*p*(cb - ca))/(ca*(cb + b2*p))))

return
end function

!--------------------------------------------------------------------------
