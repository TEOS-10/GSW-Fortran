!==========================================================================
elemental function gsw_alpha (sa, ct, p)  
!==========================================================================
!
!  Calculates the thermal expansion coefficient of seawater with respect to 
!  Conservative Temperature using the computationally-efficient 48-term 
!  expression for density in terms of SA, CT and p (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_alpha : thermal expansion coefficient of seawater (48 term equation)
!--------------------------------------------------------------------------

use gsw_mod_rho_coefficients

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: sa, ct, p  

real (r8) :: gsw_alpha

real (r8) :: sqrtsa, v_hat_denominator, v_hat_numerator
real (r8) :: dvhatden_dct, dvhatnum_dct

sqrtsa = sqrt(sa)

v_hat_denominator = v01 + ct*(v02 + ct*(v03 + v04*ct))  &
             + sa*(v05 + ct*(v06 + v07*ct) &
         + sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct)))) &
              + p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct) &
              + p*(v17 + ct*(v18 + v19*ct) + v20*sa))

v_hat_numerator = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct))) &
           + sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa &
       + sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))  &
            + p*(v37 + ct*(v38 + ct*(v39 + v40*ct))  &
           + sa*(v41 + v42*ct) &
            + p*(v43 + ct*(v44 + v45*ct + v46*sa) &
            + p*(v47 + v48*ct)))
       
dvhatden_dct = a01 + ct*(a02 + a03*ct) &
        + sa*(a04 + a05*ct &
    + sqrtsa*(a06 + ct*(a07 + a08*ct))) &
         + p*(a09 + a10*ct + a11*sa &
         + p*(a12 + a13*ct))

dvhatnum_dct = a14 + ct*(a15 + ct*(a16 + a17*ct)) &
        + sa*(a18 + ct*(a19 + ct*(a20 + a21*ct)) &
    + sqrtsa*(a22 + ct*(a23 + ct*(a24 + a25*ct)))) &
         + p*(a26 + ct*(a27 + a28*ct) + a29*sa &
         + p*(a30 + a31*ct + a32*sa + a33*p))
 
gsw_alpha = (v_hat_denominator*dvhatnum_dct - v_hat_numerator*dvhatden_dct)/ &
            (v_hat_numerator*v_hat_denominator)

return
end function

!--------------------------------------------------------------------------
