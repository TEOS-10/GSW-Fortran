!==========================================================================
elemental function gsw_kappa (sa, ct, p)  
!==========================================================================
!
!  Calculates isentropic compressibility of seawater using the computationally
!  efficient 48-term expression for density in terms of SA, CT and p
!  (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_kappa  :  isentropic compressibility (48 term equation)
!--------------------------------------------------------------------------

use gsw_mod_rho_coefficients

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: sa, ct, p  

real (r8) :: gsw_kappa

real (r8) :: v_hat_denominator, v_hat_numerator
real (r8) :: sqrtsa, dvden_dp, dvnum_dp

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

dvden_dp = c01 + ct*(c02 + c03*ct) &
    + sa*(c04 + c05*ct) &
    + p*(c06 + ct*(c07 + c08*ct) + c09*sa)

dvnum_dp = c10 + ct*(c11 + ct*(c12 + c13*ct)) &
    + sa*(c14 + c15*ct) &
    + p*(c16 + ct*(c17 + c18*ct + c19*sa) &
    + p*(c20 + c21*ct))

gsw_kappa = 1e-4_r8*(dvden_dp/v_hat_denominator  -  dvnum_dp/v_hat_numerator)

return
end function

!--------------------------------------------------------------------------
