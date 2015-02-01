!==========================================================================
elemental function gsw_thermobaric (sa, ct, p)  
!==========================================================================
!
!  Calculates the thermobaric coefficient of seawater with respect to
!  Conservative Temperature.  This routine calculates rho from the 
!  computationally-efficient 48-term expression for density in terms of
!  SA, CT and p (IOC et al., 2010).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_thermobaric  : thermobaric coefficient with          [ 1/(K Pa) ] 
!                    respect to Conservative Temperature (48 term equation)
!--------------------------------------------------------------------------

use gsw_mod_rho_coefficients

use gsw_mod_teos10_constants, only : rec_db2pa

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: sa, ct, p  

real (r8) :: gsw_thermobaric

real (r8) :: sqrtsa, v_hat_denominator, v_hat_numerator
real (r8) :: dvhatden_dct, dvhatnum_dct, dvhatden_dp, dvhatnum_dp
real (r8) :: dvhatden_dsa, dvhatnum_dsa, dvhatden_dpdct, dvhatnum_dpdct
real (r8) :: dvhatden_dpdsa, dvhatnum_dpdsa 
real (r8) :: p1a, p1b, p1c, p1d, p1e, part1, factor2
real (r8) :: p2a, p2b, p2c, p2d, p2e, part2  

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

dvhatden_dsa = b01 + ct*(b02 + b03*ct) &
     + sqrtsa*(b04 + ct*(b05 + ct*(b06 + b07*ct))) &
          + p*(b08 + b09*ct + b10*p) 

dvhatnum_dsa = b11 + ct*(b12 + ct*(b13 + ct*(b14 + b15*ct))) &
     + sqrtsa*(b16 + ct*(b17 + ct*(b18 + ct*(b19 + b20*ct)))) + b21*sa &
          + p*(b22 + ct*(b23 + b24*p))

dvhatden_dp = c01 + ct*(c02 + c03*ct) &
        + sa*(c04 + c05*ct) &
         + p*(c06 + ct*(c07 + c08*ct) + c09*sa)

dvhatnum_dp = c10 + ct*(c11 + ct*(c12 + c13*ct)) &
        + sa*(c14 + c15*ct) &
         + p*(c16 + ct*(c17 + c18*ct + c19*sa) &
         + p*(c20 + c21*ct))

dvhatden_dpdct = c02 + 2.0_r8*c03*ct + c05*sa &
            + p*(c07 + 2.0_r8*c08*ct)
       
dvhatnum_dpdct = c11 + ct*(2.0_r8*c12 + 3.0_r8*c13*ct) + c15*sa &
            + p*(c17 + ct*2.0_r8*c18 + c19*sa + c21*p)

dvhatden_dpdsa = c04 + c05*ct + c09*p

dvhatnum_dpdsa = c14 + c15*ct + c19*ct*p

p1a = dvhatnum_dpdct/v_hat_numerator
p1b = (dvhatnum_dct*dvhatden_dp)/(v_hat_numerator*v_hat_denominator)
p1c = (dvhatnum_dp*dvhatden_dct)/(v_hat_numerator*v_hat_denominator)
p1d = (dvhatden_dp*dvhatden_dct)/(v_hat_denominator*v_hat_denominator)
p1e = dvhatden_dpdct/v_hat_denominator

part1 =  p1a - p1b - p1c + 2.0_r8*p1d - p1e

factor2 = (v_hat_denominator*dvhatnum_dct - v_hat_numerator*dvhatden_dct)/ &
           (v_hat_denominator*dvhatnum_dsa - v_hat_numerator*dvhatden_dsa)

p2a = dvhatnum_dpdsa/v_hat_numerator
p2b = (dvhatnum_dsa*dvhatden_dp)/(v_hat_numerator*v_hat_denominator)
p2c = (dvhatnum_dp*dvhatden_dsa)/(v_hat_numerator*v_hat_denominator)
p2d = (dvhatden_dp*dvhatden_dsa)/(v_hat_denominator*v_hat_denominator)
p2e = dvhatden_dpdsa/v_hat_denominator

part2 =  p2a - p2b - p2c + 2.0_r8*p2d - p2e

gsw_thermobaric = (part1 - factor2*part2)*rec_db2pa

return
end function

!--------------------------------------------------------------------------
