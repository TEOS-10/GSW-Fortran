!==========================================================================
elemental function gsw_cabbeling (sa, ct, p)  
!==========================================================================
!
!  Calculates the cabbeling coefficient of seawater with respect to  
!  Conservative Temperature.  This function uses the computationally-
!  efficient 48-term expression for density in terms of SA, CT and p
!  (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_cabbeling  : cabbeling coefficient with respect to            [ 1/K^2 ]
!                  Conservative Temperature. (48 term equation)
!--------------------------------------------------------------------------

use gsw_mod_rho_coefficients

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: sa, ct, p  

real (r8) :: gsw_cabbeling

real (r8) :: sqrtsa, v_hat_denominator, v_hat_numerator
real (r8) :: dvhatden_dct, dvhatnum_dct, dvhatden_dctdct, dvhatden_dctdsa
real (r8) :: dvhatden_dsa, dvhatnum_dsa, dvhatnum_dsadsa, dvhatden_dsadsa
real (r8) :: dvhatnum_dctdct, dvhatnum_dsadct, dvhatnum_dctdsa
real (r8) :: p1a, p1b, p1c, p1d, part1, factor2a, factor2b
real (r8) :: p2a, p2b, p2c, p2d, p2e, part2 
real (r8) :: factor3a, factor3b, p3a, p3b, p3c, p3d, part3 

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

dvhatden_dctdct = a02 + 2.0_r8*a03*ct &
            + sa*(a05 + sqrtsa*(a07 + 2.0_r8*a08*ct)) &
             + p*(a10 + a13*p)
     
dvhatden_dctdsa = a04 + a05*ct &
  + (3.0_r8/2.0_r8)*sqrtsa*(a06 + ct*(a07 + a08*ct)) &
                + a11*p

dvhatden_dsa = b01 + ct*(b02 + b03*ct) &
     + sqrtsa*(b04 + ct*(b05 + ct*(b06 + b07*ct))) &
          + p*(b08 + b09*ct + b10*p) 

dvhatnum_dsa = b11 + ct*(b12 + ct*(b13 + ct*(b14 + b15*ct))) &
     + sqrtsa*(b16 + ct*(b17 + ct*(b18 + ct*(b19 + b20*ct)))) + b21*sa &
          + p*(b22 + ct*(b23 + b24*p))

dvhatnum_dctdct = a15 + ct*(2.0_r8*a16 + 3.0_r8*a17*ct) &
           + sa*(a19 + ct*(2.0_r8*a20 + 3.0_r8*a21*ct) &
        +sqrtsa*(a23 + ct*(2.0_r8*a24 + 3.0_r8*a25*ct))) &
           + p*(a27 + 2.0_r8*a28*ct + a31*p)

dvhatnum_dctdsa = a18 + ct*(a19 + ct*(a20 + a21*ct)) &
 + (3.0_r8/2.0_r8)*sqrtsa*(a22 + ct*(a23 + ct*(a24 + a25*ct))) &
            + p*(a29 + p*a32)

dvhatnum_dsadsa = (1.0_r8/(2.0_r8*sqrtsa)) &
                  *(b16 + ct*(b17 + ct*(b18 + ct*(b19 + b20*ct)))) + b21

dvhatden_dsadsa = (1.0_r8/(2.0_r8*sqrtsa))*(b04 + ct*(b05 + ct*(b06 + b07*ct)))

p1a = dvhatnum_dctdct/v_hat_numerator
p1b = (dvhatnum_dct*dvhatden_dct)/(v_hat_numerator*v_hat_denominator)
p1c = dvhatden_dctdct/v_hat_denominator
p1d = dvhatden_dct/v_hat_denominator

part1 = p1a - 2.0_r8*p1b - p1c + 2.0_r8*p1d*p1d

factor2a = (v_hat_denominator*dvhatnum_dct - v_hat_numerator*dvhatden_dct)
factor2b = (v_hat_denominator*dvhatnum_dsa - v_hat_numerator*dvhatden_dsa)

p2a = dvhatnum_dctdsa/v_hat_numerator
p2b = (dvhatnum_dct*dvhatden_dsa)/(v_hat_numerator*v_hat_denominator)
p2c = (dvhatnum_dsa*dvhatden_dct)/(v_hat_numerator*v_hat_denominator)
p2d = dvhatden_dctdsa/v_hat_denominator
p2e = (dvhatden_dct*dvhatden_dsa)/(v_hat_denominator*v_hat_denominator)

part2 = p2a - p2b - p2c - p2d + 2.0_r8*p2e

factor3a = (v_hat_denominator*dvhatnum_dct - v_hat_numerator*dvhatden_dct)
factor3b = (v_hat_denominator*dvhatnum_dsa - v_hat_numerator*dvhatden_dsa)

p3a = dvhatnum_dsadsa/v_hat_numerator
p3b = (dvhatnum_dsa*dvhatden_dsa)/(v_hat_numerator*v_hat_denominator)
p3c = (v_hat_numerator*dvhatden_dsadsa)/(v_hat_numerator*v_hat_denominator)
p3d = dvhatden_dsa/v_hat_denominator

part3 = p3a - 2.0_r8*p3b - p3c + 2.0_r8*p3d*p3d

gsw_cabbeling = part1 - 2.0_r8*(factor2a/factor2b)*part2 &
                + (factor3a/factor3b)*(factor3a/factor3b)*part3

return
end function

!--------------------------------------------------------------------------
