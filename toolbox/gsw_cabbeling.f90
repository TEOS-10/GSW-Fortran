!==========================================================================
function gsw_cabbeling(sa,ct,p)  
!==========================================================================

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

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: v01 =  9.998420897506056d2,   v02 =  2.839940833161907d0
real (r14), parameter :: v03 = -3.147759265588511d-2,  v04 =  1.181805545074306d-3
real (r14), parameter :: v05 = -6.698001071123802d0,   v06 = -2.986498947203215d-2
real (r14), parameter :: v07 =  2.327859407479162d-4,  v08 = -3.988822378968490d-2
real (r14), parameter :: v09 =  5.095422573880500d-4,  v10 = -1.426984671633621d-5
real (r14), parameter :: v11 =  1.645039373682922d-7,  v12 = -2.233269627352527d-2
real (r14), parameter :: v13 = -3.436090079851880d-4,  v14 =  3.726050720345733d-6
real (r14), parameter :: v15 = -1.806789763745328d-4,  v16 =  6.876837219536232d-7
real (r14), parameter :: v17 = -3.087032500374211d-7,  v18 = -1.988366587925593d-8
real (r14), parameter :: v19 = -1.061519070296458d-11, v20 =  1.550932729220080d-10
real (r14), parameter :: v21 =  1.0d0, v22 = 2.775927747785646d-3, v23 = -2.349607444135925d-5
real (r14), parameter :: v24 =  1.119513357486743d-6,  v25 =  6.743689325042773d-10
real (r14), parameter :: v26 = -7.521448093615448d-3,  v27 = -2.764306979894411d-5
real (r14), parameter :: v28 =  1.262937315098546d-7,  v29 =  9.527875081696435d-10
real (r14), parameter :: v30 = -1.811147201949891d-11, v31 = -3.303308871386421d-5
real (r14), parameter :: v32 =  3.801564588876298d-7,  v33 = -7.672876869259043d-9
real (r14), parameter :: v34 = -4.634182341116144d-11, v35 =  2.681097235569143d-12
real (r14), parameter :: v36 =  5.419326551148740d-6,  v37 = -2.742185394906099d-5
real (r14), parameter :: v38 = -3.212746477974189d-7,  v39 =  3.191413910561627d-9
real (r14), parameter :: v40 = -1.931012931541776d-12, v41 = -1.105097577149576d-7
real (r14), parameter :: v42 =  6.211426728363857d-10, v43 = -1.119011592875110d-10
real (r14), parameter :: v44 = -1.941660213148725d-11, v45 = -1.864826425365600d-14
real (r14), parameter :: v46 =  1.119522344879478d-14, v47 = -1.200507748551599d-15
real (r14), parameter :: v48 =  6.057902487546866d-17 
real (r14), parameter :: a01 =  2.839940833161907d0, a02 = -6.295518531177023d-2
real (r14), parameter :: a03 =  3.545416635222918d-3, a04 = -2.986498947203215d-2
real (r14), parameter :: a05 =  4.655718814958324d-4, a06 =  5.095422573880500d-4
real (r14), parameter :: a07 = -2.853969343267241d-5, a08 =  4.935118121048767d-7
real (r14), parameter :: a09 = -3.436090079851880d-4, a10 =  7.452101440691467d-6
real (r14), parameter :: a11 =  6.876837219536232d-7, a12 = -1.988366587925593d-8
real (r14), parameter :: a13 = -2.123038140592916d-11, a14 =  2.775927747785646d-3
real (r14), parameter :: a15 = -4.699214888271850d-5, a16 =  3.358540072460230d-6
real (r14), parameter :: a17 =  2.697475730017109d-9, a18 = -2.764306979894411d-5
real (r14), parameter :: a19 =  2.525874630197091d-7, a20 =  2.858362524508931d-9
real (r14), parameter :: a21 = -7.244588807799565d-11, a22 =  3.801564588876298d-7
real (r14), parameter :: a23 = -1.534575373851809d-8, a24 = -1.390254702334843d-10
real (r14), parameter :: a25 =  1.072438894227657d-11, a26 = -3.212746477974189d-7
real (r14), parameter :: a27 =  6.382827821123254d-9, a28 = -5.793038794625329d-12
real (r14), parameter :: a29 =  6.211426728363857d-10, a30 = -1.941660213148725d-11
real (r14), parameter :: a31 = -3.729652850731201d-14, a32 =  1.119522344879478d-14
real (r14), parameter :: a33 =  6.057902487546866d-17
real (r14), parameter :: b01 = -6.698001071123802d0, b02 = -2.986498947203215d-2
real (r14), parameter :: b03 =  2.327859407479162d-4, b04 = -5.983233568452735d-2
real (r14), parameter :: b05 =  7.643133860820750d-4, b06 = -2.140477007450431d-5
real (r14), parameter :: b07 =  2.467559060524383d-7, b08 = -1.806789763745328d-4
real (r14), parameter :: b09 =  6.876837219536232d-7, b10 =  1.550932729220080d-10
real (r14), parameter :: b11 = -7.521448093615448d-3, b12 = -2.764306979894411d-5
real (r14), parameter :: b13 =  1.262937315098546d-7, b14 =  9.527875081696435d-10
real (r14), parameter :: b15 = -1.811147201949891d-11, b16 = -4.954963307079632d-5
real (r14), parameter :: b17 =  5.702346883314446d-7, b18 = -1.150931530388857d-8
real (r14), parameter :: b19 = -6.951273511674217d-11, b20 =  4.021645853353715d-12
real (r14), parameter :: b21 =  1.083865310229748d-5, b22 = -1.105097577149576d-7
real (r14), parameter :: b23 =  6.211426728363857d-10, b24 =  1.119522344879478d-14

real (r14) :: sa, ct, p, gsw_cabbeling, sqrtsa, v_hat_denominator, v_hat_numerator
real (r14) :: dvhatden_dct, dvhatnum_dct, dvhatden_dctdct, dvhatden_dctdsa
real (r14) :: dvhatden_dsa, dvhatnum_dsa, dvhatnum_dsadsa, dvhatden_dsadsa
real (r14) :: dvhatnum_dctdct, dvhatnum_dsadct, dvhatnum_dctdsa
real (r14) :: p1a, p1b, p1c, p1d, part1, factor2a, factor2b
real (r14) :: p2a, p2b, p2c, p2d, p2e, part2 
real (r14) :: factor3a, factor3b, p3a, p3b, p3c, p3d, part3 

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

dvhatden_dctdct = a02 + 2d0*a03*ct &
            + sa*(a05 + sqrtsa*(a07 + 2d0*a08*ct)) &
             + p*(a10 + a13*p)
     
dvhatden_dctdsa = a04 + a05*ct &
  + (3d0/2d0)*sqrtsa*(a06 + ct*(a07 + a08*ct)) &
                + a11*p

dvhatden_dsa = b01 + ct*(b02 + b03*ct) &
     + sqrtsa*(b04 + ct*(b05 + ct*(b06 + b07*ct))) &
          + p*(b08 + b09*ct + b10*p) 

dvhatnum_dsa = b11 + ct*(b12 + ct*(b13 + ct*(b14 + b15*ct))) &
     + sqrtsa*(b16 + ct*(b17 + ct*(b18 + ct*(b19 + b20*ct)))) + b21*sa &
          + p*(b22 + ct*(b23 + b24*p))

dvhatnum_dctdct = a15 + ct*(2*a16 + 3*a17*ct) &
           + sa*(a19 + ct*(2*a20 + 3*a21*ct) &
        +sqrtsa*(a23 + ct*(2*a24 + 3*a25*ct))) &
           + p*(a27 + 2*a28*ct + a31*p)

dvhatnum_dctdsa = a18 + ct*(a19 + ct*(a20 + a21*ct)) &
 + (3d0/2d0)*sqrtsa*(a22 + ct*(a23 + ct*(a24 + a25*ct))) &
            + p*(a29 + p*a32)

dvhatnum_dsadsa = (1d0/(2d0*sqrtsa))*(b16 + ct*(b17 + ct*(b18 + ct*(b19 + b20*ct)))) + b21

dvhatden_dsadsa = (1d0/(2d0*sqrtsa))*(b04 + ct*(b05 + ct*(b06 + b07*ct)))

p1a = dvhatnum_dctdct/v_hat_numerator
p1b = (dvhatnum_dct*dvhatden_dct)/(v_hat_numerator*v_hat_denominator)
p1c = dvhatden_dctdct/v_hat_denominator
p1d = dvhatden_dct/v_hat_denominator

part1 = p1a - 2d0*p1b - p1c + 2d0*p1d*p1d

factor2a = (v_hat_denominator*dvhatnum_dct - v_hat_numerator*dvhatden_dct)
factor2b = (v_hat_denominator*dvhatnum_dsa - v_hat_numerator*dvhatden_dsa)

p2a = dvhatnum_dctdsa/v_hat_numerator
p2b = (dvhatnum_dct*dvhatden_dsa)/(v_hat_numerator*v_hat_denominator)
p2c = (dvhatnum_dsa*dvhatden_dct)/(v_hat_numerator*v_hat_denominator)
p2d = dvhatden_dctdsa/v_hat_denominator
p2e = (dvhatden_dct*dvhatden_dsa)/(v_hat_denominator*v_hat_denominator)

part2 = p2a - p2b - p2c - p2d + 2d0*p2e

factor3a = (v_hat_denominator*dvhatnum_dct - v_hat_numerator*dvhatden_dct)
factor3b = (v_hat_denominator*dvhatnum_dsa - v_hat_numerator*dvhatden_dsa)

p3a = dvhatnum_dsadsa/v_hat_numerator
p3b = (dvhatnum_dsa*dvhatden_dsa)/(v_hat_numerator*v_hat_denominator)
p3c = (v_hat_numerator*dvhatden_dsadsa)/(v_hat_numerator*v_hat_denominator)
p3d = dvhatden_dsa/v_hat_denominator

part3 = p3a - 2d0*p3b - p3c + 2d0*p3d*p3d

gsw_cabbeling = part1 - 2d0*(factor2a/factor2b)*part2 + (factor3a/factor3b)*(factor3a/factor3b)*part3

return
end function

!--------------------------------------------------------------------------

