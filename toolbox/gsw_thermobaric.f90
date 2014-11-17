!==========================================================================
function gsw_thermobaric(sa,ct,p)  
!==========================================================================

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
real (r14), parameter :: c01 = -2.233269627352527d-2,  c02 = -3.436090079851880d-4
real (r14), parameter :: c03 =  3.726050720345733d-6,  c04 = -1.806789763745328d-4
real (r14), parameter :: c05 =  6.876837219536232d-7,  c06 = -6.174065000748422d-7
real (r14), parameter :: c07 = -3.976733175851186d-8,  c08 = -2.123038140592916d-11
real (r14), parameter :: c09 =  3.101865458440160d-10, c10 = -2.742185394906099d-5
real (r14), parameter :: c11 = -3.212746477974189d-7,  c12 =  3.191413910561627d-9
real (r14), parameter :: c13 = -1.931012931541776d-12, c14 = -1.105097577149576d-7
real (r14), parameter :: c15 =  6.211426728363857d-10, c16 = -2.238023185750219d-10
real (r14), parameter :: c17 = -3.883320426297450d-11, c18 = -3.729652850731201d-14
real (r14), parameter :: c19 =  2.239044689758956d-14, c20 = -3.601523245654798d-15
real (r14), parameter :: c21 =  1.817370746264060d-16, db2pa = 1d-4

real (r14) :: sa, ct, p, gsw_thermobaric, sqrtsa, v_hat_denominator, v_hat_numerator
real (r14) :: dvhatden_dct, dvhatnum_dct, dvhatden_dp, dvhatnum_dp
real (r14) :: dvhatden_dsa, dvhatnum_dsa, dvhatden_dpdct, dvhatnum_dpdct
real (r14) :: dvhatden_dpdsa, dvhatnum_dpdsa 
real (r14) :: p1a, p1b, p1c, p1d, p1e, part1, factor2
real (r14) :: p2a, p2b, p2c, p2d, p2e, part2  

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

dvhatden_dpdct = c02 + 2d0*c03*ct + c05*sa &
           + p*(c07 + 2d0*c08*ct)
       
dvhatnum_dpdct = c11 + ct*(2d0*c12 + 3d0*c13*ct) + c15*sa &
           + p*(c17 + ct*2d0*c18 + c19*sa + c21*p)

dvhatden_dpdsa = c04 + c05*ct + c09*p

dvhatnum_dpdsa = c14 + c15*ct + c19*ct*p

p1a = dvhatnum_dpdct/v_hat_numerator
p1b = (dvhatnum_dct*dvhatden_dp)/(v_hat_numerator*v_hat_denominator)
p1c = (dvhatnum_dp*dvhatden_dct)/(v_hat_numerator*v_hat_denominator)
p1d = (dvhatden_dp*dvhatden_dct)/(v_hat_denominator*v_hat_denominator)
p1e = dvhatden_dpdct/v_hat_denominator

part1 =  p1a - p1b - p1c + 2d0*p1d - p1e

factor2 = (v_hat_denominator*dvhatnum_dct - v_hat_numerator*dvhatden_dct)/ &
           (v_hat_denominator*dvhatnum_dsa - v_hat_numerator*dvhatden_dsa)

p2a = dvhatnum_dpdsa/v_hat_numerator
p2b = (dvhatnum_dsa*dvhatden_dp)/(v_hat_numerator*v_hat_denominator)
p2c = (dvhatnum_dp*dvhatden_dsa)/(v_hat_numerator*v_hat_denominator)
p2d = (dvhatden_dp*dvhatden_dsa)/(v_hat_denominator*v_hat_denominator)
p2e = dvhatden_dpdsa/v_hat_denominator

part2 =  p2a - p2b - p2c + 2d0*p2d - p2e

gsw_thermobaric = (part1 - factor2*part2)*db2pa

return
end function

!--------------------------------------------------------------------------

