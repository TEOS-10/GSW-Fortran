!==========================================================================
module gsw_mod_rho_coefficients
!==========================================================================

use gsw_mod_kinds

implicit none

real (r8), parameter :: v01 =  9.998420897506056e2_r8
real (r8), parameter :: v02 =  2.839940833161907_r8
real (r8), parameter :: v03 = -3.147759265588511e-2_r8
real (r8), parameter :: v04 =  1.181805545074306e-3_r8
real (r8), parameter :: v05 = -6.698001071123802_r8
real (r8), parameter :: v06 = -2.986498947203215e-2_r8
real (r8), parameter :: v07 =  2.327859407479162e-4_r8
real (r8), parameter :: v08 = -3.988822378968490e-2_r8
real (r8), parameter :: v09 =  5.095422573880500e-4_r8
real (r8), parameter :: v10 = -1.426984671633621e-5_r8
real (r8), parameter :: v11 =  1.645039373682922e-7_r8
real (r8), parameter :: v12 = -2.233269627352527e-2_r8
real (r8), parameter :: v13 = -3.436090079851880e-4_r8
real (r8), parameter :: v14 =  3.726050720345733e-6_r8
real (r8), parameter :: v15 = -1.806789763745328e-4_r8
real (r8), parameter :: v16 =  6.876837219536232e-7_r8
real (r8), parameter :: v17 = -3.087032500374211e-7_r8
real (r8), parameter :: v18 = -1.988366587925593e-8_r8
real (r8), parameter :: v19 = -1.061519070296458e-11_r8
real (r8), parameter :: v20 =  1.550932729220080e-10_r8
real (r8), parameter :: v21 =  1.0_r8
real (r8), parameter :: v22 =  2.775927747785646e-3_r8
real (r8), parameter :: v23 = -2.349607444135925e-5_r8
real (r8), parameter :: v24 =  1.119513357486743e-6_r8
real (r8), parameter :: v25 =  6.743689325042773e-10_r8
real (r8), parameter :: v26 = -7.521448093615448e-3_r8
real (r8), parameter :: v27 = -2.764306979894411e-5_r8
real (r8), parameter :: v28 =  1.262937315098546e-7_r8
real (r8), parameter :: v29 =  9.527875081696435e-10_r8
real (r8), parameter :: v30 = -1.811147201949891e-11_r8
real (r8), parameter :: v31 = -3.303308871386421e-5_r8
real (r8), parameter :: v32 =  3.801564588876298e-7_r8
real (r8), parameter :: v33 = -7.672876869259043e-9_r8
real (r8), parameter :: v34 = -4.634182341116144e-11_r8
real (r8), parameter :: v35 =  2.681097235569143e-12_r8
real (r8), parameter :: v36 =  5.419326551148740e-6_r8
real (r8), parameter :: v37 = -2.742185394906099e-5_r8
real (r8), parameter :: v38 = -3.212746477974189e-7_r8
real (r8), parameter :: v39 =  3.191413910561627e-9_r8
real (r8), parameter :: v40 = -1.931012931541776e-12_r8
real (r8), parameter :: v41 = -1.105097577149576e-7_r8
real (r8), parameter :: v42 =  6.211426728363857e-10_r8
real (r8), parameter :: v43 = -1.119011592875110e-10_r8
real (r8), parameter :: v44 = -1.941660213148725e-11_r8
real (r8), parameter :: v45 = -1.864826425365600e-14_r8
real (r8), parameter :: v46 =  1.119522344879478e-14_r8
real (r8), parameter :: v47 = -1.200507748551599e-15_r8
real (r8), parameter :: v48 =  6.057902487546866e-17_r8

real (r8), parameter :: a01 =  2.839940833161907_r8
real (r8), parameter :: a02 = -6.295518531177023e-2_r8
real (r8), parameter :: a03 =  3.545416635222918e-3_r8
real (r8), parameter :: a04 = -2.986498947203215e-2_r8
real (r8), parameter :: a05 =  4.655718814958324e-4_r8
real (r8), parameter :: a06 =  5.095422573880500e-4_r8
real (r8), parameter :: a07 = -2.853969343267241e-5_r8
real (r8), parameter :: a08 =  4.935118121048767e-7_r8
real (r8), parameter :: a09 = -3.436090079851880e-4_r8
real (r8), parameter :: a10 =  7.452101440691467e-6_r8
real (r8), parameter :: a11 =  6.876837219536232e-7_r8
real (r8), parameter :: a12 = -1.988366587925593e-8_r8
real (r8), parameter :: a13 = -2.123038140592916e-11_r8
real (r8), parameter :: a14 =  2.775927747785646e-3_r8
real (r8), parameter :: a15 = -4.699214888271850e-5_r8
real (r8), parameter :: a16 =  3.358540072460230e-6_r8
real (r8), parameter :: a17 =  2.697475730017109e-9_r8
real (r8), parameter :: a18 = -2.764306979894411e-5_r8
real (r8), parameter :: a19 =  2.525874630197091e-7_r8
real (r8), parameter :: a20 =  2.858362524508931e-9_r8
real (r8), parameter :: a21 = -7.244588807799565e-11_r8
real (r8), parameter :: a22 =  3.801564588876298e-7_r8
real (r8), parameter :: a23 = -1.534575373851809e-8_r8
real (r8), parameter :: a24 = -1.390254702334843e-10_r8
real (r8), parameter :: a25 =  1.072438894227657e-11_r8
real (r8), parameter :: a26 = -3.212746477974189e-7_r8
real (r8), parameter :: a27 =  6.382827821123254e-9_r8
real (r8), parameter :: a28 = -5.793038794625329e-12_r8
real (r8), parameter :: a29 =  6.211426728363857e-10_r8
real (r8), parameter :: a30 = -1.941660213148725e-11_r8
real (r8), parameter :: a31 = -3.729652850731201e-14_r8
real (r8), parameter :: a32 =  1.119522344879478e-14_r8
real (r8), parameter :: a33 =  6.057902487546866e-17_r8

real (r8), parameter :: b01 = -6.698001071123802_r8
real (r8), parameter :: b02 = -2.986498947203215e-2_r8
real (r8), parameter :: b03 =  2.327859407479162e-4_r8
real (r8), parameter :: b04 = -5.983233568452735e-2_r8
real (r8), parameter :: b05 =  7.643133860820750e-4_r8
real (r8), parameter :: b06 = -2.140477007450431e-5_r8
real (r8), parameter :: b07 =  2.467559060524383e-7_r8
real (r8), parameter :: b08 = -1.806789763745328e-4_r8
real (r8), parameter :: b09 =  6.876837219536232e-7_r8
real (r8), parameter :: b10 =  1.550932729220080e-10_r8
real (r8), parameter :: b11 = -7.521448093615448e-3_r8
real (r8), parameter :: b12 = -2.764306979894411e-5_r8
real (r8), parameter :: b13 =  1.262937315098546e-7_r8
real (r8), parameter :: b14 =  9.527875081696435e-10_r8
real (r8), parameter :: b15 = -1.811147201949891e-11_r8
real (r8), parameter :: b16 = -4.954963307079632e-5_r8
real (r8), parameter :: b17 =  5.702346883314446e-7_r8
real (r8), parameter :: b18 = -1.150931530388857e-8_r8
real (r8), parameter :: b19 = -6.951273511674217e-11_r8
real (r8), parameter :: b20 =  4.021645853353715e-12_r8
real (r8), parameter :: b21 =  1.083865310229748e-5_r8
real (r8), parameter :: b22 = -1.105097577149576e-7_r8
real (r8), parameter :: b23 =  6.211426728363857e-10_r8
real (r8), parameter :: b24 =  1.119522344879478e-14_r8

real (r8), parameter :: c01 = -2.233269627352527e-2_r8
real (r8), parameter :: c02 = -3.436090079851880e-4_r8
real (r8), parameter :: c03 =  3.726050720345733e-6_r8
real (r8), parameter :: c04 = -1.806789763745328e-4_r8
real (r8), parameter :: c05 =  6.876837219536232e-7_r8
real (r8), parameter :: c06 = -6.174065000748422e-7_r8
real (r8), parameter :: c07 = -3.976733175851186e-8_r8
real (r8), parameter :: c08 = -2.123038140592916e-11_r8
real (r8), parameter :: c09 =  3.101865458440160e-10_r8
real (r8), parameter :: c10 = -2.742185394906099e-5_r8
real (r8), parameter :: c11 = -3.212746477974189e-7_r8
real (r8), parameter :: c12 =  3.191413910561627e-9_r8
real (r8), parameter :: c13 = -1.931012931541776e-12_r8
real (r8), parameter :: c14 = -1.105097577149576e-7_r8
real (r8), parameter :: c15 =  6.211426728363857e-10_r8
real (r8), parameter :: c16 = -2.238023185750219e-10_r8
real (r8), parameter :: c17 = -3.883320426297450e-11_r8
real (r8), parameter :: c18 = -3.729652850731201e-14_r8
real (r8), parameter :: c19 =  2.239044689758956e-14_r8
real (r8), parameter :: c20 = -3.601523245654798e-15_r8
real (r8), parameter :: c21 =  1.817370746264060e-16_r8

end module

!--------------------------------------------------------------------------
