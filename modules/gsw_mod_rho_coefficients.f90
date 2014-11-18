!==========================================================================
module gsw_mod_rho_coefficients
!==========================================================================

implicit none

integer, parameter :: grc_r14 = selected_real_kind(14,30)

real (grc_r14), parameter :: v01 =  9.998420897506056d2
real (grc_r14), parameter :: v02 =  2.839940833161907d0
real (grc_r14), parameter :: v03 = -3.147759265588511d-2
real (grc_r14), parameter :: v04 =  1.181805545074306d-3
real (grc_r14), parameter :: v05 = -6.698001071123802d0
real (grc_r14), parameter :: v06 = -2.986498947203215d-2
real (grc_r14), parameter :: v07 =  2.327859407479162d-4
real (grc_r14), parameter :: v08 = -3.988822378968490d-2
real (grc_r14), parameter :: v09 =  5.095422573880500d-4
real (grc_r14), parameter :: v10 = -1.426984671633621d-5
real (grc_r14), parameter :: v11 =  1.645039373682922d-7
real (grc_r14), parameter :: v12 = -2.233269627352527d-2
real (grc_r14), parameter :: v13 = -3.436090079851880d-4
real (grc_r14), parameter :: v14 =  3.726050720345733d-6
real (grc_r14), parameter :: v15 = -1.806789763745328d-4
real (grc_r14), parameter :: v16 =  6.876837219536232d-7
real (grc_r14), parameter :: v17 = -3.087032500374211d-7
real (grc_r14), parameter :: v18 = -1.988366587925593d-8
real (grc_r14), parameter :: v19 = -1.061519070296458d-11
real (grc_r14), parameter :: v20 =  1.550932729220080d-10
real (grc_r14), parameter :: v21 =  1.0d0
real (grc_r14), parameter :: v22 =  2.775927747785646d-3
real (grc_r14), parameter :: v23 = -2.349607444135925d-5
real (grc_r14), parameter :: v24 =  1.119513357486743d-6
real (grc_r14), parameter :: v25 =  6.743689325042773d-10
real (grc_r14), parameter :: v26 = -7.521448093615448d-3
real (grc_r14), parameter :: v27 = -2.764306979894411d-5
real (grc_r14), parameter :: v28 =  1.262937315098546d-7
real (grc_r14), parameter :: v29 =  9.527875081696435d-10
real (grc_r14), parameter :: v30 = -1.811147201949891d-11
real (grc_r14), parameter :: v31 = -3.303308871386421d-5
real (grc_r14), parameter :: v32 =  3.801564588876298d-7
real (grc_r14), parameter :: v33 = -7.672876869259043d-9
real (grc_r14), parameter :: v34 = -4.634182341116144d-11
real (grc_r14), parameter :: v35 =  2.681097235569143d-12
real (grc_r14), parameter :: v36 =  5.419326551148740d-6
real (grc_r14), parameter :: v37 = -2.742185394906099d-5
real (grc_r14), parameter :: v38 = -3.212746477974189d-7
real (grc_r14), parameter :: v39 =  3.191413910561627d-9
real (grc_r14), parameter :: v40 = -1.931012931541776d-12
real (grc_r14), parameter :: v41 = -1.105097577149576d-7
real (grc_r14), parameter :: v42 =  6.211426728363857d-10
real (grc_r14), parameter :: v43 = -1.119011592875110d-10
real (grc_r14), parameter :: v44 = -1.941660213148725d-11
real (grc_r14), parameter :: v45 = -1.864826425365600d-14
real (grc_r14), parameter :: v46 =  1.119522344879478d-14
real (grc_r14), parameter :: v47 = -1.200507748551599d-15
real (grc_r14), parameter :: v48 =  6.057902487546866d-17 

real (grc_r14), parameter :: a01 =  2.839940833161907d0
real (grc_r14), parameter :: a02 = -6.295518531177023d-2
real (grc_r14), parameter :: a03 =  3.545416635222918d-3
real (grc_r14), parameter :: a04 = -2.986498947203215d-2
real (grc_r14), parameter :: a05 =  4.655718814958324d-4
real (grc_r14), parameter :: a06 =  5.095422573880500d-4
real (grc_r14), parameter :: a07 = -2.853969343267241d-5
real (grc_r14), parameter :: a08 =  4.935118121048767d-7
real (grc_r14), parameter :: a09 = -3.436090079851880d-4
real (grc_r14), parameter :: a10 =  7.452101440691467d-6
real (grc_r14), parameter :: a11 =  6.876837219536232d-7
real (grc_r14), parameter :: a12 = -1.988366587925593d-8
real (grc_r14), parameter :: a13 = -2.123038140592916d-11
real (grc_r14), parameter :: a14 =  2.775927747785646d-3
real (grc_r14), parameter :: a15 = -4.699214888271850d-5
real (grc_r14), parameter :: a16 =  3.358540072460230d-6
real (grc_r14), parameter :: a17 =  2.697475730017109d-9
real (grc_r14), parameter :: a18 = -2.764306979894411d-5
real (grc_r14), parameter :: a19 =  2.525874630197091d-7
real (grc_r14), parameter :: a20 =  2.858362524508931d-9
real (grc_r14), parameter :: a21 = -7.244588807799565d-11
real (grc_r14), parameter :: a22 =  3.801564588876298d-7
real (grc_r14), parameter :: a23 = -1.534575373851809d-8
real (grc_r14), parameter :: a24 = -1.390254702334843d-10
real (grc_r14), parameter :: a25 =  1.072438894227657d-11
real (grc_r14), parameter :: a26 = -3.212746477974189d-7
real (grc_r14), parameter :: a27 =  6.382827821123254d-9
real (grc_r14), parameter :: a28 = -5.793038794625329d-12
real (grc_r14), parameter :: a29 =  6.211426728363857d-10
real (grc_r14), parameter :: a30 = -1.941660213148725d-11
real (grc_r14), parameter :: a31 = -3.729652850731201d-14
real (grc_r14), parameter :: a32 =  1.119522344879478d-14
real (grc_r14), parameter :: a33 =  6.057902487546866d-17

real (grc_r14), parameter :: b01 = -6.698001071123802d0
real (grc_r14), parameter :: b02 = -2.986498947203215d-2
real (grc_r14), parameter :: b03 =  2.327859407479162d-4
real (grc_r14), parameter :: b04 = -5.983233568452735d-2
real (grc_r14), parameter :: b05 =  7.643133860820750d-4
real (grc_r14), parameter :: b06 = -2.140477007450431d-5
real (grc_r14), parameter :: b07 =  2.467559060524383d-7
real (grc_r14), parameter :: b08 = -1.806789763745328d-4
real (grc_r14), parameter :: b09 =  6.876837219536232d-7
real (grc_r14), parameter :: b10 =  1.550932729220080d-10
real (grc_r14), parameter :: b11 = -7.521448093615448d-3
real (grc_r14), parameter :: b12 = -2.764306979894411d-5
real (grc_r14), parameter :: b13 =  1.262937315098546d-7
real (grc_r14), parameter :: b14 =  9.527875081696435d-10
real (grc_r14), parameter :: b15 = -1.811147201949891d-11
real (grc_r14), parameter :: b16 = -4.954963307079632d-5
real (grc_r14), parameter :: b17 =  5.702346883314446d-7
real (grc_r14), parameter :: b18 = -1.150931530388857d-8
real (grc_r14), parameter :: b19 = -6.951273511674217d-11
real (grc_r14), parameter :: b20 =  4.021645853353715d-12
real (grc_r14), parameter :: b21 =  1.083865310229748d-5
real (grc_r14), parameter :: b22 = -1.105097577149576d-7
real (grc_r14), parameter :: b23 =  6.211426728363857d-10
real (grc_r14), parameter :: b24 =  1.119522344879478d-14

real (grc_r14), parameter :: c01 = -2.233269627352527d-2
real (grc_r14), parameter :: c02 = -3.436090079851880d-4
real (grc_r14), parameter :: c03 =  3.726050720345733d-6
real (grc_r14), parameter :: c04 = -1.806789763745328d-4
real (grc_r14), parameter :: c05 =  6.876837219536232d-7
real (grc_r14), parameter :: c06 = -6.174065000748422d-7
real (grc_r14), parameter :: c07 = -3.976733175851186d-8
real (grc_r14), parameter :: c08 = -2.123038140592916d-11
real (grc_r14), parameter :: c09 =  3.101865458440160d-10
real (grc_r14), parameter :: c10 = -2.742185394906099d-5
real (grc_r14), parameter :: c11 = -3.212746477974189d-7
real (grc_r14), parameter :: c12 =  3.191413910561627d-9
real (grc_r14), parameter :: c13 = -1.931012931541776d-12
real (grc_r14), parameter :: c14 = -1.105097577149576d-7
real (grc_r14), parameter :: c15 =  6.211426728363857d-10
real (grc_r14), parameter :: c16 = -2.238023185750219d-10
real (grc_r14), parameter :: c17 = -3.883320426297450d-11
real (grc_r14), parameter :: c18 = -3.729652850731201d-14
real (grc_r14), parameter :: c19 =  2.239044689758956d-14
real (grc_r14), parameter :: c20 = -3.601523245654798d-15
real (grc_r14), parameter :: c21 =  1.817370746264060d-16

end module

!--------------------------------------------------------------------------
