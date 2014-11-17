!==========================================================================
function gsw_sound_speed(sa,ct,p) 
!==========================================================================

!  Calculates sound speed of seawater using the computationally
!  efficient 48-term expression for density in terms of SA, CT and p
!  (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_sound_speed  :  sound speed of seawater (48 term equation)

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
real (r14), parameter :: c21 =  1.817370746264060d-16

real (r14) :: sa, ct, p, gsw_sound_speed, v_hat_denominator, v_hat_numerator
real (r14) :: sqrtsa, dvhatden_dp, dvhatnum_dp, dp_drho

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

dvhatden_dp = c01 + ct*(c02 + c03*ct) &
    + sa*(c04 + c05*ct) &
    + p*(c06 + ct*(c07 + c08*ct) + c09*sa)

dvhatnum_dp = c10 + ct*(c11 + ct*(c12 + c13*ct)) &
    + sa*(c14 + c15*ct) &
    + p*(c16 + ct*(c17 + c18*ct + c19*sa) &
    + p*(c20 + c21*ct))

dp_drho = (v_hat_numerator*v_hat_numerator)/(dvhatden_dp*v_hat_numerator - dvhatnum_dp*v_hat_denominator)
    
gsw_sound_speed = 100d0*sqrt(dp_drho)

return
end function

!--------------------------------------------------------------------------

