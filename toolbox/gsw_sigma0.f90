!==========================================================================
function gsw_sigma0(sa,ct) 
!==========================================================================

!  Calculates potential density anomaly with reference pressure of 0 dbar,
!  this being this particular potential density minus 1000 kg/m^3.  This
!  function has inputs of Absolute Salinity and Conservative Temperature.
!  This function uses the computationally-efficient 48-term expression for 
!  density in terms of SA, CT and p (IOC et al., 2010).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! 
! gsw_sigma0  : potential density anomaly with reference pressure of 0
!                                                      (48 term equation)

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: v01 =  9.998420897506056d2, v02 =  2.839940833161907d0
real (r14), parameter :: v03 = -3.147759265588511d-2, v04 =  1.181805545074306d-3
real (r14), parameter :: v05 = -6.698001071123802d0, v06 = -2.986498947203215d-2
real (r14), parameter :: v07 =  2.327859407479162d-4, v08 = -3.988822378968490d-2
real (r14), parameter :: v09 =  5.095422573880500d-4, v10 = -1.426984671633621d-5
real (r14), parameter :: v11 =  1.645039373682922d-7
real (r14), parameter :: v21 =  1.0d0, v22 =  2.775927747785646d-3, v23 = -2.349607444135925d-5
real (r14), parameter :: v24 =  1.119513357486743d-6, v25 =  6.743689325042773d-10
real (r14), parameter :: v26 = -7.521448093615448d-3, v27 = -2.764306979894411d-5
real (r14), parameter :: v28 =  1.262937315098546d-7, v29 =  9.527875081696435d-10
real (r14), parameter :: v30 = -1.811147201949891d-11, v31 = -3.303308871386421d-5
real (r14), parameter :: v32 =  3.801564588876298d-7, v33 = -7.672876869259043d-9
real (r14), parameter :: v34 = -4.634182341116144d-11, v35 =  2.681097235569143d-12
real (r14), parameter :: v36 =  5.419326551148740d-6

real (r14) :: sa, ct, gsw_sigma0, v_hat_denominator, v_hat_numerator, sqrtsa

sqrtsa = sqrt(sa)

v_hat_denominator = v01 + ct*(v02 + ct*(v03 + v04*ct))  &
             + sa*(v05 + ct*(v06 + v07*ct) &
         + sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))))

v_hat_numerator = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct))) &
           + sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa &
       + sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))

gsw_sigma0 = v_hat_denominator/v_hat_numerator  - 1000d0

return
end function

!--------------------------------------------------------------------------

