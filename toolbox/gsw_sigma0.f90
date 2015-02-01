!==========================================================================
elemental function gsw_sigma0 (sa, ct) 
!==========================================================================
!
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
!--------------------------------------------------------------------------

use gsw_mod_rho_coefficients

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: sa, ct 

real (r8) :: gsw_sigma0

real (r8) :: v_hat_denominator, v_hat_numerator, sqrtsa

sqrtsa = sqrt(sa)

v_hat_denominator = v01 + ct*(v02 + ct*(v03 + v04*ct))  &
              + sa*(v05 + ct*(v06 + v07*ct) &
          + sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))))

v_hat_numerator = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct))) &
            + sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa &
        + sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))

gsw_sigma0 = v_hat_denominator/v_hat_numerator - 1e3_r8

return
end function

!--------------------------------------------------------------------------
