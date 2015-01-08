!==========================================================================
elemental subroutine gsw_rho_first_derivatives (sa, ct, p, drho_dsa, &
                                                drho_dct, drho_dp)
!==========================================================================
!
!  Calculates the three (3) partial derivatives of in situ density with 
!  respect to Absolute Salinity, Conservative Temperature and pressure.  
!  Note that the pressure derivative is done with respect to pressure in 
!  Pa, not dbar.  This function uses the computationally-efficient 48-term 
!  expression for density in terms of SA, CT and p.
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
!
! drho_dsa  : partial derivatives of density             [ kg^2/(g m^3) ]
!              with respect to Absolute Salinity
! drho_dct  : partial derivatives of density               [ kg/(K m^3) ]
!              with respect to Conservative Temperature
! drho_dp   : partial derivatives of density              [ kg/(Pa m^3) ]
!              with respect to pressure in Pa
!--------------------------------------------------------------------------

use gsw_mod_teos10_constants, only : pa2db

use gsw_mod_rho_coefficients

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, ct, p
real (r14), intent(out), optional :: drho_dsa, drho_dct, drho_dp

real (r14) :: sqrtsa, v_hat_denominator, v_hat_numerator
real (r14) :: dvhatden_dct, dvhatnum_dct, dvhatden_dsa, dvhatnum_dsa
real (r14) :: dvhatden_dp, dvhatnum_dp, rho, rec_num

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
       
rec_num = 1d0/v_hat_numerator
       
rho = rec_num*v_hat_denominator

if (present(drho_dsa)) then

    dvhatden_dsa = b01 + ct*(b02 + b03*ct) &
         + sqrtsa*(b04 + ct*(b05 + ct*(b06 + b07*ct))) &
              + p*(b08 + b09*ct + b10*p) 

    dvhatnum_dsa = b11 + ct*(b12 + ct*(b13 + ct*(b14 + b15*ct))) &
         + sqrtsa*(b16 + ct*(b17 + ct*(b18 + ct*(b19 + b20*ct)))) + b21*sa &
              + p*(b22 + ct*(b23 + b24*p))
		
    drho_dsa = (dvhatden_dsa - dvhatnum_dsa*rho)*rec_num

end if

if (present(drho_dct)) then

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

    drho_dct = (dvhatden_dct - dvhatnum_dct*rho)*rec_num

end if

if (present(drho_dp)) then

    dvhatden_dp = c01 + ct*(c02 + c03*ct) &
        + sa*(c04 + c05*ct) &
        + p*(c06 + ct*(c07 + c08*ct) + c09*sa)

    dvhatnum_dp = c10 + ct*(c11 + ct*(c12 + c13*ct)) &
        + sa*(c14 + c15*ct) &
        + p*(c16 + ct*(c17 + c18*ct + c19*sa) &
        + p*(c20 + c21*ct))

    drho_dp = pa2db*(dvhatden_dp - dvhatnum_dp*rho)*rec_num

end if

return
end subroutine

!--------------------------------------------------------------------------
