!==========================================================================
elemental function gsw_ct_freezing_derivative_poly (sa, p, &
                                                    saturation_fraction)
!==========================================================================
!
!  Calculates the first derivative of the Conservative Temperature at
!  which seawater freezes, with respect to Absolute Salinity SA.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!  saturation_fraction = the saturation fraction of dissolved air in 
!                        seawater
!
!  CTfreezing_SA = the derivative of the Conservative Temperature at
!                  freezing (ITS-90) with respect to Absolute Salinity at
!                  fixed pressure              [ K/(g/kg) ] i.e. [ K kg/g ]
!--------------------------------------------------------------------------

use gsw_mod_teos10_constants, only : gsw_sso

use gsw_mod_freezing_poly_coefficients

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, p, saturation_fraction

real (r14) :: gsw_ct_freezing_derivative_poly

real (r14) :: x, sa_r, p_r

sa_r = sa*1d-2
x = sqrt(sa_r)
p_r = p*1d-4

gsw_ct_freezing_derivative_poly = 2*c1 &
    + x*(3*c2 + x*(4*c3 + x*(5*c4 + x*(6*c5 + 7*c6*x)))) &
    + p_r*(2*c10 + p_r*(2*c12 + p_r*(2*c15 + 4*c21*x*x)) &
    + x*x*(4*c13 + 4*c17*p_r + 6*c19*x*x) &
    + x*(3*c11 + 3*p_r*(c14 + c18*p_r) &
    + x*x*(5*c16 + 5*c20*p_r + 7*c22*x*x)))

gsw_ct_freezing_derivative_poly = 0.5d0*1d-2*gsw_ct_freezing_derivative_poly &
    - saturation_fraction*(1d-3)*(-a*(1d0 + b*(1d0 - sa/gsw_sso)) &
    - b*(2.4d0 - a*sa)/gsw_sso)

return
end function

!--------------------------------------------------------------------------
