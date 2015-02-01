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

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: sa, p, saturation_fraction

real (r8) :: gsw_ct_freezing_derivative_poly

real (r8) :: x, sa_r, p_r

sa_r = sa*1e-2_r8
x = sqrt(sa_r)
p_r = p*1e-4_r8

gsw_ct_freezing_derivative_poly = 2.0_r8*c1 &
  + x*(3.0_r8*c2 + x*(4.0_r8*c3 + x*(5.0_r8*c4 + x*(6.0_r8*c5 + 7.0_r8*c6*x))))&
  + p_r*(2.0_r8*c10 + p_r*(2.0_r8*c12 + p_r*(2.0_r8*c15 + 4.0_r8*c21*x*x)) &
  + x*x*(4.0_r8*c13 + 4.0_r8*c17*p_r + 6.0_r8*c19*x*x) &
  + x*(3.0_r8*c11 + 3.0_r8*p_r*(c14 + c18*p_r) &
  + x*x*(5.0_r8*c16 + 5.0_r8*c20*p_r + 7.0_r8*c22*x*x)))

gsw_ct_freezing_derivative_poly = 0.5e-2_r8*gsw_ct_freezing_derivative_poly &
    - saturation_fraction*(1e-3_r8)*(-a*(1.0_r8 + b*(1.0_r8 - sa/gsw_sso)) &
    - b*(2.4_r8 - a*sa)/gsw_sso)

return
end function

!--------------------------------------------------------------------------
