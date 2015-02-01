!==========================================================================
elemental function gsw_t_freezing_derivative_poly (sa, p, saturation_fraction)
!==========================================================================
!
!  Calculates the first derivative of the in-situ temperature at which 
!  seawater freezes with respect to Absolute Salinity SA.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!  saturation_fraction = the saturation fraction of dissolved air in 
!                        seawater
!
!  tfreezing_SA = the derivative of the in-situ freezing temperature 
!                 (ITS-90) with respect to Absolute Salinity at fixed    
!                 pressure                     [ K/(g/kg) ] i.e. [ K kg/g ] 
!--------------------------------------------------------------------------

use gsw_mod_teos10_constants, only : gsw_sso

use gsw_mod_freezing_poly_coefficients

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: sa, p, saturation_fraction

real (r8) :: gsw_t_freezing_derivative_poly

real (r8) :: x, sa_r, p_r

sa_r = sa*1e-2_r8
x = sqrt(sa_r)
p_r = p*1e-4_r8

gsw_t_freezing_derivative_poly = 2.0_r8*t1 + x* &
     (3.0_r8*t2 + x*(4.0_r8*t3 + x*(5.0_r8*t4 + x*(6.0_r8*t5 + 7.0_r8*t6*x)))) &
     + p_r*(2.0_r8*t10 + p_r*(2.0_r8*t12 + p_r*(2.0_r8*t15 + 4.0_r8*t21*sa_r)) &
     + sa_r*(4.0_r8*t13 + 4.0_r8*t17*p_r + 6.0_r8*t19*sa_r) &
     + x*(3.0_r8*t11 + 3.0_r8*p_r*(t14 + t18*p_r) &
     + sa_r*(5.0_r8*t16 + 5.0_r8*t20*p_r + 7.0_r8*t22*sa_r)))

gsw_t_freezing_derivative_poly = 0.5e-2_r8*gsw_t_freezing_derivative_poly + &
                                 saturation_fraction*(1e-3_r8)/(2.0_r8*gsw_sso)

return
end function

!--------------------------------------------------------------------------
