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

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, p, saturation_fraction

real (r14) :: gsw_t_freezing_derivative_poly

real (r14) :: x, sa_r, p_r

sa_r = sa*1d-2
x = sqrt(sa_r)
p_r = p*1d-4

gsw_t_freezing_derivative_poly = 2d0*t1 &
     + x*(3d0*t2 + x*(4d0*t3 + x*(5d0*t4 + x*(6d0*t5 + 7d0*t6*x)))) &
     + p_r*(2d0*t10 + p_r*(2d0*t12 + p_r*(2d0*t15 + 4d0*t21*sa_r)) &
     + sa_r*(4d0*t13 + 4d0*t17*p_r + 6d0*t19*sa_r) &
     + x*(3d0*t11 + 3d0*p_r*(t14 + t18*p_r) &
     + sa_r*(5d0*t16 + 5d0*t20*p_r + 7d0*t22*sa_r)))

gsw_t_freezing_derivative_poly = 0.5d0*1d-2*gsw_t_freezing_derivative_poly &
                                 + saturation_fraction*(1d-3)/(2d0*gsw_sso)

return
end function

!--------------------------------------------------------------------------
