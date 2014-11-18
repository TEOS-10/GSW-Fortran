!==========================================================================
elemental function gsw_pt0_cold_ice_poly (pot_enthalpy_ice)
!==========================================================================
!
!  Calculates an initial estimate of pt0_ice when it is less than about
!  -100 deg C. 
!
!  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
!
!  pt0_cold_ice_poly  =  initial estimate of potential temperatur 
!                        of very cold ice in dgress C (not K)     [ deg C ] 
!--------------------------------------------------------------------------

use gsw_mod_teos10_constants, only : gsw_t0

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: pot_enthalpy_ice

real (r14) :: gsw_pt0_cold_ice_poly

real (r14) :: log_abs_theta0, log_h_diff

! h00 = gsw_enthalpy_ice(-gsw_t0,0)
real (r14), parameter :: h00 = -6.320202333358860d5

real (r14), parameter :: s0 =  1.493103204647916d0
real (r14), parameter :: s1 =  2.372788609320607d-1
real (r14), parameter :: s2 = -2.014996002119374d-3
real (r14), parameter :: s3 =  2.640600197732682d-6
real (r14), parameter :: s4 =  3.134706016844293d-5
real (r14), parameter :: s5 =  2.733592344937913d-6
real (r14), parameter :: s6 =  4.726828010223258d-8
real (r14), parameter :: s7 = -2.735193883189589d-9
real (r14), parameter :: s8 = -8.547714991377670d-11

log_h_diff = log(pot_enthalpy_ice - h00)

log_abs_theta0 = s0 + log_h_diff*(s1 + log_h_diff*(s2 + log_h_diff*(s3 &
                + log_h_diff*(s4 + log_h_diff*(s5 + log_h_diff*(s6 &
                + log_h_diff*(s7 + log_h_diff*s8)))))))

gsw_pt0_cold_ice_poly = exp(log_abs_theta0) - gsw_t0
    
return
end function

!--------------------------------------------------------------------------
