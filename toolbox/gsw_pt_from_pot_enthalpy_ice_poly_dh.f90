!==========================================================================
elemental function gsw_pt_from_pot_enthalpy_ice_poly_dh (pot_enthalpy_ice)
!==========================================================================
!
!  Calculates the derivative of potential temperature of ice with respect 
!  to potential enthalpy.  This is based on the compuationally-efficient 
!  polynomial fit to the potential enthalpy of ice. 
!
!  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
!
!  dpt0_ice_dh  =  derivative of potential temperature of ice 
!                  with respect to potential enthalpy             [ deg C ]
!--------------------------------------------------------------------------

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: pot_enthalpy_ice

real (r14) :: gsw_pt_from_pot_enthalpy_ice_poly_dh

real (r14), parameter :: q1 = 2.594351081876611d-3
real (r14), parameter :: r2 = 3.530155620427630d-8
real (r14), parameter :: r3 = 2.330421169287162d-13
real (r14), parameter :: r4 = 8.139369017110120d-19
real (r14), parameter :: r5 = 1.610007265856420d-24
real (r14), parameter :: r6 = 1.707103685781641d-30
real (r14), parameter :: r7 = 7.658041152250651d-37

gsw_pt_from_pot_enthalpy_ice_poly_dh = q1 &
    + pot_enthalpy_ice*(r2 + pot_enthalpy_ice*(r3 &
    + pot_enthalpy_ice*(r4 + pot_enthalpy_ice*(r5 + pot_enthalpy_ice*(r6 &
    + pot_enthalpy_ice*r7)))))

return
end function

!--------------------------------------------------------------------------
