!==========================================================================
elemental function gsw_pot_enthalpy_from_pt_ice_poly (pt0_ice)
!==========================================================================
!
!  Calculates the potential enthalpy of ice from potential temperature of
!  ice (whose reference sea pressure is zero dbar).  This is a
!  compuationally efficient polynomial fit to the potential enthalpy of
!  ice.
!   
!  pt0_ice  =  potential temperature of ice (ITS-90)              [ deg C ]
!
!  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_pt_from_pot_enthalpy_ice_poly
use gsw_mod_toolbox, only : gsw_pt_from_pot_enthalpy_ice_poly_dh

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: pt0_ice

real (r14) :: gsw_pot_enthalpy_from_pt_ice_poly

integer :: iteration
real (r14) :: df_dt, f, pot_enthalpy_ice
real (r14) :: pot_enthalpy_ice_mid, pot_enthalpy_ice_old

real (r14), parameter :: p0 = -3.333601570157700d5
real (r14), parameter :: p1 =  2.096693916810367d3
real (r14), parameter :: p2 =  3.687110754043292d0
real (r14), parameter :: p3 =  4.559401565980682d-4
real (r14), parameter :: p4 = -2.516011957758120d-6
real (r14), parameter :: p5 = -1.040364574632784d-8
real (r14), parameter :: p6 = -1.701786588412454d-10
real (r14), parameter :: p7 = -7.667191301635057d-13
    
! initial estimate of the potential enthalpy.
pot_enthalpy_ice = p0 + pt0_ice*(p1 + pt0_ice*(p2 + pt0_ice*(p3 &
                   + pt0_ice*(p4 + pt0_ice*(p5 + pt0_ice*(p6 &
                   + pt0_ice*p7))))))
               
df_dt = gsw_pt_from_pot_enthalpy_ice_poly_dh(pot_enthalpy_ice)

do iteration = 1, 5
    pot_enthalpy_ice_old = pot_enthalpy_ice
    f = gsw_pt_from_pot_enthalpy_ice_poly(pot_enthalpy_ice_old) - pt0_ice
    pot_enthalpy_ice = pot_enthalpy_ice_old - f/df_dt
    pot_enthalpy_ice_mid = 0.5d0*(pot_enthalpy_ice + pot_enthalpy_ice_old)
    df_dt = gsw_pt_from_pot_enthalpy_ice_poly_dh(pot_enthalpy_ice_mid)
    pot_enthalpy_ice = pot_enthalpy_ice_old - f/df_dt
end do
               
! The error of this fit ranges between -6e-3 and 6e-3 J/kg over the potential 
! temperature range of -100 to 2 deg C, or the potential enthalpy range of 
! -5.7 x 10^5 to -3.3 x 10^5 J/kg. 
       
gsw_pot_enthalpy_from_pt_ice_poly = pot_enthalpy_ice

return
end function

!--------------------------------------------------------------------------
