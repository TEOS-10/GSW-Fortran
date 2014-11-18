!==========================================================================
elemental function gsw_pot_rho_t_exact (sa, t, p, p_ref)  
!==========================================================================
!
! Calculates the potential density of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! p_ref  : reference sea pressure                          [dbar]
! 
! gsw_pot_rho_t_exact : potential density                  [kg/m^3]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_pt_from_t, gsw_rho_t_exact

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, t, p, p_ref  

real (r14) :: gsw_pot_rho_t_exact

real (r14) :: pt

pt = gsw_pt_from_t(sa,t,p,p_ref)

gsw_pot_rho_t_exact = gsw_rho_t_exact(sa,pt,p_ref)

return
end function

!--------------------------------------------------------------------------
