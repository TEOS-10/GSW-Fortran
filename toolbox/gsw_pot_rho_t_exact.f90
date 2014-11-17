!==========================================================================
function gsw_pot_rho_t_exact(sa,t,p,p_ref)  
!==========================================================================

! Calculates the potential density of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! p_ref  : reference sea pressure                          [dbar]
! 
! gsw_pot_rho_t_exact : potential density                  [kg/m^3]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, t, p, p_ref, gsw_pot_rho_t_exact, pt
real (r14) :: gsw_pt_from_t, gsw_rho_t_exact

pt = gsw_pt_from_t(sa,t,p,p_ref)

gsw_pot_rho_t_exact = gsw_rho_t_exact(sa,pt,p_ref)

return
end

!--------------------------------------------------------------------------

