!==========================================================================
elemental function gsw_rho_t_exact (sa, t, p) 
!==========================================================================
!
! Calculates in-situ density of seawater from Absolute Salinity and 
! in-situ temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_rho_t_exact : in-situ density                        [kg/m^3]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_gibbs

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, t, p 

real (r14) :: gsw_rho_t_exact

integer, parameter :: n0=0, n1=1

gsw_rho_t_exact = 1d0/gsw_gibbs(n0,n0,n1,sa,t,p)

return
end function

!--------------------------------------------------------------------------
