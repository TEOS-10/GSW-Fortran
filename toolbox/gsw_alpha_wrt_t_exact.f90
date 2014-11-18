!==========================================================================
elemental function gsw_alpha_wrt_t_exact (sa, t, p) 
!==========================================================================
!
! Calculates thermal expansion coefficient of seawater with respect to 
! in-situ temperature
!
! sa     : Absolute Salinity                               [g/kg]
! t      : insitu temperature                              [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_alpha_wrt_t_exact : thermal expansion coefficient    [1/K]
!                         wrt (in-situ) temperature
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_gibbs

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, t, p 

real (r14) :: gsw_alpha_wrt_t_exact

integer, parameter :: n0=0, n1=1

gsw_alpha_wrt_t_exact = gsw_gibbs(n0,n1,n1,sa,t,p) / gsw_gibbs(n0,n0,n1,sa,t,p)

return
end function
