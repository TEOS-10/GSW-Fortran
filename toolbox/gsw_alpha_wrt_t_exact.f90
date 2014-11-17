!==========================================================================
function gsw_alpha_wrt_t_exact(sa,t,p) 
!==========================================================================

! Calculates thermal expansion coefficient of seawater with respect to 
! in-situ temperature
!
! sa     : Absolute Salinity                               [g/kg]
! t      : insitu temperature                              [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_alpha_wrt_t_exact : thermal expansion coefficient    [1/K]
!                         wrt (in-situ) temperature

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n0, n1
real (r14) :: sa, t, p, gsw_gibbs, gsw_alpha_wrt_t_exact

n0 = 0
n1 = 1

gsw_alpha_wrt_t_exact = gsw_gibbs(n0,n1,n1,sa,t,p)/gsw_gibbs(n0,n0,n1,sa,t,p)

end

!--------------------------------------------------------------------------

