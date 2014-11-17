!==========================================================================
function gsw_beta_const_t_exact(sa,t,p)  
!==========================================================================

! Calculates saline (haline) contraction coefficient of seawater at 
! constant in-situ temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_beta_const_t_exact : haline contraction coefficient  [kg/g]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n0, n1
real (r14) :: sa, t, p, gsw_beta_const_t_exact, gsw_gibbs, uPS

n0 = 0
n1 = 1 

gsw_beta_const_t_exact = -gsw_gibbs(n1,n0,n1,sa,t,p)/gsw_gibbs(n0,n0,n1,sa,t,p)

return
end function

!--------------------------------------------------------------------------

