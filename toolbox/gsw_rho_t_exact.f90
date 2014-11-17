!==========================================================================
function gsw_rho_t_exact(sa,t,p) 
!==========================================================================

! Calculates in-situ density of seawater from Absolute Salinity and 
! in-situ temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_rho_t_exact : in-situ density                        [kg/m^3]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n0, n1
real (r14) :: sa, t, p, gsw_rho_t_exact, gsw_gibbs

n0 = 0
n1 = 1

gsw_rho_t_exact = 1.d0/gsw_gibbs(n0,n0,n1,sa,t,p)

return
end

!--------------------------------------------------------------------------

