!==========================================================================
function gsw_enthalpy_t_exact(sa,t,p) 
!==========================================================================

! Calculates the specific enthalpy of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_enthalpy_t_exact : specific enthalpy                 [J/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n0, n1
real (r14) :: sa, t, p, gsw_enthalpy_t_exact, gsw_gibbs

n0 = 0
n1 = 1

gsw_enthalpy_t_exact = gsw_gibbs(n0,n0,n0,sa,t,p) - (t+273.15d0)*gsw_gibbs(n0,n1,n0,sa,t,p)

return
end

!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Library functions of the GSW toolbox
!--------------------------------------------------------------------------

