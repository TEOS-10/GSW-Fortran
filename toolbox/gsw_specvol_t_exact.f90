!==========================================================================
function gsw_specvol_t_exact(sa,t,p)  
!==========================================================================

! Calulates the specific volume of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_specvol_t_exact : specific volume                    [kg/m^3]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n0, n1
real (r14) :: sa, t, p, gsw_specvol_t_exact, gsw_gibbs

n0 = 0
n1 = 1

gsw_specvol_t_exact = gsw_gibbs(n0,n0,n1,sa,t,p)

return
end

!--------------------------------------------------------------------------

