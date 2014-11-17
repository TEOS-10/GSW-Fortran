!==========================================================================
function gsw_entropy_from_t(sa,t,p) 
!==========================================================================

! Calculates the specific entropy of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_entropy_from_t : specific entropy                    [J/(kg K)]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n0, n1
real (r14) :: sa, t, p, gsw_entropy_from_t, gsw_gibbs

n0 = 0
n1 = 1

gsw_entropy_from_t = -gsw_gibbs(n0,n1,n0,sa,t,p)

return
end

!--------------------------------------------------------------------------

