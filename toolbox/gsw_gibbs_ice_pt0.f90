! =========================================================================
elemental function gsw_gibbs_ice_pt0 (pt0)
! =========================================================================
!
!  Part of the the first temperature derivative of Gibbs energy of ice
!  that is the outout is "gibbs_ice(1,0,pt0,0) + s0"
!
!  pt0  =  potential temperature with reference sea pressure of zero dbar
!                                                                 [ deg C ]
!
!  gsw_gibbs_ice_pt0 = part of temperature derivative     [ J kg^-1 K^-1 ]
!--------------------------------------------------------------------------

use gsw_mod_teos10_constants, only : gsw_t0

use gsw_mod_gibbs_ice_coefficients

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: pt0

real (r14) :: gsw_gibbs_ice_pt0

real (r14) :: tau
complex (r14) :: g, tau_t1, tau_t2

tau = (pt0 + gsw_t0)*rec_tt

tau_t1 = tau/t1
tau_t2 = tau/t2

g = r1*(log((1 + tau_t1)/(1 - tau_t1)) - 2*tau_t1) &
    + r20*(log((1 + tau_t2)/(1 - tau_t2)) - 2*tau_t2)

gsw_gibbs_ice_pt0 = real(g)

return
end function

!--------------------------------------------------------------------------
