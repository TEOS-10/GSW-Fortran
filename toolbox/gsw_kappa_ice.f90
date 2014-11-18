!==========================================================================
elemental function gsw_kappa_ice (t, p)
!==========================================================================
!
!  Calculates the isentropic compressibility of ice. 
!  
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  kappa_ice  =  isentropic compressibility                        [ 1/Pa ]
!   Note. The output units are 1/Pa not 1/dbar.
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_gibbs_ice

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: t, p

real (r14) :: gsw_kappa_ice

real (r14) :: gi_tp, gi_tt

gi_tt = gsw_gibbs_ice(2,0,t,p) 
gi_tp = gsw_gibbs_ice(1,1,t,p)

gsw_kappa_ice = (gi_tp*gi_tp - gi_tt*gsw_gibbs_ice(0,2,t,p))/ &
                  (gsw_gibbs_ice(0,1,t,p)*gi_tt)

return
end function

!--------------------------------------------------------------------------
