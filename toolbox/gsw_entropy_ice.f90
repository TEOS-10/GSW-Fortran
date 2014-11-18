!==========================================================================
elemental function gsw_entropy_ice (t, p)
!==========================================================================
!
!  Calculates specific entropy of ice. 
!
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!
!  ice_entropy  =  specific entropy of ice                 [ J kg^-1 K^-1 ]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_gibbs_ice

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: t, p

real (r14) :: gsw_entropy_ice

gsw_entropy_ice = -gsw_gibbs_ice(1,0,t,p)

return
end function

!--------------------------------------------------------------------------
