!==========================================================================
elemental function gsw_specvol_ice (t, p)
!==========================================================================
!
!  Calculates the specific volume of ice. 
! 
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  specvol_ice  =  specific volume                               [ m^3/kg ]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_gibbs_ice

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: t, p

real (r14) :: gsw_specvol_ice

gsw_specvol_ice = gsw_gibbs_ice(0,1,t,p)

return
end function

!--------------------------------------------------------------------------
