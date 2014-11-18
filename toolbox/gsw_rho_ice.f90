!==========================================================================
elemental function gsw_rho_ice (t, p)
!==========================================================================
! 
!  Calculates in-situ density of ice from in-situ temperature and pressure.
!  Note that the output, rho_ice, is density, not density anomaly;  that 
!  is, 1000 kg/m^3 is not subracted from it.  
!
!  t   =  in-situ temperature (ITS-90)                            [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  rho_ice  =  in-situ density of ice (not density anomaly)      [ kg/m^3 ]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_gibbs_ice

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: t, p

real (r14) :: gsw_rho_ice

gsw_rho_ice = 1.d0/gsw_gibbs_ice(0,1,t,p)

return
end function

!--------------------------------------------------------------------------
