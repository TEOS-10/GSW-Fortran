!==========================================================================
elemental function gsw_adiabatic_lapse_rate_ice (t, p)
!==========================================================================
!
!  Calculates the adiabatic lapse rate of ice.
!
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!
!    Note.  The output is in unit of degress Celsius per Pa,
!      (or equivilently K/Pa) not in units of K/dbar. 
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_gibbs_ice

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: t, p

real (r14) :: gsw_adiabatic_lapse_rate_ice

gsw_adiabatic_lapse_rate_ice = -gsw_gibbs_ice(1,1,t,p)/(gsw_gibbs_ice(2,0,t,p))

return
end function

!--------------------------------------------------------------------------
