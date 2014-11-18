!==========================================================================
elemental function gsw_chem_potential_water_ice (t, p)
!==========================================================================
! 
!  Calculates the chemical potential of water in ice from in-situ
!  temperature and pressure.
!
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  chem_potential_water_ice  =  chemical potential of ice          [ J/kg ]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_gibbs_ice

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: t, p

real (r14) :: gsw_chem_potential_water_ice

gsw_chem_potential_water_ice = gsw_gibbs_ice(0,0,t,p)

return
end function

!--------------------------------------------------------------------------
