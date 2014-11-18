!==========================================================================
elemental function gsw_specvol_t_exact (sa, t, p)  
!==========================================================================
!
! Calulates the specific volume of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_specvol_t_exact : specific volume                    [kg/m^3]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_gibbs

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, t, p  

real (r14) :: gsw_specvol_t_exact

integer, parameter :: n0=0, n1=1

gsw_specvol_t_exact = gsw_gibbs(n0,n0,n1,sa,t,p)

return
end function

!--------------------------------------------------------------------------
