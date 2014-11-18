!==========================================================================
elemental function gsw_sp_from_sstar (sstar, p, long, lat)  
!==========================================================================
!
! Calculates Practical Salinity, SP, from Preformed Salinity, Sstar. 
!
! sstar  : Preformed Salinity                              [g/kg]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sp_from_Sstar : Preformed Salinity                   [g/kg]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_saar, gsw_sp_from_sa_baltic

use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit

use gsw_mod_teos10_constants, only : gsw_ups

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sstar, p, long, lat  

real (r14) :: gsw_sp_from_sstar

real (r14) :: saar, sp_baltic

character (*), parameter :: func_name = "gsw_sp_from_sstar"

! In the Baltic Sea, SA = Sstar.
sp_baltic = gsw_sp_from_sa_baltic(sstar,long,lat)

if (sp_baltic .lt. 1d10) then

   gsw_sp_from_sstar = sp_baltic

else

   saar = gsw_saar(p,long,lat)
   if (saar .gt. gsw_error_limit) then
      gsw_sp_from_sstar = gsw_error_code(1,func_name,saar)
   else
      gsw_sp_from_sstar = (sstar/gsw_ups)/(1 - 0.35d0*saar)
   end if

end if

return
end function

!--------------------------------------------------------------------------
