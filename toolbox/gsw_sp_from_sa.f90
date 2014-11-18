!==========================================================================
elemental function gsw_sp_from_sa (sa, p, long, lat) 
!==========================================================================
!
! Calculates Practical salinity, sp, from Absolute salinity, sa  
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [DEG E]     
! lat    : latitude                                        [DEG N]
!
! gsw_sp_from_sa      : Practical Salinity                 [unitless]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_saar, gsw_sp_from_sa_baltic

use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit

use gsw_mod_teos10_constants, only : gsw_ups

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, p, long, lat 

real (r14) :: gsw_sp_from_sa

real (r14) :: saar, sp_baltic

character (*), parameter :: func_name = "gsw_sp_from_sa"

sp_baltic = gsw_sp_from_sa_baltic(sa,long,lat)

if (sp_baltic .lt. 1d10) then

   gsw_sp_from_sa = sp_baltic

else

   saar = gsw_saar(p,long,lat)
   if (saar .gt. gsw_error_limit) then
      gsw_sp_from_sa = gsw_error_code(1,func_name,saar)
   else
      gsw_sp_from_sa = (sa/gsw_ups)/(1d0 + saar)
   end if

end if

return
end function

!--------------------------------------------------------------------------
