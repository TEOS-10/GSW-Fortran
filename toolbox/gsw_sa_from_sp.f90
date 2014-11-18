!==========================================================================
elemental function gsw_sa_from_sp (sp, p, long, lat)       
!==========================================================================
!
! Calculates Absolute Salinity, SA, from Practical Salinity, SP
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [DEG E]     
! lat    : latitude                                        [DEG N]
!
! gsw_sa_from_sp   : Absolute Salinity                     [g/kg]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_saar, gsw_sa_from_sp_baltic

use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit

use gsw_mod_teos10_constants, only : gsw_ups

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sp, p, long, lat       

real (r14) :: gsw_sa_from_sp

real (r14) :: saar, sa_baltic

character (*), parameter :: func_name = "gsw_sa_from_sp"

sa_baltic = gsw_sa_from_sp_baltic(sp,long,lat)

if (sa_baltic .lt. 1d10) then

   gsw_sa_from_sp = sa_baltic

else

   saar = gsw_saar(p,long,lat)
   if (saar .gt. gsw_error_limit) then
      gsw_sa_from_sp = gsw_error_code(1,func_name,saar)
   else
      gsw_sa_from_sp = gsw_ups*sp*(1d0 + saar)
   end if

end if

return
end function

!--------------------------------------------------------------------------
