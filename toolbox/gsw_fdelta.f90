!==========================================================================
elemental function gsw_fdelta (p, long, lat)
!==========================================================================
!
! Calculates fdelta. 
!
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_fdelta : Absolute Salinty Anomaly                    [unitless]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_saar

use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: p, long, lat

real (r14) :: gsw_fdelta

real (r14) ::  saar

character (*), parameter :: func_name = "gsw_fdelta"

saar = gsw_saar(p,long,lat)

if (saar.gt.gsw_error_limit) then
   gsw_fdelta = gsw_error_code(1,func_name)
else
   gsw_fdelta = ((1d0 + 0.35d0)*saar)/(1d0 - 0.35d0*saar)
end if

return
end function

!--------------------------------------------------------------------------
