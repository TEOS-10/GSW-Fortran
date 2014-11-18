!==========================================================================
elemental function gsw_sstar_from_sa (sa, p, long, lat) 
!==========================================================================
!
! Calculates Preformed Salinity, Sstar, from Absolute Salinity, SA. 
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sstar_from_sa : Preformed Salinity                   [g/kg]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_saar

use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, p, long, lat 

real (r14) :: gsw_sstar_from_sa

real (r14) :: saar

character (*), parameter :: func_name = "gsw_sstar_from_sa"

! In the Baltic Sea, Sstar = sa, and note that gsw_saar returns zero
! for saar in the Baltic.

saar = gsw_saar(p,long,lat)

if (saar.gt.gsw_error_limit) then
   gsw_sstar_from_sa = gsw_error_code(1,func_name,saar)
else
   gsw_sstar_from_sa = sa*(1d0 - 0.35d0*saar)/(1d0 + saar)
end if

return
end function

!--------------------------------------------------------------------------
