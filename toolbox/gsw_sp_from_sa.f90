!==========================================================================
function gsw_sp_from_sa(sa,p,long,lat) 
!==========================================================================

! Calculates Practical salinity, sp, from Absolute salinity, sa  
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [DEG E]     
! lat    : latitude                                        [DEG N]
!
! gsw_sp_from_sa      : Practical Salinity                 [unitless]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, long, lat, p, gsw_sp_from_sa, gsw_saar, saar
real (r14) :: gsw_sp_baltic, gsw_sp_from_sa_baltic

saar = gsw_saar(p,long,lat)

gsw_sp_from_sa = (35.d0/35.16504d0)*sa/(1d0 + saar)

gsw_sp_baltic = gsw_sp_from_sa_baltic(sa,long,lat);

if (gsw_sp_baltic.lt.1d10) then
   gsw_sp_from_sa = gsw_sp_baltic
end if

if (saar.eq.9d15) then
   gsw_sp_from_sa = 9d15
end if

return
end function

!--------------------------------------------------------------------------

