!==========================================================================
function gsw_sa_from_sp(sp,p,long,lat)       
!==========================================================================

! Calculates Absolute Salinity, SA, from Practical Salinity, SP
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [DEG E]     
! lat    : latitude                                        [DEG N]
!
! gsw_sa_from_sp   : Absolute Salinity                     [g/kg]

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sp, long, lat, p, gsw_sa_from_sp, gsw_saar, saar
real (r14) :: gsw_sa_baltic, gsw_sa_from_sp_baltic

saar = gsw_saar(p,long,lat)

gsw_sa_from_sp = (35.16504d0/35.d0)*sp*(1.d0 + saar)

gsw_sa_baltic = gsw_sa_from_sp_baltic(sp,long,lat)

if (gsw_sa_baltic.lt.1d10) then
   gsw_sa_from_sp = gsw_sa_baltic
end if

if (saar.eq.9d15) then
   gsw_sa_from_sp = 9d15
end if

return
end function

!--------------------------------------------------------------------------

