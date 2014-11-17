!==========================================================================
function gsw_fdelta(p,long,lat)
!==========================================================================

! Calculates fdelta. 
!
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_fdelta : Absolute Salinty Anomaly                    [unitless]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) ::  long, lat, p, gsw_saar, saar, gsw_fdelta

saar = gsw_saar(p,long,lat)

gsw_fdelta = ((1d0 + 0.35d0)*saar)/(1d0 - 0.35d0*saar);

if (saar.gt.1d10) then
    gsw_fdelta = 9d15
end if

return
end function

!--------------------------------------------------------------------------

