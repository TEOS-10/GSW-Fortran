!==========================================================================
function gsw_sa_from_sstar(sstar,p,long,lat)  
!==========================================================================

! Calculates Absolute Salinity, SA, from Preformed Salinity, Sstar.
!
! Sstar  : Preformed Salinity                              [g/kg]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sa_from_sstar   : Absolute Salinity                  [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, long, lat, p, gsw_saar, gsw_sp_from_sa_baltic
real (r14) :: saar, gsw_sa_from_sstar, sstar

saar = gsw_saar(p,long,lat)

gsw_sa_from_sstar = sstar*(1d0 + saar)/(1d0 - 0.35d0*saar)

! In the Baltic Sea, Sstar = SA, and note that gsw_saar returns zero
! for SAAR in the Baltic.

if (saar.eq.9d15) then
    gsw_sa_from_sstar = 9d15
end if

return
end function

!--------------------------------------------------------------------------

