!==========================================================================
function gsw_sstar_from_sa(sa,p,long,lat) 
!==========================================================================

! Calculates Preformed Salinity, Sstar, from Absolute Salinity, SA. 
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sstar_from_sa : Preformed Salinity                   [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, long, lat, p, gsw_saar, gsw_sp_from_sa_baltic
real (r14) :: saar, gsw_sstar_from_sa

saar = gsw_saar(p,long,lat)

gsw_sstar_from_sa = sa*(1d0 - 0.35d0*saar)/(1d0 + saar)

! In the Baltic Sea, Sstar = sa, and note that gsw_saar returns zero
! for saar in the Baltic.

if (saar.eq.9d15) then
    gsw_sstar_from_sa = 9d15
end if

return
end function

!--------------------------------------------------------------------------

