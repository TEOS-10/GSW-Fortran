!==========================================================================
function gsw_sstar_from_sp(sp,p,long,lat) 
!==========================================================================

! Calculates Preformed Salinity, Sstar, from Practical Salinity, SP. 
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sstar_from_sp  : Preformed Salinity                  [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sp, long, lat, p, gsw_saar, gsw_sa_from_sp_baltic
real (r14) :: saar, gsw_sstar_from_sp, sstar_baltic

saar = gsw_saar(p,long,lat)

gsw_sstar_from_sp = (35.16504d0/35.d0)*sp*(1 - 0.35d0*saar);

!In the Baltic Sea, Sstar = SA.
sstar_baltic = gsw_sa_from_sp_baltic(sp,long,lat);

if (sstar_baltic.lt.1d10) then
    gsw_sstar_from_sp = sstar_baltic;
end if

if (saar.eq.9d15) then
    gsw_sstar_from_sp = 9d15
end if

return
end function

!--------------------------------------------------------------------------

