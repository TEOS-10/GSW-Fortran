!==========================================================================
function gsw_sp_from_sstar(sstar,p,long,lat)  
!==========================================================================

! Calculates Practical Salinity, SP, from Preformed Salinity, Sstar. 
!
! sstar  : Preformed Salinity                              [g/kg]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sp_from_Sstar : Preformed Salinity                   [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: long, lat, p, gsw_saar, gsw_sp_from_sa_baltic
real (r14) :: saar, gsw_sp_from_sstar, sp_baltic, Sstar

saar = gsw_saar(p,long,lat)

gsw_sp_from_sstar = (35.d0/35.16504d0)*Sstar/(1 - 0.35d0*saar);

!In the Baltic Sea, SA = Sstar.
sp_baltic = gsw_sp_from_sa_baltic(sstar,long,lat);

if (sp_baltic.lt.1d10) then
    gsw_sp_from_sstar = sp_baltic;
end if

if (saar.eq.9d15) then
    gsw_sp_from_sstar = 9d15
end if

return
end function

!--------------------------------------------------------------------------

