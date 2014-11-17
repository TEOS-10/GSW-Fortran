!==========================================================================
function gsw_specvol(sa,ct,p) 
!==========================================================================

!  Calculates specific volume of seawater using the computationally-
!  efficient 48-term expression for density in terms of SA, CT and p
!  (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_specvol  :  specific volume of seawater (48 term equation)

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, ct, p, gsw_specvol, gsw_rho

gsw_specvol = 1d0/gsw_rho(sa,ct,p)

return
end function

!--------------------------------------------------------------------------

