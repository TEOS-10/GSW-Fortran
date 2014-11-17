!==========================================================================
function gsw_internal_energy(sa,ct,p)  
!==========================================================================

!  Calculates internal energy of seawater using the computationally
!  efficient 48-term expression for density in terms of SA, CT and p
!  (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_internal_energy  :  internal_energy of seawater (48 term equation)

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: p0 = 101325d0, db2pa = 1d4

real (r14) :: sa, ct, p, gsw_internal_energy, gsw_enthalpy, gsw_specvol

gsw_internal_energy = gsw_enthalpy(sa,ct,p) - (p0 + db2pa*p)*gsw_specvol(sa,ct,p)

return
end function

!--------------------------------------------------------------------------

