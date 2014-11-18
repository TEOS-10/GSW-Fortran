!==========================================================================
elemental function gsw_ct_from_enthalpy (sa, h, p)
!==========================================================================
!
!  Calculates the Conservative Temperature of seawater, given the Absolute 
!  Salinity, specific enthalpy, h, and pressure p.  The specific enthalpy 
!  input is the one calculated from the computationally-efficient 48-term
!  expression for specific volume in terms of SA, CT and p (IOC et al., 
!  2010).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  h   =  specific enthalpy                                        [ J/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325d0 dbar ) 
!
!  CT  =  Conservative Temperature ( ITS-90)                      [ deg C ]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_ct_freezing, gsw_enthalpy
use gsw_mod_toolbox, only : gsw_enthalpy_first_derivatives

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, h, p

real (r14) :: gsw_ct_from_enthalpy

real (r14) :: ct, ct_freezing, ct_mean, ct_old, f, h_freezing
real (r14) :: h_ct, h_40

ct_freezing = gsw_ct_freezing(sa,p,0d0)
h_freezing = gsw_enthalpy(sa,ct_freezing,p)

h_40 = gsw_enthalpy(sa,40d0,p)

! first guess of ct
ct = ct_freezing + (40d0 - ct_freezing)*(h - h_freezing)/(h_40 - h_freezing)
call gsw_enthalpy_first_derivatives(sa,ct,p,h_ct=h_ct)

!--------------------------------------------------------------------------
! Begin the modified Newton-Raphson iterative procedure 
!--------------------------------------------------------------------------

ct_old = ct
f = gsw_enthalpy(sa,ct_old,p) - h
ct = ct_old - f/h_ct
ct_mean = 0.5d0*(ct + ct_old)
call gsw_enthalpy_first_derivatives(sa,ct_mean,p,h_ct=h_ct)
ct = ct_old - f/h_ct

ct_old = ct
f = gsw_enthalpy(sa,ct_old,p) - h
ct = ct_old - f/h_ct

! After 1.5d0 iterations of this modified Newton-Raphson iteration,
! the error in CT is no larger than 4x10^-13 degrees C, which 
! is machine precision for this calculation. 

gsw_ct_from_enthalpy = ct

return
end function

!--------------------------------------------------------------------------
