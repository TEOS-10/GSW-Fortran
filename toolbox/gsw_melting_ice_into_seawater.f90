!==========================================================================
elemental subroutine gsw_melting_ice_into_seawater (sa, ct, p, &
                       saturation_fraction, w_ih, t_ih, sa_final, ct_final)
!==========================================================================
!
!  Calculates the Absolute Salinity and Conservative Temperature that 
!  results when a given mass of ice melts and is mixed into a known mass of
!  seawater (whose properties are (SA,CT,p)).  
!
!  SA  =  Absolute Salinity of seawater                            [ g/kg ]
!  CT  =  Conservative Temperature of seawater (ITS-90)           [ deg C ]
!  p   =  sea pressure at which the melting occurs                 [ dbar ]
!         ( i.e. absolute pressure - 10.1325d0 dbar ) 
!  saturation_fraction = the saturation fraction of dissolved air in 
!               seawater.  The saturation_fraction must be between 0 and 1.
!  w_Ih  =  mass fraction of ice, that is the mass of ice divided by the
!           sum of the masses of ice and seawater.  That is, the mass of 
!           ice divided by the mass of the final mixed fluid.  
!           w_Ih must be between 0 and 1.                      [ unitless ]
!  t_Ih  =  the in-situ temperature of the ice (ITS-90)           [ deg C ]
!
!  SA_final  =  Absolute Salinity of the mixture of the melted ice and the
!               orignal seawater                                   [ g/kg ]
!  CT_final  =  Conservative Temperature of the mixture of the melted ice 
!               and the orignal seawater                          [ deg C ] 
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_ct_freezing, gsw_enthalpy, gsw_t_freezing
use gsw_mod_toolbox, only : gsw_enthalpy_ice, gsw_ct_from_enthalpy

use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, ct, p, saturation_fraction, w_ih, t_ih
real (r14), intent(out) :: sa_final, ct_final

real (r14) :: ctf, h, h_final, h_ih, tf_ih

character (*), parameter :: func_name = "gsw_melting_ice_into_seawater"

ctf = gsw_ct_freezing(sa,p,saturation_fraction)
if (ct .lt. ctf) then
    ! The seawater ct input is below the freezing temp
    sa_final = gsw_error_code(1,func_name)
    ct_final = sa_final
    return
end if

tf_ih = gsw_t_freezing(0d0,p,saturation_fraction)
if (t_ih .gt. tf_ih) then
    ! t_ih input exceeds the freezing temp
    sa_final = gsw_error_code(2,func_name)
    ct_final = sa_final
    return
end if

h = gsw_enthalpy(sa,ct,p)
h_ih = gsw_enthalpy_ice(t_ih,p)

h_final = h - w_ih*(h - h_ih)

sa_final = sa*(1d0 - w_ih)

ctf = gsw_ct_freezing(sa_final,p,saturation_fraction)

if (h_final .lt. gsw_enthalpy(sa_final,ctf,p)) then
    ! Melting this much ice is not possible as it would result in
    ! frozen seawater
    sa_final = gsw_error_code(3,func_name)
    ct_final = sa_final
    return
end if

ct_final = gsw_ct_from_enthalpy(sa_final,h_final,p)

return
end subroutine

!--------------------------------------------------------------------------
