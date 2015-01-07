!==========================================================================
elemental subroutine gsw_enthalpy_second_derivatives (sa, ct, p, &
                                                 h_sa_sa, h_sa_ct, h_ct_ct)
!==========================================================================
!
!  Calculates three second-order derivatives of specific enthalpy (h),
!  using the computationally-efficient 48-term expression for 
!  density in terms of sa, ct and p (IOC et al., 2010).
!
!  sa  =  Absolute Salinity                                        [ g/kg ]
!  ct  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  Sea Pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  h_sa_sa  =  The second derivative of specific enthalpy with respect to 
!              Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
!  h_sa_ct  =  The second derivative of specific enthalpy with respect to 
!              SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
!  h_ct_ct  =  The second derivative of specific enthalpy with respect to 
!              CT at constant SA and p.                      [ J/(kg K^2) ]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_dynamic_enthalpy
use gsw_mod_toolbox, only : gsw_enthalpy_first_derivatives
use gsw_mod_toolbox, only : gsw_enthalpy_second_derivatives_ct_exact

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, ct, p
real (r14), intent(out), optional :: h_sa_sa, h_sa_ct, h_ct_ct

real (r14) :: ct_l, ct_u, hct, hct_l, hct_u, hsa_l
real (r14) :: h_sa_sa_48, h_sa_sa_exact, hsa_u, sa_l, sa_taper, sa_u

real (r14), parameter :: dsa = 5d-2, dct = 5d-3

if (present(h_sa_sa)) then

   sa_l = max(sa - dsa, 0d0)
   sa_u = sa + dsa  

   call gsw_enthalpy_first_derivatives(sa_l,ct,p,h_sa=hsa_l)
   call gsw_enthalpy_first_derivatives(sa_u,ct,p,h_sa=hsa_u)

   h_sa_sa = (hsa_u - hsa_l)/(sa_u - sa_l)

   if (sa .lt. 1d0) then

      call gsw_enthalpy_second_derivatives_ct_exact(sa,ct,p, &
                                                    h_sa_sa=h_sa_sa_exact)
      if (sa .le. 0.5d0) then

         h_sa_sa = h_sa_sa_exact

      else

         h_sa_sa_48 = h_sa_sa
         sa_taper = 2d0 - 2d0*sa
         h_sa_sa = sa_taper*h_sa_sa_exact + (1d0-sa_taper)*h_sa_sa_48

      end if
   
   end if

end if

ct_l = ct - dct
ct_u = ct + dct

if (present(h_sa_ct)) then

   call gsw_enthalpy_first_derivatives(sa,ct_l,p,h_sa=hsa_l)
   call gsw_enthalpy_first_derivatives(sa,ct_u,p,h_sa=hsa_u)

   h_sa_ct  = (hsa_u - hsa_l)/(ct_u - ct_l)

end if

if (present(h_ct_ct)) then

   hct_l = gsw_dynamic_enthalpy(sa,ct_l,p)
   hct = gsw_dynamic_enthalpy(sa,ct,p)
   hct_u = gsw_dynamic_enthalpy(sa,ct_u,p)

   h_ct_ct = 4d0*(hct_u - 2d0*hct + hct_l)/((ct_u - ct_l)**2)

end if

return
end subroutine

!--------------------------------------------------------------------------
