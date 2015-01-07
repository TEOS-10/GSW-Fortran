!==========================================================================
elemental subroutine gsw_enthalpy_first_derivatives (sa, ct, p, h_sa, h_ct)
!==========================================================================
!
!  Calculates the following two derivatives of specific enthalpy (h) of
!  seawater using the computationally-efficient 48-term expression for 
!  density in terms of SA, CT and p (IOC et al., 2010).  
!   (1) h_SA, the derivative with respect to Absolute Salinity at 
!       constant CT and p, and
!   (2) h_CT, derivative with respect to CT at constant SA and p. 
!  Note that h_P is specific volume (1/rho) it can be calulated by calling
!  gsw_specvol(SA,CT,p).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325d0 dbar )
!
!  h_SA  =  The first derivative of specific enthalpy with respect to 
!           Absolute Salinity at constant CT and p.     
!                                            [ J/(kg (g/kg))]  i.e. [ J/g ]
!  h_CT  =  The first derivative of specific enthalpy with respect to 
!           CT at constant SA and p.                           [ J/(kg K) ]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_dynamic_enthalpy

use gsw_mod_teos10_constants, only : gsw_cp0

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, ct, p
real (r14), intent(out), optional :: h_sa, h_ct

real (r14) :: ct_l, ct_u, dct, dsa, hct_l, hct_u, hsa_l, hsa_u, sa_l, sa_u

real (r14), parameter :: dsa = 0.01d0, dct = 0.005d0

if (present(h_sa)) then
    sa_l = max(sa - dsa, 0.d0)
    sa_u = sa + dsa

    hsa_l = gsw_dynamic_enthalpy(sa_l,ct,p)
    hsa_u = gsw_dynamic_enthalpy(sa_u,ct,p)

    h_sa = (hsa_u - hsa_l)/(sa_u - sa_l)
endif

if (present(h_ct)) then
    ct_l = ct - dct
    ct_u = ct + dct

    hct_l = gsw_dynamic_enthalpy(sa,ct_l,p)
    hct_u = gsw_dynamic_enthalpy(sa,ct_u,p)

    h_ct = gsw_cp0 + (hct_u - hct_l)/(ct_u - ct_l)
endif

return
end subroutine

!--------------------------------------------------------------------------
