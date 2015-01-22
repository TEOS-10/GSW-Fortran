!==========================================================================
elemental subroutine gsw_enthalpy_first_derivatives_ct_exact (sa, ct, p, &
                                                              h_sa, h_ct)
!==========================================================================
!
!  Calculates the following two derivatives of specific enthalpy (h)
!   (1) h_SA, the derivative with respect to Absolute Salinity at 
!       constant CT and p, and
!   (2) h_CT, derivative with respect to CT at constant SA and p. 
!  Note that h_P is specific volume (1/rho) it can be calulated by calling
!  gsw_specvol_CT_exact(SA,CT,p).
!
!  Note that this function uses the full Gibbs function.  There is an 
!  alternative to calling gsw_enthalpy_first_derivatives(SA,CT,p)
!  which uses the computationally
!  efficient 48-term expression for density in terms of SA, CT and p 
!  (IOC et al., 2010).   
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  h_SA  =  The first derivative of specific enthalpy with respect to 
!           Absolute Salinity at constant CT and p.     
!                                            [ J/(kg (g/kg))]  i.e. [ J/g ]
!  h_CT  =  The first derivative of specific enthalpy with respect to 
!           CT at constant SA and p.                           [ J/(kg K) ]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_pt_from_ct, gsw_t_from_ct

use gsw_mod_teos10_constants, only : gsw_cp0, gsw_sfac, gsw_t0, rec_db2pa

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, ct, p
real (r14), intent(out), optional :: h_sa, h_ct

real (r14) :: g_sa_mod_pt, g_sa_mod_t, pt0, t, temp_ratio
real (r14) :: x, y, y_pt, z

t = gsw_t_from_ct(sa,ct,p)
pt0 = gsw_pt_from_ct(sa,ct)  

temp_ratio = (gsw_t0 + t)/(gsw_t0 + pt0)

if (present(h_ct)) h_ct = gsw_cp0*temp_ratio

if (.not. present(h_sa)) return

x = sqrt(gsw_sfac*sa)
y = 0.025d0*t
z = rec_db2pa*p !note.the input pressure (p) is sea pressure in units of dbar.

g_sa_mod_t = 8645.36753595126d0 + z*(-6620.98308089678d0 + &
       z*(769.588305957198d0 + z*(-193.0648640214916d0 + (31.6816345533648d0 - 5.24960313181984d0*z)*z))) + &
       x*(-7296.43987145382d0 + x*(8103.20462414788d0 + &
       y*(2175.341332000392d0 + y*(-274.2290036817964d0 + &
       y*(197.4670779425016d0 + y*(-68.5590309679152d0 + 9.98788038278032d0*y))) - 90.6734234051316d0*z) + &
       x*(-5458.34205214835d0 - 980.14153344888d0*y + &
       x*(2247.60742726704d0 - 340.1237483177863d0*x + 220.542973797483d0*y) + 180.142097805543d0*z) + &
       z*(-219.1676534131548d0 + (-16.32775915649044d0 - 120.7020447884644d0*z)*z)) + &
       z*(598.378809221703d0 + z*(-156.8822727844005d0 + (204.1334828179377d0 - 10.23755797323846d0*z)*z)) + &
       y*(-1480.222530425046d0 + z*(-525.876123559641d0 + (249.57717834054571d0 - 88.449193048287d0*z)*z) + &
       y*(-129.1994027934126d0 + z*(1149.174198007428d0 + z*(-162.5751787551336d0 + 76.9195462169742d0*z)) + &
       y*(-30.0682112585625d0 - 1380.9597954037708d0*z + y*(2.626801985426835d0 + 703.695562834065d0*z))))) + &
       y*(1187.3715515697959d0 + z*(1458.233059470092d0 + &
       z*(-687.913805923122d0 + z*(249.375342232496d0 + z*(-63.313928772146d0 + 14.09317606630898d0*z)))) + &
       y*(1760.062705994408d0 + y*(-450.535298526802d0 + &
       y*(182.8520895502518d0 + y*(-43.3206481750622d0 + 4.26033941694366d0*y) + &
       z*(-595.457483974374d0 + (149.452282277512d0 - 72.9745838003176d0*z)*z)) + &
       z*(1388.489628266536d0 + z*(-409.779283929806d0 + (227.123395681188d0 - 22.2565468652826d0*z)*z))) + &
       z*(-1721.528607567954d0 + z*(674.819060538734d0 + &
       z*(-356.629112415276d0 + (88.4080716616d0 - 15.84003094423364d0*z)*z)))))
  
g_sa_mod_t = 0.5d0*gsw_sfac*g_sa_mod_t   
    
y_pt = 0.025d0*pt0

g_sa_mod_pt = 8645.36753595126d0 + &
        x*(-7296.43987145382d0 + x*(8103.20462414788d0 + &
        y_pt*(2175.341332000392d0 + y_pt*(-274.2290036817964d0 + &
        y_pt*(197.4670779425016d0 + y_pt*(-68.5590309679152d0 + 9.98788038278032d0*y_pt)))) + &
        x*(-5458.34205214835d0 - 980.14153344888d0*y_pt + &
        x*(2247.60742726704d0 - 340.1237483177863d0*x + 220.542973797483d0*y_pt))) + &
        y_pt*(-1480.222530425046d0 + y_pt*(-129.1994027934126d0 + &
        y_pt*(-30.0682112585625d0 + y_pt*2.626801985426835d0)))) + &
        y_pt*(1187.3715515697959d0 + y_pt*(1760.062705994408d0 + y_pt*(-450.535298526802d0 + &
        y_pt*(182.8520895502518d0 + y_pt*(-43.3206481750622d0 + 4.26033941694366d0*y_pt)))))
    
g_sa_mod_pt = 0.5d0*gsw_sfac*g_sa_mod_pt   

h_sa = g_sa_mod_t - temp_ratio*g_sa_mod_pt

return
end subroutine

!--------------------------------------------------------------------------
