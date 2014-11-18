!==========================================================================
elemental subroutine gsw_ct_first_derivatives_wrt_t_exact (sa, t, p, &
				       ct_sa_wrt_t, ct_t_wrt_t, ct_p_wrt_t)
!==========================================================================
!
!  Calculates the following three derivatives of Conservative Temperature.
!  These derivatives are done with respect to in-situ temperature t (in the
!  case of CT_T_wrt_t) or at constant in-situ tempertature (in the cases of
!  CT_SA_wrt_t and CT_P_wrt_t).  
!   (1) CT_SA_wrt_t, the derivative of CT with respect to Absolute Salinity 
!       at constant t and p, and
!   (2) CT_T_wrt_t, derivative of CT with respect to in-situ temperature t 
!       at constant SA and p. 
!   (3) CT_P_wrt_t, derivative of CT with respect to pressure P (in Pa) at  
!       constant SA and t.    
!
!  This function uses the full Gibbs function. Note that this function
!  avoids the NaN that would exist in CT_SA_wrt_t at SA = 0 if it were
!  evaluated in the straightforward way from the derivatives of the Gibbs 
!  function function.
!   
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  t   =  in-situ temperature (ITS-90)                            [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar)
!
!  CT_SA_wrt_t  =  The first derivative of Conservative Temperature with 
!                  respect to Absolute Salinity at constant t and p.     
!                                              [ K/(g/kg)]  i.e. [ K kg/g ]
!  CT_T_wrt_t  =  The first derivative of Conservative Temperature with 
!                 respect to in-situ temperature, t, at constant SA and p.     
!                                                              [ unitless ]
!  CT_P_wrt_t  =  The first derivative of Conservative Temperature with 
!                 respect to pressure P (in Pa) at constant SA and t. 
!                                                                  [ K/Pa ]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_gibbs, gsw_pt0_from_t

use gsw_mod_teos10_constants, only : gsw_cp0, gsw_sfac, gsw_t0, rec_db2pa

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, t, p
real (r14), intent(out) :: ct_p_wrt_t, ct_sa_wrt_t, ct_t_wrt_t

real (r14) :: g_sa_mod, g_sa_t_mod, pt0, x, y, y_pt, z

pt0 = gsw_pt0_from_t(sa,t,p)

!  Note that the following call is an alternative to the code from lines
!  121 - 154, however, at SA = 0 the following will return a NaN, whereas 
!  the code from 121 - 154 will not.
!  CT_SA_wrt_t = (gsw_gibbs(1,0,0,SA,pt0,zeros(size(SA))) ...
!                - (273.15+pt0).*gsw_gibbs(1,1,0,SA,t,p))...
!                   ./gsw_cp0;

x = sqrt(gsw_sfac*sa)
y = 0.025d0*t
y_pt = 0.025d0*pt0
z = rec_db2pa*p !note. the input pressure (p) is sea pressure in units of dbar.

g_sa_t_mod = 1187.3715515697959d0 + z*(1458.233059470092d0 + &
        z*(-687.913805923122d0 + z*(249.375342232496d0 + &
	z*(-63.313928772146d0 + 14.09317606630898d0*z)))) + &
        x*(-1480.222530425046d0 + x*(2175.341332000392d0 + &
	x*(-980.14153344888d0 + 220.542973797483d0*x) + &
        y*(-548.4580073635929d0 + y*(592.4012338275047d0 + &
	y*(-274.2361238716608d0 + 49.9394019139016d0*y))) - &
        90.6734234051316d0*z) + z*(-525.876123559641d0 + &
	(249.57717834054571d0 - 88.449193048287d0*z)*z) + &
        y*(-258.3988055868252d0 + z*(2298.348396014856d0 + &
	z*(-325.1503575102672d0 + 153.8390924339484d0*z)) + &
        y*(-90.2046337756875d0 - 4142.8793862113125d0*z + &
	y*(10.50720794170734d0 + 2814.78225133626d0*z)))) + &
        y*(3520.125411988816d0 + y*(-1351.605895580406d0 + &
        y*(731.4083582010072d0 + y*(-216.60324087531103d0 + &
	25.56203650166196d0*y) + z*(-2381.829935897496d0 + &
	(597.809129110048d0 - 291.8983352012704d0*z)*z)) + &
        z*(4165.4688847996085d0 + z*(-1229.337851789418d0 + &
	(681.370187043564d0 - 66.7696405958478d0*z)*z))) + &
        z*(-3443.057215135908d0 + z*(1349.638121077468d0 + &
        z*(-713.258224830552d0 + &
	(176.8161433232d0 - 31.68006188846728d0*z)*z))))
g_sa_t_mod = 0.5d0*gsw_sfac*0.025d0*g_sa_t_mod
   
g_sa_mod = 8645.36753595126d0 + &
        x*(-7296.43987145382d0 + x*(8103.20462414788d0 + &
        y_pt*(2175.341332000392d0 + y_pt*(-274.2290036817964d0 + &
        y_pt*(197.4670779425016d0 + y_pt*(-68.5590309679152d0 + &
	9.98788038278032d0*y_pt)))) + &
        x*(-5458.34205214835d0 - 980.14153344888d0*y_pt + &
        x*(2247.60742726704d0 - 340.1237483177863d0*x + &
	220.542973797483d0*y_pt))) + &
        y_pt*(-1480.222530425046d0 + &
        y_pt*(-129.1994027934126d0 + &
        y_pt*(-30.0682112585625d0 + y_pt*(2.626801985426835d0 ))))) + &
        y_pt*(1187.3715515697959d0 + &
        y_pt*(1760.062705994408d0 + y_pt*(-450.535298526802d0 + &
        y_pt*(182.8520895502518d0 + y_pt*(-43.3206481750622d0 + &
	4.26033941694366d0*y_pt)))))
g_sa_mod = 0.5d0*gsw_sfac*g_sa_mod   

ct_sa_wrt_t = (g_sa_mod - (gsw_t0+pt0)*g_sa_t_mod)/gsw_cp0

ct_t_wrt_t = -(gsw_t0+pt0)*gsw_gibbs(0,2,0,sa,t,p)/gsw_cp0

ct_p_wrt_t = -(gsw_t0+pt0)*gsw_gibbs(0,1,1,sa,t,p)/gsw_cp0
                       
return
end subroutine

!--------------------------------------------------------------------------
