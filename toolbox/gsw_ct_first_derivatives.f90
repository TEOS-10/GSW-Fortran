!==========================================================================
elemental subroutine gsw_ct_first_derivatives (sa, pt, ct_sa, ct_pt)
!==========================================================================
!
!  Calculates the following two derivatives of Conservative Temperature
!  (1) CT_SA, the derivative with respect to Absolute Salinity at 
!      constant potential temperature (with pr = 0 dbar), and
!   2) CT_pt, the derivative with respect to potential temperature
!      (the regular potential temperature which is referenced to 0 dbar)
!      at constant Absolute Salinity.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  pt  =  potential temperature (ITS-90)                          [ deg C ]   
!         (whose reference pressure is 0 dbar)
!
!  CT_SA  =  The derivative of Conservative Temperature with respect to 
!            Absolute Salinity at constant potential temperature 
!            (the regular potential temperature which has reference 
!            sea pressure of 0 dbar).    
!            The CT_SA output has units of:                     [ K/(g/kg)]
!  CT_pt  =  The derivative of Conservative Temperature with respect to 
!            potential temperature (the regular one with pr = 0 dbar)
!            at constant SA. CT_pt is dimensionless.           [ unitless ]
!--------------------------------------------------------------------------

use gsw_mod_teos10_constants, only : gsw_cp0, gsw_sfac, gsw_t0

use gsw_mod_toolbox, only : gsw_gibbs_pt0_pt0

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, pt
real (r14), intent(out), optional :: ct_sa, ct_pt

real (r14) :: abs_pt, g_sa_mod, g_sa_t_mod, x, y_pt

abs_pt = gsw_t0 + pt 

if (present(ct_pt)) ct_pt = -(abs_pt*gsw_gibbs_pt0_pt0(sa,pt))/gsw_cp0

if (.not. present(ct_sa)) return

!--------------------------------------------------------------------------

x = sqrt(gsw_sfac*sa)
y_pt = 0.025d0*pt

g_sa_t_mod = 1187.3715515697959d0 + &
        x*(-1480.222530425046d0 + x*(2175.341332000392d0 + x*(-980.14153344888d0 + 220.542973797483d0*x) + &
        y_pt*(-548.4580073635929d0 + y_pt*(592.4012338275047d0 + y_pt*(-274.2361238716608d0 + 49.9394019139016d0*y_pt)))) + &
        y_pt*(-258.3988055868252d0 + y_pt*(-90.2046337756875d0 + y_pt*10.50720794170734d0))) + &
        y_pt*(3520.125411988816d0  + y_pt*(-1351.605895580406d0 + &
        y_pt*(731.4083582010072d0  + y_pt*(-216.60324087531103d0 + 25.56203650166196d0*y_pt))))
g_sa_t_mod = 0.5d0*gsw_sfac*0.025d0*g_sa_t_mod
   
g_sa_mod = 8645.36753595126d0 + &
        x*(-7296.43987145382d0 + x*(8103.20462414788d0 + &
        y_pt*(2175.341332000392d0 + y_pt*(-274.2290036817964d0 + &
        y_pt*(197.4670779425016d0 + y_pt*(-68.5590309679152d0 + 9.98788038278032d0*y_pt)))) + &
        x*(-5458.34205214835d0 - 980.14153344888d0*y_pt + &
        x*(2247.60742726704d0 - 340.1237483177863d0*x + 220.542973797483d0*y_pt))) + &
        y_pt*(-1480.222530425046d0 + &
        y_pt*(-129.1994027934126d0 + &
        y_pt*(-30.0682112585625d0 + y_pt*(2.626801985426835d0 ))))) + &
        y_pt*(1187.3715515697959d0 + &
        y_pt*(1760.062705994408d0 + y_pt*(-450.535298526802d0 + &
        y_pt*(182.8520895502518d0 + y_pt*(-43.3206481750622d0 + 4.26033941694366d0*y_pt)))))
g_sa_mod = 0.5d0*gsw_sfac*g_sa_mod   

ct_sa = (g_sa_mod - abs_pt*g_sa_t_mod)/gsw_cp0

return
end subroutine

!--------------------------------------------------------------------------
