!==========================================================================
elemental function gsw_ct_from_pt (sa, pt) 
!==========================================================================
!
! Calculates Conservative Temperature from potential temperature of seawater  
!
! sa      : Absolute Salinity                              [g/kg]
! pt      : potential temperature with                     [deg C]
!           reference pressure of 0 dbar
!
! gsw_ct_from_pt : Conservative Temperature                [deg C]
!--------------------------------------------------------------------------

use gsw_mod_teos10_constants, only : gsw_cp0, gsw_sfac

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, pt 

real (r14) :: gsw_ct_from_pt

real (r14) :: pot_enthalpy, x2, x, y

x2 = gsw_sfac*sa
x = sqrt(x2)
y = pt*0.025d0        ! normalize for F03 and F08

pot_enthalpy =  61.01362420681071d0 + y*(168776.46138048015d0 + &
               y*(-2735.2785605119625d0 + y*(2574.2164453821433d0 + &
               y*(-1536.6644434977543d0 + y*(545.7340497931629d0 + &
               (-50.91091728474331d0 - 18.30489878927802d0*y)*y))))) + &
               x2*(268.5520265845071d0 + y*(-12019.028203559312d0 + &
               y*(3734.858026725145d0 + y*(-2046.7671145057618d0 + &
               y*(465.28655623826234d0 + (-0.6370820302376359d0 - &
               10.650848542359153d0*y)*y)))) + &
               x*(937.2099110620707d0 + y*(588.1802812170108d0 + &
               y*(248.39476522971285d0 + (-3.871557904936333d0 - &
               2.6268019854268356d0*y)*y)) + &
               x*(-1687.914374187449d0 + x*(246.9598888781377d0 + &
               x*(123.59576582457964d0 - 48.5891069025409d0*x)) + &
               y*(936.3206544460336d0 + &
               y*(-942.7827304544439d0 + y*(369.4389437509002d0 + &
               (-33.83664947895248d0 - 9.987880382780322d0*y)*y))))))

gsw_ct_from_pt = pot_enthalpy/gsw_cp0

return
end function

!--------------------------------------------------------------------------
