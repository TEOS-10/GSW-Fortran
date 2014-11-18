!==========================================================================
elemental function gsw_latentheat_evap_ct (sa, ct) 
!==========================================================================
!
! Calculates latent heat, or enthalpy, of evaporation.
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! 
! gsw_latentheat_evaporation : latent heat of evaporation  [J/kg]
!--------------------------------------------------------------------------

use gsw_mod_teos10_constants, only : gsw_sfac

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, ct 

real (r14) :: gsw_latentheat_evap_ct

real (r14), parameter :: c0 =   2.499065844825125d6, c1 =  -1.544590633515099d-1
real (r14), parameter :: c2 =  -9.096800915831875d4, c3 =   1.665513670736000d2
real (r14), parameter :: c4 =   4.589984751248335d1, c5 =   1.894281502222415d1
real (r14), parameter :: c6 =   1.192559661490269d3, c7 =  -6.631757848479068d3
real (r14), parameter :: c8 =  -1.104989199195898d2, c9 =  -1.207006482532330d3
real (r14), parameter :: c10 = -3.148710097513822d3, c11 =  7.437431482069087d2
real (r14), parameter :: c12 =  2.519335841663499d3, c13 =  1.186568375570869d1
real (r14), parameter :: c14 =  5.731307337366114d2, c15 =  1.213387273240204d3
real (r14), parameter :: c16 =  1.062383995581363d3, c17 = -6.399956483223386d2
real (r14), parameter :: c18 = -1.541083032068263d3, c19 =  8.460780175632090d1
real (r14), parameter :: c20 = -3.233571307223379d2, c21 = -2.031538422351553d2
real (r14), parameter :: c22 =  4.351585544019463d1, c23 = -8.062279018001309d2
real (r14), parameter :: c24 =  7.510134932437941d2, c25 =  1.797443329095446d2
real (r14), parameter :: c26 = -2.389853928747630d1, c27 =  1.021046205356775d2

real (r14) :: x, y

x = sqrt(gsw_sfac*sa)
y = ct/40d0

gsw_latentheat_evap_ct = c0 + x*(c1 + c4*y + x*(c3   &
    + y*(c7 + c12*y) + x*(c6 + y*(c11 + y*(c17 + c24*y)) &
    + x*(c10 + y*(c16 + c23*y) + x*(c15 + c22*y + c21*x)))))  &
    + y*(c2 + y*(c5 + c8*x + y*(c9 + x*(c13 + c18*x)  &
    + y*(c14 + x*(c19 + c25*x) + y*(c20 + c26*x + c27*y)))))

return
end function

!--------------------------------------------------------------------------
