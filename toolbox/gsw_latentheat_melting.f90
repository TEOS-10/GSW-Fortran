!==========================================================================
elemental function gsw_latentheat_melting (sa, p)  
!==========================================================================
!
! Calculates latent heat, or enthalpy, of melting.
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! 
! gsw_latentheat_melting : latent heat of melting          [kg/m^3]
!--------------------------------------------------------------------------

use gsw_mod_teos10_constants, only : gsw_sfac

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: sa, p  

real (r8) :: gsw_latentheat_melting

real (r8), parameter :: c0 =  3.334265169240710e5_r8
real (r8), parameter :: c1 = -2.789444646733159_r8
real (r8), parameter :: c2 = -1.822150156453350e4_r8
real (r8), parameter :: c3 = -4.984585692734338e3_r8
real (r8), parameter :: c4 = -7.371966528571920e1_r8
real (r8), parameter :: c5 = -7.605802553358546e3_r8
real (r8), parameter :: c6 =  1.195857305019339e3_r8
real (r8), parameter :: c7 =  1.233720336206392e3_r8
real (r8), parameter :: c8 =  2.294798676591890e2_r8
real (r8), parameter :: c9 =  9.655751370889338e2_r8
real (r8), parameter :: c10 = -5.792068522727968e2_r8
real (r8), parameter :: c11 = -1.649446955902331e3_r8
real (r8), parameter :: c12 = -1.029021448430547e3_r8
real (r8), parameter :: c13 = -3.171558017172501e2_r8
real (r8), parameter :: c14 = -1.751401389905041e2_r8
real (r8), parameter :: c15 =  6.836527214265952e2_r8
real (r8), parameter :: c16 =  1.078283734113611e3_r8
real (r8), parameter :: c17 =  5.613896351265648e2_r8
real (r8), parameter :: c18 =  6.968934948667265e2_r8
real (r8), parameter :: c19 =  1.793032021946783e2_r8
real (r8), parameter :: c20 =  8.692558481134256e1_r8
real (r8), parameter :: c21 = -2.371103254714944e2_r8
real (r8), parameter :: c22 = -5.775033277201674e2_r8
real (r8), parameter :: c23 = -3.019749254648732e2_r8
real (r8), parameter :: c24 = -6.420420579160927e2_r8
real (r8), parameter :: c25 = -2.657570848596042e2_r8
real (r8), parameter :: c26 = -1.646738151143109e1_r8
real (r8), parameter :: c27 =  4.618228988300871_r8

real (r8) :: x, y

x = sqrt(gsw_sfac*sa)
y = p*1e-4_r8

gsw_latentheat_melting = c0 + x*(c1 + c4*y + x*(c3   &
    + y*(c7 + c12*y) + x*(c6 + y*(c11 + y*(c17 + c24*y))  &
    + x*(c10  + y*(c16 + c23*y) + x*(c15 + c22*y + c21*x)))))  &
    + y*(c2 + y*(c5 + c8*x + y*(c9 + x*(c13 + c18*x)  &
    + y*(c14 + x*(c19 + c25*x) + y*(c20 + c26*x + c27*y)))))

return
end function

!--------------------------------------------------------------------------
