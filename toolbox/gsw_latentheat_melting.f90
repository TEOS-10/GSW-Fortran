!==========================================================================
function gsw_latentheat_melting(sa,p)  
!==========================================================================

! Calculates latent heat, or enthalpy, of melting.
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! 
! gsw_latentheat_melting : latent heat of melting          [kg/m^3]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: c0 =  3.334265169240710d5, c1 = -2.789444646733159d0
real (r14), parameter :: c2 = -1.822150156453350d4, c3 = -4.984585692734338d3
real (r14), parameter :: c4 = -7.371966528571920d1, c5 = -7.605802553358546d3
real (r14), parameter :: c6 =  1.195857305019339d3, c7 =  1.233720336206392d3
real (r14), parameter :: c8 =  2.294798676591890d2, c9 =  9.655751370889338d2
real (r14), parameter :: c10 = -5.792068522727968d2, c11 = -1.649446955902331d3
real (r14), parameter :: c12 = -1.029021448430547d3, c13 = -3.171558017172501d2
real (r14), parameter :: c14 = -1.751401389905041d2, c15 =  6.836527214265952d2
real (r14), parameter :: c16 =  1.078283734113611d3, c17 =  5.613896351265648d2
real (r14), parameter :: c18 =  6.968934948667265d2, c19 =  1.793032021946783d2
real (r14), parameter :: c20 =  8.692558481134256d1, c21 = -2.371103254714944d2
real (r14), parameter :: c22 = -5.775033277201674d2, c23 = -3.019749254648732d2
real (r14), parameter :: c24 = -6.420420579160927d2, c25 = -2.657570848596042d2
real (r14), parameter :: c26 = -1.646738151143109d1, c27 =  4.618228988300871d0

real (r14) :: sa, p, s_u, x, y, gsw_latentheat_melting

s_u = 40d0*(35.16504d0/35d0)
x = sqrt(sa/s_u)
y = p*1d-4

gsw_latentheat_melting = c0 + x*(c1 + c4*y + x*(c3   &
    + y*(c7 + c12*y) + x*(c6 + y*(c11 + y*(c17 + c24*y))  &
    + x*(c10  + y*(c16 + c23*y) + x*(c15 + c22*y + c21*x)))))  &
    + y*(c2 + y*(c5 + c8*x + y*(c9 + x*(c13 + c18*x)  &
    + y*(c14 + x*(c19 + c25*x) + y*(c20 + c26*x + c27*y)))))

return
end

!--------------------------------------------------------------------------

