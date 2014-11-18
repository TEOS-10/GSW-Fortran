!==========================================================================
module gsw_mod_gibbs_ice_coefficients
!==========================================================================

implicit none

integer, parameter :: gr14 = selected_real_kind(14,30)

complex (gr14), parameter :: t1 = (3.68017112855051d-2, 5.10878114959572d-2)
complex (gr14), parameter :: t2 = (3.37315741065416d-1, 3.35449415919309d-1)

complex (gr14), parameter :: r1 = (4.47050716285388d1, 6.56876847463481d1)
complex (gr14), parameter :: r20 = (-7.25974574329220d1, -7.81008427112870d1)
complex (gr14), parameter :: r21 = (-5.57107698030123d-5, 4.64578634580806d-5)
complex (gr14), parameter :: r22 = (2.34801409215913d-11, -2.85651142904972d-11)

! 1./Pt, where Pt = 611.657;  Experimental triple-point pressure in Pa.
real (gr14), parameter :: rec_pt = 1.634903221903779d-3
real (gr14), parameter :: tt = 273.16d0 ! Triple-point temperature, kelvin (K).
real (gr14), parameter :: rec_tt = 3.660858105139845d-3   ! = 1/tt

real (gr14), parameter :: g00 = -6.32020233335886d5
real (gr14), parameter :: g01 =  6.55022213658955d-1
real (gr14), parameter :: g02 = -1.89369929326131d-8
real (gr14), parameter :: g03 =  3.3974612327105304d-15
real (gr14), parameter :: g04 = -5.564648690589909d-22

end module

!--------------------------------------------------------------------------
