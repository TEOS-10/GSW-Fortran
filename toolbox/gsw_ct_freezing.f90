!==========================================================================
function gsw_ct_freezing(sa,p,saturation_fraction)
!==========================================================================

! Calculates the Conservative Temperature at which of seawater freezes 
! from Absolute Salinity and pressure.
!
! sa     : Absolute Salinity                                 [g/kg]
! p      : sea pressure                                      [dbar]
! saturation_fraction : saturation fraction
!
! gsw_ct_freezing : Conservative Temperature freezing point  [deg C]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: c0 = 0.017947064327968736d0, c1 = -6.076099099929818d0
real (r14), parameter :: c2 = 4.883198653547851d0, c3 = -11.88081601230542d0
real (r14), parameter :: c4 = 13.34658511480257d0, c5 = -8.722761043208607d0
real (r14), parameter :: c6 = 2.082038908808201d0, c7 = -7.389420998107497d0
real (r14), parameter :: c8 = -2.110913185058476d0, c9 = 0.2295491578006229d0
real (r14), parameter :: c10 = -0.9891538123307282d0, c11 = -0.08987150128406496d0
real (r14), parameter :: c12 = 0.3831132432071728d0, c13 = 1.054318231187074d0
real (r14), parameter :: c14 = 1.065556599652796d0, c15 = -0.7997496801694032d0
real (r14), parameter :: c16 = 0.3850133554097069d0, c17 = -2.078616693017569d0
real (r14), parameter :: c18 = 0.8756340772729538d0, c19 = -2.079022768390933d0
real (r14), parameter :: c20 = 1.596435439942262d0, c21 = 0.1338002171109174d0
real (r14), parameter :: c22 = 1.242891021876471d0

real (r14) :: sa, p, saturation_fraction, ct_freezing, gsw_ct_freezing
real (r14) :: sa_r, x, p_r, a, b

sa_r = sa*1d-2
x = sqrt(sa_r)
p_r = p*1d-4

ct_freezing = c0 &
 + sa_r*(c1 + x*(c2 + x*(c3 + x*(c4 + x*(c5 + c6*x))))) &
 + p_r*(c7 + p_r*(c8 + c9*p_r)) &
 + sa_r*p_r*(c10 + p_r*(c12 + p_r*(c15 + c21*sa_r)) + sa_r*(c13 + c17*p_r + c19*sa_r) &
 + x*(c11 + p_r*(c14 + c18*p_r)  + sa_r*(c16 + c20*p_r + c22*sa_r)))

! Adjust for the effects of dissolved air 
a = 0.014289763856964d0             ! Note that a = 0.502500117621/35.16504.
b = 0.057000649899720d0
ct_freezing = ct_freezing - saturation_fraction*(1d-3)*(2.4d0 - a*sa)*(1d0 + b*(1d0 - sa/35.16504d0))

if (p.gt.10000d0 .or. sa.gt.120d0 .or. (p+sa*71.428571428571402d0).gt.13571.42857142857d0) then
    ct_freezing = 9d15
end if

gsw_ct_freezing = ct_freezing

return
end function

!--------------------------------------------------------------------------
