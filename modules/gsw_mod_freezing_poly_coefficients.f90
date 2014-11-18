!==========================================================================
module gsw_mod_freezing_poly_coefficients
!==========================================================================

implicit none
integer, parameter :: fpc_r14 = selected_real_kind(14,30)

real (fpc_r14), parameter :: c0  =  0.017947064327968736d0
real (fpc_r14), parameter :: c1 =  -6.076099099929818d0
real (fpc_r14), parameter :: c2 =   4.883198653547851d0
real (fpc_r14), parameter :: c3 =  -11.88081601230542d0
real (fpc_r14), parameter :: c4 =   13.34658511480257d0
real (fpc_r14), parameter :: c5 =  -8.722761043208607d0
real (fpc_r14), parameter :: c6 =   2.082038908808201d0
real (fpc_r14), parameter :: c7 =  -7.389420998107497d0
real (fpc_r14), parameter :: c8 =  -2.110913185058476d0
real (fpc_r14), parameter :: c9 =   0.2295491578006229d0 
real (fpc_r14), parameter :: c10 = -0.9891538123307282d0
real (fpc_r14), parameter :: c11 = -0.08987150128406496d0
real (fpc_r14), parameter :: c12 =  0.3831132432071728d0
real (fpc_r14), parameter :: c13 =  1.054318231187074d0
real (fpc_r14), parameter :: c14 =  1.065556599652796d0
real (fpc_r14), parameter :: c15 = -0.7997496801694032d0
real (fpc_r14), parameter :: c16 =  0.3850133554097069d0
real (fpc_r14), parameter :: c17 = -2.078616693017569d0
real (fpc_r14), parameter :: c18 =  0.8756340772729538d0
real (fpc_r14), parameter :: c19 = -2.079022768390933d0
real (fpc_r14), parameter :: c20 =  1.596435439942262d0
real (fpc_r14), parameter :: c21 =  0.1338002171109174d0
real (fpc_r14), parameter :: c22 =  1.242891021876471d0

! Note that a = 0.502500117621d0/gsw_sso
real (fpc_r14), parameter :: a = 0.014289763856964d0
real (fpc_r14), parameter :: b = 0.057000649899720d0

real (fpc_r14), parameter :: t0 = 0.002519d0
real (fpc_r14), parameter :: t1 = -5.946302841607319d0
real (fpc_r14), parameter :: t2 =  4.136051661346983d0
real (fpc_r14), parameter :: t3 = -1.115150523403847d1
real (fpc_r14), parameter :: t4 =  1.476878746184548d1
real (fpc_r14), parameter :: t5 = -1.088873263630961d1
real (fpc_r14), parameter :: t6 =  2.961018839640730d0
real (fpc_r14), parameter :: t7 = -7.433320943962606d0
real (fpc_r14), parameter :: t8 = -1.561578562479883d0
real (fpc_r14), parameter :: t9 =  4.073774363480365d-2
real (fpc_r14), parameter :: t10 =  1.158414435887717d-2
real (fpc_r14), parameter :: t11 = -4.122639292422863d-1
real (fpc_r14), parameter :: t12 = -1.123186915628260d-1
real (fpc_r14), parameter :: t13 =  5.715012685553502d-1
real (fpc_r14), parameter :: t14 =  2.021682115652684d-1
real (fpc_r14), parameter :: t15 =  4.140574258089767d-2
real (fpc_r14), parameter :: t16 = -6.034228641903586d-1
real (fpc_r14), parameter :: t17 = -1.205825928146808d-2
real (fpc_r14), parameter :: t18 = -2.812172968619369d-1
real (fpc_r14), parameter :: t19 =  1.877244474023750d-2
real (fpc_r14), parameter :: t20 = -1.204395563789007d-1
real (fpc_r14), parameter :: t21 =  2.349147739749606d-1
real (fpc_r14), parameter :: t22 =  2.748444541144219d-3

end module

!--------------------------------------------------------------------------
