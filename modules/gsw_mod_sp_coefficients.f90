!==========================================================================
module gsw_mod_sp_coefficients
!==========================================================================
!
!--------------------------------------------------------------------------

implicit none
integer, parameter :: gsc_r14 = selected_real_kind(14,30)

real (gsc_r14), parameter :: a0 =  0.0080d0
real (gsc_r14), parameter :: a1 = -0.1692d0
real (gsc_r14), parameter :: a2 = 25.3851d0
real (gsc_r14), parameter :: a3 = 14.0941d0
real (gsc_r14), parameter :: a4 = -7.0261d0
real (gsc_r14), parameter :: a5 =  2.7081d0

real (gsc_r14), parameter :: b0 =  0.0005d0
real (gsc_r14), parameter :: b1 = -0.0056d0
real (gsc_r14), parameter :: b2 = -0.0066d0
real (gsc_r14), parameter :: b3 = -0.0375d0
real (gsc_r14), parameter :: b4 =  0.0636d0
real (gsc_r14), parameter :: b5 = -0.0144d0

real (gsc_r14), parameter :: c0 =  0.6766097d0
real (gsc_r14), parameter :: c1 =  2.00564d-2
real (gsc_r14), parameter :: c2 =  1.104259d-4
real (gsc_r14), parameter :: c3 = -6.9698d-7
real (gsc_r14), parameter :: c4 =  1.0031d-9

real (gsc_r14), parameter :: d1 =  3.426d-2
real (gsc_r14), parameter :: d2 =  4.464d-4
real (gsc_r14), parameter :: d3 =  4.215d-1
real (gsc_r14), parameter :: d4 = -3.107d-3

real (gsc_r14), parameter :: e1 =  2.070d-5
real (gsc_r14), parameter :: e2 = -6.370d-10
real (gsc_r14), parameter :: e3 =  3.989d-15

real (gsc_r14), parameter :: k  =  0.0162d0

end module

!--------------------------------------------------------------------------
