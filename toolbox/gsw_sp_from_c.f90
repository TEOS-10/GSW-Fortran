!==========================================================================
function gsw_sp_from_c(c,t,p)       
!==========================================================================

!  Calculates Practical Salinity, SP, from conductivity, C, primarily using
!  the PSS-78 algorithm.  Note that the PSS-78 algorithm for Practical 
!  Salinity is only valid in the range 2 < SP < 42.  If the PSS-78 
!  algorithm produces a Practical Salinity that is less than 2 then the 
!  Practical Salinity is recalculated with a modified form of the Hill et 
!  al. (1986) formula.  The modification of the Hill et al. (1986)
!  expression is to ensure that it is exactly consistent with PSS-78 
!  at SP = 2.  Note that the input values of conductivity need to be in 
!  units of mS/cm (not S/m). 
!
! c      : conductivity                                     [ mS/cm ]
! t      : in-situ temperature [ITS-90]                     [deg C]
! p      : sea pressure                                     [dbar]
!
! sp     : Practical Salinity                               [unitless]

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: a0 = 0.0080d0, a1 = -0.1692d0, a2 = 25.3851d0
real (r14), parameter :: a3 = 14.0941d0, a4 = -7.0261d0, a5 = 2.7081d0
real (r14), parameter :: b0 = 0.0005d0, b1 = -0.0056d0, b2 = -0.0066d0
real (r14), parameter :: b3 = -0.0375d0, b4 = 0.0636d0, b5 = -0.0144d0
real (r14), parameter :: c0 = 0.6766097d0, c1 = 2.00564d-2
real (r14), parameter :: c2 = 1.104259d-4, c3 = -6.9698d-7, c4 = 1.0031d-9
real (r14), parameter :: d1 = 3.426d-2, d2 = 4.464d-4, d3 =  4.215d-1
real (r14), parameter :: d4 = -3.107d-3, e1 = 2.070d-5, e2 = -6.370d-10
real (r14), parameter :: e3 = 3.989d-15, k  = 0.0162

real (r14) :: c, t, p, gsw_sp_from_c, sp, t68, ft68, r, rt_lc, rp, rt, rtx
real (r14) :: hill_ratio, gsw_hill_ratio_at_sp2, x, sqrty, part1, part2
real (r14) :: sp_hill_raw

t68 = t*1.00024d0
ft68 = (t68 - 15d0)/(1d0 + k*(t68 - 15d0))

! The dimensionless conductivity ratio, R, is the conductivity input, C,
! divided by the present estimate of C(SP=35, t_68=15, p=0) which is 
! 42.9140 mS/cm (=4.29140 S/m), (Culkin and Smith, 1980). 

r = 0.023302418791070513d0*c          !   0.023302418791070513 = 1./42.9140

! rt_lc corresponds to rt as defined in the UNESCO 44 (1983) routines.  
rt_lc = c0 + (c1 + (c2 + (c3 + c4*t68)*t68)*t68)*t68
rp = 1d0 + (p*(e1 + e2*p + e3*p*p))/(1d0 + d1*t68 + d2*t68*t68 + (d3 + d4*t68)*r)
rt = r/(rp*rt_lc)  

if (rt.lt.0) then
  rt = 9d15
endif

rtx = sqrt(rt)

sp = a0 + (a1 + (a2 + (a3 + (a4 + a5*rtx)*rtx)*rtx)*rtx)*rtx + &
    ft68*(b0 + (b1 + (b2 + (b3 + (b4 + b5*rtx)*rtx)*rtx)*rtx)*rtx)

! The following section of the code is designed for SP < 2 based on the
! Hill et al. (1986) algorithm.  This algorithm is adjusted so that it is
! exactly equal to the PSS-78 algorithm at SP = 2.

if (sp < 2) then
    hill_ratio = gsw_hill_ratio_at_sp2(t);
    x = 400d0*rt
    sqrty = 10d0*rtx
    part1 = 1d0 + x*(1.5d0 + x)
    part2 = 1d0 + sqrty*(1d0 + sqrty*(1d0 + sqrty))
    sp_hill_raw = sp - a0/part1 - b0*ft68/part2
    sp = hill_ratio*sp_hill_raw
endif

! This line ensures that SP is non-negative.
if (sp.lt.0) then
   sp = 9d15
end if

gsw_sp_from_c = sp

return
end function

!--------------------------------------------------------------------------

