!==========================================================================
function  gsw_hill_ratio_at_sp2(t)
!==========================================================================

!  Calculates the Hill ratio, which is the adjustment needed to apply for
!  Practical Salinities smaller than 2.  This ratio is defined at a 
!  Practical Salinity = 2 and in-situ temperature, t using PSS-78. The Hill
!  ratio is the ratio of 2 to the output of the Hill et al. (1986) formula
!  for Practical Salinity at the conductivity ratio, Rt, at which Practical
!  Salinity on the PSS-78 scale is exactly 2.

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), parameter :: a0 =  0.0080d0, a1 = -0.1692d0, a2 = 25.3851d0
real (r14), parameter :: a3 = 14.0941d0, a4 = -7.0261d0, a5 =  2.7081d0
real (r14), parameter :: b0 =  0.0005d0, b1 = -0.0056d0, b2 = -0.0066d0
real (r14), parameter :: b3 = -0.0375d0, b4 =  0.0636d0, b5 = -0.0144d0
real (r14), parameter :: g0 = 2.641463563366498d-1, g1 = 2.007883247811176d-4
real (r14), parameter :: g2 = -4.107694432853053d-6, g3 = 8.401670882091225d-8
real (r14), parameter :: g4 = -1.711392021989210d-9, g5 = 3.374193893377380d-11
real (r14), parameter :: g6 = -5.923731174730784d-13, g7 = 8.057771569962299d-15
real (r14), parameter :: g8 = -7.054313817447962d-17, g9 = 2.859992717347235d-19
real (r14), parameter :: rk  =  0.0162d0, sp2 = 2d0

real (r14) :: t, t68, ft68, rtx0, dsp_drtx, sp_est, rtx, rtxm, x, part1, part2
real (r14) :: sqrty, sp_hill_raw_at_sp2, gsw_hill_ratio_at_sp2

t68 = t*1.00024d0
ft68 = (t68 - 15d0)/(1d0 + rk*(t68 - 15d0))

!--------------------------------------------------------------------------
! Find the initial estimates of Rtx (Rtx0) and of the derivative dSP_dRtx
! at SP = 2. 
!--------------------------------------------------------------------------
rtx0 = g0 + t68*(g1 + t68*(g2 + t68*(g3 + t68*(g4 + t68*(g5 &
         + t68*(g6 + t68*(g7 + t68*(g8 + t68*g9))))))))
     
dsp_drtx =  a1 + (2*a2 + (3*a3 + (4*a4 + 5*a5*rtx0)*rtx0)*rtx0)*rtx0 + &
    ft68*(b1 + (2*b2 + (3*b3 + (4*b4 + 5*b5*rtx0)*rtx0)*rtx0)*rtx0)    

!--------------------------------------------------------------------------
! Begin a single modified Newton-Raphson iteration to find Rt at SP = 2.
!--------------------------------------------------------------------------
sp_est = a0 + (a1 + (a2 + (a3 + (a4 + a5*rtx0)*rtx0)*rtx0)*rtx0)*rtx0 &
        + ft68*(b0 + (b1 + (b2+ (b3 + (b4 + b5*rtx0)*rtx0)*rtx0)*rtx0)*rtx0)
rtx = rtx0 - (sp_est - sp2)/dsp_drtx
rtxm = 0.5d0*(rtx + rtx0)
dsp_drtx =  a1 + (2*a2 + (3*a3 + (4*a4 + 5*a5*rtxm)*rtxm)*rtxm)*rtxm &
        + ft68*(b1 + (2*b2 + (3*b3 + (4*b4 + 5*b5*rtxm)*rtxm)*rtxm)*rtxm)
rtx = rtx0 - (sp_est - sp2)/dsp_drtx

! This is the end of one full iteration of the modified Newton-Raphson 
! iterative equation solver.  The error in Rtx at this point is equivalent 
! to an error in SP of 9e-16 psu.  
                                
x = 400d0*rtx*rtx
sqrty = 10d0*rtx
part1 = 1 + x*(1.5d0 + x) 
part2 = 1 + sqrty*(1 + sqrty*(1 + sqrty))
sp_hill_raw_at_sp2 = SP2 - a0/part1 - b0*ft68/part2

gsw_hill_ratio_at_sp2 = 2/sp_hill_raw_at_sp2

return
end function

!==========================================================================
