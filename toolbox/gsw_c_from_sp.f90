!==========================================================================
elemental function gsw_c_from_sp (sp, t, p)       
!==========================================================================
!
!  Calculates conductivity, C, from (SP,t,p) using PSS-78 in the range 
!  2 < SP < 42.  If the input Practical Salinity is less than 2 then a 
!  modified form of the Hill et al. (1986) fomula is used for Practical 
!  Salinity.  The modification of the Hill et al. (1986) expression is to
!  ensure that it is exactly consistent with PSS-78 at SP = 2.
!
!  The conductivity ratio returned by this function is consistent with the
!  input value of Practical Salinity, SP, to 2x10^-14 psu over the full 
!  range of input parameters (from pure fresh water up to SP = 42 psu).  
!  This error of 2x10^-14 psu is machine precision at typical seawater 
!  salinities.  This accuracy is achieved by having four different 
!  polynomials for the starting value of Rtx (the square root of Rt) in 
!  four different ranges of SP, and by using one and a half iterations of 
!  a computationally efficient modified Newton-Raphson technique (McDougall 
!  and Wotherspoon, 2012) to find the root of the equation.  
!
!  Note that strictly speaking PSS-78 (Unesco, 1983) defines Practical
!  Salinity in terms of the conductivity ratio, R, without actually
!  specifying the value of C(35,15,0) (which we currently take to be
!  42.9140 mS/cm).
!
! sp     : Practical Salinity                               [unitless]
! t      : in-situ temperature [ITS-90]                     [deg C]
! p      : sea pressure                                     [dbar]
!
! c      : conductivity                                     [ mS/cm ]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_hill_ratio_at_sp2

use gsw_mod_teos10_constants, only : gsw_c3515

use gsw_mod_sp_coefficients

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sp, t, p       

real (r14) :: gsw_c_from_sp

real (r14) :: t68, ft68, x, rtx, dsp_drtx, sqrty
real (r14) :: part1, part2, hill_ratio, sp_hill_raw, sp_est
real (r14) :: rtx_old, rt, aa, bb, cc, dd, ee, ra,r, rt_lc, rtxm

real (r14), parameter :: p0 = 4.577801212923119d-3
real (r14), parameter :: p1 = 1.924049429136640d-1
real (r14), parameter :: p2 = 2.183871685127932d-5
real (r14), parameter :: p3 = -7.292156330457999d-3
real (r14), parameter :: p4 = 1.568129536470258d-4
real (r14), parameter :: p5 = -1.478995271680869d-6
real (r14), parameter :: p6 = 9.086442524716395d-4
real (r14), parameter :: p7 = -1.949560839540487d-5
real (r14), parameter :: p8 = -3.223058111118377d-6
real (r14), parameter :: p9 = 1.175871639741131d-7
real (r14), parameter :: p10 = -7.522895856600089d-5
real (r14), parameter :: p11 = -2.254458513439107d-6
real (r14), parameter :: p12 = 6.179992190192848d-7
real (r14), parameter :: p13 = 1.005054226996868d-8
real (r14), parameter :: p14 = -1.923745566122602d-9
real (r14), parameter :: p15 = 2.259550611212616d-6
real (r14), parameter :: p16 = 1.631749165091437d-7
real (r14), parameter :: p17 = -5.931857989915256d-9
real (r14), parameter :: p18 = -4.693392029005252d-9
real (r14), parameter :: p19 = 2.571854839274148d-10
real (r14), parameter :: p20 = 4.198786822861038d-12

real (r14), parameter :: q0 = 5.540896868127855d-5
real (r14), parameter :: q1 = 2.015419291097848d-1
real (r14), parameter :: q2 = -1.445310045430192d-5 
real (r14), parameter :: q3 = -1.567047628411722d-2
real (r14), parameter :: q4 = 2.464756294660119d-4
real (r14), parameter :: q5 = -2.575458304732166d-7
real (r14), parameter :: q6 = 5.071449842454419d-3
real (r14), parameter :: q7 = 9.081985795339206d-5
real (r14), parameter :: q8 = -3.635420818812898d-6
real (r14), parameter :: q9 = 2.249490528450555d-8
real (r14), parameter :: q10 = -1.143810377431888d-3
real (r14), parameter :: q11 = 2.066112484281530d-5
real (r14), parameter :: q12 = 7.482907137737503d-7
real (r14), parameter :: q13 = 4.019321577844724d-8
real (r14), parameter :: q14 = -5.755568141370501d-10
real (r14), parameter :: q15 = 1.120748754429459e-4
real (r14), parameter :: q16 = -2.420274029674485d-6
real (r14), parameter :: q17 = -4.774829347564670d-8
real (r14), parameter :: q18 = -4.279037686797859d-9
real (r14), parameter :: q19 = -2.045829202713288d-10
real (r14), parameter :: q20 = 5.025109163112005d-12

real (r14), parameter :: s0 = 3.432285006604888d-3
real (r14), parameter :: s1 = 1.672940491817403d-1
real (r14), parameter :: s2 = 2.640304401023995d-5
real (r14), parameter :: s3 = 1.082267090441036d-1
real (r14), parameter :: s4 = -6.296778883666940d-5
real (r14), parameter :: s5 = -4.542775152303671d-7
real (r14), parameter :: s6 = -1.859711038699727d-1
real (r14), parameter :: s7 = 7.659006320303959d-4
real (r14), parameter :: s8 = -4.794661268817618d-7
real (r14), parameter :: s9 = 8.093368602891911d-9
real (r14), parameter :: s10 = 1.001140606840692d-1 
real (r14), parameter :: s11 = -1.038712945546608d-3
real (r14), parameter :: s12 = -6.227915160991074d-6
real (r14), parameter :: s13 = 2.798564479737090d-8
real (r14), parameter :: s14 = -1.343623657549961d-10
real (r14), parameter :: s15 = 1.024345179842964d-2
real (r14), parameter :: s16 = 4.981135430579384d-4
real (r14), parameter :: s17 = 4.466087528793912d-6
real (r14), parameter :: s18 = 1.960872795577774d-8
real (r14), parameter :: s19 = -2.723159418888634d-10
real (r14), parameter :: s20 = 1.122200786423241d-12

real (r14), parameter :: u0 = 5.180529787390576d-3
real (r14), parameter :: u1 = 1.052097167201052d-3
real (r14), parameter :: u2 = 3.666193708310848d-5
real (r14), parameter :: u3 = 7.112223828976632d0
real (r14), parameter :: u4 = -3.631366777096209d-4
real (r14), parameter :: u5 = -7.336295318742821d-7
real (r14), parameter :: u6 = -1.576886793288888d+2
real (r14), parameter :: u7 = -1.840239113483083d-3
real (r14), parameter :: u8 = 8.624279120240952d-6
real (r14), parameter :: u9 = 1.233529799729501d-8
real (r14), parameter :: u10 = 1.826482800939545d+3
real (r14), parameter :: u11 = 1.633903983457674d-1
real (r14), parameter :: u12 = -9.201096427222349d-5
real (r14), parameter :: u13 = -9.187900959754842d-8
real (r14), parameter :: u14 = -1.442010369809705d-10
real (r14), parameter :: u15 = -8.542357182595853d+3
real (r14), parameter :: u16 = -1.408635241899082d0
real (r14), parameter :: u17 = 1.660164829963661d-4
real (r14), parameter :: u18 = 6.797409608973845d-7
real (r14), parameter :: u19 = 3.345074990451475d-10
real (r14), parameter :: u20 = 8.285687652694768d-13

t68 = t*1.00024d0
ft68 = (t68 - 15d0)/(1d0 + k*(t68 - 15d0))

x = sqrt(sp)

!--------------------------------------------------------------------------
! Finding the starting value of Rtx, the square root of Rt, using four 
! different polynomials of SP and t68.  
!--------------------------------------------------------------------------

if (sp.ge.9d0) then

    rtx = p0 + x*(p1 + p4*t68 + x*(p3 + p7*t68 + x*(p6  &
        + p11*t68 + x*(p10 + p16*t68 + x*p15))))  &
        + t68*(p2+ t68*(p5 + x*x*(p12 + x*p17) + p8*x  &
        + t68*(p9 + x*(p13 + x*p18)+ t68*(p14 + p19*x + p20*t68))))

else if (sp.ge.0.25d0.and.sp.lt.9d0) then

    rtx = q0 + x*(q1 + q4*t68 + x*(q3 + q7*t68 + x*(q6  &
        + q11*t68 + x*(q10 + q16*t68 + x*q15))))  &
        + t68*(q2+ t68*(q5 + x*x*(q12 + x*q17) + q8*x  &
        + t68*(q9 + x*(q13 + x*q18)+ t68*(q14 + q19*x + q20*t68))))

else if (sp.ge.0.003d0.and.sp.lt.0.25d0) then

    rtx = s0 + x*(s1 + s4*t68 + x*(s3 + s7*t68 + x*(s6  &
        + s11*t68 + x*(s10 + s16*t68 + x*s15))))  &
        + t68*(s2+ t68*(s5 + x*x*(s12 + x*s17) + s8*x  &
        + t68*(s9 + x*(s13 + x*s18)+ t68*(s14 + s19*x + s20*t68))))

else if (sp.lt.0.003d0) then

    rtx = u0 + x*(u1 + u4*t68 + x*(u3 + u7*t68 + x*(u6  &
        + u11*t68 + x*(u10 + u16*t68 + x*u15))))  &
        + t68*(u2+ t68*(u5 + x*x*(u12 + x*u17) + u8*x  &
        + t68*(u9 + x*(u13 + x*u18)+ t68*(u14 + u19*x + u20*t68))))

end if

!--------------------------------------------------------------------------
! Finding the starting value of dSP_dRtx, the derivative of SP with respect
! to Rtx.  
!--------------------------------------------------------------------------
dsp_drtx =  a1 + (2d0*a2 + (3d0*a3 + (4d0*a4 + 5d0*a5*rtx)*rtx)*rtx)*rtx  &
    + ft68*(b1 + (2d0*b2 + (3d0*b3 + (4d0*b4 + 5d0*b5*rtx)*rtx)*rtx)*rtx)

if (sp.lt.2d0) then
    x = 400d0*(rtx*rtx)
    sqrty = 10*rtx
    part1 = 1d0 + x*(1.5d0 + x) 
    part2 = 1d0 + sqrty*(1d0 + sqrty*(1d0 + sqrty))
    hill_ratio = gsw_hill_ratio_at_sp2(t)
    dsp_drtx = dsp_drtx  &
        + a0*800d0*Rtx*(1.5d0 + 2d0*x)/(part1*part1)  &
        + b0*ft68*(10d0 + sqrty*(20d0 + 30d0*sqrty))/(part2*part2)
    dsp_drtx = hill_ratio*dsp_drtx
end if

!--------------------------------------------------------------------------
! One iteration through the modified Newton-Raphson method (McDougall and 
! Wotherspoon, 2012) achieves an error in Practical Salinity of about 
! 10^-12 for all combinations of the inputs.  One and a half iterations of 
! the modified Newton-Raphson method achevies a maximum error in terms of 
! Practical Salinity of better than 2x10^-14 everywhere. 
!
! We recommend one and a half iterations of the modified Newton-Raphson
! method. 
!
! Begin the modified Newton-Raphson method.  
!--------------------------------------------------------------------------
    sp_est = a0 + (a1 + (a2 + (a3 + (a4 + a5*rtx)*rtx)*rtx)*rtx)*rtx &
        + ft68*(b0 + (b1 + (b2+ (b3 + (b4 + b5*rtx)*rtx)*rtx)*rtx)*rtx)
    if (sp_est .lt. 2) then
        x = 400d0*(rtx*rtx)
        sqrty = 10d0*rtx
        part1 = 1d0 + x*(1.5d0 + x) 
        part2 = 1d0 + sqrty*(1d0 + sqrty*(1d0 + sqrty))
        sp_hill_raw = sp_est - a0/part1 - b0*ft68/part2
        hill_ratio = gsw_hill_ratio_at_sp2(t)
        sp_est = hill_ratio*sp_hill_raw
    end if
 
    rtx_old = rtx
    rtx = rtx_old - (sp_est - sp)/dsp_drtx
    
    rtxm = 0.5d0*(rtx + rtx_old)      ! This mean value of Rtx, Rtxm, is the  
!                 value of Rtx at which the derivative dSP_dRtx is evaluated.
    
    dsp_drtx = a1 + (2d0*a2 + (3d0*a3 + (4d0*a4 + 5d0*a5*rtxm)*rtxm)*rtxm)*rtxm&
       + ft68*(b1 + (2d0*b2 + (3d0*b3 + (4d0*b4 + 5d0*b5*rtxm)*rtxm)*rtxm)*rtxm)
    if (sp_est .lt. 2) then
        x = 400d0*(rtxm*rtxm)
        sqrty = 10d0*rtxm
        part1 = 1d0 + x*(1.5d0 + x) 
        part2 = 1d0 + sqrty*(1d0 + sqrty*(1d0 + sqrty))
        dsp_drtx = dsp_drtx  &
            + a0*800d0*rtxm*(1.5d0 + 2d0*x)/(part1*part1)  &
            + b0*ft68*(10d0 + sqrty*(20d0 + 30d0*sqrty))/(part2*part2)
        hill_ratio = gsw_hill_ratio_at_sp2(t)
        dsp_drtx = hill_ratio*dsp_drtx
    end if

!--------------------------------------------------------------------------
! The line below is where Rtx is updated at the end of the one full 
! iteration of the modified Newton-Raphson technique.
!--------------------------------------------------------------------------
    rtx = rtx_old - (sp_est - sp)/dsp_drtx
!--------------------------------------------------------------------------
! Now we do another half iteration of the modified Newton-Raphson  
! technique, making a total of one and a half modified N-R iterations.
!-------------------------------------------------------------------------- 
    sp_est = a0 + (a1 + (a2 + (a3 + (a4 + a5*rtx)*rtx)*rtx)*rtx)*rtx  &
        + ft68*(b0 + (b1 + (b2+ (b3 + (b4 + b5*rtx)*rtx)*rtx)*rtx)*rtx)
    if (sp_est .lt. 2) then
        x = 400d0*(rtx*rtx)
        sqrty = 10d0*rtx
        part1 = 1d0 + x*(1.5d0 + x) 
        part2 = 1d0 + sqrty*(1d0 + sqrty*(1d0 + sqrty))
        sp_hill_raw = sp_est - a0/part1 - b0*ft68/part2
        hill_ratio = gsw_hill_ratio_at_sp2(t)
        sp_est = hill_ratio*sp_hill_raw
    end if
    rtx = rtx - (sp_est - sp)/dsp_drtx

!--------------------------------------------------------------------------
! Now go from Rtx to Rt and then to the conductivity ratio R at pressure p.
!--------------------------------------------------------------------------
rt = rtx*rtx
aa  = d3 + d4*t68
bb  = 1d0 + t68*(d1 + d2*t68)
cc  = p*(e1 + p*(e2 + e3*p))
! rt_lc (i.e. rt_lower_case) corresponds to rt as defined in 
! the UNESCO 44 (1983) routines.
rt_lc = c0 + (c1 + (c2 + (c3 + c4*t68)*t68)*t68)*t68

dd  = bb - aa*rt_lc*rt
ee  = rt_lc*rt*aa*(bb + cc)
ra = sqrt(dd*dd + 4d0*ee) - dd
r  = 0.5d0*ra/aa

! The dimensionless conductivity ratio, R, is the conductivity input, C,
! divided by the present estimate of C(SP=35, t_68=15, p=0) which is 
! 42.9140 mS/cm (=4.29140 S/m^). 

gsw_c_from_sp = gsw_c3515*r      

return
end function

!--------------------------------------------------------------------------
