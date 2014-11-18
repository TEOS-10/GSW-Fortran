!==========================================================================
elemental function gsw_brinesa_estimate (p, saturation_fraction, ct, t)
!==========================================================================
!
! Form an estimate of brineSA_t from a polynomial in CT and p 
!
!--------------------------------------------------------------------------

use gsw_mod_teos10_constants, only : gsw_sso

use gsw_mod_toolbox, only : gsw_ct_from_t

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: p, saturation_fraction
real (r14), intent(in), optional :: ct, t

real (r14) :: gsw_brinesa_estimate

real (r14) :: ctx, ctsat, sa

! note that aa = 0.502500117621d0/35.16504
real (r14), parameter :: aa = 0.014289763856964d0
real (r14), parameter :: bb = 0.057000649899720d0

real (r14), parameter :: p0  =  2.570124672768757d-1
real (r14), parameter :: p1  = -1.917742353032266d1
real (r14), parameter :: p2  = -1.413382858617969d-2
real (r14), parameter :: p3  = -5.427484830917552d-1
real (r14), parameter :: p4  = -4.126621135193472d-4
real (r14), parameter :: p5  = -4.176407833276121d-7
real (r14), parameter :: p6  =  4.688217641883641d-5
real (r14), parameter :: p7  = -3.039808885885726d-8
real (r14), parameter :: p8  = -4.990118091261456d-11
real (r14), parameter :: p9  = -9.733920711119464d-9
real (r14), parameter :: p10 = -7.723324202726337d-12
real (r14), parameter :: p11 =  7.121854166249257d-16
real (r14), parameter :: p12 =  1.256474634100811d-12
real (r14), parameter :: p13 =  2.105103897918125d-15
real (r14), parameter :: p14 =  8.663811778227171d-19

! a rough estimate to get the saturated ct
if (present(ct)) then
    sa = max(-(ct + 9d-4*p)/0.06d0, 0d0)
    ctx = ct
else if (present(t)) then
    sa = max(-(t + 9d-4*p)/0.06d0, 0d0)
    ctx = gsw_ct_from_t(sa,t,p)
else
    gsw_brinesa_estimate = 0d0
    return
end if

! CTsat is the estimated value of CT if the seawater were saturated with
! dissolved air, recognizing that it actually has the air fraction
! saturation_fraction; see McDougall, Barker and Feistel, 2014).  

ctsat = ctx - &
        (1d0-saturation_fraction)*(1d-3)*(2.4d0-aa*sa)*(1d0+bb*(1d0-sa/gsw_sso))

gsw_brinesa_estimate = p0 &
        + p*(p2 + p4*ctsat + p*(p5 + ctsat*(p7 + p9*ctsat) &
        + p*(p8  + ctsat*(p10 + p12*ctsat) + p*(p11 + p13*ctsat + p14*p)))) &
        + ctsat*(p1 + ctsat*(p3 + p6*p))

return
end function

!--------------------------------------------------------------------------
