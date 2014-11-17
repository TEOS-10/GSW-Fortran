!==========================================================================
function gsw_pt_from_t(sa,t,p,p_ref) 
!==========================================================================
   
! Calculates potential temperature of seawater from in-situ temperature 
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! p_ref  : reference sea pressure                          [dbar]
!
! gsw_pt_from_t : potential temperature                    [deg C]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer n0, n2, n, no_iter
real (r14) :: sa, t, p, p_ref, s1, gsw_entropy_t_exact, gsw_gibbs, gsw_pt_from_t
real (r14) :: pt, pt_old, de_dt, dentropy, dentropy_dt, sso, cp0
real (r14) :: gsw_entropy_part, true_entropy_part, ptm

n0 = 0
n2 = 2

cp0 = 3991.86795711963d0
sso = 35.16504d0

s1 = sa*35d0/sso

pt = t + (p-p_ref)*( 8.65483913395442d-6  - &
               s1 *  1.41636299744881d-6  - &
        (p+p_ref) *  7.38286467135737d-9  + &
               t  *(-8.38241357039698d-6  + &
               s1 *  2.83933368585534d-8  + &
               t  *  1.77803965218656d-8  + &
        (p+p_ref) *  1.71155619208233d-10))

dentropy_dt = cp0/((273.15d0 + pt)*(1d0-0.05d0*(1d0 - sa/sso)));

true_entropy_part = gsw_entropy_part(sa,t,p)

do no_iter = 1,2
    pt_old = pt
    dentropy = gsw_entropy_part(sa,pt_old,p_ref) - true_entropy_part
    pt = pt_old - dentropy/dentropy_dt 
    ptm = 0.5d0*(pt + pt_old)
    dentropy_dt = -gsw_gibbs(n0,n2,n0,sa,ptm,p_ref)
    pt = pt_old - dentropy/dentropy_dt
end do

gsw_pt_from_t = pt

end

!--------------------------------------------------------------------------

