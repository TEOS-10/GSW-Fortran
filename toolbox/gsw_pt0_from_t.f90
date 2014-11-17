!==========================================================================
function gsw_pt0_from_t(sa,t,p) 
!==========================================================================
   
! Calculates potential temperature with reference pressure, p_ref = 0 dbar. 
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_pt0_from_t : potential temperature, p_ref = 0        [deg C]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer n0, n2, n, no_iter
real (r14) :: sa, t, p, s1, gsw_entropy_t_exact, gsw_gibbs_pt0_pt0, gsw_pt0_from_t
real (r14) :: pt0, pt0_old, de_dt, dentropy, dentropy_dt, sso, cp0
real (r14) :: gsw_entropy_part_zerop, true_entropy_part, pt0m, gsw_entropy_part

n0 = 0
n2 = 2

cp0 = 3991.86795711963d0
sso = 35.16504d0

s1 = sa*35d0/sso

pt0 = t + p*( 8.65483913395442d-6  - &
        s1 *  1.41636299744881d-6  - &
         p *  7.38286467135737d-9  + &
         t *(-8.38241357039698d-6  + &
        s1 *  2.83933368585534d-8  + &
         t *  1.77803965218656d-8  + &
         p *  1.71155619208233d-10))

dentropy_dt = cp0/((273.15d0 + pt0)*(1d0-0.05d0*(1d0 - sa/sso)))

true_entropy_part = gsw_entropy_part(sa,t,p)

do no_iter = 1,2
    pt0_old = pt0
    dentropy = gsw_entropy_part_zerop(sa,pt0_old) - true_entropy_part
    pt0 = pt0_old - dentropy/dentropy_dt 
    pt0m = 0.5d0*(pt0 + pt0_old)
    dentropy_dt = -gsw_gibbs_pt0_pt0(SA,pt0m)
    pt0 = pt0_old - dentropy/dentropy_dt
end do

gsw_pt0_from_t = pt0

end

!--------------------------------------------------------------------------

