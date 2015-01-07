!==========================================================================
elemental function gsw_sa_from_rho (rho, ct, p)
!==========================================================================
!
!  Calculates the Absolute Salinity of a seawater sample, for given values
!  of its density, Conservative Temperature and sea pressure (in dbar). 
!  This function uses the computationally-efficient 48-term expression for 
!  density in terms of SA, CT and p (IOC et al., 2010).

!  rho =  density of a seawater sample (e.g. 1026 kg/m^3).       [ kg/m^3 ]
!   Note. This input has not had 1000 kg/m^3 subtracted from it. 
!     That is, it is 'density', not 'density anomaly'.
!  ct  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!
!  sa  =  Absolute Salinity                                          [g/kg]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_rho_alpha_beta, gsw_specvol

use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: rho, ct, p

real (r14) :: gsw_sa_from_rho

integer no_iter

real (r14) :: sa, v_lab, v_0, v_50, v_sa, sa_old, delta_v, sa_mean
real (r14) :: beta_mean, rho_mean

character (*), parameter :: func_name = "gsw_sa_from_rho"

v_lab = 1d0/rho
v_0 = gsw_specvol(0d0,ct,p)
v_50 = gsw_specvol(50d0,ct,p)

sa = 50d0*(v_lab - v_0)/(v_50 - v_0)

if (sa.lt.0d0.or.sa.gt.50d0) then
    gsw_sa_from_rho = gsw_error_code(1,func_name)
    return
end if

v_sa = (v_50 - v_0)/50d0

do no_iter = 1, 2 
    sa_old = sa
    delta_v = gsw_specvol(sa_old,ct,p) - v_lab
    sa = sa_old - delta_v/v_sa 
    sa_mean = 0.5d0*(sa + sa_old)
    call gsw_rho_alpha_beta(sa_mean,ct,p,rho=rho_mean,beta=beta_mean)
    v_sa = -beta_mean/rho_mean
    sa = sa_old - delta_v/v_sa
    if (sa.lt.0d0.or.sa.gt.50d0) then
        gsw_sa_from_rho = gsw_error_code(no_iter+1,func_name)
        return
    end if
end do

gsw_sa_from_rho = sa

return
end function

!--------------------------------------------------------------------------
