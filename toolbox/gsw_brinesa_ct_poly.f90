!==========================================================================
elemental function gsw_brinesa_ct_poly (ct, p, saturation_fraction)
!==========================================================================
!
!  Calculates the Absolute Salinity of seawater at the freezing temperature.  
!  That is, the output is the Absolute Salinity of seawater, with the 
!  fraction saturation_fraction of dissolved air, that is in equilibrium 
!  with ice at Conservative Temperature CT and pressure p.  If the input 
!  values are such that there is no positive value of Absolute Salinity for
!  which seawater is frozen, the output, brineSA_CT, is put equal to Nan.
!
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!  saturation_fraction  =  the saturation fraction of dissolved air in 
!                          seawater
!
!  brineSA_CT  =  Absolute Salinity of seawater when it freezes, for 
!                 given input values of Conservative Temperature
!                 pressure and air saturation fraction.            [ g/kg ]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_ct_freezing_poly, gsw_sa_p_inrange
use gsw_mod_toolbox, only : gsw_ct_freezing_derivative_poly
use gsw_mod_toolbox, only : gsw_brinesa_estimate

use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: ct, p, saturation_fraction

real (r14) :: gsw_brinesa_ct_poly

integer :: i_iter
real (r14) :: ct_freezing, ct_freezing_zero_sa, dct_dsa
real (r14) :: sa, sa_cut_off, sa_old, sa_r, sa_mean

integer, parameter :: number_of_iterations = 2

character (*), parameter :: func_name = "gsw_brinesa_ct_poly"

! Find CT > CT_freezing_zero_SA.  If this is the case, the input values
! represent seawater that is not frozen (at any positive SA). 
ct_freezing_zero_sa = gsw_ct_freezing_poly(0d0,p,saturation_fraction)
if (ct .gt. ct_freezing_zero_sa) then
    gsw_brinesa_ct_poly = gsw_error_code(1,func_name)
    return
end if

! Form the first estimate of brineSA_CT from a polynomial in CT and p 
sa = gsw_brinesa_estimate(p,saturation_fraction,ct=ct)

! Find -SA_cut_off < SA < SA_cut_off, replace the above estimate of SA  
! with one based on (CT_freezing_zero_SA - CT).
sa_cut_off = 2.5d0 ! this is the band of sa within +- 2.5d0 g/kg of sa = 0, 
!                   which we treat differently in calculating the initial
!                   values of both SA and dCT_dSA. 
if (sa .lt. -sa_cut_off) then
    gsw_brinesa_ct_poly = gsw_error_code(2,func_name)
    return
end if

!--------------------------------------------------------------------------
! Form the first estimate of dCT_dSA, the derivative of CT with respect 
! to SA at fixed p.  
!--------------------------------------------------------------------------
sa = max(sa,0d0)
dct_dsa = gsw_ct_freezing_derivative_poly(sa,p,saturation_fraction)

if (abs(sa) .lt. sa_cut_off) sa = (ct - ct_freezing_zero_sa)/dct_dsa

!--------------------------------------------------------------------------
! Begin the modified Newton-Raphson method to solve the root of 
! CT_freezing = CT for SA. 
!--------------------------------------------------------------------------
do i_iter = 1, number_of_iterations
    sa_old = sa
    ct_freezing = gsw_ct_freezing_poly(sa_old,p,saturation_fraction)
    sa = sa_old - (ct_freezing - ct)/dct_dsa
    sa_mean = 0.5d0*(sa + sa_old)
    dct_dsa = gsw_ct_freezing_derivative_poly(sa_mean,p,saturation_fraction)
    sa = sa_old - (ct_freezing - ct)/dct_dsa
end do

if (gsw_sa_p_inrange(sa,p)) then
    gsw_brinesa_ct_poly = sa
else
    gsw_brinesa_ct_poly = gsw_error_code(3,func_name)
end if

return
end function

!--------------------------------------------------------------------------
