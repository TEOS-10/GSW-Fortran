!==========================================================================
subroutine gsw_turner_rsubrho (sa, ct, p, nz, tu, rsubrho, p_mid)
!==========================================================================
!
!  Calculates the Turner angle and the Rsubrho as a function of pressure 
!  down a vertical water column.  These quantities express the relative 
!  contributions of the vertical gradients of Conservative Temperature 
!  and Absolute Salinity to the vertical stability (the square of the 
!  Brunt-Vaisala Frequency squared, N^2).  Tu and Rsubrho are evaluated at 
!  the mid pressure between the individual data points in the vertical.  
!  This function uses computationally-efficient 48-term expression for 
!  density in terms of SA, CT and p (IOC et al., 2010).  Note that 
!  in the double-diffusive literature, papers concerned with the 
!  "diffusive" form of double-diffusive convection often define the 
!  stability ratio as the reciprocal of what is defined here as the 
!  stability ratio.  
!
!  Note. The 48-term equation has been fitted in a restricted range of parameter
!  space, and is most accurate inside the "oceanographic funnel" described 
!  in IOC et al. (2010).  The GSW library function "gsw_infunnel(SA,CT,p)" 
!  is available to be used if one wants to test if some of one's data lies
!  outside this "funnel".  

! sa      : Absolute Salinity         (a profile (length nz))     [g/kg]
! ct      : Conservative Temperature  (a profile (length nz))     [deg C]
! p       : sea pressure              (a profile (length nz))     [dbar]                            
! nz      : number of bottles                             
! tu      : Turner angle, on the same (nz-1) grid as p_mid.
!           Turner angle has units of:           [ degrees of rotation ]
! rsubrho : Stability Ratio, on the same (nz-1) grid as p_mid.
!           Rsubrho is dimensionless.                       [ unitless ]
! p_mid   : Mid pressure between p grid  (length nz-1)           [dbar]
!--------------------------------------------------------------------------

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)
integer :: nz, k

real (r14), parameter :: pi = 3.141592653589793d0
real (r14) :: dsa, sa_mid, dct, ct_mid, dp
real (r14) :: alpha_mid, gsw_alpha, beta_mid, gsw_beta
real (r14), dimension(nz) :: sa, ct, p
real (r14), dimension(nz-1) :: tu, rsubrho, p_mid

do k = 1, nz-1
	dsa = (sa(k) - sa(k+1))
	sa_mid = 0.5d0*(sa(k) + sa(k+1))
	dct = (ct(k) - ct(k+1))
	ct_mid = 0.5d0*(ct(k) + ct(k+1))
	dp = (p(k) - p(k+1))
	p_mid(k) = 0.5d0*(p(k) + p(k+1))

	alpha_mid = gsw_alpha(sa_mid,ct_mid,p_mid(k))
	beta_mid = gsw_beta(sa_mid,ct_mid,p_mid(k))

	tu(k) = (180d0/pi)*atan2((alpha_mid*dct + beta_mid*dsa), &
	                         (alpha_mid*dct - beta_mid*dsa))

	if (dsa.eq.0d0) then
		rsubrho(k) = 9d15
	else 
		rsubrho(k) = (alpha_mid*dct)/(beta_mid*dsa)
	end if
end do

return
end subroutine

!--------------------------------------------------------------------------
