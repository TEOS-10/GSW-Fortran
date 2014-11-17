!==========================================================================
subroutine gsw_ipv_vs_fnsquared_ratio (sa, ct, p, nz, ipv_vs_fnsquared_ratio, p_mid)
!==========================================================================
!
!  Calculates the ratio of the vertical gradient of potential density to 
!  the vertical gradient of locally-referenced potential density.  This 
!  ratio is also the ratio of the planetary Isopycnal Potential Vorticity
!  (IPV) to f times N^2, hence the name for this variable,
!  IPV_vs_fNsquared_ratio (see Eqn. (3.20.5) of IOC et al. (2010)). 
!  The reference sea pressure, p_ref, of the potential density surface must
!  have a constant value.
!
!  IPV_vs_fNsquared_ratio is evaluated at the mid pressure between the 
!  individual data points in the vertical.  This function uses the 
!  computationally-efficient 48-term expression for density in terms of 
!  SA, CT and p (IOC et al., 2010). 
!  Note. The 48-term equation has been fitted in a restricted range of parameter
!  space, and is most accurate inside the "oceanographic funnel" described 
!  in IOC et al. (2010).  The GSW library function "gsw_infunnel(SA,CT,p)" 
!  is avaialble to be used if one wants to test if some of one's data lies
!  outside this "funnel".  

! sa      : Absolute Salinity         (a profile (length nz))     [g/kg]
! ct      : Conservative Temperature  (a profile (length nz))     [deg C]
! p       : sea pressure              (a profile (length nz))     [dbar]
! nz      : number of bottles                             
! IPV_vs_fNsquared_ratio
!         : The ratio of the vertical gradient of potential density
!           referenced to p_ref, to the vertical gradient of locally-
!           referenced potential density.  It is ouput on the same
!           vertical (M-1)xN grid as p_mid. 
!           IPV_vs_fNsquared_ratio is dimensionless.          [ unitless ]
! p_mid   : Mid pressure between p grid  (length nz-1)           [dbar]
!--------------------------------------------------------------------------

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)
integer :: nz, k

real (r14) :: dsa, sa_mid, dct, ct_mid, dp, p_ref
real (r14) :: alpha_mid, gsw_alpha, beta_mid, gsw_beta
real (r14) :: alpha_pref, beta_pref, numerator, denominator
real (r14), dimension(nz) :: sa, ct, p
real (r14), dimension(nz-1) :: ipv_vs_fnsquared_ratio, p_mid

do k = 1, nz-1
	dsa = (sa(k+1) - sa(k))
	sa_mid = 0.5d0*(sa(k) + sa(k+1))
	dct = (ct(k+1) - ct(k))
	ct_mid = 0.5d0*(ct(k) + ct(k+1))
	dp = (p(k+1) - p(k))
	p_mid(k) = 0.5d0*(p(k) + p(k+1))

	alpha_mid = gsw_alpha(sa_mid,ct_mid,p_mid(k))
	beta_mid = gsw_beta(sa_mid,ct_mid,p_mid(k))
	alpha_pref = gsw_alpha(sa_mid,ct_mid,p_ref)
	beta_pref = gsw_beta(sa_mid,ct_mid,p_ref)

	numerator = dct*alpha_pref - dsa*beta_pref
	denominator = dct*alpha_mid - dsa*beta_mid

	if (denominator.eq.0d0) then
		ipv_vs_fnsquared_ratio(k) = 9d15
	else
		ipv_vs_fnsquared_ratio(k) = numerator/denominator
	end if
end do

return
end

!--------------------------------------------------------------------------
