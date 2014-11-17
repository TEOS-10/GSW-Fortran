!==========================================================================
subroutine gsw_nsquared (sa, ct, p, lat, nz, n2, p_mid)
!==========================================================================
!
!  Calculates the buoyancy frequency squared (N^2)(i.e. the Brunt-Vaisala 
!  frequency squared) at the mid pressure from the equation,
!
!
!           2      2             beta.d(SA) - alpha.d(CT)
!         N   =  g  .rho_local. -------------------------
!                                          dP
!
!  The pressure increment, dP, in the above formula is in Pa, so that it is
!  10^4 times the pressure increment dp in dbar. 
!
!  Note. This routine uses rho from "gsw_rho", which is the computationally
!  efficient 48-term expression for density in terms of SA, CT and p.  The    
!  48-term equation has been fitted in a restricted range of parameter 
!  space, and is most accurate inside the "oceanographic funnel" described 
!  in IOC et al. (2010).  The GSW library function "gsw_infunnel(SA,CT,p)" 
!  is avaialble to be used if one wants to test if some of one's data lies
!  outside this "funnel".  

! sa     : Absolute Salinity         (a profile (length nz))     [g/kg]
! ct     : Conservative Temperature  (a profile (length nz))     [deg C]
! p      : sea pressure              (a profile (length nz))     [dbar]
! lat    : latitude                  (a profile (length nz))     [deg N]                            
! nz     : number of bottles                             
! n2     : Brunt-Vaisala Frequency squared  (length nz-1)        [s^-2]
! p_mid  : Mid pressure between p grid      (length nz-1)        [dbar]
!--------------------------------------------------------------------------

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)
integer :: nz, k

real (r14), parameter :: db2pa = 1d4
real (r14) :: grav_local, dsa, sa_mid, dct, ct_mid, dp, rho_mid, gsw_rho
real (r14) :: alpha_mid, gsw_alpha, beta_mid, gsw_beta, gsw_grav
real (r14), dimension(nz) :: sa, ct, p, lat
real (r14), dimension(nz-1) :: p_mid, n2

do k = 1, nz-1
	grav_local = 0.5*(gsw_grav(lat(k),p(k)) + gsw_grav(lat(k+1),p(k+1)))

	dsa = (sa(k+1) - sa(k))
	sa_mid = 0.5*(sa(k) + sa(k+1))
	dct = (ct(k+1) - ct(k))
	ct_mid = 0.5*(ct(k) + ct(k+1))
	dp = (p(k+1) - p(k))
	p_mid(k) = 0.5*(p(k) + p(k+1))

	rho_mid = gsw_rho(sa_mid,ct_mid,p_mid(k))
	alpha_mid = gsw_alpha(sa_mid,ct_mid,p_mid(k))
	beta_mid = gsw_beta(sa_mid,ct_mid,p_mid(k))

	n2(k) = (grav_local*grav_local) * (rho_mid/(db2pa*dp)) &
	        * (beta_mid*dsa - alpha_mid*dct)
end do

return
end subroutine

!--------------------------------------------------------------------------
