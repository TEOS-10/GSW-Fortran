!==========================================================================
pure subroutine gsw_nsquared_min (sa, ct, p, lat, n2, n2_p, &
                               n2_specvol, n2_alpha, n2_beta, dsa, dct, dp)
!==========================================================================
! 
!  Calculates the minimum buoyancy frequency squared (N^2) (i.e. the 
!  Brunt-Vaisala frequency squared) between two bottles from the equation,
!
!           2      2     beta.dSA - alpha.dCT
!         N   =  g  . -------------------------
!                         specvol_local.dP
!
!  The pressure increment, dP, in the above formula is in Pa, so that it is
!  10^4 times the pressure increment dp in dbar. 
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!  lat =  latitude in decimal degrees north                 [ -90 ... +90 ]
!  Note: If lat is outside this range, a default gravitational
!        acceleration of 9.7963 m/s^2 (Griffies, 2004) will be applied.
!
!  N2         =  minimum Brunt-Vaisala Frequency squared          [ 1/s^2 ]
!  N2_p       =  pressure of minimum N2                            [ dbar ]
!  N2_specvol =  specific volume at the minimum N2                [ kg/m3 ]
!  N2_alpha   =  thermal expansion coefficient with respect         [ 1/K ]
!                to Conservative Temperature at the minimum N2
!  N2_beta    =  saline contraction coefficient at constant        [ kg/g ]
!                Conservative Temperature at the minimum N2
!  dSA        =  difference in salinity between bottles            [ g/kg ]
!  dCT        =  difference in Conservative Temperature between   [ deg C ]
!                bottles
!  dp         =  difference in pressure between bottles            [ dbar ]
!
!==========================================================================

use gsw_mod_toolbox, only : gsw_grav, gsw_specvol_alpha_beta

use gsw_mod_kinds

use gsw_mod_teos10_constants, only : db2pa

implicit none

real (r8), intent(in) :: sa(:), ct(:), p(:), lat
real (r8), intent(out) :: n2(:), n2_p(:), n2_specvol(:), n2_alpha(:)
real (r8), intent(out) :: n2_beta(:), dsa(:), dct(:), dp(:)

integer :: i, ideep, ishallow, mp

real (r8) :: n2_deep, n2_shallow
real (r8), allocatable :: alpha(:), beta(:), specvol(:), grav2(:)

mp = size(sa)
allocate(grav2(mp),specvol(mp),alpha(mp),beta(mp))

if (lat .lt. -90.0_r8 .or. lat .gt. +90.0_r8) then
    grav2 = (/ (9.7963_r8**2, i=1,mp) /)
else
    grav2 = gsw_grav(lat,p)**2
end if

dp  =  p(2:mp) -  p(1:mp-1)
dsa = sa(2:mp) - sa(1:mp-1)
dct = ct(2:mp) - ct(1:mp-1)

call gsw_specvol_alpha_beta(sa,ct,p,specvol,alpha,beta)

ishallow = 1
ideep = 2
do i = 1, mp-1
    n2_shallow =  grav2(ishallow)/(specvol(ishallow)*db2pa*dp(i))* &
                 (beta(ishallow)*dsa(i) - alpha(ishallow)*dct(i))
    n2_deep =  grav2(ideep)/(specvol(ideep)*db2pa*dp(i))* &
	      (beta(ideep)*dsa(i) - alpha(ideep)*dct(i))
    if (n2_shallow .lt. n2_deep) then
        n2(i) = n2_shallow
        n2_p(i) = p(ishallow)
        n2_specvol(i) = specvol(ishallow)
        n2_alpha(i) = alpha(ishallow)
        n2_beta(i) = beta(ishallow)
    else
	n2(i) = n2_deep
        n2_p(i) = p(ideep)
        n2_specvol(i) = specvol(ideep)
        n2_alpha(i) = alpha(ideep)
        n2_beta(i) = beta(ideep)
    end if
    ishallow = ishallow + 1
    ideep = ideep + 1
end do

return
end subroutine

!--------------------------------------------------------------------------
