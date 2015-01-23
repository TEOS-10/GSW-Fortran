!==========================================================================
elemental function gsw_t_freezing_poly (sa, p, saturation_fraction, polynomial)
!==========================================================================
!
!  Calculates the in-situ temperature at which seawater freezes from a 
!  computationally efficient polynomial.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!  saturation_fraction = the saturation fraction of dissolved air in 
!                        seawater
!
!  t_freezing = in-situ temperature at which seawater freezes.    [ deg C ]
!               (ITS-90)                
!--------------------------------------------------------------------------

use gsw_mod_teos10_constants, only : gsw_sso

use gsw_mod_freezing_poly_coefficients

use gsw_mod_toolbox, only : gsw_ct_freezing_poly, gsw_t_from_ct

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, p
real (r14), intent(in), optional :: saturation_fraction
logical, intent(in), optional :: polynomial

real (r14) :: gsw_t_freezing_poly

real (r14) :: p_r, sa_r, x, ctf, sfrac
logical :: direct_poly

if (present(polynomial)) then
   direct_poly = polynomial
else
   direct_poly = .false.
end if

if (.not. direct_poly) then

   if (present(saturation_fraction)) then
      sfrac = saturation_fraction
   else
      sfrac = 1d0
   end if

   ctf = gsw_ct_freezing_poly(sa,p,sfrac)
   gsw_t_freezing_poly = gsw_t_from_ct(sa,ctf,p)

else

   ! Alternative calculation ...
   sa_r = sa*1d-2
   x = sqrt(sa_r)
   p_r = p*1d-4

   gsw_t_freezing_poly = t0 &
       + sa_r*(t1 + x*(t2 + x*(t3 + x*(t4 + x*(t5 + t6*x))))) &
       + p_r*(t7 + p_r*(t8 + t9*p_r)) &
       + sa_r*p_r*(t10 + p_r*(t12 + p_r*(t15 + t21*sa_r)) &
       + sa_r*(t13 + t17*p_r + t19*sa_r) &
       + x*(t11 + p_r*(t14 + t18*p_r) + sa_r*(t16 + t20*p_r + t22*sa_r)))

   if (.not. present(saturation_fraction)) return

   ! Adjust for the effects of dissolved air
   gsw_t_freezing_poly = gsw_t_freezing_poly - &
                         saturation_fraction*(1d-3)*(2.4d0 - sa/(2d0*gsw_sso))
end if

return
end function

!--------------------------------------------------------------------------
