!==========================================================================
elemental function gsw_t_freezing (sa, p, saturation_fraction, exact)
!==========================================================================
!
!  Calculates the in-situ temperature at which seawater freezes.
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

use gsw_mod_toolbox, only : gsw_t_freezing_exact, gsw_t_freezing_poly

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, p, saturation_fraction
logical, intent(in), optional :: exact

real (r14) :: gsw_t_freezing

logical :: do_exact

if (present(exact)) then
	do_exact = exact
else
	do_exact = .false.
end if

if (do_exact) then
	gsw_t_freezing = gsw_t_freezing_exact(sa,p,saturation_fraction)
else
	gsw_t_freezing = gsw_t_freezing_poly(sa,p,saturation_fraction)
end if

return
end function

!--------------------------------------------------------------------------
