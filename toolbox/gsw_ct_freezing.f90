!==========================================================================
elemental function gsw_ct_freezing (sa, p, saturation_fraction, exact)
!==========================================================================
!
!  Calculates the Conservative Temperature at which seawater freezes.
!  The error of this fit ranges between -5e-4 K and 6e-4 K when compared 
!  with the Conservative Temperature calculated from the exact in-situ 
!  freezing temperature which is found by a Newton-Raphson iteration of the 
!  equality of the chemical potentials of water in seawater and in ice.  
!  Note that the Conservative temperature freezing temperature can be found
!  by this exact method using the function gsw_CT_freezing.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!  saturation_fraction = the saturation fraction of dissolved air in 
!                        seawater
!
!  CT_freezing = Conservative Temperature at freezing of seawater [ deg C ]
!                That is, the freezing temperature expressed in
!                terms of Conservative Temperature (ITS-90).                
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_ct_freezing_exact, gsw_ct_freezing_poly

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, p, saturation_fraction
logical, intent(in), optional :: exact

real (r14) :: gsw_ct_freezing

logical :: do_exact

if (present(exact)) then
	do_exact = exact
else
	do_exact = .false.
end if

if (do_exact) then
	gsw_ct_freezing = gsw_ct_freezing_exact(sa, p, saturation_fraction)
else
	gsw_ct_freezing = gsw_ct_freezing_poly(sa, p, saturation_fraction)
end if

return
end function

!--------------------------------------------------------------------------
