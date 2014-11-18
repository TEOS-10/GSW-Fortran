! =========================================================================
elemental function gsw_t_from_pt0_ice (pt0_ice, p)
! =========================================================================
!
!  Calculates in-situ temperature from the potential temperature of ice Ih 
!  with reference pressure, p_ref, of 0 dbar (the surface), and the 
!  in-situ pressure.
!
!  pt0_ice  =  potential temperature of ice Ih with reference pressure of 
!              zero dbar (ITS-90)                                 [ deg C ]
!  p        =  sea pressure                                        [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_pt_from_t_ice

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: pt0_ice, p

real (r14) :: gsw_t_from_pt0_ice

real (r14), parameter :: p0 = 0d0

gsw_t_from_pt0_ice = gsw_pt_from_t_ice(pt0_ice,p0,p)

return
end function

!--------------------------------------------------------------------------
