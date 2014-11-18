!==========================================================================
elemental function gsw_t_from_ct (sa, ct, p)  
!==========================================================================
!
! Calculates in-situ temperature from Conservative Temperature of seawater  
!
! sa      : Absolute Salinity                              [g/kg]
! ct      : Conservative Temperature                       [deg C]
!
! gsw_t_from_ct : in-situ temperature                      [deg C]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_pt_from_ct, gsw_pt_from_t

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, ct, p  

real (r14) :: gsw_t_from_ct

real (r14) :: pt0

real (r14), parameter :: p0 = 0d0

pt0 = gsw_pt_from_ct(sa,ct)
gsw_t_from_ct = gsw_pt_from_t(sa,pt0,p0,p)

return
end function

!--------------------------------------------------------------------------
