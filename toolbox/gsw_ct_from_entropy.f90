!=========================================================================
elemental function gsw_ct_from_entropy (sa, entropy)
!=========================================================================
!
!  Calculates Conservative Temperature with entropy as an input variable.  
!
!  SA       =  Absolute Salinity                                   [ g/kg ]
!  entropy  =  specific entropy                                   [ deg C ]
!
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_ct_from_pt, gsw_pt_from_entropy

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, entropy

real (r14) :: gsw_ct_from_entropy

real (r14) :: pt

pt = gsw_pt_from_entropy(sa,entropy)
gsw_ct_from_entropy = gsw_ct_from_pt(sa,pt)

return
end function

!--------------------------------------------------------------------------
