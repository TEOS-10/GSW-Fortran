!==========================================================================
elemental function gsw_sp_from_sa_baltic (sa, long, lat)
!==========================================================================
!
! For the Baltic Sea, calculates Practical Salinity with a value
! computed analytically from Absolute Salinity
!
! sa     : Absolute Salinity                               [g/kg]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sp_from_sa_baltic  : Practical Salinity              [unitless]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_util_xinterp1

use gsw_mod_baltic_data

use gsw_mod_teos10_constants, only : gsw_sso

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), intent(in) :: sa, long, lat

real (r14) :: gsw_sp_from_sa_baltic

real (r14) :: xx_left, xx_right

if (xb_left(2).lt.long .and. long.lt.xb_right(1) .and. &
    yb_left(1).lt.lat  .and.  lat.lt.yb_left(3)) then
  
    xx_left = gsw_util_xinterp1(yb_left, xb_left, 3, lat)
    
    xx_right = gsw_util_xinterp1(yb_right, xb_right, 2, lat)
    
    if(xx_left.le.long .and. long.le.xx_right) then
        gsw_sp_from_sa_baltic = (35.d0/(gsw_sso - 0.087d0))*(sa - 0.087d0)
    else
        gsw_sp_from_sa_baltic = 9d15
    end if
     
else
    gsw_sp_from_sa_baltic = 9d15
end if

return
end function

!--------------------------------------------------------------------------
