!==========================================================================
function gsw_sp_from_sa_baltic(sa,long,lat)
!==========================================================================

! For the Baltic Sea, calculates Practical Salinity with a value
! computed analytically from Absolute Salinity
!
! sa     : Absolute Salinity                               [g/kg]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sp_from_sa_baltic  : Practical Salinity              [unitless]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), dimension(2) :: xb_right, yb_right
real (r14), dimension(3) :: xb_left, yb_left
real (r14) :: sa, long, lat, gsw_sp_from_sa_baltic, xinterp1, xx_left, xx_right

data xb_left/12.6d0, 7.d0, 26.d0/, yb_left/50.d0, 59.d0, 69.d0/
data xb_right/45.d0, 26.d0/, yb_right/50.d0, 69.d0/

if(xb_left(2).lt.long .and. long.lt.xb_right(1) .and. yb_left(1).lt.lat .and. lat.lt.yb_left(3)) then
  
    xx_left = xinterp1(yb_left, xb_left, 3, lat)
    
    xx_right = xinterp1(yb_right, xb_right, 2, lat)
    
    if(xx_left.le.long .and. long.le.xx_right) then
        gsw_sp_from_sa_baltic = (35.d0/(35.16504d0 - 0.087d0))*(sa - 0.087d0)
    else
        gsw_sp_from_sa_baltic = 9d15
    end if
     
else
    gsw_sp_from_sa_baltic = 9d15
end if

return
end

!--------------------------------------------------------------------------

