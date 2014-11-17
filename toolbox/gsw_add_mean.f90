!==========================================================================
subroutine gsw_add_mean(data_in,long,lat,data_out)
!==========================================================================

! Replaces NaN's with non-nan mean of the 4 adjacent neighbours
!
! data_in   : data set of the 4 adjacent neighbours   
! long      : longitude
! lat       : latitude
!
! data_out : non-nan mean of the 4 adjacent neighbours     [unitless]

implicit none

integer, parameter :: int9 = selected_int_kind(9)
integer, parameter :: r14 = selected_real_kind(14,30)

integer :: k, nmean

real (r14), dimension(4) :: data_in, data_out 
real (r14) :: data_mean, long, lat

nmean = 0
data_mean = 0.d0

do k = 1,4
   if (abs(data_in(k)).le.100d0) then
      nmean = nmean+1
      data_mean = data_mean+data_in(k)
   end if
end do

if(nmean.eq.0)then
   data_mean = 0d0    !errorreturn
else
   data_mean = data_mean/nmean
endif

do k = 1,4
   if(abs(data_in(k)).ge.100d0) then
      data_out(k) = data_mean
   else
      data_out(k) = data_in(k)
   end if
end do

return
end subroutine

!--------------------------------------------------------------------------

