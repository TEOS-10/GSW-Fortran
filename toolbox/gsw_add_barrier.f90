!==========================================================================
subroutine gsw_add_barrier(input_data,long,lat,long_grid,lat_grid,dlong_grid,dlat_grid,output_data)
!==========================================================================

!  Adds a barrier through Central America (Panama) and then averages
!  over the appropriate side of the barrier
! 
!  data_in      :  data                                                     [unitless]
!  long         :  Longitudes of data in decimal degrees east               [ 0 ... +360 ]
!  lat          :  Latitudes of data in decimal degrees north               [ -90 ... +90 ]
!  longs_grid   :  Longitudes of regular grid in decimal degrees east       [ 0 ... +360 ]
!  lats_grid    :  Latitudes of regular grid in decimal degrees north       [ -90 ... +90 ]
!  dlongs_grid  :  Longitude difference of regular grid in decimal degrees  [ deg longitude ]
!  dlats_grid   :  Latitude difference of regular grid in decimal degrees   [ deg latitude ]
!
! gsw_add_barrier  : average of data depending on which side of the 
!                    Panama cannal it is on                                 [unitless]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer, dimension(4) :: above_line
integer k, nmean, above_line0, kk
real (r14), dimension(4) :: input_data, output_data
real (r14), dimension(6) :: longs_pan, lats_pan
real (r14) :: long, lat, r, lats_line, long_grid, lat_grid
real (r14) :: dlong_grid, dlat_grid, data_mean

data longs_pan/260.0000d0, 272.5900d0, 276.5000d0, 278.6500d0, 280.7300d0, 292.000d0/ 
data  lats_pan/ 19.5500d0,  13.9700d0,   9.6000d0,   8.1000d0,   9.3300d0,   3.400d0/ 

call indx(longs_pan,6,long,k)                            !   the long/lat point
r = (long-longs_pan(k))/(longs_pan(k+1)-longs_pan(k))
lats_line = lats_pan(k) + r*(lats_pan(k+1)-lats_pan(k))

if(lats_line.le.lat) then
   above_line0 = 1
else
   above_line0 = 0
end if

call indx(longs_pan,6,long_grid,k)                                     !  the 1 and 4 long/lat points 
r = (long_grid-longs_pan(k))/(longs_pan(k+1)-longs_pan(k))
lats_line = lats_pan(k) + r*(lats_pan(k+1)-lats_pan(k))

if(lats_line.le.lat_grid) then
   above_line(1) = 1
else
   above_line(1) = 0
end if

if(lats_line.le.lat_grid+dlat_grid) then
   above_line(4) = 1
else
   above_line(4) = 0
end if

call indx(longs_pan,6,long_grid+dlong_grid,k)                              !  the 2 and 3 long/lat points 
r = (long_grid+dlong_grid-longs_pan(k))/(longs_pan(k+1)-longs_pan(k))
lats_line = lats_pan(k) + r*(lats_pan(k+1)-lats_pan(k))

if(lats_line.le.lat_grid) then
   above_line(2) = 1
else
   above_line(2) = 0
end if

if(lats_line.le.lat_grid+dlat_grid) then
   above_line(3) = 1
else
   above_line(3) = 0
end if

nmean = 0 
data_mean = 0.d0

do kk = 1,4
   if ((abs(input_data(kk)).le.100d0).and.above_line0.eq.above_line(kk)) then
      nmean = nmean+1
      data_mean = data_mean+input_data(kk)
   end if
end do

if(nmean .eq. 0d0)then
   data_mean = 0d0    !errorreturn
else
   data_mean = data_mean/nmean
endif

do kk = 1,4
   if((abs(input_data(kk)).ge.1d10).or.above_line0.ne.above_line(kk)) then
      output_data(kk) = data_mean
   else
      output_data(kk) = input_data(kk)
   end if
end do

return
end subroutine

!--------------------------------------------------------------------------
