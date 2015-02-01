!==========================================================================
elemental function gsw_deltasa_atlas (p, long, lat)
!==========================================================================
!
! Calculates the Absolute Salinity Anomaly atlas value, deltaSA_atlas.
!
! p      : sea pressure                                    [dbar]
! long   : longiture                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_deltasa_atlas : Absolute Salinity Anomaly atlas value    [g/kg]
!--------------------------------------------------------------------------

use gsw_mod_toolbox, only : gsw_add_barrier, gsw_add_mean, gsw_util_indx

use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit

use gsw_mod_saar_data

use gsw_mod_kinds

implicit none

real (r8), intent(in) :: p, long, lat

real (r8) :: gsw_deltasa_atlas

integer :: indx0, indy0, indz0, i, j, k, nmean

real (r8), dimension(4) :: dsar, dsar_old
real (r8) :: dlong, dlat
real (r8) :: p0_original, lon0_in, sa_upper, sa_lower 
real (r8) :: r1, s1, t1, dsar_mean, ndepth_max, p_tmp, long_tmp

character (*), parameter :: func_name = "gsw_deltasa_atlas"

if (.not. delta_loaded) then
   gsw_deltasa_atlas = gsw_error_code(1,func_name)
   return
end if

if (lat .lt. -86.0_r8 .or. lat .gt. 90.0_r8) then
   gsw_deltasa_atlas = gsw_error_code(2,func_name)
   return
end if

long_tmp = long
if (long_tmp.lt.0.0_r8) long_tmp = long_tmp + 360.0_r8

dlong = longs_ref(2) - longs_ref(1)
dlat = lats_ref(2) - lats_ref(1)

indx0 = floor(1.0_r8 + (nx-1)*(long_tmp-longs_ref(1))/(longs_ref(nx)-longs_ref(1)))
if (indx0.eq.nx) indx0 = nx - 1

indy0 = floor(1.0_r8 + (ny-1.0_r8)*(lat-lats_ref(1))/(lats_ref(ny)-lats_ref(1)))
if (indy0.eq.ny) indy0 = ny - 1

ndepth_max = -1
do k = 1, 4
   if (ndepth_ref(indy0+delj(k),indx0+deli(k)).gt.0.0_r8) &
      ndepth_max = max(ndepth_max,ndepth_ref(indy0+delj(k),indx0+deli(k)))
end do

if (ndepth_max.eq.-1.0_r8) then
   gsw_deltasa_atlas = 0.0_r8 
   return
end if 

p0_original = p
p_tmp = p
if (p_tmp.gt.p_ref(int(ndepth_max))) p_tmp = p_ref(int(ndepth_max))
call gsw_util_indx(p_ref,nz,p_tmp,indz0)
    
r1 = (long_tmp-longs_ref(indx0))/(longs_ref(indx0+1)-longs_ref(indx0))
s1 = (lat-lats_ref(indy0))/(lats_ref(indy0+1)-lats_ref(indy0))
t1 = (p_tmp-p_ref(indz0))/(p_ref(indz0+1)-p_ref(indz0))

do k = 1, 4
   dsar(k) = delta_sa_ref(indz0,indy0+delj(k),indx0+deli(k))
end do

if ( longs_pan(1).le.long_tmp .and. long_tmp.le.longs_pan(npan)-0.001_r8 .and. &
   lats_pan(npan).le.lat      .and.      lat.le.lats_pan(1)) then
   dsar_old = dsar
   call gsw_add_barrier(dsar_old,long_tmp,lat,longs_ref(indx0), &
   			lats_ref(indy0),dlong,dlat,dsar)
else if (abs(sum(dsar)).ge.1e10_r8) then 
   dsar_old = dsar
   call gsw_add_mean(dsar_old,dsar)
end if

sa_upper = (1.0_r8-s1)*(dsar(1) + r1*(dsar(2)-dsar(1))) + s1*(dsar(4) + &
		r1*(dsar(3)-dsar(4)))

do k = 1,4
   dsar(k) = delta_sa_ref(indz0+1,indy0+delj(k),indx0+deli(k))
end do

if ( longs_pan(1).le.long_tmp .and. long_tmp.le.longs_pan(npan)-0.001_r8 .and. &
   lats_pan(npan).le.lat      .and.      lat.le.lats_pan(1)) then
   dsar_old = dsar
   call gsw_add_barrier(dsar_old,long_tmp,lat,longs_ref(indx0), &
   			lats_ref(indy0),dlong,dlat,dsar)
else if (abs(sum(dsar)).ge.1e10_r8) then 
   dsar_old = dsar
   call gsw_add_mean(dsar_old,dsar)
end if

sa_lower = (1.0_r8-s1)*(dsar(1) + r1*(dsar(2)-dsar(1))) + s1*(dsar(4) + &
		r1*(dsar(3)-dsar(4)))
if (abs(sa_lower).ge.1e10_r8) sa_lower = sa_upper
gsw_deltasa_atlas = sa_upper + t1*(sa_lower-sa_upper)

if (abs(gsw_deltasa_atlas).ge.1e10_r8) &
	gsw_deltasa_atlas = gsw_error_code(3,func_name)
  
return
end function

!--------------------------------------------------------------------------
