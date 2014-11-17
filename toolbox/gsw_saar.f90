!==========================================================================
function gsw_saar(p,long,lat)
!==========================================================================

! Calculates the Absolute Salinity Anomaly Ratio, SAAR.
!
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_saar : Absolute Salinity Anomaly Ratio               [unitless]

implicit none

!integer, parameter :: int9 = selected_int_kind(9) 
integer, parameter :: r14 = selected_real_kind(14,30)

integer, parameter :: nx=91, ny=45, nz=45

integer :: indx0, indy0, indz0, i, j, icalled, k
integer :: nmean, flag_saar
integer, dimension(4) :: deli, delj

real (r14), dimension(4) :: saar, saar_old
real (r14), dimension(nz) :: p_ref
real (r14), dimension(ny) :: lats_ref
real (r14), dimension(nx) :: longs_ref
real (r14), dimension(ny,nx) :: ndepth_ref 
real (r14), dimension(nz,ny,nx) :: saar_ref
!real (r14), dimension(nz,ny,nx) :: delta_sa_ref
real (r14) :: p, lat, long, dlong, dlat
real (r14) :: gsw_saar, p0_original, lon0_in, sa_upper, sa_lower 
real (r14) :: r1, s1, t1, saar_mean, ndepth_max, p_tmp, long_tmp

data deli/0,1,1,0/, delj/0,0,1,1/

data icalled/0/, flag_saar/0/

save icalled, flag_saar, longs_ref, lats_ref, p_ref, ndepth_ref, saar_ref

gsw_saar = 9d15

if(lat .lt. -86d0 .or. lat .gt. 90d0) then
 gsw_saar = 9d15
 return
end if

long_tmp = long
if(long_tmp.lt.0d0) then
 long_tmp = long_tmp + 360d0
end if

if(icalled.eq.0d0) then
  icalled = 1
   open(10,file='gsw_data_v3_0.dat',status='old',err=1)
   flag_saar = 1
   read(10,*) (longs_ref(i), i=1,nx)
   read(10,*) (lats_ref(i), i=1,ny)
   read(10,*) (p_ref(i), i=1,nz)
   read(10,*) ((ndepth_ref(j,i), j=1,ny), i=1,nx)
   read(10,*) (((saar_ref(k,j,i), k=1,nz), j=1,ny), i=1,nx)
   !read(10,*) (((delta_sa_ref(k,j,i), k=1,nz), j=1,ny), i=1,nx)
   close(10)
   go to 2
1  saar_ref = 9d15
   flag_saar = 0
2  continue
end if

if (flag_saar.eq.0d0) then
   write(*,*) "*** gsw_data_v3_0.dat is missing !!! ***"
   write(*,*) "Set the full path of gsw_data_v3_0.dat in gsw_saar"
end if

!Set gsw_saar = 9d15 and return if there is no data file present
if(flag_saar.eq.0d0) then
 gsw_saar = 9d15
 return
endif

dlong = longs_ref(2)-longs_ref(1)
dlat = lats_ref(2)-lats_ref(1)

indx0 = floor(1d0 + (nx-1d0)*(long_tmp-longs_ref(1))/(longs_ref(nx)-longs_ref(1)))
if(indx0.eq.nx) then
   indx0 = nx-1
end if

indy0 = floor(1d0 + (ny-1d0)*(lat-lats_ref(1))/(lats_ref(ny)-lats_ref(1)))
if(indy0.eq.ny) then
   indy0 = ny-1d0
end if

ndepth_max = -1
do k = 1,4
   if(ndepth_ref(indy0+delj(k),indx0+deli(k)).gt.0.d0) then
      ndepth_max = max(ndepth_max,ndepth_ref(indy0+delj(k),indx0+deli(k)))
   end if
end do

if(ndepth_max.eq.-1d0) then
  gsw_saar = 0d0 
   return
end if 

p0_original = p
p_tmp = p
if(p_tmp.gt.p_ref(int(ndepth_max))) then
 p_tmp = p_ref(int(ndepth_max))
end if
call indx(p_ref,nz,p_tmp,indz0)
    
r1 = (long_tmp-longs_ref(indx0))/(longs_ref(indx0+1)-longs_ref(indx0));
s1 = (lat-lats_ref(indy0))/(lats_ref(indy0+1)-lats_ref(indy0));
t1 = (p_tmp-p_ref(indz0))/(p_ref(indz0+1)-p_ref(indz0));

do k = 1,4
   saar(k) = saar_ref(indz0,indy0+delj(k),indx0+deli(k))
end do

if(260d0.le.long_tmp .and.long_tmp.le.291.999d0.and.3.4d0.le.lat.and.lat.le.19.55d0) then
  saar_old = saar
  call gsw_add_barrier(saar_old,long_tmp,lat,longs_ref(indx0),lats_ref(indy0),dlong,dlat,saar)
else if(abs(sum(saar)).ge.1d10) then
  saar_old = saar
  call gsw_add_mean(saar_old,long_tmp,lat,saar)
end if

sa_upper = (1d0-s1)*(saar(1) + r1*(saar(2)-saar(1))) + s1*(saar(4) + r1*(saar(3)-saar(4)))

do k = 1,4
   saar(k) = saar_ref(indz0+1,indy0+delj(k),indx0+deli(k))
end do

if(260d0.le.long_tmp.and.long_tmp.le.291.999d0.and.3.4d0.le.lat.and.lat.le.19.55d0) then
   saar_old = saar
   call gsw_add_barrier(saar_old,long_tmp,lat,longs_ref(indx0),lats_ref(indy0),dlong,dlat,saar)
else if(abs(sum(saar)).ge.1d10) then 
   saar_old = saar
   call gsw_add_mean(saar_old,long_tmp,lat,saar)
end if

sa_lower = (1d0-s1)*(saar(1) + r1*(saar(2)-saar(1))) + s1*(saar(4) + r1*(saar(3)-saar(4)))
if(abs(sa_lower).ge.1d10) then
  sa_lower = sa_upper
end if
gsw_saar = sa_upper + t1*(sa_lower-sa_upper)

if(abs(gsw_saar).ge.1d10) then
   gsw_saar = 9d15
endif
  
return
end function

!--------------------------------------------------------------------------

