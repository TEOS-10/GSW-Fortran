!==========================================================================
function gsw_deltasa_atlas(p,long,lat)
!==========================================================================

! Calculates the Absolute Salinity Anomaly atlas value, deltaSA_atlas.
!
! p      : sea pressure                                    [dbar]
! long   : longiture                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_deltasa_atlas : Absolute Salinity Anomaly atlas value    [g/kg]

implicit none

!integer, parameter :: int9 = selected_int_kind(9) 
integer, parameter :: r14 = selected_real_kind(14,30)

integer, parameter :: nx=91, ny=45, nz=45

integer :: indx0, indy0, indz0, i, j, icalled2, k
integer :: nmean, flag_dsar
integer, dimension(4) :: deli, delj

real (r14), dimension(4) :: dsar, dsar_old
real (r14), dimension(nz) :: p_ref
real (r14), dimension(ny) :: lats_ref
real (r14), dimension(nx) :: longs_ref
real (r14), dimension(ny,nx) :: ndepth_ref 
real (r14), dimension(nz,ny,nx) :: saar_ref, delta_sa_ref
real (r14) :: p, lat, long, dlong, dlat
real (r14) :: gsw_deltasa_atlas, p0_original, lon0_in, sa_upper, sa_lower 
real (r14) :: r1, s1, t1, dsar_mean, ndepth_max, p_tmp, long_tmp

data deli/0,1,1,0/, delj/0,0,1,1/

data icalled2/0/, flag_dsar/0/

save icalled2, flag_dsar, longs_ref, lats_ref, p_ref, ndepth_ref, delta_sa_ref

gsw_deltasa_atlas = 9d15

if(lat .lt. -86d0 .or. lat .gt. 90d0) then
 gsw_deltasa_atlas = 9d15
 return
end if

long_tmp = long
if(long_tmp.lt.0d0) then
 long_tmp = long_tmp + 360
end if

if(icalled2.eq.0d0) then
   icalled2 = 1
   open(10,file='gsw_data_v3_0.dat',status='old',err=1)
   flag_dsar = 1
   read(10,*) (longs_ref(i), i=1,nx)
   read(10,*) (lats_ref(i), i=1,ny)
   read(10,*) (p_ref(i), i=1,nz)
   read(10,*) ((ndepth_ref(j,i), j=1,ny), i=1,nx)
   read(10,*) (((saar_ref(k,j,i), k=1,nz), j=1,ny), i=1,nx)
   read(10,*) (((delta_sa_ref(k,j,i), k=1,nz), j=1,ny), i=1,nx)
   close(10)
   go to 2
1  delta_sa_ref = 9d15
   flag_dsar = 0
2  continue
end if

if (flag_dsar.eq.0d0) then
   write(*,*) "*** gsw_data_v3_0.dat is missing !!! ***"
   write(*,*) "Set the full path of gsw_data_v3_0.dat in gsw_deltasa_atlas"
end if

!Set gsw_deltasa_atlas = 9d15 and return if there is no data set present
if(flag_dsar.eq.0d0) then
 gsw_deltasa_atlas = 9d15
 return
endif

dlong = longs_ref(2)-longs_ref(1)
dlat = lats_ref(2)-lats_ref(1)

indx0 = floor(1d0 + (nx-1)*(long_tmp-longs_ref(1))/(longs_ref(nx)-longs_ref(1)))
if(indx0.eq.nx) then
   indx0 = nx-1d0
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
  gsw_deltasa_atlas = 0d0 
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
   dsar(k) = delta_sa_ref(indz0,indy0+delj(k),indx0+deli(k))
end do

if(260d0.le.long_tmp.and.long_tmp.le.291.999d0.and.3.4d0.le.lat.and.lat.le.19.55d0) then
   dsar_old = dsar
   call gsw_add_barrier(dsar_old,long_tmp,lat,longs_ref(indx0),lats_ref(indy0),dlong,dlat,dsar)
else if(abs(sum(dsar)).ge.1d10) then 
   dsar_old = dsar
   call gsw_add_mean(dsar_old,long_tmp,lat,dsar)
end if

sa_upper = (1d0-s1)*(dsar(1) + r1*(dsar(2)-dsar(1))) + s1*(dsar(4) + r1*(dsar(3)-dsar(4)))

do k = 1,4
   dsar(k) = delta_sa_ref(indz0+1,indy0+delj(k),indx0+deli(k))
end do

if(260d0.le.long_tmp.and.long_tmp.le.291.999d0.and.3.4d0.le.lat.and.lat.le.19.55d0) then
   dsar_old = dsar
   call gsw_add_barrier(dsar_old,long_tmp,lat,longs_ref(indx0),lats_ref(indy0),dlong,dlat,dsar)
else if(abs(sum(dsar)).ge.1d10) then 
   dsar_old = dsar
   call gsw_add_mean(dsar_old,long_tmp,lat,dsar)
end if

sa_lower = (1d0-s1)*(dsar(1) + r1*(dsar(2)-dsar(1))) + s1*(dsar(4) + r1*(dsar(3)-dsar(4)))
if(abs(sa_lower).ge.1d10) then
  sa_lower = sa_upper
end if
gsw_deltasa_atlas = sa_upper + t1*(sa_lower-sa_upper)

if(abs(gsw_deltasa_atlas).ge.1d10) then
   gsw_deltasa_atlas = 9d15
endif
  
return
end function

!--------------------------------------------------------------------------
