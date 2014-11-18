!==========================================================================
subroutine gsw_saar_init (load_delta)
!==========================================================================
!
! Loads gsw_data_v3_0.dat into memory
!
! load_delta  :  optionally loads delta_sa_ref data
!--------------------------------------------------------------------------

use gsw_mod_saar_data

implicit none

logical, intent(in) :: load_delta

integer :: i, j, k, iostat

open(10,file='gsw_data_v3_0.dat',status='old',iostat=iostat)

if (iostat .eq. 0) then

   read(10,*) (longs_ref(i), i=1,nx)
   read(10,*) (lats_ref(i), i=1,ny)
   read(10,*) (p_ref(i), i=1,nz)
   read(10,*) ((ndepth_ref(j,i), j=1,ny), i=1,nx)
   read(10,*) (((saar_ref(k,j,i), k=1,nz), j=1,ny), i=1,nx)
   if (load_delta) then
      read(10,*) (((delta_sa_ref(k,j,i), k=1,nz), j=1,ny), i=1,nx)
      delta_loaded = .true.
   endif
   saar_loaded = .true.
   close(10)

else

   write(*,*) "*** gsw_data_v3_0.dat is missing !!! ***"
   write(*,*) "Set the full path of gsw_data_v3_0.dat in gsw_saar_init"
   stop

end if

return
end subroutine

!--------------------------------------------------------------------------



