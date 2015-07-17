!==========================================================================
subroutine gsw_saar_init (load_delta)
!==========================================================================
!
! Loads gsw_data_v3_0.nc contents into memory
!
! load_delta  :  optionally loads delta_sa_ref data
!--------------------------------------------------------------------------

use gsw_mod_netcdf

use gsw_mod_saar_data

implicit none

logical, intent(in) :: load_delta

call ncdf_open('gsw_data_v3_0.nc')

nx = ncdf_get_dim("nx")
ny = ncdf_get_dim("ny")
nz = ncdf_get_dim("nz")

allocate(p_ref(nz))
allocate(lats_ref(ny))
allocate(longs_ref(nx))
allocate(ndepth_ref(ny,nx))
allocate(saar_ref(nz,ny,nx))

call ncdf_get_var("p_ref", var1=p_ref)
call ncdf_get_var("lats_ref", var1=lats_ref)
call ncdf_get_var("longs_ref", var1=longs_ref)
call ncdf_get_var("ndepth_ref", var2=ndepth_ref)
call ncdf_get_var("SAAR_ref", var3=saar_ref)
if (load_delta) then
    allocate(delta_sa_ref(nz,ny,nx))
    call ncdf_get_var("deltaSA_ref", var3=delta_sa_ref)
    delta_loaded = .true.
endif
saar_loaded = .true.

call ncdf_close()
return

end subroutine

!--------------------------------------------------------------------------
