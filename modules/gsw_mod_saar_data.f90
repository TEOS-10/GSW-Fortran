!==========================================================================
module gsw_mod_saar_data
!==========================================================================

implicit none

integer, parameter :: gsd_r14 = selected_real_kind(14,30)

logical :: saar_loaded = .false.
logical :: delta_loaded = .false.

integer, dimension(4) :: deli = (/0,1,1,0/), delj = (/0,0,1,1/)

integer, parameter :: nx = 91, ny = 45, nz = 45

real (gsd_r14), save, dimension(nz) :: p_ref
real (gsd_r14), save, dimension(ny) :: lats_ref
real (gsd_r14), save, dimension(nx) :: longs_ref
real (gsd_r14), save, dimension(ny,nx) :: ndepth_ref 
real (gsd_r14), save, dimension(nz,ny,nx) :: saar_ref
real (gsd_r14), save, dimension(nz,ny,nx) :: delta_sa_ref

integer, parameter :: npan = 6
real (gsd_r14), dimension(npan) :: longs_pan, lats_pan

data longs_pan /260.00d0, 272.59d0, 276.50d0, 278.65d0, 280.73d0, 292.0d0/
data  lats_pan / 19.55d0,  13.97d0,   9.60d0,   8.10d0,   9.33d0,   3.4d0/

end module

!--------------------------------------------------------------------------



