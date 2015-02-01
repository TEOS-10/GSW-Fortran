!==========================================================================
module gsw_mod_saar_data
!==========================================================================

use gsw_mod_kinds

implicit none

logical, save :: saar_loaded = .false.
logical, save :: delta_loaded = .false.

integer, dimension(4) :: deli = (/0,1,1,0/), delj = (/0,0,1,1/)

integer, save :: nx, ny, nz

real (r8), save, dimension(:), allocatable :: p_ref, lats_ref, longs_ref
real (r8), save, dimension(:,:), allocatable :: ndepth_ref 
real (r8), save, dimension(:,:,:), allocatable :: saar_ref, delta_sa_ref

integer, parameter :: npan = 6
real (r8), dimension(npan) :: longs_pan, lats_pan

data longs_pan /260.00_r8, 272.59_r8, 276.50_r8, 278.65_r8, 280.73_r8, 292.0_r8/
data  lats_pan / 19.55_r8,  13.97_r8,   9.60_r8,   8.10_r8,   9.33_r8,   3.4_r8/

end module

!--------------------------------------------------------------------------



