program gsw_poly_check

use gsw_mod_kinds
use gsw_mod_netcdf
use gsw_mod_toolbox

implicit none

integer :: gsw_error_flag = 0

integer :: cast_m, cast_n, cast_mpres_m, cast_mpres_n, cast_ice_m
integer :: cast_ice_n, i

real (r8) :: saturation_fraction, pref

real (r8), dimension(:,:), allocatable :: ct, rt, sa, sk, sp, t, p
real (r8), dimension(:,:), allocatable :: p_shallow, p_deep, lat, long

real (r8), dimension(:,:), allocatable :: ct_arctic, sa_arctic, t_arctic
real (r8), dimension(:,:), allocatable :: p_arctic, t_ice, w_ice
real (r8), dimension(:,:), allocatable :: sa_seaice, t_seaice, w_seaice
real (r8), dimension(:,:), allocatable :: lat_ice, long_ice

real (r8), dimension(:), allocatable :: lat_cast, long_cast

real (r8), dimension(:,:), allocatable :: value, check_value
real (r8), dimension(:,:), allocatable :: val1, val2, val3, val4, val5, val6

real (r8), dimension(:,:), allocatable :: c, sr, sstar, pt, entropy
real (r8), dimension(:,:), allocatable :: h, ctf, tf, rho, diff
real (r8), dimension(:,:), allocatable :: ctf_poly, tf_poly, pt0

call gsw_saar_init (.true.)

call ncdf_open('gsw_data_v3_0.nc')

cast_m = ncdf_get_dim("test_cast_length")
cast_n = ncdf_get_dim("test_cast_number")

cast_ice_m = ncdf_get_dim("Arctic_test_cast_length")
cast_ice_n = ncdf_get_dim("Arctic_test_cast_number")

cast_mpres_m = ncdf_get_dim("test_cast_midpressure_length")
cast_mpres_n = ncdf_get_dim("test_cast_midpressure_number")

allocate(ct(cast_m,cast_n))
allocate(rt(cast_m,cast_n))
allocate(sa(cast_m,cast_n))
allocate(sk(cast_m,cast_n))
allocate(sp(cast_m,cast_n))
allocate(t(cast_m,cast_n))
allocate(p(cast_m,cast_n))
allocate(p_shallow(cast_m,cast_n))
allocate(p_deep(cast_m,cast_n))
allocate(lat(cast_m,cast_n))
allocate(long(cast_m,cast_n))
allocate(lat_cast(cast_n))
allocate(long_cast(cast_n))

allocate(ct_arctic(cast_ice_m,cast_ice_n))
allocate(sa_arctic(cast_ice_m,cast_ice_n))
allocate(t_arctic(cast_ice_m,cast_ice_n))
allocate(p_arctic(cast_ice_m,cast_ice_n))
allocate(sa_seaice(cast_ice_m,cast_ice_n))
allocate(t_seaice(cast_ice_m,cast_ice_n))
allocate(w_seaice(cast_ice_m,cast_ice_n))
allocate(t_ice(cast_ice_m,cast_ice_n))
allocate(w_ice(cast_ice_m,cast_ice_n))
allocate(lat_ice(cast_ice_m,cast_ice_n))
allocate(long_ice(cast_ice_m,cast_ice_n))
allocate(pt0(cast_ice_m,cast_ice_n))

allocate(value(cast_m,cast_n))
allocate(check_value(cast_m,cast_n))
allocate(val1(cast_m,cast_n))
allocate(val2(cast_m,cast_n))
allocate(val3(cast_m,cast_n))
allocate(val4(cast_m,cast_n))
allocate(val5(cast_m,cast_n))
allocate(val6(cast_m,cast_n))

allocate(c(cast_m,cast_n))
allocate(sr(cast_m,cast_n))
allocate(sstar(cast_m,cast_n))
allocate(pt(cast_m,cast_n))
allocate(entropy(cast_m,cast_n))
allocate(ctf(cast_m,cast_n))
allocate(rho(cast_m,cast_n))
allocate(tf(cast_m,cast_n))
allocate(ctf_poly(cast_m,cast_n))
allocate(tf_poly(cast_m,cast_n))
allocate(h(cast_m,cast_n))
allocate(diff(cast_m,cast_n))

call ncdf_get_var("CT_chck_cast", var2=ct)
call ncdf_get_var("Rt_chck_cast", var2=rt)
call ncdf_get_var("SA_chck_cast", var2=sa)
call ncdf_get_var("SK_chck_cast", var2=sk)
call ncdf_get_var("SP_chck_cast", var2=sp)
call ncdf_get_var("t_chck_cast", var2=t)
call ncdf_get_var("p_chck_cast", var2=p)
call ncdf_get_var("p_chck_cast_shallow", var2=p_shallow)
call ncdf_get_var("p_chck_cast_deep", var2=p_deep)

call ncdf_get_var("lat_chck_cast", var1=lat_cast)
call ncdf_get_var("long_chck_cast", var1=long_cast)
do i = 1, cast_n
    lat(:,i) = lat_cast(i)
    long(:,i) = long_cast(i)
end do

call ncdf_get_var("CT_Arctic", var2=ct_arctic)
call ncdf_get_var("SA_Arctic", var2=sa_arctic)
call ncdf_get_var("t_Arctic", var2=t_arctic)
call ncdf_get_var("p_Arctic", var2=p_arctic)
call ncdf_get_var("SA_seaice", var2=sa_seaice)
call ncdf_get_var("t_seaice", var2=t_seaice)
call ncdf_get_var("w_seaice", var2=w_seaice)
call ncdf_get_var("t_ice", var2=t_ice)
call ncdf_get_var("w_ice", var2=w_ice)

call ncdf_get_var("pr", var0=pref)

!------------------------------------------------------------------------------
call section_title('Freezing temperatures')

saturation_fraction = 0.5_r8

ctf = gsw_ct_freezing_exact(sa,p,saturation_fraction)
ctf_poly = gsw_ct_freezing_poly(sa,p,saturation_fraction)
call check_accuracy('CT_freezing',ctf,ctf_poly)

tf = gsw_t_freezing_exact(sa,p,saturation_fraction)
tf_poly = gsw_t_freezing_poly(sa,p,saturation_fraction)
call check_accuracy('t_freezing',tf,tf_poly)

val1 = gsw_sa_freezing_from_ct(ctf,p,saturation_fraction)
val2 = gsw_sa_freezing_from_ct_poly(ctf,p,saturation_fraction)
call check_accuracy('SA_freezing_from_CT',val1,val2)

val1 = gsw_sa_freezing_from_t(tf,p,saturation_fraction)
val2 = gsw_sa_freezing_from_t_poly(tf,p,saturation_fraction)
call check_accuracy('SA_freezing_from_t',val1,val2)

call gsw_ct_freezing_first_derivatives(sa,p,saturation_fraction,val1,val2)
call gsw_ct_freezing_first_derivatives_poly(sa,p,saturation_fraction,val3,val4)
call check_accuracy('CT_freezing_first_derivatives (ctf_sa)',val1,val3)
call check_accuracy('CT_freezing_first_derivatives (ctf_p)',val2,val4)

call gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,val1,val2)
call gsw_t_freezing_first_derivatives_poly(sa,p,saturation_fraction,val3,val4)
call check_accuracy('t_freezing_first_derivatives (tf_sa)',val1,val3)
call check_accuracy('t_freezing_first_derivatives (tf_p)',val2,val4)

!------------------------------------------------------------------------------
call section_title('Themodynamic properties of ice Ih')

deallocate(h,val1,val2,val3,val4,val5,val6)
allocate(h(cast_ice_m,cast_ice_n))
allocate(val1(cast_ice_m,cast_ice_n))
allocate(val2(cast_ice_m,cast_ice_n))
allocate(val3(cast_ice_m,cast_ice_n))
allocate(val4(cast_ice_m,cast_ice_n))
allocate(val5(cast_ice_m,cast_ice_n))
allocate(val6(cast_ice_m,cast_ice_n))

pt0 = gsw_pt0_from_t_ice(t_seaice,p_arctic)
h = gsw_pot_enthalpy_from_pt_ice(pt0)

val1 = gsw_pot_enthalpy_from_pt_ice(pt0)
val2 = gsw_pot_enthalpy_from_pt_ice_poly(pt0)
call check_accuracy('pot_enthalpy_from_pt_ice',val1,val2)

val1 = gsw_pt_from_pot_enthalpy_ice(h)
val2 = gsw_pt_from_pot_enthalpy_ice_poly(h)
call check_accuracy('pt_from_pot_enthalpy_ice',val1,val2)

!------------------------------------------------------------------------------
call section_title('Thermodynamic interaction between ice and seawater')

saturation_fraction = 0.0_r8

val1 = gsw_melting_ice_sa_ct_ratio(sa_arctic,ct_arctic,p_arctic,t_ice)
val2 = gsw_melting_ice_sa_ct_ratio_poly(sa_arctic,ct_arctic,p_arctic,t_ice)
call check_accuracy('melting_ice_SA_CT_ratio',val1,val2)

val1 = gsw_melting_ice_equilibrium_sa_ct_ratio(sa_arctic,p_arctic)
val2 = gsw_melting_ice_equilibrium_sa_ct_ratio_poly(sa_arctic,p_arctic)
call check_accuracy('melting_ice_equilibrium_SA_CT_ratio',val1,val2)

call gsw_frazil_ratios_adiabatic(sa_arctic,p_arctic,w_ice,val1,val2,val3)
call gsw_frazil_ratios_adiabatic_poly(sa_arctic,p_arctic,w_ice,val4,val5,val6)
call check_accuracy('frazil_ratios_adiabatic (dsa_dct)',val1,val4)
call check_accuracy('frazil_ratios_adiabatic (dsa_dp)',val2,val5)
call check_accuracy('frazil_ratios_adiabatic (dct_dp)',val3,val6)

!------------------------------------------------------------------------------
call section_title('Thermodynamic interaction between seaice and seawater')

value = gsw_melting_seaice_sa_ct_ratio(sa_arctic,ct_arctic,p_arctic, &
                                       sa_seaice,t_seaice)
!call check_accuracy('melting_seaice_SA_CT_ratio',value)

value = gsw_melting_seaice_equilibrium_sa_ct_ratio(sa_arctic,p_arctic)
!call check_accuracy('melting_seaice_equilibrium_SA_CT_ratio',value)

call gsw_melting_seaice_into_seawater(sa_arctic,ct_arctic,p_arctic, &
                                      w_seaice,sa_seaice,t_seaice,val1,val2)
!call check_accuracy('melting_seaice_into_seawater',val1, &
!                    'melting_seaice_into_seawater_SA_final')
!call check_accuracy('melting_seaice_into_seawater',val2, &
!                    'melting_seaice_into_seawater_CT_final')

call gsw_seaice_fraction_to_freeze_seawater(sa_arctic,ct_arctic,p_arctic, &
                                            sa_seaice,t_seaice,val1,val2,val3)
!call check_accuracy('seaice_fraction_to_freeze_seawater',val1, &
!                    'seaice_fraction_to_freeze_seawater_SA_freeze')
!call check_accuracy('seaice_fraction_to_freeze_seawater',val2, &
!                    'seaice_fraction_to_freeze_seawater_CT_freeze')
!call check_accuracy('seaice_fraction_to_freeze_seawater',val3, &
!                    'seaice_fraction_to_freeze_seawater_w_Ih')

!------------------------------------------------------------------------------
if (gsw_error_flag.eq.1) then
  print*
  print*; print*, 'Your installation of the Gibbs SeaWater (GSW) Oceanographic Toolbox has errors!'
else  
  print*
  print*; print*, 'Well done! The gsw_check_fuctions confirms that the Gibbs'
  print*; print*, 'SeaWater (GSW) Oceanographic Toolbox is installed correctly.'
  print*
endif

contains

    !--------------------------------------------------------------------------

    subroutine section_title (title)

    character (*), intent(in) :: title

    print *
    print *, "----------------------------------------------------------------------------"
    print *, title
    print *

    return
    end subroutine section_title

    !--------------------------------------------------------------------------

    subroutine check_accuracy (func_name, fvalue1, fvalue2, vprint)

    use gsw_mod_error_functions, only : gsw_error_limit

    implicit none

    character (*), intent(in) :: func_name
    real (r8), intent(in) :: fvalue1(:,:), fvalue2(:,:)
    logical, intent(in), optional :: vprint

    integer :: ndots, i, j, k, ik, jk
    real (r8) :: check_limit, dmax, drel
    real (r8) :: diff(size(fvalue1,1),size(fvalue1,2))
    character (len(func_name)+3) :: message
    character (4) :: errflg

    character (*), parameter :: att_name = 'computation_accuracy'
    character (*), parameter :: &
        dots = ' .............................................................'

    message = func_name

    diff = fvalue1 - fvalue2

    if (present(vprint)) then
        if (vprint) then
            print '(i3,3ES24.15)', ((i,fvalue1(i,j),fvalue2(i,j),diff(i,j),&
	            i=1,size(fvalue1,1)), j=1,size(fvalue1,2))
            print *
	end if
    end if

    if (any(fvalue1 .gt. gsw_error_limit)) then
        where (fvalue1 .gt. gsw_error_limit) diff = 0.0_r8
	errflg = ' (*)'
    else
	errflg = '    '
    end if
    ndots = 50 - len(trim(message))

    print *
    print '(2a,2es12.3)', trim(message), dots(:ndots), minval(diff), maxval(diff)

    return
    end subroutine check_accuracy

end

!--------------------------------------------------------------------------
