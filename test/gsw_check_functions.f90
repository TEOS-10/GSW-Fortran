program gsw_check_functions

use gsw_mod_kinds
use gsw_mod_netcdf
use gsw_mod_toolbox

use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit

implicit none

integer :: gsw_error_flag = 0

integer :: cast_m, cast_n, cast_mpres_m, cast_mpres_n, cast_ice_m
integer :: cast_ice_n, i, m, n

real (r8) :: saturation_fraction, pref

real (r8), dimension(:,:), allocatable :: ct, rt, sa, sk, sp, t, p, delta_p
real (r8), dimension(:,:), allocatable :: p_shallow, p_deep, lat, long

real (r8), dimension(:,:), allocatable :: ct_arctic, sa_arctic, t_arctic
real (r8), dimension(:,:), allocatable :: p_arctic, t_ice, w_ice
real (r8), dimension(:,:), allocatable :: sa_seaice, t_seaice, w_seaice
real (r8), dimension(:,:), allocatable :: lat_ice, long_ice

real (r8), dimension(:), allocatable :: lat_cast, long_cast

real (r8), dimension(:,:), allocatable :: value, check_value
real (r8), dimension(:,:), allocatable :: val1, val2, val3, val4, val5

real (r8), dimension(:,:), allocatable :: c, sr, sstar, pt, entropy
real (r8), dimension(:,:), allocatable :: h, ctf, tf, rho, diff
real (r8), dimension(:,:), allocatable :: ctf_poly, tf_poly, pt0
real (r8), dimension(:,:), allocatable :: sa_bulk, h_pot_bulk, h_bulk

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
allocate(delta_p(cast_m,cast_n))
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
allocate(sa_bulk(cast_ice_m,cast_ice_n))
allocate(h_pot_bulk(cast_ice_m,cast_ice_n))
allocate(h_bulk(cast_ice_m,cast_ice_n))

allocate(value(cast_m,cast_n))
allocate(check_value(cast_m,cast_n))
allocate(val1(cast_m,cast_n))
allocate(val2(cast_m,cast_n))
allocate(val3(cast_m,cast_n))
allocate(val4(cast_m,cast_n))
allocate(val5(cast_m,cast_n))

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
call ncdf_get_var("delta_p_chck_cast", var2=delta_p)
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

call ncdf_get_var("SA_bulk", var2=sa_bulk)
call ncdf_get_var("h_pot_bulk", var2=h_pot_bulk)
call ncdf_get_var("h_bulk", var2=h_bulk)

call ncdf_get_var("pr", var0=pref)

print*
print*,'============================================================================'
print*
print*,' Gibbs SeaWater (GSW) Oceanographic Toolbox of TEOS-10 (Fortran)'
print*
print*,'============================================================================'
print*
print*,' These are the check values for the subset of functions that have been '
print*
print*,' converted into FORTRAN 95 from the Gibbs SeaWater (GSW) Oceanographic '
print*
print*,' Toolbox of TEOS-10.'
print*

!------------------------------------------------------------------------------
call section_title('Practical Salinity, PSS-78')

c = gsw_c_from_sp(sp,t,p)
call check_accuracy('C_from_SP',c)

value = gsw_sp_from_c(c,t,p)
call check_accuracy('SP_from_C',value)

value = gsw_sp_from_sk(sk)
call check_accuracy('SP_from_SK',value)

!------------------------------------------------------------------------------
call section_title('Absolute Salinity, Preformed Salinity and Conservative Temperature')

value = gsw_sa_from_sp(sp,p,long,lat)
call check_accuracy('SA_from_SP',value)

value = gsw_sstar_from_sp(sp,p,long,lat)
call check_accuracy('Sstar_from_SP',value)

value = gsw_ct_from_t(sa,t,p)
call check_accuracy('CT_from_t',value)

!------------------------------------------------------------------------------
call section_title('Other conversions between Temperatures, Salinities, Entropy, Pressure and Height')

value = gsw_deltasa_from_sp(sp,p,long,lat)
call check_accuracy('deltaSA_from_SP',value)

sr = gsw_sr_from_sp(sp)
call check_accuracy('SR_from_SP',sr)

value = gsw_sp_from_sr(sr)
call check_accuracy('SP_from_SR',value)

value = gsw_sp_from_sa(sa,p,long,lat)
call check_accuracy('SP_from_SA',value)

sstar = gsw_sstar_from_sa(sa,p,long,lat)
call check_accuracy('Sstar_from_SA',sstar)

value = gsw_sa_from_sstar(sstar,p,long,lat)
call check_accuracy('SA_from_Sstar',value)

value = gsw_sp_from_sstar(sstar,p,long,lat)
call check_accuracy('SP_from_Sstar',value)

pt = gsw_pt_from_ct(sa,ct)
call check_accuracy('pt_from_CT',pt)

value = gsw_t_from_ct(sa,ct,p)
call check_accuracy('t_from_CT',value)

value = gsw_ct_from_pt(sa,pt)
call check_accuracy('CT_from_pt',value)

value = gsw_pt0_from_t(sa,t,p)
call check_accuracy('pt0_from_t',value)

value = gsw_pt_from_t(sa,t,p,pref)
call check_accuracy('pt_from_t',value)

value = gsw_z_from_p(p,lat)
call check_accuracy('z_from_p',value)

entropy = gsw_entropy_from_pt(sa,pt)
call check_accuracy('entropy_from_pt',entropy)

value = gsw_pt_from_entropy(sa,entropy)
call check_accuracy('pt_from_entropy',value)

value = gsw_CT_from_entropy(sa,entropy)
call check_accuracy('CT_from_entropy',value)

value = gsw_entropy_from_t(sa,t,p)
call check_accuracy('entropy_from_t',value)

value = gsw_adiabatic_lapse_rate_from_ct(sa,ct,p)
call check_accuracy('adiabatic_lapse_rate_from_CT',value)

!------------------------------------------------------------------------------
call section_title('Specific Volume, Density and Enthalpy')

value = gsw_specvol(sa,ct,p)
call check_accuracy('specvol',value)

value = gsw_alpha(sa,ct,p)
call check_accuracy('alpha',value)

value = gsw_beta(sa,ct,p)
call check_accuracy('beta',value)

value = gsw_alpha_on_beta(sa,ct,p)
call check_accuracy('alpha_on_beta',value)

call gsw_specvol_alpha_beta(sa,ct,p,val1,val2,val3)
call check_accuracy('specvol_alpha_beta',val1,'v_vab')
call check_accuracy('specvol_alpha_beta',val2,'alpha_vab')
call check_accuracy('specvol_alpha_beta',val3,'beta_vab')

call gsw_specvol_first_derivatives(sa,ct,p,val1,val2,val3)
call check_accuracy('specvol_first_derivatives',val1,'v_SA')
call check_accuracy('specvol_first_derivatives',val2,'v_CT')
call check_accuracy('specvol_first_derivatives',val3,'v_P')

call gsw_specvol_second_derivatives(sa,ct,p,val1,val2,val3,val4,val5)
call check_accuracy('specvol_second_derivatives',val1,'v_SA_SA')
call check_accuracy('specvol_second_derivatives',val2,'v_SA_CT')
call check_accuracy('specvol_second_derivatives',val3,'v_CT_CT')
call check_accuracy('specvol_second_derivatives',val4,'v_SA_P')
call check_accuracy('specvol_second_derivatives',val5,'v_CT_P')

call gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,val1,val2)
call check_accuracy('specvol_first_derivatives_wrt_enthalpy',val1,'v_SA_wrt_h')
call check_accuracy('specvol_first_derivatives_wrt_enthalpy',val2,'v_h')

call gsw_specvol_second_derivatives_wrt_enthalpy(sa,ct,p,val1,val2,val3)
call check_accuracy('specvol_second_derivatives_wrt_enthalpy',val1,'v_SA_SA_wrt_h')
call check_accuracy('specvol_second_derivatives_wrt_enthalpy',val2,'v_SA_h')
call check_accuracy('specvol_second_derivatives_wrt_enthalpy',val3,'v_h_h')

value = gsw_specvol_anom_standard(sa,ct,p)
call check_accuracy('specvol_anom_standard',value)

rho = gsw_rho(sa,ct,p)
call check_accuracy('rho',rho)

call gsw_rho_alpha_beta(sa,ct,p,val1,val2,val3)
call check_accuracy('rho_alpha_beta',val1,'rho_rab')
call check_accuracy('rho_alpha_beta',val2,'alpha_rab')
call check_accuracy('rho_alpha_beta',val3,'beta_rab')

call gsw_rho_first_derivatives(sa,ct,p,val1,val2,val3)
call check_accuracy('rho_first_derivatives',val1,'rho_SA')
call check_accuracy('rho_first_derivatives',val2,'rho_CT')
call check_accuracy('rho_first_derivatives',val3,'rho_P')

call gsw_rho_second_derivatives(sa,ct,p,val1,val2,val3,val4,val5)
call check_accuracy('rho_second_derivatives',val1,'rho_SA_SA')
call check_accuracy('rho_second_derivatives',val2,'rho_SA_CT')
call check_accuracy('rho_second_derivatives',val3,'rho_CT_CT')
call check_accuracy('rho_second_derivatives',val4,'rho_SA_P')
call check_accuracy('rho_second_derivatives',val5,'rho_CT_P')

call gsw_rho_first_derivatives_wrt_enthalpy(sa,ct,p,val1,val2)
call check_accuracy('rho_first_derivatives_wrt_enthalpy',val1,'rho_SA_wrt_h')
call check_accuracy('rho_first_derivatives_wrt_enthalpy',val2,'rho_h')

call gsw_rho_second_derivatives_wrt_enthalpy(sa,ct,p,val1,val2,val3)
call check_accuracy('rho_second_derivatives_wrt_enthalpy',val1,'rho_SA_SA_wrt_h')
call check_accuracy('rho_second_derivatives_wrt_enthalpy',val2,'rho_SA_h')
call check_accuracy('rho_second_derivatives_wrt_enthalpy',val3,'rho_h_h')

value = gsw_sigma0(sa,ct)
call check_accuracy('sigma0',value)

value = gsw_sigma1(sa,ct)
call check_accuracy('sigma1',value)

value = gsw_sigma2(sa,ct)
call check_accuracy('sigma2',value)

value = gsw_sigma3(sa,ct)
call check_accuracy('sigma3',value)

value = gsw_sigma4(sa,ct)
call check_accuracy('sigma4',value)

value = gsw_sound_speed(sa,ct,p)
call check_accuracy('sound_speed',value)

value = gsw_kappa(sa,ct,p)
call check_accuracy('kappa',value)

value = gsw_cabbeling(sa,ct,p)
call check_accuracy('cabbeling',value)

value = gsw_thermobaric(sa,ct,p)
call check_accuracy('thermobaric',value)

value = gsw_sa_from_rho(rho,ct,p)
call check_accuracy('SA_from_rho',value)

call gsw_ct_from_rho(rho,sa,p,value)
call check_accuracy('CT_from_rho',value)

value = gsw_ct_maxdensity(sa,p)
call check_accuracy('CT_maxdensity',value)

value = gsw_internal_energy(sa,ct,p)
call check_accuracy('internal_energy',value)

h = gsw_enthalpy(sa,ct,p)
call check_accuracy('enthalpy',h)

value = gsw_enthalpy_diff(sa,ct,p_shallow,p_deep)
call check_accuracy('enthalpy_diff',value)

value = gsw_ct_from_enthalpy(sa,h,p)
call check_accuracy('CT_from_enthalpy',value)

value = gsw_dynamic_enthalpy(sa,ct,p)
call check_accuracy('dynamic_enthalpy',value)

call gsw_enthalpy_first_derivatives(sa,ct,p,val1,val2)
call check_accuracy('enthalpy_first_derivatives',val1,'h_SA')
call check_accuracy('enthalpy_first_derivatives',val2,'h_CT')

call gsw_enthalpy_second_derivatives(sa,ct,p,val1,val2,val3)
call check_accuracy('enthalpy_second_derivatives',val1,'h_SA_SA')
call check_accuracy('enthalpy_second_derivatives',val2,'h_SA_CT')
call check_accuracy('enthalpy_second_derivatives',val3,'h_CT_CT')

!------------------------------------------------------------------------------
call section_title('Derivatives of entropy, CT and pt')

call gsw_ct_first_derivatives(sa,pt,val1,val2)
call check_accuracy('CT_first_derivatives',val1,'CT_SA')
call check_accuracy('CT_first_derivatives',val2,'CT_pt')

call gsw_ct_second_derivatives(sa,pt,val1,val2,val3)
call check_accuracy('CT_second_derivatives',val1,'CT_SA_SA')
call check_accuracy('CT_second_derivatives',val2,'CT_SA_pt')
call check_accuracy('CT_second_derivatives',val3,'CT_pt_pt')

call gsw_entropy_first_derivatives(sa,ct,val1,val2)
call check_accuracy('entropy_first_derivatives',val1,'eta_SA')
call check_accuracy('entropy_first_derivatives',val2,'eta_CT')

call gsw_entropy_second_derivatives(sa,ct,val1,val2,val3)
call check_accuracy('entropy_second_derivatives',val1,'eta_SA_SA')
call check_accuracy('entropy_second_derivatives',val2,'eta_SA_CT')
call check_accuracy('entropy_second_derivatives',val3,'eta_CT_CT')

call gsw_pt_first_derivatives(sa,ct,val1,val2)
call check_accuracy('pt_first_derivatives',val1,'pt_SA')
call check_accuracy('pt_first_derivatives',val2,'pt_CT')

call gsw_pt_second_derivatives(sa,ct,val1,val2,val3)
call check_accuracy('pt_second_derivatives',val1,'pt_SA_SA')
call check_accuracy('pt_second_derivatives',val2,'pt_SA_CT')
call check_accuracy('pt_second_derivatives',val3,'pt_CT_CT')

!------------------------------------------------------------------------------
call section_title('Freezing temperatures')

saturation_fraction = 0.5_r8

ctf = gsw_ct_freezing(sa,p,saturation_fraction)
call check_accuracy('CT_freezing',ctf)

ctf_poly = gsw_ct_freezing_poly(sa,p,saturation_fraction)
call check_accuracy('CT_freezing_poly',ctf_poly)

tf = gsw_t_freezing(sa,p,saturation_fraction)
call check_accuracy('t_freezing',tf)

tf_poly = gsw_t_freezing_poly(sa,p,saturation_fraction)
call check_accuracy('t_freezing_poly',tf_poly)

value = gsw_pot_enthalpy_ice_freezing(sa,p)
call check_accuracy('pot_enthalpy_ice_freezing',value)

value = gsw_pot_enthalpy_ice_freezing_poly(sa,p)
call check_accuracy('pot_enthalpy_ice_freezing_poly',value)

value = gsw_sa_freezing_from_ct(ctf,p,saturation_fraction)
call check_accuracy('SA_freezing_from_CT',value)

value = gsw_sa_freezing_from_ct_poly(ctf_poly,p,saturation_fraction)
call check_accuracy('SA_freezing_from_CT_poly',value)

value = gsw_sa_freezing_from_t(tf,p,saturation_fraction)
call check_accuracy('SA_freezing_from_t',value)

value = gsw_sa_freezing_from_t_poly(tf_poly,p,saturation_fraction)
call check_accuracy('SA_freezing_from_t_poly',value)

call gsw_ct_freezing_first_derivatives(sa,p,saturation_fraction,val1,val2)
call check_accuracy('CT_freezing_first_derivatives',val1,'CTfreezing_SA')
call check_accuracy('CT_freezing_first_derivatives',val2,'CTfreezing_P')

call gsw_ct_freezing_first_derivatives_poly(sa,p,saturation_fraction,val4,val5)
call check_accuracy('CT_freezing_first_derivatives_poly',val4,'CTfreezing_SA_poly')
call check_accuracy('CT_freezing_first_derivatives_poly',val5,'CTfreezing_P_poly')

call gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,val1,val2)
call check_accuracy('t_freezing_first_derivatives',val1,'tfreezing_SA')
call check_accuracy('t_freezing_first_derivatives',val2,'tfreezing_P')

call gsw_t_freezing_first_derivatives_poly(sa,p,saturation_fraction,val4,val5)
call check_accuracy('t_freezing_first_derivatives_poly',val4,'tfreezing_SA_poly')
call check_accuracy('t_freezing_first_derivatives_poly',val5,'tfreezing_P_poly')

call gsw_pot_enthalpy_ice_freezing_first_derivatives(sa,p,val1,val2)
call check_accuracy('pot_enthalpy_ice_freezing_first_derivatives',val1, &
                    'pot_enthalpy_ice_freezing_SA')
call check_accuracy('pot_enthalpy_ice_freezing_first_derivatives',val2, &
                    'pot_enthalpy_ice_freezing_P')

call gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(sa,p,val1,val2)
call check_accuracy('pot_enthalpy_ice_freezing_first_derivatives_poly',val1, &
                    'pot_enthalpy_ice_freezing_SA_poly')
call check_accuracy('pot_enthalpy_ice_freezing_first_derivatives_poly',val2, &
                    'pot_enthalpy_ice_freezing_P_poly')

!------------------------------------------------------------------------------
call section_title('Isobaric Melting Enthalpy and Isobaric Evaporation Enthalpy')

value = gsw_latentheat_melting(sa,p)
call check_accuracy('latentheat_melting',value)

value = gsw_latentheat_evap_ct(sa,ct)
call check_accuracy('latentheat_evap_CT',value)

value = gsw_latentheat_evap_t(sa,t)
call check_accuracy('latentheat_evap_t',value)

!------------------------------------------------------------------------------
call section_title('Planet Earth properties')

value = gsw_grav(lat,p)
call check_accuracy('grav',value)

!------------------------------------------------------------------------------
call section_title('Density and enthalpy in terms of CT, derived from the exact Gibbs function')

value = gsw_enthalpy_ct_exact(sa,ct,p)
call check_accuracy('enthalpy_CT_exact',value)

call gsw_enthalpy_first_derivatives_ct_exact(sa,ct,p,val1,val2)
call check_accuracy('enthalpy_first_derivatives_CT_exact',val1,'h_SA_CT_exact')
call check_accuracy('enthalpy_first_derivatives_CT_exact',val2,'h_CT_CT_exact')

call gsw_enthalpy_second_derivatives_ct_exact(sa,ct,p,val1,val2,val3)
call check_accuracy('enthalpy_second_derivatives_CT_exact',val1, &
                    'h_SA_SA_CT_exact')
call check_accuracy('enthalpy_second_derivatives_CT_exact',val2, &
                    'h_SA_CT_CT_exact')
call check_accuracy('enthalpy_second_derivatives_CT_exact',val3, &
                    'h_CT_CT_CT_exact')

!------------------------------------------------------------------------------
call section_title('Basic thermodynamic properties in terms of in-situ t, based on the exact Gibbs function')

value = gsw_rho_t_exact(sa,t,p)
call check_accuracy('rho_t_exact',value)

value = gsw_pot_rho_t_exact(sa,t,p,pref)
call check_accuracy('pot_rho_t_exact',value)

value = gsw_alpha_wrt_t_exact(sa,t,p)
call check_accuracy('alpha_wrt_t_exact',value)

value = gsw_beta_const_t_exact(sa,t,p)
call check_accuracy('beta_const_t_exact',value)

value = gsw_specvol_t_exact(sa,t,p)
call check_accuracy('specvol_t_exact',value)

value = gsw_sound_speed_t_exact(sa,t,p)
call check_accuracy('sound_speed_t_exact',value)

value = gsw_kappa_t_exact(sa,t,p)
call check_accuracy('kappa_t_exact',value)

value = gsw_enthalpy_t_exact(sa,t,p)
call check_accuracy('enthalpy_t_exact',value)

call gsw_ct_first_derivatives_wrt_t_exact(sa,t,p,val1,val2,val3)
call check_accuracy('CT_first_derivatives_wrt_t_exact',val1,'CT_SA_wrt_t')
call check_accuracy('CT_first_derivatives_wrt_t_exact',val2,'CT_T_wrt_t')
call check_accuracy('CT_first_derivatives_wrt_t_exact',val3,'CT_P_wrt_t')

value = gsw_chem_potential_water_t_exact(sa,t,p)
call check_accuracy('chem_potential_water_t_exact',value)

value = gsw_t_deriv_chem_potential_water_t_exact(sa,t,p)
call check_accuracy('t_deriv_chem_potential_water_t_exact',value)

value = gsw_dilution_coefficient_t_exact(sa,t,p)
call check_accuracy('dilution_coefficient_t_exact',value)

!------------------------------------------------------------------------------
call section_title('Library functions of the GSW Toolbox')

value = gsw_deltasa_atlas(p,long,lat)
call check_accuracy('deltaSA_atlas',value)

value = gsw_fdelta(p,long,lat)
call check_accuracy('Fdelta',value)

!------------------------------------------------------------------------------
call section_title('Water column properties, based on the 75-term polynomial for specific volume')

deallocate(check_value,val1,val2,val3)
allocate(check_value(cast_mpres_m,cast_mpres_n))
allocate(val1(cast_mpres_m,cast_mpres_n))
allocate(val2(cast_mpres_m,cast_mpres_n))
allocate(val3(cast_mpres_m,cast_mpres_n))

do i = 1, cast_mpres_n
    call gsw_nsquared(sa(:,i),ct(:,i),p(:,i),lat(:,i),val1(:,i),val2(:,i))
end do
call check_accuracy('Nsquared',val1,'n2')
call check_accuracy('Nsquared',val2,'p_mid_n2')

do i = 1, cast_mpres_n
    call gsw_turner_rsubrho(sa(:,i),ct(:,i),p(:,i),val1(:,i),val2(:,i), &
                            val3(:,i))
end do
call check_accuracy('Turner_Rsubrho',val1,'Tu')
call check_accuracy('Turner_Rsubrho',val2,'Rsubrho')
call check_accuracy('Turner_Rsubrho',val3,'p_mid_TuRsr')

do i = 1, cast_mpres_n
    call gsw_ipv_vs_fnsquared_ratio(sa(:,i),ct(:,i),p(:,i),pref,val1(:,i), &
                                    val2(:,i))
end do
call check_accuracy('IPV_vs_fNsquared_ratio',val1,'IPVfN2')
call check_accuracy('IPV_vs_fNsquared_ratio',val2,'p_mid_IPVfN2')

do i = 1, cast_mpres_n
    n = count(sa(:,i) .eq. sa(:,i))  ! check for NaN's
    val1(:,i) = gsw_geo_strf_dyn_height(sa(:n,i),ct(:n,i),p(:n,i),pref)
    if (n .lt. cast_mpres_m) val1(n+1:,i) = sa(n+1:cast_mpres_m,i)
end do
call check_accuracy('geo_strf_dyn_height',val1)

do i = 1, cast_mpres_n
    call gsw_geo_strf_dyn_height_pc(sa(:,i),ct(:,i),delta_p(:,i), &
                                    val1(:,i),val2(:,i))
end do
call check_accuracy('geo_strf_dyn_height_pc',val1,'geo_strf_dyn_height_pc')
call check_accuracy('geo_strf_dyn_height_pc',val2,'geo_strf_dyn_height_pc_p_mid')

!------------------------------------------------------------------------------
call section_title('Thermodynamic properties of ice Ih')

deallocate(value,check_value,h)
allocate(value(cast_ice_m,cast_ice_n))
allocate(check_value(cast_ice_m,cast_ice_n))
allocate(h(cast_ice_m,cast_ice_n))

value = gsw_rho_ice(t_seaice,p_arctic)
call check_accuracy('rho_ice',value)

value = gsw_alpha_wrt_t_ice(t_seaice,p_arctic)
call check_accuracy('alpha_wrt_t_ice',value)

value = gsw_specvol_ice(t_seaice,p_arctic)
call check_accuracy('specvol_ice',value)

value = gsw_pressure_coefficient_ice(t_seaice,p_arctic)
call check_accuracy('pressure_coefficient_ice',value)

value = gsw_sound_speed_ice(t_seaice,p_arctic)
call check_accuracy('sound_speed_ice',value)

value = gsw_kappa_ice(t_seaice,p_arctic)
call check_accuracy('kappa_ice',value)

value = gsw_kappa_const_t_ice(t_seaice,p_arctic)
call check_accuracy('kappa_const_t_ice',value)

value = gsw_internal_energy_ice(t_seaice,p_arctic)
call check_accuracy('internal_energy_ice',value)

value = gsw_enthalpy_ice(t_seaice,p_arctic)
call check_accuracy('enthalpy_ice',value)

value = gsw_entropy_ice(t_seaice,p_arctic)
call check_accuracy('entropy_ice',value)

value = gsw_cp_ice(t_seaice,p_arctic)
call check_accuracy('cp_ice',value)

value = gsw_chem_potential_water_ice(t_seaice,p_arctic)
call check_accuracy('chem_potential_water_ice',value)

value = gsw_helmholtz_energy_ice(t_seaice,p_arctic)
call check_accuracy('Helmholtz_energy_ice',value)

value = gsw_adiabatic_lapse_rate_ice(t_seaice,p_arctic)
call check_accuracy('adiabatic_lapse_rate_ice',value)

pt0 = gsw_pt0_from_t_ice(t_seaice,p_arctic)
call check_accuracy('pt0_from_t_ice',pt0)

value = gsw_pt_from_t_ice(t_seaice,p_arctic,pref)
call check_accuracy('pt_from_t_ice',value)

value = gsw_t_from_pt0_ice(pt0,p_arctic)
call check_accuracy('t_from_pt0_ice',value)

h = gsw_pot_enthalpy_from_pt_ice(pt0)
call check_accuracy('pot_enthalpy_from_pt_ice',h)

value = gsw_pt_from_pot_enthalpy_ice(h)
call check_accuracy('pt_from_pot_enthalpy_ice',value)

h = gsw_pot_enthalpy_from_pt_ice_poly(pt0)
call check_accuracy('pot_enthalpy_from_pt_ice_poly',h)

value = gsw_pt_from_pot_enthalpy_ice_poly(h)
call check_accuracy('pt_from_pot_enthalpy_ice_poly',value)

saturation_fraction = 0.5_r8

value = gsw_pressure_freezing_ct(sa_arctic,ct_arctic-1.0_r8,saturation_fraction)
call check_accuracy('pressure_freezing_CT',value)

!------------------------------------------------------------------------------
call section_title('Thermodynamic interaction between ice and seawater')

deallocate(val1,val2,val3)
allocate(val1(cast_ice_m,cast_ice_n))
allocate(val2(cast_ice_m,cast_ice_n))
allocate(val3(cast_ice_m,cast_ice_n))

value = gsw_melting_ice_sa_ct_ratio(sa_arctic,ct_arctic,p_arctic,t_ice)
call check_accuracy('melting_ice_SA_CT_ratio',value)

value = gsw_melting_ice_sa_ct_ratio_poly(sa_arctic,ct_arctic,p_arctic,t_ice)
call check_accuracy('melting_ice_SA_CT_ratio_poly',value)

value = gsw_melting_ice_equilibrium_sa_ct_ratio(sa_arctic,p_arctic)
call check_accuracy('melting_ice_equilibrium_SA_CT_ratio',value)

value = gsw_melting_ice_equilibrium_sa_ct_ratio_poly(sa_arctic,p_arctic)
call check_accuracy('melting_ice_equilibrium_SA_CT_ratio_poly',value)

call gsw_melting_ice_into_seawater(sa_arctic,ct_arctic+0.1_r8,p_arctic,w_ice,t_ice, &
                                   val1,val2,val3)
call check_accuracy('melting_ice_into_seawater',val1, &
                    'melting_ice_into_seawater_SA_final')
call check_accuracy('melting_ice_into_seawater',val2, &
                    'melting_ice_into_seawater_CT_final')
!call check_accuracy('melting_ice_into_seawater',val3, &
!                    'melting_ice_into_seawater_w_Ih')

call gsw_ice_fraction_to_freeze_seawater(sa_arctic,ct_arctic,p_arctic,t_ice, &
                                         val1,val2,val3)
call check_accuracy('ice_fraction_to_freeze_seawater',val1, &
                    'ice_fraction_to_freeze_seawater_SA_freeze')
call check_accuracy('ice_fraction_to_freeze_seawater',val2, &
                    'ice_fraction_to_freeze_seawater_CT_freeze')
call check_accuracy('ice_fraction_to_freeze_seawater',val3, &
                    'ice_fraction_to_freeze_seawater_w_Ih')

call gsw_frazil_ratios_adiabatic(sa_arctic,p_arctic,w_ice,val1,val2,val3)
call check_accuracy('frazil_ratios_adiabatic',val1,'dSA_dCT_frazil')
call check_accuracy('frazil_ratios_adiabatic',val2,'dSA_dP_frazil')
call check_accuracy('frazil_ratios_adiabatic',val3,'dCT_dP_frazil')

call gsw_frazil_ratios_adiabatic_poly(sa_arctic,p_arctic,w_ice,val1,val2,val3)
call check_accuracy('frazil_ratios_adiabatic_poly',val1,'dSA_dCT_frazil_poly')
call check_accuracy('frazil_ratios_adiabatic_poly',val2,'dSA_dP_frazil_poly')
call check_accuracy('frazil_ratios_adiabatic_poly',val3,'dCT_dP_frazil_poly')

call gsw_frazil_properties_potential(sa_bulk,h_pot_bulk,p_arctic,val1,val2,val3)
call check_accuracy('frazil_properties_potential',val1, &
                    'frazil_properties_potential_SA_final')
call check_accuracy('frazil_properties_potential',val2, &
                    'frazil_properties_potential_CT_final')
call check_accuracy('frazil_properties_potential',val3, &
                    'frazil_properties_potential_w_Ih_final')

call gsw_frazil_properties_potential_poly(sa_bulk,h_pot_bulk,p_arctic,val1, &
                                          val2,val3)
call check_accuracy('frazil_properties_potential_poly',val1, &
                    'frazil_properties_potential_poly_SA_final')
call check_accuracy('frazil_properties_potential_poly',val2, &
                    'frazil_properties_potential_poly_CT_final')
call check_accuracy('frazil_properties_potential_poly',val3, &
                    'frazil_properties_potential_poly_w_Ih_final')

call gsw_frazil_properties(sa_bulk,h_bulk,p_arctic,val1,val2,val3)
call check_accuracy('frazil_properties',val1,'frazil_properties_SA_final')
call check_accuracy('frazil_properties',val2,'frazil_properties_CT_final')
call check_accuracy('frazil_properties',val3,'frazil_properties_w_Ih_final')

!------------------------------------------------------------------------------
call section_title('Thermodynamic interaction between seaice and seawater')

value = gsw_melting_seaice_sa_ct_ratio(sa_arctic,ct_arctic,p_arctic, &
                                       sa_seaice,t_seaice)
call check_accuracy('melting_seaice_SA_CT_ratio',value)

value = gsw_melting_seaice_sa_ct_ratio_poly(sa_arctic,ct_arctic,p_arctic, &
                                            sa_seaice,t_seaice)
call check_accuracy('melting_seaice_SA_CT_ratio_poly',value)

value = gsw_melting_seaice_equilibrium_sa_ct_ratio(sa_arctic,p_arctic)
call check_accuracy('melting_seaice_equilibrium_SA_CT_ratio',value)

value = gsw_melting_seaice_equilibrium_sa_ct_ratio_poly(sa_arctic,p_arctic)
call check_accuracy('melting_seaice_equilibrium_SA_CT_ratio_poly',value)

call gsw_melting_seaice_into_seawater(sa_arctic,ct_arctic,p_arctic, &
                                      w_seaice,sa_seaice,t_seaice,val1,val2)
call check_accuracy('melting_seaice_into_seawater',val1, &
                    'melting_seaice_into_seawater_SA_final')
call check_accuracy('melting_seaice_into_seawater',val2, &
                    'melting_seaice_into_seawater_CT_final')

call gsw_seaice_fraction_to_freeze_seawater(sa_arctic,ct_arctic,p_arctic, &
                                            sa_seaice,t_seaice,val1,val2,val3)
call check_accuracy('seaice_fraction_to_freeze_seawater',val1, &
                    'seaice_fraction_to_freeze_seawater_SA_freeze')
call check_accuracy('seaice_fraction_to_freeze_seawater',val2, &
                    'seaice_fraction_to_freeze_seawater_CT_freeze')
call check_accuracy('seaice_fraction_to_freeze_seawater',val3, &
                    'seaice_fraction_to_freeze_seawater_w_Ih')

!------------------------------------------------------------------------------
if (gsw_error_flag.eq.1) then
  print*
  print*
  print*, 'Your installation of the Gibbs SeaWater (GSW) Oceanographic Toolbox has errors!'
else  
  print*
  print*
  print*, 'Well done! The gsw_check_fuctions confirms that the'
  print*
  print*, 'Gibbs SeaWater (GSW) Oceanographic Toolbox is installed correctly.'
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

    subroutine check_accuracy (func_name, fvalue, var_name, vprint)

    use gsw_mod_error_functions, only : gsw_error_limit

    implicit none

    character (*), intent(in) :: func_name
    real (r8), intent(in) :: fvalue(:,:)
    character (*), intent(in), optional :: var_name
    logical, intent(in), optional :: vprint

    integer :: ndots, i, j, k, ik, jk
    real (r8) :: check_limit, dmax, drel
    real (r8) :: diff(size(fvalue,1),size(fvalue,2))
    character (80) :: message
    character (4) :: errflg

    character (*), parameter :: att_name = 'computation_accuracy'
    character (*), parameter :: &
        dots = ' .............................................................'

    if (present(var_name)) then

        call ncdf_get_var_att(var_name, check_value, att_name, check_limit)

	if (len(func_name)+len(var_name).gt.55) then
	    k = len(func_name) + len(var_name) - 55
	    message = func_name // ' (..' // var_name(k:) // ')'
	else
	    message = func_name // ' (' // var_name // ')'
	end if

    else

        call ncdf_get_var_att(func_name, check_value, att_name, check_limit)
	message = func_name

    end if

    diff = abs(fvalue - check_value)

    if (present(vprint)) then
        if (vprint) then
	    print *, "Limit =", check_limit
            print '(i3,3ES24.15)', ((i,fvalue(i,j),check_value(i,j),diff(i,j), &
	            i=1,size(fvalue,1)), j=1,size(fvalue,2))
            print *
	end if
    end if

    if (any(fvalue .gt. gsw_error_limit)) then
        where (fvalue .gt. gsw_error_limit) diff = 0.0_r8
	errflg = ' (*)'
    else
	errflg = '    '
    end if
    ndots = 65 - len(trim(message))

    if (any(diff .gt. check_limit)) then
        gsw_error_flag = 1
	dmax = 0.0_r8
	do i = 1, size(fvalue,1)
	    do j = 1, size(fvalue,2)
	         if (diff(i,j) .gt. dmax) then
		     dmax = diff(i,j)
		     ik = i
		     jk = j
		 end if
	    end do
	end do
	drel = dmax*100.0_r8/abs(fvalue(ik,jk))
        print *, trim(message), dots(:ndots-3), ' << failed >>'
	print *
	print *, "  Max. difference =", dmax, ", limit =", check_limit
	print *, "  Max. diff (rel) =", drel, ", limit =", check_limit*100.0_r8/abs(fvalue(ik,jk))
	print *
    else
        print *, trim(message), dots(:ndots), ' passed', errflg
    endif

    return
    end subroutine check_accuracy

end

!--------------------------------------------------------------------------
