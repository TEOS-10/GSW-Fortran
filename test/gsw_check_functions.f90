program gsw_check_functions

use gsw_mod_kinds
use gsw_mod_toolbox

use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit
use gsw_mod_check_data

implicit none

integer :: gsw_error_flag = 0

integer :: i, n

real (r8) :: saturation_fraction

real (r8), dimension(:,:), allocatable :: lat, long


real (r8), dimension(:,:), allocatable :: value, check_value
real (r8), dimension(:,:), allocatable :: val1, val2, val3, val4, val5
real (r8), dimension(:,:), allocatable :: val6, val7, val8


real (r8), dimension(:,:), allocatable :: c, sr, sstar, pt, entropy
real (r8), dimension(:,:), allocatable :: h, ctf, tf, diff, z
real (r8), dimension(:,:), allocatable :: ctf_poly, tf_poly, pt0


allocate(lat(cast_m,cast_n))
allocate(long(cast_m,cast_n))
allocate(pt0(cast_ice_m,cast_ice_n))

allocate(value(cast_m,cast_n))
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
allocate(tf(cast_m,cast_n))
allocate(ctf_poly(cast_m,cast_n))
allocate(tf_poly(cast_m,cast_n))
allocate(h(cast_m,cast_n))
allocate(diff(cast_m,cast_n))
allocate(z(cast_m,cast_n))

do i = 1, cast_n
    lat(:,i) = lat_cast(i)
    long(:,i) = long_cast(i)
end do


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
call check_accuracy('C_from_SP',c,C_from_SP,C_from_SP_ca)

value = gsw_sp_from_c(c,t,p)
call check_accuracy('SP_from_C',value,SP_from_C,SP_from_C_ca)

value = gsw_sp_from_sk(sk)
call check_accuracy('SP_from_SK',value,SP_from_SK,SP_from_SK_ca)

!------------------------------------------------------------------------------
call section_title('Absolute Salinity, Preformed Salinity and Conservative Temperature')

value = gsw_sa_from_sp(sp,p,long,lat)
call check_accuracy('SA_from_SP',value,SA_from_SP,SA_from_SP_ca)

value = gsw_sstar_from_sp(sp,p,long,lat)
call check_accuracy('Sstar_from_SP',value,Sstar_from_SP,Sstar_from_SP_ca)

value = gsw_ct_from_t(sa,t,p)
call check_accuracy('CT_from_t',value,CT_from_t,CT_from_t_ca)

!------------------------------------------------------------------------------
call section_title('Other conversions between Temperatures, Salinities, Entropy, Pressure and Height')

value = gsw_deltasa_from_sp(sp,p,long,lat)
call check_accuracy('deltaSA_from_SP',value,deltaSA_from_SP,deltaSA_from_SP_ca)

sr = gsw_sr_from_sp(sp)
call check_accuracy('SR_from_SP',sr,SR_from_SP,SR_from_SP_ca)

value = gsw_sp_from_sr(sr)
call check_accuracy('SP_from_SR',value,SP_from_SR,SP_from_SR_ca)

value = gsw_sp_from_sa(sa,p,long,lat)
call check_accuracy('SP_from_SA',value,SP_from_SA,SP_from_SA_ca)

sstar = gsw_sstar_from_sa(sa,p,long,lat)
call check_accuracy('Sstar_from_SA',sstar,Sstar_from_SA,Sstar_from_SA_ca)

value = gsw_sa_from_sstar(sstar,p,long,lat)
call check_accuracy('SA_from_Sstar',value,SA_from_Sstar,SA_from_Sstar_ca)

value = gsw_sp_from_sstar(sstar,p,long,lat)
call check_accuracy('SP_from_Sstar',value,SP_from_Sstar,SP_from_Sstar_ca)

pt = gsw_pt_from_ct(sa,ct)
call check_accuracy('pt_from_CT',pt,pt_from_CT,pt_from_CT_ca)

value = gsw_t_from_ct(sa,ct,p)
call check_accuracy('t_from_CT',value,t_from_CT,t_from_CT_ca)

value = gsw_ct_from_pt(sa,pt)
call check_accuracy('CT_from_pt',value,CT_from_pt,CT_from_pt_ca)

value = gsw_pt0_from_t(sa,t,p)
call check_accuracy('pt0_from_t',value,pt0_from_t,pt0_from_t_ca)

value = gsw_pt_from_t(sa,t,p,pref(1))
call check_accuracy('pt_from_t',value,pt_from_t,pt_from_t_ca)

z = gsw_z_from_p(p,lat)
call check_accuracy('z_from_p',z,z_from_p,z_from_p_ca)

value = gsw_p_from_z(z,lat)
call check_accuracy('p_from_z',value,p_from_z,p_from_z_ca)

entropy = gsw_entropy_from_pt(sa,pt)
call check_accuracy('entropy_from_pt',entropy,entropy_from_pt,entropy_from_pt_ca)

value = gsw_pt_from_entropy(sa,entropy)
call check_accuracy('pt_from_entropy',value,pt_from_entropy,pt_from_entropy_ca)

value = gsw_CT_from_entropy(sa,entropy)
call check_accuracy('CT_from_entropy',value,CT_from_entropy,CT_from_entropy_ca)

value = gsw_entropy_from_t(sa,t,p)
call check_accuracy('entropy_from_t',value,entropy_from_t,entropy_from_t_ca)

value = gsw_adiabatic_lapse_rate_from_ct(sa,ct,p)
call check_accuracy('adiabatic_lapse_rate_from_CT',value,adiabatic_lapse_rate_from_CT,adiabatic_lapse_rate_from_CT_ca)

!------------------------------------------------------------------------------
call section_title('Specific Volume, Density and Enthalpy')

value = gsw_specvol(sa,ct,p)
call check_accuracy('specvol',value,specvol,specvol_ca)

value = gsw_alpha(sa,ct,p)
call check_accuracy('alpha',value,alpha,alpha_ca)

value = gsw_beta(sa,ct,p)
call check_accuracy('beta',value,beta,beta_ca)

value = gsw_alpha_on_beta(sa,ct,p)
call check_accuracy('alpha_on_beta',value,alpha_on_beta,alpha_on_beta_ca)

call gsw_specvol_alpha_beta(sa,ct,p,val1,val2,val3)
call check_accuracy('specvol_alpha_beta',val1,v_vab,v_vab_ca,'v_vab')
call check_accuracy('specvol_alpha_beta',val2,alpha_vab,alpha_vab_ca,'alpha_vab')
call check_accuracy('specvol_alpha_beta',val3,beta_vab,beta_vab_ca,'beta_vab')

call gsw_specvol_first_derivatives(sa,ct,p,val1,val2,val3)
call check_accuracy('specvol_first_derivatives',val1,v_SA,v_SA_ca,'v_SA')
call check_accuracy('specvol_first_derivatives',val2,v_CT,v_CT_ca,'v_CT')
call check_accuracy('specvol_first_derivatives',val3,v_P,v_P_ca,'v_P')

call gsw_specvol_second_derivatives(sa,ct,p,val1,val2,val3,val4,val5)
call check_accuracy('specvol_second_derivatives',val1,v_SA_SA,v_SA_SA_ca,'v_SA_SA')
call check_accuracy('specvol_second_derivatives',val2,v_SA_CT,v_SA_CT_ca,'v_SA_CT')
call check_accuracy('specvol_second_derivatives',val3,v_CT_CT,v_CT_CT_ca,'v_CT_CT')
call check_accuracy('specvol_second_derivatives',val4,v_SA_P,v_SA_P_ca,'v_SA_P')
call check_accuracy('specvol_second_derivatives',val5,v_CT_P,v_CT_P_ca,'v_CT_P')

call gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,val1,val2)
call check_accuracy('specvol_first_derivatives_wrt_enthalpy',val1,v_SA_wrt_h,v_SA_wrt_h_ca,'v_SA_wrt_h')
call check_accuracy('specvol_first_derivatives_wrt_enthalpy',val2,v_h,v_h_ca,'v_h')

call gsw_specvol_second_derivatives_wrt_enthalpy(sa,ct,p,val1,val2,val3)
call check_accuracy('specvol_second_derivatives_wrt_enthalpy',val1,v_SA_SA_wrt_h,v_SA_SA_wrt_h_ca,'v_SA_SA_wrt_h')
call check_accuracy('specvol_second_derivatives_wrt_enthalpy',val2,v_SA_h,v_SA_h_ca,'v_SA_h')
call check_accuracy('specvol_second_derivatives_wrt_enthalpy',val3,v_h_h,v_h_h_ca,'v_h_h')

value = gsw_specvol_anom_standard(sa,ct,p)
call check_accuracy('specvol_anom_standard',value,specvol_anom_standard,specvol_anom_standard_ca)

rho = gsw_rho(sa,ct,p)
call check_accuracy('rho',rho,rho,rho_ca)

call gsw_rho_alpha_beta(sa,ct,p,val1,val2,val3)
call check_accuracy('rho_alpha_beta',val1,rho_rab,rho_rab_ca,'rho_rab')
call check_accuracy('rho_alpha_beta',val2,alpha_rab,alpha_rab_ca,'alpha_rab')
call check_accuracy('rho_alpha_beta',val3,beta_rab,beta_rab_ca,'beta_rab')

call gsw_rho_first_derivatives(sa,ct,p,val1,val2,val3)
call check_accuracy('rho_first_derivatives',val1,rho_SA,rho_SA_ca,'rho_SA')
call check_accuracy('rho_first_derivatives',val2,rho_CT,rho_CT_ca,'rho_CT')
call check_accuracy('rho_first_derivatives',val3,rho_P,rho_P_ca,'rho_P')

call gsw_rho_second_derivatives(sa,ct,p,val1,val2,val3,val4,val5)
call check_accuracy('rho_second_derivatives',val1,rho_SA_SA,rho_SA_SA_ca,'rho_SA_SA')
call check_accuracy('rho_second_derivatives',val2,rho_SA_CT,rho_SA_CT_ca,'rho_SA_CT')
call check_accuracy('rho_second_derivatives',val3,rho_CT_CT,rho_CT_CT_ca,'rho_CT_CT')
call check_accuracy('rho_second_derivatives',val4,rho_SA_P,rho_SA_P_ca,'rho_SA_P')
call check_accuracy('rho_second_derivatives',val5,rho_CT_P,rho_CT_P_ca,'rho_CT_P')

call gsw_rho_first_derivatives_wrt_enthalpy(sa,ct,p,val1,val2)
call check_accuracy('rho_first_derivatives_wrt_enthalpy',val1,rho_SA_wrt_h,rho_SA_wrt_h_ca,'rho_SA_wrt_h')
call check_accuracy('rho_first_derivatives_wrt_enthalpy',val2,rho_h,rho_h_ca,'rho_h')

call gsw_rho_second_derivatives_wrt_enthalpy(sa,ct,p,val1,val2,val3)
call check_accuracy('rho_second_derivatives_wrt_enthalpy',val1,rho_SA_SA_wrt_h,rho_SA_SA_wrt_h_ca,'rho_SA_SA_wrt_h')
call check_accuracy('rho_second_derivatives_wrt_enthalpy',val2,rho_SA_h,rho_SA_h_ca,'rho_SA_h')
call check_accuracy('rho_second_derivatives_wrt_enthalpy',val3,rho_h_h,rho_h_h_ca,'rho_h_h')

value = gsw_sigma0(sa,ct)
call check_accuracy('sigma0',value,sigma0,sigma0_ca)

value = gsw_sigma1(sa,ct)
call check_accuracy('sigma1',value,sigma1,sigma1_ca)

value = gsw_sigma2(sa,ct)
call check_accuracy('sigma2',value,sigma2,sigma2_ca)

value = gsw_sigma3(sa,ct)
call check_accuracy('sigma3',value,sigma3,sigma3_ca)

value = gsw_sigma4(sa,ct)
call check_accuracy('sigma4',value,sigma4,sigma4_ca)

value = gsw_sound_speed(sa,ct,p)
call check_accuracy('sound_speed',value,sound_speed,sound_speed_ca)

value = gsw_kappa(sa,ct,p)
call check_accuracy('kappa',value,kappa,kappa_ca)

value = gsw_cabbeling(sa,ct,p)
call check_accuracy('cabbeling',value,cabbeling,cabbeling_ca)

value = gsw_thermobaric(sa,ct,p)
call check_accuracy('thermobaric',value,thermobaric,thermobaric_ca)

value = gsw_sa_from_rho(rho,ct,p)
call check_accuracy('SA_from_rho',value,SA_from_rho,SA_from_rho_ca)

call gsw_ct_from_rho(rho,sa,p,value)
call check_accuracy('CT_from_rho',value,CT_from_rho,CT_from_rho_ca)

value = gsw_ct_maxdensity(sa,p)
call check_accuracy('CT_maxdensity',value,CT_maxdensity,CT_maxdensity_ca)

value = gsw_internal_energy(sa,ct,p)
call check_accuracy('internal_energy',value,internal_energy,internal_energy_ca)

h = gsw_enthalpy(sa,ct,p)
call check_accuracy('enthalpy',h,enthalpy,enthalpy_ca)

value = gsw_enthalpy_diff(sa,ct,p_shallow,p_deep)
call check_accuracy('enthalpy_diff',value,enthalpy_diff,enthalpy_diff_ca)

value = gsw_ct_from_enthalpy(sa,h,p)
call check_accuracy('CT_from_enthalpy',value,CT_from_enthalpy,CT_from_enthalpy_ca)

value = gsw_dynamic_enthalpy(sa,ct,p)
call check_accuracy('dynamic_enthalpy',value,dynamic_enthalpy,dynamic_enthalpy_ca)

call gsw_enthalpy_first_derivatives(sa,ct,p,val1,val2)
call check_accuracy('enthalpy_first_derivatives',val1,h_SA,h_SA_ca,'h_SA')
call check_accuracy('enthalpy_first_derivatives',val2,h_CT,h_CT_ca,'h_CT')

call gsw_enthalpy_second_derivatives(sa,ct,p,val1,val2,val3)
call check_accuracy('enthalpy_second_derivatives',val1,h_SA_SA,h_SA_SA_ca,'h_SA_SA')
call check_accuracy('enthalpy_second_derivatives',val2,h_SA_CT,h_SA_CT_ca,'h_SA_CT')
call check_accuracy('enthalpy_second_derivatives',val3,h_CT_CT,h_CT_CT_ca,'h_CT_CT')

!------------------------------------------------------------------------------
call section_title('Derivatives of entropy, CT and pt')

call gsw_ct_first_derivatives(sa,pt,val1,val2)
call check_accuracy('CT_first_derivatives',val1,CT_SA,CT_SA_ca,'CT_SA')
call check_accuracy('CT_first_derivatives',val2,CT_pt,CT_pt_ca,'CT_pt')

call gsw_ct_second_derivatives(sa,pt,val1,val2,val3)
call check_accuracy('CT_second_derivatives',val1,CT_SA_SA,CT_SA_SA_ca,'CT_SA_SA')
call check_accuracy('CT_second_derivatives',val2,CT_SA_pt,CT_SA_pt_ca,'CT_SA_pt')
call check_accuracy('CT_second_derivatives',val3,CT_pt_pt,CT_pt_pt_ca,'CT_pt_pt')

call gsw_entropy_first_derivatives(sa,ct,val1,val2)
call check_accuracy('entropy_first_derivatives',val1,eta_SA,eta_SA_ca,'eta_SA')
call check_accuracy('entropy_first_derivatives',val2,eta_CT,eta_CT_ca,'eta_CT')

call gsw_entropy_second_derivatives(sa,ct,val1,val2,val3)
call check_accuracy('entropy_second_derivatives',val1,eta_SA_SA,eta_SA_SA_ca,'eta_SA_SA')
call check_accuracy('entropy_second_derivatives',val2,eta_SA_CT,eta_SA_CT_ca,'eta_SA_CT')
call check_accuracy('entropy_second_derivatives',val3,eta_CT_CT,eta_CT_CT_ca,'eta_CT_CT')

call gsw_pt_first_derivatives(sa,ct,val1,val2)
call check_accuracy('pt_first_derivatives',val1,pt_SA,pt_SA_ca,'pt_SA')
call check_accuracy('pt_first_derivatives',val2,pt_CT,pt_CT_ca,'pt_CT')

call gsw_pt_second_derivatives(sa,ct,val1,val2,val3)
call check_accuracy('pt_second_derivatives',val1,pt_SA_SA,pt_SA_SA_ca,'pt_SA_SA')
call check_accuracy('pt_second_derivatives',val2,pt_SA_CT,pt_SA_CT_ca,'pt_SA_CT')
call check_accuracy('pt_second_derivatives',val3,pt_CT_CT,pt_CT_CT_ca,'pt_CT_CT')

!------------------------------------------------------------------------------
call section_title('Freezing temperatures')

saturation_fraction = 0.5_r8

ctf = gsw_ct_freezing(sa,p,saturation_fraction)
call check_accuracy('CT_freezing',ctf,CT_freezing,CT_freezing_ca)

ctf_poly = gsw_ct_freezing_poly(sa,p,saturation_fraction)
call check_accuracy('CT_freezing_poly',ctf_poly,CT_freezing_poly,CT_freezing_poly_ca)

tf = gsw_t_freezing(sa,p,saturation_fraction)
call check_accuracy('t_freezing',tf,t_freezing,t_freezing_ca)

tf_poly = gsw_t_freezing_poly(sa,p,saturation_fraction)
call check_accuracy('t_freezing_poly',tf_poly,t_freezing_poly,t_freezing_poly_ca)

value = gsw_pot_enthalpy_ice_freezing(sa,p)
call check_accuracy('pot_enthalpy_ice_freezing',value,pot_enthalpy_ice_freezing,pot_enthalpy_ice_freezing_ca)

value = gsw_pot_enthalpy_ice_freezing_poly(sa,p)
call check_accuracy('pot_enthalpy_ice_freezing_poly',value,pot_enthalpy_ice_freezing_poly,pot_enthalpy_ice_freezing_poly_ca)

value = gsw_sa_freezing_from_ct(ctf,p,saturation_fraction)
call check_accuracy('SA_freezing_from_CT',value,SA_freezing_from_CT,SA_freezing_from_CT_ca)

value = gsw_sa_freezing_from_ct_poly(ctf_poly,p,saturation_fraction)
call check_accuracy('SA_freezing_from_CT_poly',value,SA_freezing_from_CT_poly,SA_freezing_from_CT_poly_ca)

value = gsw_sa_freezing_from_t(tf,p,saturation_fraction)
call check_accuracy('SA_freezing_from_t',value,SA_freezing_from_t,SA_freezing_from_t_ca)

value = gsw_sa_freezing_from_t_poly(tf_poly,p,saturation_fraction)
call check_accuracy('SA_freezing_from_t_poly',value,SA_freezing_from_t_poly,SA_freezing_from_t_poly_ca)

call gsw_ct_freezing_first_derivatives(sa,p,saturation_fraction,val1,val2)
call check_accuracy('CT_freezing_first_derivatives',val1,CTfreezing_SA,CTfreezing_SA_ca,'CTfreezing_SA')
call check_accuracy('CT_freezing_first_derivatives',val2,CTfreezing_P,CTfreezing_P_ca,'CTfreezing_P')

call gsw_ct_freezing_first_derivatives_poly(sa,p,saturation_fraction,val4,val5)
call check_accuracy('CT_freezing_first_derivatives_poly',val4,CTfreezing_SA_poly,CTfreezing_SA_poly_ca,'CTfreezing_SA_poly')
call check_accuracy('CT_freezing_first_derivatives_poly',val5,CTfreezing_P_poly,CTfreezing_P_poly_ca,'CTfreezing_P_poly')

call gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,val1,val2)
call check_accuracy('t_freezing_first_derivatives',val1,tfreezing_SA,tfreezing_SA_ca,'tfreezing_SA')
call check_accuracy('t_freezing_first_derivatives',val2,tfreezing_P,tfreezing_P_ca,'tfreezing_P')

call gsw_t_freezing_first_derivatives_poly(sa,p,saturation_fraction,val4,val5)
call check_accuracy('t_freezing_first_derivatives_poly',val4,tfreezing_SA_poly,tfreezing_SA_poly_ca,'tfreezing_SA_poly')
call check_accuracy('t_freezing_first_derivatives_poly',val5,tfreezing_P_poly,tfreezing_P_poly_ca,'tfreezing_P_poly')

call gsw_pot_enthalpy_ice_freezing_first_derivatives(sa,p,val1,val2)
call check_accuracy('pot_enthalpy_ice_freezing_first_derivatives',val1, &
                    pot_enthalpy_ice_freezing_SA,pot_enthalpy_ice_freezing_SA_ca,'pot_enthalpy_ice_freezing_SA')
call check_accuracy('pot_enthalpy_ice_freezing_first_derivatives',val2, &
                    pot_enthalpy_ice_freezing_P,pot_enthalpy_ice_freezing_P_ca,'pot_enthalpy_ice_freezing_P')

call gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(sa,p,val1,val2)
call check_accuracy('pot_enthalpy_ice_freezing_first_derivatives_poly',val1, &
                    pot_enthalpy_ice_freezing_SA_poly,pot_enthalpy_ice_freezing_SA_poly_ca,'pot_enthalpy_ice_freezing_SA_poly')
call check_accuracy('pot_enthalpy_ice_freezing_first_derivatives_poly',val2, &
                    pot_enthalpy_ice_freezing_P_poly,pot_enthalpy_ice_freezing_P_poly_ca,'pot_enthalpy_ice_freezing_P_poly')

!------------------------------------------------------------------------------
call section_title('Isobaric Melting Enthalpy and Isobaric Evaporation Enthalpy')

value = gsw_latentheat_melting(sa,p)
call check_accuracy('latentheat_melting',value,latentheat_melting,latentheat_melting_ca)

value = gsw_latentheat_evap_ct(sa,ct)
call check_accuracy('latentheat_evap_CT',value,latentheat_evap_CT,latentheat_evap_CT_ca)

value = gsw_latentheat_evap_t(sa,t)
call check_accuracy('latentheat_evap_t',value,latentheat_evap_t,latentheat_evap_t_ca)

!------------------------------------------------------------------------------
call section_title('Planet Earth properties')

value = gsw_grav(lat,p)
call check_accuracy('grav',value,grav,grav_ca)

!------------------------------------------------------------------------------
call section_title('Density and enthalpy in terms of CT, derived from the exact Gibbs function')

value = gsw_enthalpy_ct_exact(sa,ct,p)
call check_accuracy('enthalpy_CT_exact',value,enthalpy_CT_exact,enthalpy_CT_exact_ca)

call gsw_enthalpy_first_derivatives_ct_exact(sa,ct,p,val1,val2)
call check_accuracy('enthalpy_first_derivatives_CT_exact',val1,h_SA_CT_exact,h_SA_CT_exact_ca,'h_SA_CT_exact')
call check_accuracy('enthalpy_first_derivatives_CT_exact',val2,h_CT_CT_exact,h_CT_CT_exact_ca,'h_CT_CT_exact')

call gsw_enthalpy_second_derivatives_ct_exact(sa,ct,p,val1,val2,val3)
call check_accuracy('enthalpy_second_derivatives_CT_exact',val1, &
                    h_SA_SA_CT_exact,h_SA_SA_CT_exact_ca,'h_SA_SA_CT_exact')
call check_accuracy('enthalpy_second_derivatives_CT_exact',val2, &
                    h_SA_CT_CT_exact,h_SA_CT_CT_exact_ca,'h_SA_CT_CT_exact')
call check_accuracy('enthalpy_second_derivatives_CT_exact',val3, &
                    h_CT_CT_CT_exact,h_CT_CT_CT_exact_ca,'h_CT_CT_CT_exact')

!------------------------------------------------------------------------------
call section_title('Basic thermodynamic properties in terms of in-situ t, based on the exact Gibbs function')

value = gsw_rho_t_exact(sa,t,p)
call check_accuracy('rho_t_exact',value,rho_t_exact,rho_t_exact_ca)

value = gsw_pot_rho_t_exact(sa,t,p,pref(1))
call check_accuracy('pot_rho_t_exact',value,pot_rho_t_exact,pot_rho_t_exact_ca)

value = gsw_alpha_wrt_t_exact(sa,t,p)
call check_accuracy('alpha_wrt_t_exact',value,alpha_wrt_t_exact,alpha_wrt_t_exact_ca)

value = gsw_beta_const_t_exact(sa,t,p)
call check_accuracy('beta_const_t_exact',value,beta_const_t_exact,beta_const_t_exact_ca)

value = gsw_specvol_t_exact(sa,t,p)
call check_accuracy('specvol_t_exact',value,specvol_t_exact,specvol_t_exact_ca)

value = gsw_sound_speed_t_exact(sa,t,p)
call check_accuracy('sound_speed_t_exact',value,sound_speed_t_exact,sound_speed_t_exact_ca)

value = gsw_kappa_t_exact(sa,t,p)
call check_accuracy('kappa_t_exact',value,kappa_t_exact,kappa_t_exact_ca)

value = gsw_enthalpy_t_exact(sa,t,p)
call check_accuracy('enthalpy_t_exact',value,enthalpy_t_exact,enthalpy_t_exact_ca)

call gsw_ct_first_derivatives_wrt_t_exact(sa,t,p,val1,val2,val3)
call check_accuracy('CT_first_derivatives_wrt_t_exact',val1,CT_SA_wrt_t,CT_SA_wrt_t_ca,'CT_SA_wrt_t')
call check_accuracy('CT_first_derivatives_wrt_t_exact',val2,CT_T_wrt_t,CT_T_wrt_t_ca,'CT_T_wrt_t')
call check_accuracy('CT_first_derivatives_wrt_t_exact',val3,CT_P_wrt_t,CT_P_wrt_t_ca,'CT_P_wrt_t')

value = gsw_chem_potential_water_t_exact(sa,t,p)
call check_accuracy('chem_potential_water_t_exact',value,chem_potential_water_t_exact,chem_potential_water_t_exact_ca)

value = gsw_t_deriv_chem_potential_water_t_exact(sa,t,p)
call check_accuracy('t_deriv_chem_potential_water_t_exact',value,t_deriv_chem_potential_water_t_exact, &
    t_deriv_chem_potential_water_t_exact_ca)

value = gsw_dilution_coefficient_t_exact(sa,t,p)
call check_accuracy('dilution_coefficient_t_exact',value,dilution_coefficient_t_exact,dilution_coefficient_t_exact_ca)

!------------------------------------------------------------------------------
call section_title('Library functions of the GSW Toolbox')

value = gsw_deltasa_atlas(p,long,lat)
call check_accuracy('deltaSA_atlas',value,deltaSA_atlas,deltaSA_atlas_ca)

value = gsw_fdelta(p,long,lat)
call check_accuracy('Fdelta',value,Fdelta,Fdelta_ca)

!------------------------------------------------------------------------------
call section_title('Vertical stability')

deallocate(val1,val2,val3,val4,val5)
allocate(val1(cast_mpres_m,cast_mpres_n))
allocate(val2(cast_mpres_m,cast_mpres_n))
allocate(val3(cast_mpres_m,cast_mpres_n))
allocate(val4(cast_mpres_m,cast_mpres_n))
allocate(val5(cast_mpres_m,cast_mpres_n))
allocate(val6(cast_mpres_m,cast_mpres_n))
allocate(val7(cast_mpres_m,cast_mpres_n))
allocate(val8(cast_mpres_m,cast_mpres_n))

do i = 1, cast_mpres_n
    call gsw_turner_rsubrho(sa(:,i),ct(:,i),p(:,i),val1(:,i),val2(:,i), &
                            val3(:,i))
end do
call check_accuracy('Turner_Rsubrho',val1,Tu,Tu_ca,'Tu')
call check_accuracy('Turner_Rsubrho',val2,Rsubrho,Rsubrho_ca,'Rsubrho')
call check_accuracy('Turner_Rsubrho',val3,p_mid_TuRsr,p_mid_TuRsr_ca,'p_mid_TuRsr')

do i = 1, cast_mpres_n
    call gsw_nsquared(sa(:,i),ct(:,i),p(:,i),lat(:,i),val1(:,i),val2(:,i))
end do
call check_accuracy('Nsquared',val1,n2,n2_ca,'n2')
call check_accuracy('Nsquared',val2,p_mid_n2,p_mid_n2_ca,'p_mid_n2')

do i = 1, cast_mpres_n
    call gsw_nsquared_min(sa(:,i),ct(:,i),p(:,i),lat(:,i),val1(:,i),val2(:,i), &
                val3(:,i),val4(:,i),val5(:,i),val6(:,i),val7(:,i),val8(:,i))
end do
call check_accuracy('Nsquared_min',val1,n2min,n2min_ca,'n2min')
call check_accuracy('Nsquared_min',val2,n2min_pmid,n2min_pmid_ca,'n2min_pmid')
call check_accuracy('Nsquared_min',val3,n2min_specvol,n2min_specvol_ca,'n2min_specvol')
call check_accuracy('Nsquared_min',val4,n2min_alpha,n2min_alpha_ca,'n2min_alpha')
call check_accuracy('Nsquared_min',val5,n2min_beta,n2min_beta_ca,'n2min_beta')
call check_accuracy('Nsquared_min',val6,n2min_dsa,n2min_dsa_ca,'n2min_dsa')
call check_accuracy('Nsquared_min',val7,n2min_dct,n2min_dct_ca,'n2min_dct')
call check_accuracy('Nsquared_min',val8,n2min_dp,n2min_dp_ca,'n2min_dp')

do i = 1, cast_mpres_n
    call gsw_ipv_vs_fnsquared_ratio(sa(:,i),ct(:,i),p(:,i),pref(1),val1(:,i), &
                                    val2(:,i))
end do
call check_accuracy('IPV_vs_fNsquared_ratio',val1,IPVfN2,IPVfN2_ca,'IPVfN2')
call check_accuracy('IPV_vs_fNsquared_ratio',val2,p_mid_IPVfN2,p_mid_IPVfN2_ca,'p_mid_IPVfN2')
  
value = gsw_nsquared_lowerlimit(p,long,lat)
call check_accuracy('n2_lowerlimit',value,n2_lowerlimit,n2_lowerlimit_ca)

deallocate(value)
allocate(value(cast_mpres_n,1))
do i = 1, cast_mpres_n
    value(i,1) = gsw_mlp(sa(:,i),ct(:,i),p(:,i))
end do
allocate(check_value(cast_mpres_n,1))
check_value(:,1) = mlp
call check_accuracy('mlp',value,check_value,mlp_ca)

!------------------------------------------------------------------------------
call section_title('Geostrophic streamfunctions and acoustic travel time')

do i = 1, cast_mpres_n
    n = count(sa(:,i) .eq. sa(:,i))  ! check for NaN's
    val1(:,i) = gsw_geo_strf_dyn_height(sa(:n,i),ct(:n,i),p(:n,i),pref(1))
    if (n .lt. cast_mpres_m) val1(n+1:,i) = sa(n+1:cast_mpres_m,i)
end do
call check_accuracy('geo_strf_dyn_height',val1,geo_strf_dyn_height,geo_strf_dyn_height_ca)

do i = 1, cast_mpres_n
    call gsw_geo_strf_dyn_height_pc(sa(:,i),ct(:,i),delta_p(:,i), &
                                    val1(:,i),val2(:,i))
end do
call check_accuracy('geo_strf_dyn_height_pc',val1,geo_strf_dyn_height_pc,geo_strf_dyn_height_pc_ca, &
    'geo_strf_dyn_height_pc')
call check_accuracy('geo_strf_dyn_height_pc',val2,geo_strf_dyn_height_pc_p_mid,geo_strf_dyn_height_pc_p_mid_ca, &
    'geo_strf_dyn_height_pc_p_mid')

!------------------------------------------------------------------------------
call section_title('Thermodynamic properties of ice Ih')

deallocate(value)
allocate(value(cast_ice_m,cast_ice_n))

value = gsw_rho_ice(t_seaice,p_arctic)
call check_accuracy('rho_ice',value,rho_ice,rho_ice_ca)

value = gsw_alpha_wrt_t_ice(t_seaice,p_arctic)
call check_accuracy('alpha_wrt_t_ice',value,alpha_wrt_t_ice,alpha_wrt_t_ice_ca)

value = gsw_specvol_ice(t_seaice,p_arctic)
call check_accuracy('specvol_ice',value,specvol_ice,specvol_ice_ca)

value = gsw_pressure_coefficient_ice(t_seaice,p_arctic)
call check_accuracy('pressure_coefficient_ice',value,pressure_coefficient_ice,pressure_coefficient_ice_ca)

value = gsw_sound_speed_ice(t_seaice,p_arctic)
call check_accuracy('sound_speed_ice',value,sound_speed_ice,sound_speed_ice_ca)

value = gsw_kappa_ice(t_seaice,p_arctic)
call check_accuracy('kappa_ice',value,kappa_ice,kappa_ice_ca)

value = gsw_kappa_const_t_ice(t_seaice,p_arctic)
call check_accuracy('kappa_const_t_ice',value,kappa_const_t_ice,kappa_const_t_ice_ca)

value = gsw_internal_energy_ice(t_seaice,p_arctic)
call check_accuracy('internal_energy_ice',value,internal_energy_ice,internal_energy_ice_ca)

value = gsw_enthalpy_ice(t_seaice,p_arctic)
call check_accuracy('enthalpy_ice',value,enthalpy_ice,enthalpy_ice_ca)

value = gsw_entropy_ice(t_seaice,p_arctic)
call check_accuracy('entropy_ice',value,entropy_ice,entropy_ice_ca)

value = gsw_cp_ice(t_seaice,p_arctic)
call check_accuracy('cp_ice',value,cp_ice,cp_ice_ca)

value = gsw_chem_potential_water_ice(t_seaice,p_arctic)
call check_accuracy('chem_potential_water_ice',value,chem_potential_water_ice,chem_potential_water_ice_ca)

value = gsw_helmholtz_energy_ice(t_seaice,p_arctic)
call check_accuracy('Helmholtz_energy_ice',value,Helmholtz_energy_ice,Helmholtz_energy_ice_ca)

value = gsw_adiabatic_lapse_rate_ice(t_seaice,p_arctic)
call check_accuracy('adiabatic_lapse_rate_ice',value,adiabatic_lapse_rate_ice,adiabatic_lapse_rate_ice_ca)

pt0 = gsw_pt0_from_t_ice(t_seaice,p_arctic)
call check_accuracy('pt0_from_t_ice',pt0,pt0_from_t_ice,pt0_from_t_ice_ca)

value = gsw_pt_from_t_ice(t_seaice,p_arctic,pref(1))
call check_accuracy('pt_from_t_ice',value,pt_from_t_ice,pt_from_t_ice_ca)

value = gsw_t_from_pt0_ice(pt0,p_arctic)
call check_accuracy('t_from_pt0_ice',value,t_from_pt0_ice,t_from_pt0_ice_ca)

h = gsw_pot_enthalpy_from_pt_ice(pt0)
call check_accuracy('pot_enthalpy_from_pt_ice',h,pot_enthalpy_from_pt_ice,pot_enthalpy_from_pt_ice_ca)

value = gsw_pt_from_pot_enthalpy_ice(h)
call check_accuracy('pt_from_pot_enthalpy_ice',value,pt_from_pot_enthalpy_ice,pt_from_pot_enthalpy_ice_ca)

h = gsw_pot_enthalpy_from_pt_ice_poly(pt0)
call check_accuracy('pot_enthalpy_from_pt_ice_poly',h,pot_enthalpy_from_pt_ice_poly,pot_enthalpy_from_pt_ice_poly_ca)

value = gsw_pt_from_pot_enthalpy_ice_poly(h)
call check_accuracy('pt_from_pot_enthalpy_ice_poly',value,pt_from_pot_enthalpy_ice_poly,pt_from_pot_enthalpy_ice_poly_ca)

saturation_fraction = 0.5_r8

value = gsw_pressure_freezing_ct(sa_arctic,ct_arctic-1.0_r8,saturation_fraction)
call check_accuracy('pressure_freezing_CT',value,pressure_freezing_CT,pressure_freezing_CT_ca)

!------------------------------------------------------------------------------
call section_title('Thermodynamic interaction between ice and seawater')

deallocate(val1,val2,val3)
allocate(val1(cast_ice_m,cast_ice_n))
allocate(val2(cast_ice_m,cast_ice_n))
allocate(val3(cast_ice_m,cast_ice_n))

value = gsw_melting_ice_sa_ct_ratio(sa_arctic,ct_arctic,p_arctic,t_ice)
call check_accuracy('melting_ice_SA_CT_ratio',value,melting_ice_SA_CT_ratio,melting_ice_SA_CT_ratio_ca)

value = gsw_melting_ice_sa_ct_ratio_poly(sa_arctic,ct_arctic,p_arctic,t_ice)
call check_accuracy('melting_ice_SA_CT_ratio_poly',value,melting_ice_SA_CT_ratio_poly,melting_ice_SA_CT_ratio_poly_ca)

value = gsw_melting_ice_equilibrium_sa_ct_ratio(sa_arctic,p_arctic)
call check_accuracy('melting_ice_equilibrium_SA_CT_ratio',value,melting_ice_equilibrium_SA_CT_ratio, &
    melting_ice_equilibrium_SA_CT_ratio_ca)

value = gsw_melting_ice_equilibrium_sa_ct_ratio_poly(sa_arctic,p_arctic)
call check_accuracy('melting_ice_equilibrium_SA_CT_ratio_poly',value,melting_ice_equilibrium_SA_CT_ratio_poly, &
    melting_ice_equilibrium_SA_CT_ratio_poly_ca)

call gsw_melting_ice_into_seawater(sa_arctic,ct_arctic+0.1_r8,p_arctic,w_ice,t_ice, &
                                   val1,val2,val3)
call check_accuracy('melting_ice_into_seawater',val1, &
                    melting_ice_into_seawater_SA_final,melting_ice_into_seawater_SA_final_ca,'melting_ice_into_seawater_SA_final')
call check_accuracy('melting_ice_into_seawater',val2, &
                    melting_ice_into_seawater_CT_final,melting_ice_into_seawater_CT_final_ca,'melting_ice_into_seawater_CT_final')
!call check_accuracy('melting_ice_into_seawater',val3, &
!                    'melting_ice_into_seawater_w_Ih')

call gsw_ice_fraction_to_freeze_seawater(sa_arctic,ct_arctic,p_arctic,t_ice, &
                                         val1,val2,val3)
call check_accuracy('ice_fraction_to_freeze_seawater',val1, &
                    ice_fraction_to_freeze_seawater_SA_freeze,ice_fraction_to_freeze_seawater_SA_freeze_ca, &
                    'ice_fraction_to_freeze_seawater_SA_freeze')
call check_accuracy('ice_fraction_to_freeze_seawater',val2, &
                    ice_fraction_to_freeze_seawater_CT_freeze,ice_fraction_to_freeze_seawater_CT_freeze_ca, &
                    'ice_fraction_to_freeze_seawater_CT_freeze')
call check_accuracy('ice_fraction_to_freeze_seawater',val3, &
                    ice_fraction_to_freeze_seawater_w_Ih,ice_fraction_to_freeze_seawater_w_Ih_ca, &
                    'ice_fraction_to_freeze_seawater_w_Ih')

call gsw_frazil_ratios_adiabatic(sa_arctic,p_arctic,w_ice,val1,val2,val3)
call check_accuracy('frazil_ratios_adiabatic',val1,dSA_dCT_frazil,dSA_dCT_frazil_ca,'dSA_dCT_frazil')
call check_accuracy('frazil_ratios_adiabatic',val2,dSA_dP_frazil,dSA_dP_frazil_ca,'dSA_dP_frazil')
call check_accuracy('frazil_ratios_adiabatic',val3,dCT_dP_frazil,dCT_dP_frazil_ca,'dCT_dP_frazil')

call gsw_frazil_ratios_adiabatic_poly(sa_arctic,p_arctic,w_ice,val1,val2,val3)
call check_accuracy('frazil_ratios_adiabatic_poly',val1,dSA_dCT_frazil_poly,dSA_dCT_frazil_poly_ca,'dSA_dCT_frazil_poly')
call check_accuracy('frazil_ratios_adiabatic_poly',val2,dSA_dP_frazil_poly,dSA_dP_frazil_poly_ca,'dSA_dP_frazil_poly')
call check_accuracy('frazil_ratios_adiabatic_poly',val3,dCT_dP_frazil_poly,dCT_dP_frazil_poly_ca,'dCT_dP_frazil_poly')

call gsw_frazil_properties_potential(sa_bulk,h_pot_bulk,p_arctic,val1,val2,val3)
call check_accuracy('frazil_properties_potential',val1, &
                    frazil_properties_potential_SA_final,frazil_properties_potential_SA_final_ca, &
                    'frazil_properties_potential_SA_final')
call check_accuracy('frazil_properties_potential',val2, &
                    frazil_properties_potential_CT_final,frazil_properties_potential_CT_final_ca, &
                    'frazil_properties_potential_CT_final')
call check_accuracy('frazil_properties_potential',val3, &
                    frazil_properties_potential_w_Ih_final,frazil_properties_potential_w_Ih_final_ca, &
                    'frazil_properties_potential_w_Ih_final')

call gsw_frazil_properties_potential_poly(sa_bulk,h_pot_bulk,p_arctic,val1, &
                                          val2,val3)
call check_accuracy('frazil_properties_potential_poly',val1, &
                    frazil_properties_potential_poly_SA_final,frazil_properties_potential_poly_SA_final_ca, &
                    'frazil_properties_potential_poly_SA_final')
call check_accuracy('frazil_properties_potential_poly',val2, &
                    frazil_properties_potential_poly_CT_final,frazil_properties_potential_poly_CT_final_ca, &
                    'frazil_properties_potential_poly_CT_final')
call check_accuracy('frazil_properties_potential_poly',val3, &
                    frazil_properties_potential_poly_w_Ih_final,frazil_properties_potential_poly_w_Ih_final_ca, &
                    'frazil_properties_potential_poly_w_Ih_final')

call gsw_frazil_properties(sa_bulk,h_bulk,p_arctic,val1,val2,val3)
call check_accuracy('frazil_properties',val1,frazil_properties_SA_final,frazil_properties_SA_final_ca, &
    'frazil_properties_SA_final')
call check_accuracy('frazil_properties',val2,frazil_properties_CT_final,frazil_properties_CT_final_ca, &
    'frazil_properties_CT_final')
call check_accuracy('frazil_properties',val3,frazil_properties_w_Ih_final,frazil_properties_w_Ih_final_ca, &
    'frazil_properties_w_Ih_final')

!------------------------------------------------------------------------------
call section_title('Thermodynamic interaction between seaice and seawater')

value = gsw_melting_seaice_sa_ct_ratio(sa_arctic,ct_arctic,p_arctic, &
                                       sa_seaice,t_seaice)
call check_accuracy('melting_seaice_SA_CT_ratio',value,melting_seaice_SA_CT_ratio,melting_seaice_SA_CT_ratio_ca)

value = gsw_melting_seaice_sa_ct_ratio_poly(sa_arctic,ct_arctic,p_arctic, &
                                            sa_seaice,t_seaice)
call check_accuracy('melting_seaice_SA_CT_ratio_poly',value,melting_seaice_SA_CT_ratio_poly, &
    melting_seaice_SA_CT_ratio_poly_ca)

value = gsw_melting_seaice_equilibrium_sa_ct_ratio(sa_arctic,p_arctic)
call check_accuracy('melting_seaice_equilibrium_SA_CT_ratio',value,melting_seaice_equilibrium_SA_CT_ratio, &
    melting_seaice_equilibrium_SA_CT_ratio_ca)

value = gsw_melting_seaice_equilibrium_sa_ct_ratio_poly(sa_arctic,p_arctic)
call check_accuracy('melting_seaice_equilibrium_SA_CT_ratio_poly',value,melting_seaice_equilibrium_SA_CT_ratio_poly, &
    melting_seaice_equilibrium_SA_CT_ratio_poly_ca)

call gsw_melting_seaice_into_seawater(sa_arctic,ct_arctic,p_arctic, &
                                      w_seaice,sa_seaice,t_seaice,val1,val2)
call check_accuracy('melting_seaice_into_seawater',val1, &
                    melting_seaice_into_seawater_SA_final,melting_seaice_into_seawater_SA_final_ca, &
                    'melting_seaice_into_seawater_SA_final')
call check_accuracy('melting_seaice_into_seawater',val2, &
                    melting_seaice_into_seawater_CT_final,melting_seaice_into_seawater_CT_final_ca, &
                    'melting_seaice_into_seawater_CT_final')

call gsw_seaice_fraction_to_freeze_seawater(sa_arctic,ct_arctic,p_arctic, &
                                            sa_seaice,t_seaice,val1,val2,val3)
call check_accuracy('seaice_fraction_to_freeze_seawater',val1, &
                    seaice_fraction_to_freeze_seawater_SA_freeze,seaice_fraction_to_freeze_seawater_SA_freeze_ca, &
                    'seaice_fraction_to_freeze_seawater_SA_freeze')
call check_accuracy('seaice_fraction_to_freeze_seawater',val2, &
                    seaice_fraction_to_freeze_seawater_CT_freeze,seaice_fraction_to_freeze_seawater_CT_freeze_ca, &
                    'seaice_fraction_to_freeze_seawater_CT_freeze')
call check_accuracy('seaice_fraction_to_freeze_seawater',val3, &
                    seaice_fraction_to_freeze_seawater_w_Ih,seaice_fraction_to_freeze_seawater_w_Ih_ca, &
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

    subroutine check_accuracy (func_name, fvalue, check_value, check_limit, var_name, vprint)

    use gsw_mod_error_functions, only : gsw_error_limit

    implicit none

    character (*), intent(in) :: func_name
    real (r8), intent(in) :: fvalue(:,:)
    real (r8), intent(in) :: check_value(:,:)
    real (r8), intent(in) :: check_limit
    character (*), intent(in), optional :: var_name
    logical, intent(in), optional :: vprint

    integer :: ndots, i, j, k, ik=1, jk=1
    real (r8) :: dmax, drel
    real (r8) :: diff(size(fvalue,1),size(fvalue,2))
    character (80) :: message
    character (4) :: errflg

    character (*), parameter :: att_name = 'computation_accuracy'
    character (*), parameter :: &
        dots = ' .............................................................'

    if (present(var_name)) then

        if (len(func_name)+len(var_name).gt.55) then
            k = len(func_name) + len(var_name) - 55
            message = func_name // ' (..' // var_name(k:) // ')'
        else
            message = func_name // ' (' // var_name // ')'
        end if

    else

        message = func_name

    end if

    diff = abs(fvalue - check_value)
    where (check_value .eq. 9e90_r8) diff = 0.0_r8

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
