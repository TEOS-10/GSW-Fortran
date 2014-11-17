program gsw_check_functions

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)
integer  :: gsw_error_flag, nz = 3

real (r14), parameter :: sp = 35.5d0, sa = 35.7d0, sstar = 35.5d0
real (r14), parameter :: sr = 35.5d0, t = 15d0, ct = 20d0, pt = 15d0
real (r14), parameter :: p = 300d0, p_bs = 50d0, p_ref = 100d0, long = 260d0
real (r14), parameter :: long_bs = 20d0, lat = 16d0, lat_bs = 60d0, c = 43.6d0
real (r14), parameter :: saturation_fraction = 0.5d0, c_from_sp_ca = 6.163816124171717d-010
real (r14), parameter :: sp_from_c_ca = 1.297193463756230d-010, r_from_sp_ca = 1.436317731418058d-011
real (r14), parameter :: sp_from_r_ca = 2.681299626772216d-012, sp_salinometer_ca = 1.297131291266851d-010
real (r14), parameter :: sa_from_sp_ca = 1.300080043620255d-010, sstar_from_sp_ca = 1.300008989346679d-010
real (r14), parameter :: ct_from_t_ca = 6.261107188265669d-010, deltasa_from_sp_ca = 6.963318810448982d-013
real (r14), parameter :: sa_sa_sstar_from_sp_ca = 1.300008989346679d-010, sstar_sa_sstar_from_sp_ca = 1.300008989346679d-010
real (r14), parameter :: sr_from_sp_ca = 1.303233077010191d-010, sp_from_sr_ca = 1.297122409482654d-010
real (r14), parameter :: sp_from_sa_ca = 1.297113527698457d-010, sstar_from_sa_ca = 1.300008989346679d-010
real (r14), parameter :: sa_from_sstar_ca = 1.300222152167407d-010, sp_from_sstar_ca = 1.297122409482654d-010
real (r14), parameter :: t_from_CT_ca = 6.000142604989378d-010, pt_from_ct_ca = 6.054037271496782d-010
real (r14), parameter :: ct_from_pt_ca = 6.261107188265669d-010, pot_enthalpy_from_pt_ca = 2.499356924090534d-006
real (r14), parameter :: pt0_from_t_ca = 6.054037271496782d-010, pt_from_t_ca = 6.054037271496782d-010
real (r14), parameter :: entropy_from_ct_ca = 9.028163105995191d-009, ct_from_entropy_ca = 6.261107188265669d-010
real (r14), parameter :: entropy_from_pt_ca = 9.028163105995191d-009, pt_from_entropy_ca = 6.054072798633570d-010
real (r14), parameter :: entropy_from_t_ca = 9.028163105995191d-009
real (r14), parameter :: adiabatic_lapse_rate_from_t_ca = 5.699738675609156d-019
real (r14), parameter :: adiabatic_lapse_rate_from_ct_ca = 5.699743845487985d-019
real (r14), parameter :: t90_from_t68_ca = 5.998579410970706d-010, t90_from_t48_ca = 5.997407015456702d-010
real (r14), parameter :: z_from_p_ca = 2.287223921371151d-008, p_from_z_ca = 2.301931090187281d-008
real (r14), parameter :: depth_from_z_ca = 2.287223921371151d-008, z_from_depth_ca = 2.287223921371151d-008
real (r14), parameter :: molality_from_sa_ca = 4.446665258228677d-012, ionic_strength_from_sa_ca = 2.768674178810215d-012
real (r14), parameter :: rho_ca = 2.945625965367071d-010, alpha_ca = 8.264713918662917d-015
real (r14), parameter :: beta_ca = 1.846179459308317d-015, alpha_on_beta_ca = 1.052907760978883d-11
real (r14), parameter :: rho_rab_ca = 2.944489096989855d-010, drho_dsa_ca = 1.943112337698949d-12
real (r14), parameter :: drho_dct_ca = 8.315403920988729d-12, drho_dp_ca = 1.782157321023048d-18
real (r14), parameter :: alpha_rab_ca = 8.264713918662917d-015, beta_rab_ca = 1.846179459308317d-015
real (r14), parameter :: specvol_ca = 2.821094052807283d-016, specvol_anom_ca = 2.810252031082428d-016
real (r14), parameter :: sigma0_ca = 2.933120413217694d-010, sigma1_ca = 2.999058779096231d-010
real (r14), parameter :: sigma2_ca = 3.060449671465904d-010, sigma3_ca = 3.119566827081144d-010
real (r14), parameter :: sigma4_ca = 3.180957719450817d-010, sound_speed_ca = 2.596152626210824d-009
real (r14), parameter :: kappa_ca = 1.717743939542931d-21, cabbeling_ca = 1.634722766223964d-16
real (r14), parameter :: thermobaric_ca = 4.890907320111513d-23
real (r14), parameter :: internal_energy_ca = 2.499342372175306d-006, enthalpy_ca = 2.499356924090534d-006
real (r14), parameter :: enthalpy_diff_ca = 1.154347728515859d-010, dynamic_enthalpy_ca = 2.288754734930485d-007
real (r14), parameter :: sa_from_rho_ca = 1.308677610722953d-010, ct_maxdensity_ca = 6.688338771709823d-011
real (r14), parameter :: ct_from_rho_ca = 6.298552790440226d-010, n2_ca = 1.578186366313350d-014
real (r14), parameter :: p_mid_n2_ca = 2.300021151313558d-008, tu_ca = 2.190718007000214d-008
real (r14), parameter :: rsubrho_ca = 1.709803143512545d-008, p_mid_tursr_ca = 2.300021151313558d-008
real (r14), parameter :: ipvfn2_ca = 3.474816878679121d-009, p_mid_ipvfn2_ca = 2.300021151313558d-008
real (r14), parameter :: geo_strf_dyn_height_ca = 9.969444363377988d-007, geo_strf_dyn_height_pc_ca = 4.210555459849275d-008
real (r14), parameter :: geo_strf_dyn_height_pc_p_mid_ca = 1.000000000000000d-015, geo_strf_isopycnal_ca = 4.297791695861974d-007
real (r14), parameter :: geo_strf_isopycnal_pc_ca = 8.540336438045415d-009, geo_strf_isopycnal_pc_p_mid_ca = 1.000000000000000d-015
real (r14), parameter :: geo_strf_montgomery_ca = 9.915690188933013d-007, geo_strf_cunningham_ca = 9.926097845891491d-007
real (r14), parameter :: geo_strf_velocity_ca = 8.024344404222727d-009, geo_strf_velocity_mid_lat_ca = 1.000000000000000d-015
real (r14), parameter :: geo_strf_velocity_mid_long_ca = 1.000000000000000d-015
real (r14), parameter :: isopycnal_slope_ratio_ca = 3.781384094736495d-010
real (r14), parameter :: g_ct_ca = 2.854610769986721d-009, p_mid_g_ct_ca = 2.300021151313558d-008
real (r14), parameter :: ntpptct_ca = 5.024869409453459d-013, ct_sa_ca = 1.006231122729906d-012
real (r14), parameter :: ct_pt_ca = 2.964295475749168d-013, ct_sa_sa_ca = 1.431146867680866d-014
real (r14), parameter :: ct_sa_pt_ca = 1.457167719820518d-014, ct_pt_pt_ca = 5.551115123125783d-014
real (r14), parameter :: h_sa_ca = 2.371365326325758d-010, h_ct_ca = 3.160494088660926d-010
real (r14), parameter :: h_p_ca = 2.818925648462312d-016, h_sa_sa_ca = 7.264189250122399d-013
real (r14), parameter :: h_sa_ct_ca = 2.188554892867956d-012, h_ct_ct_ca = 1.135647131889073d-011
real (r14), parameter :: eta_sa_ca = 4.653527563291959d-012, eta_ct_ca = 3.137579085432662d-011
real (r14), parameter :: eta_sa_sa_ca = 6.995931611797346d-013, eta_sa_ct_ca = 4.981922535098049d-014
real (r14), parameter :: eta_ct_ct_ca = 2.381358998881922d-013, pt_sa_ca = 9.670458878119348d-013
real (r14), parameter :: pt_ct_ca = 2.733369086627135d-013, pt_sa_sa_ca = 2.081668171172169d-014
real (r14), parameter :: pt_sa_ct_ca = 1.199127602768968d-014, pt_ct_ct_ca = 4.440892098500626d-014
real (r14), parameter :: ct_freezing_ca = 2.257127817983928d-011, t_freezing_ca = 2.157829470661454d-011
real (r14), parameter :: brinesa_ct_ca = 1.300080043620255d-010, brinesa_t_ca = 1.300293206440983d-010
real (r14), parameter :: latentheat_melting_ca = 6.286427378654480d-008, latentheat_evap_ct_ca = 1.455657184123993d-006
real (r14), parameter :: latentheat_evap_t_ca = 1.443084329366684d-006, f_ca = 1.000000000000000d-015
real (r14), parameter :: grav_ca = 5.329070518200751d-014, distance_ca = 4.470348358154297d-008
real (r14), parameter :: steric_height_ca = 1.017674460257467d-007, abs_pressure_from_p_ca = 2.300031483173370d-004
real (r14), parameter :: p_from_abs_pressure_ca = 2.300066626048647d-008, rho_ct_exact_ca = 2.944489096989855d-010
real (r14), parameter :: alpha_ct_exact_ca = 8.257102480598889d-015, beta_ct_exact_ca = 1.847372081698051d-015
real (r14), parameter :: rho_ct_exact_rab_ca = 2.944489096989855d-010, alpha_ct_exact_rab_ca = 8.257102480598889d-015
real (r14), parameter :: beta_ct_exact_rab_ca = 1.847372081698051d-015, specvol_ct_exact_ca = 2.818925648462312d-016
real (r14), parameter :: specvol_anom_ct_exact_ca = 2.805915222392486d-016, sigma0_ct_exact_ca = 2.929709808086045d-010
real (r14), parameter :: sigma1_ct_exact_ca = 2.994511305587366d-010, sigma2_ct_exact_ca = 3.055902197957039d-010
real (r14), parameter :: sigma3_ct_exact_ca = 3.117293090326712d-010, sigma4_ct_exact_ca = 3.176410245941952d-010
real (r14), parameter :: sound_speed_ct_exact_ca = 2.590240910649300d-009, internal_energy_ct_exact_ca = 2.499335096217692d-006
real (r14), parameter :: enthalpy_ct_exact_ca = 2.499349648132920d-006, enthalpy_diff_ct_exact_ca = 7.275957614183426d-011
real (r14), parameter :: dynamic_enthalpy_ct_exact_ca = 2.288797986693680d-007, sa_from_rho_ct_exact_ca = 1.306119656874216d-010
real (r14), parameter :: ct_maxdensity_exact_ca = 5.839062566792563d-011, ct_from_rho_exact_ca = 6.280096442878858d-010
real (r14), parameter :: rho_t_exact_ca = 2.944489096989855d-010, pot_rho_t_exact_ca = 2.929709808086045d-010
real (r14), parameter :: sigma0_pt0_exact_ca = 2.929709808086045d-010, specvol_t_exact_ca = 2.818925648462312d-016
real (r14), parameter :: specvol_anom_t_exact_ca = 2.805915222392486d-016, alpha_wrt_ct_t_exact_ca = 8.257094010269417d-015
real (r14), parameter :: alpha_wrt_pt_t_exact_ca = 8.599068316130290d-015, alpha_wrt_t_exact_ca = 8.594856868316542d-015
real (r14), parameter :: beta_const_ct_t_exact_ca = 1.847372081698051d-015, beta_const_pt_t_exact_ca = 1.805738718274608d-015
real (r14), parameter :: beta_const_t_exact_ca = 1.804871356536619d-015, enthalpy_t_exact_ca = 2.499349648132920d-006
real (r14), parameter :: internal_energy_t_exact_ca = 2.499335096217692d-006
real (r14), parameter :: dynamic_enthalpy_t_exact_ca = 2.288943505845964d-007
real (r14), parameter :: isochoric_heat_cap_t_exact_ca = 1.614353095646948d-009, chem_potential_t_exact_ca = 1.317864928296331d-009
real (r14), parameter :: chem_potential_water_t_exact_ca = 4.811790859093890d-007
real (r14), parameter :: chem_potential_salt_t_exact_ca = 4.805315256817266d-007
real (r14), parameter :: helmholtz_energy_t_exact_ca = 2.440137905068696d-007, sound_speed_t_exact_ca = 2.590240910649300d-009
real (r14), parameter :: kappa_t_exact_ca = 1.712677458291044d-021, kappa_const_t_exact_ca = 1.697064424229105d-021
real (r14), parameter :: osmotic_coefficient_t_exact_ca = 3.583799923490005d-013
real (r14), parameter :: osmotic_pressure_t_exact_ca = 1.023465756588848d-009, t_maxdensity_exact_ca = 6.274447628129565d-011
real (r14), parameter :: sa_from_rho_t_exact_ca = 1.304769625676272d-010, t_from_rho_exact_ca = 6.032974120273593d-010
real (r14), parameter :: fdelta_ca = 2.702939055302528d-014, deltasa_atlas_ca = 6.945514042372425d-013


!real (r14)  :: sp, sa, sstar, sr, t, ct, pt, p, p_bs, p_ref 
!real (r14)  :: long, long_bs, lat, lat_bs, saturation_fraction
real (r14)  :: gsw_sa_from_sp, gsw_sp_from_sa, gsw_fdelta, gsw_ct_from_t
real (r14)  :: gsw_sstar_from_sa, gsw_sstar_from_sp, gsw_deltasa_atlas
real (r14)  :: gsw_sa_from_sp_baltic, gsw_sp_from_sa_baltic
real (r14)  :: gsw_sr_from_sp, gsw_sp_from_sr, gsw_deltasa_from_sp
real (r14)  :: gsw_sa_from_sstar, gsw_sp_from_sstar, gsw_ct_freezing
real (r14)  :: gsw_alpha_wrt_t_exact, gsw_beta_const_t_exact 
real (r14)  :: gsw_ct_from_pt, gsw_rho_t_exact, gsw_enthalpy_t_exact 
real (r14)  :: gsw_entropy_from_t, gsw_kappa_t_exact, gsw_t_freezing
real (r14)  :: gsw_pot_rho_t_exact, gsw_pt_from_t, gsw_pt_from_ct, gsw_t_from_ct
real (r14)  :: gsw_specvol_t_exact, gsw_sound_speed_t_exact, gsw_internal_energy
real (r14)  :: gsw_rho, gsw_alpha, gsw_beta, gsw_pt0_from_t, gsw_dynamic_enthalpy  
real (r14)  :: gsw_specvol, gsw_specvol_anom, gsw_enthalpy, gsw_sound_speed
real (r14)  :: gsw_latentheat_melting, gsw_latentheat_evap_ct, gsw_latentheat_evap_t
real (r14)  :: gsw_sa_from_rho, gsw_kappa, gsw_sp_from_c, gsw_c_from_sp
real (r14)  :: gsw_adiabatic_lapse_rate_from_ct, gsw_z_from_p, gsw_grav
real (r14)  :: gsw_cabbeling, gsw_thermobaric, gsw_alpha_on_beta
real (r14)  :: gsw_sigma0, gsw_sigma1, gsw_sigma2, gsw_sigma3, gsw_sigma4
real (r14)  :: gsw_rho_first_derivaties

real (r14)  :: sa_from_sp, sp_from_sa, fdelta, ct_from_t
real (r14)  :: sstar_from_sa, sstar_from_sp, deltasa_atlas
real (r14)  :: sa_from_sp_baltic, sp_from_sa_baltic
real (r14)  :: sr_from_sp, sp_from_sr, deltasa_from_sp
real (r14)  :: sa_from_sstar, sp_from_sstar, ct_freezing
real (r14)  :: alpha_wrt_t_exact, beta_const_t_exact, cp_t_exact 
real (r14)  :: ct_from_pt, rho_t_exact, enthalpy_t_exact 
real (r14)  :: entropy_from_t, kappa_t_exact, t_freezing
real (r14)  :: pot_rho_t_exact, pt_from_t, pt_from_ct, t_from_ct
real (r14)  :: specvol_t_exact, sound_speed_t_exact, internal_energy
real (r14)  :: rho, alpha, beta, pt0_from_t, dynamic_enthalpy  
real (r14)  :: specvol, specvol_anom, enthalpy, sound_speed
real (r14)  :: latentheat_melting, latentheat_evap_ct, latentheat_evap_t
real (r14)  :: sa_from_rho, kappa, sp_from_c, c_from_sp
real (r14)  :: adiabatic_lapse_rate_from_ct, z_from_p, grav
real (r14)  :: cabbeling, thermobaric, alpha_on_beta
real (r14)  :: sigma0, sigma1, sigma2, sigma3, sigma4
real (r14)  :: drho_dsa, drho_dct, drho_dp, drho_dsa_error, drho_dct_error
real (r14)  :: drho_dp_error

real (r14), dimension(2) :: n2, p_mid_n2, n2_error, p_mid_n2_error 
real (r14), dimension(2) :: ipvfn2,p_mid_ipvfn2
real (r14), dimension(2) :: ipvfn2_error, p_mid_ipvfn2_error
real (r14), dimension(2) :: tu, rsubrho, p_mid_tursr
real (r14), dimension(2) :: tu_error, rsubrho_error, p_mid_tursr_error
real (r14), dimension(3) :: sa_profile, ct_profile, p_profile, lat_profile

sa_profile(1) = 35.5d0
sa_profile(2) = 35.7d0
sa_profile(3) = 35.6d0
ct_profile(1) = 12.5d0
ct_profile(2) = 15d0
ct_profile(3) = 10d0
p_profile(1) = 00d0
p_profile(2) = 50d0
p_profile(3) = 100d0
lat_profile(1) = 10d0
lat_profile(2) = 10d0
lat_profile(3) = 10d0

gsw_error_flag = 0;

print*; print *, '============================================================================'
print*; print *, ' Gibbs SeaWater (GSW) Oceanographic Toolbox of TEOS-10 version 3.01 (Fortran)'
print*; print *, '============================================================================'
print*; print *, ' '
print*; print *, ' These are the check values for the subset of functions that have been ' 
print*; print *, ' converted into FORTRAN 90 from the Gibbs SeaWater (GSW) Oceanographic Toolbox '
print*; print *, ' of TEOS-10 (version 3.03).'
print*; print *, ' '
print*; print *, '-------------------------------------------------------------------------------------------'

print*; print *, 'Practical Salinity, PSS-78'

sp_from_c = abs(gsw_sp_from_c(c,t,p) - 35.500961780774482d0)
if (sp_from_c.lt.sp_from_c_ca) then
  print*; print *, 'gsw_sp_from_c  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sp_from_c  ........  failed'
  print*;
endif

c_from_sp = abs(gsw_c_from_sp(sp,t,p) - 43.598945605280484d0)
if (c_from_sp.lt.c_from_sp_ca) then
  print*; print *, 'gsw_c_from_sp  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_c_from_sp  ........  failed'
  print*;
endif

print*; print *, ' '
print*; print *, 'Absolute Salinity, Preformed Salinity and Conservative Temperature'

sa_from_sp = abs(gsw_sa_from_sp(sp,p,long,lat) - 35.671358392019094d0)
if (sa_from_sp.lt.sa_from_sp_ca) then
  print*; print *, 'gsw_sa_from_sp  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sa_from_sp  ........  failed'
  print*;
endif

sstar_from_sp = abs(gsw_sstar_from_sp(sa,p,long,lat) - 35.866946753006239d0)
if (sstar_from_sp.lt.sstar_from_sp_ca) then
  print*; print *, 'gsw_sstar_from_sp  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sstar_from_sp  ........  failed'
  print*;
endif
    
ct_from_t = abs(gsw_ct_from_t(sa,t,p) - 14.930280459895560d0)
if (ct_from_t.lt.ct_from_t_ca) then
  print*; print *, 'gsw_ct_from_t  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_ct_from_t  ........  failed'
  print*;
endif

print*; print *, ' '
print*; print *, 'other conversions between temperatures, salinities, entropy, pressure and height'

deltasa_from_sp = abs(gsw_deltasa_from_sp(sp,p,long,lat) - 3.96067773336028495d-3)
if (deltasa_from_sp.lt.deltasa_from_sp_ca) then
  print*; print *, 'gsw_deltasa_from_sp  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_deltasa_from_sp  ........  failed'
  print*;
endif

sr_from_sp = abs(gsw_sr_from_sp(sp) - 35.667397714285734d0)
if (sr_from_sp.lt.sr_from_sp_ca) then
  print*; print *, 'gsw_sr_from_sp  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sr_from_sp  ........  failed'
  print*;
endif

sp_from_sr = abs(gsw_sp_from_sr(sr) - 35.333387933015295d0)
if (sp_from_sr.lt.sp_from_sr_ca) then
  print*; print *, 'gsw_sp_from_sr  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sp_from_sr  ........  failed'
  print*;
endif

sp_from_sa = abs(gsw_sp_from_sa(sa,p,long,lat) - 35.528504019167094d0)
if (sp_from_sa.lt.sp_from_sa_ca) then
  print*; print *, 'gsw_sp_from_sa  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sp_from_sa  ........  failed'
  print*;
endif

sstar_from_sa = abs(gsw_sstar_from_sa(sa,p,long,lat) - 35.694648791860907d0)
if (sstar_from_sa.lt.sstar_from_sa_ca) then
  print*; print *, 'gsw_sstar_from_sa  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sstar_from_sa  ........  failed'
  print*;
endif

sp_from_sstar = abs(gsw_sp_from_sstar(sstar,p,long,lat) - 35.334761242083573d0)
if (sp_from_sstar.lt.sp_from_sstar_ca) then
  print*; print *, 'gsw_sp_from_sstar  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sp_from_sstar  ........  failed'
  print*;
endif

sa_from_sstar = abs(gsw_sa_from_sstar(sstar,p,long,lat) - 35.505322027120805d0)
if (sa_from_sstar.lt.sa_from_sstar_ca) then
  print*; print *, 'gsw_sa_from_sstar  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sa_from_sstar  ........  failed'
  print*;
endif

pt_from_ct = abs(gsw_pt_from_ct(sa,ct) - 20.023899375975017d0)
if (pt_from_ct.lt.pt_from_ct_ca) then
  print*; print *, 'gsw_pt_from_ct  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_pt_from_ct  ........  failed'
  print*;
endif

t_from_ct = abs(gsw_t_from_ct(sa,ct,p) - 20.079820359223014d0)
if (t_from_ct.lt.t_from_ct_ca) then
  print*; print *, 'gsw_t_from_ct  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_t_from_ct  ........  failed'
  print*;
endif

ct_from_pt = abs(gsw_ct_from_pt(sa,pt) - 14.976021403957613d0)
if (ct_from_pt.lt.ct_from_pt_ca) then
  print*; print *, 'gsw_ct_from_pt  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_ct_from_pt  ........  failed'
  print*;
endif

pt0_from_t = abs(gsw_pt0_from_t(sa,t,p) - 14.954241363902305d0)
if (pt0_from_t.lt.pt0_from_t_ca) then
  print*; print *, 'gsw_pt0_from_t  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_pt0_from_t  ........  failed'
  print*;
endif

pt_from_t = abs(gsw_pt_from_t(sa,t,p,p_ref) - 14.969381237883740d0)
if (pt_from_t.lt.pt_from_t_ca) then
  print*; print *, 'gsw_pt_from_t  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_pt_from_t  ........  failed'
  print*;
endif

z_from_p = abs(gsw_z_from_p(p,lat) + 2.980161553316402d2)
if (z_from_p.lt.z_from_p_ca) then
  print*; print *, 'gsw_z_from_p  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_z_from_p  ........  failed'
  print*;
endif

entropy_from_t = abs(gsw_entropy_from_t(sa,t,p) - 212.30166821093002d0)
if (entropy_from_t.lt.entropy_from_t_ca) then
  print*; print *, 'gsw_entropy_from_t  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_entropy_from_t  ........  failed'
  print*;
endif

adiabatic_lapse_rate_from_ct = abs(gsw_adiabatic_lapse_rate_from_ct(sa,ct,p) - 1.877941744191155d-8)
if (adiabatic_lapse_rate_from_ct.lt.adiabatic_lapse_rate_from_ct_ca) then
  print*; print *, 'gsw_adiabatic_lapse_rate_from_ct  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_adiabatic_lapse_rate_from_ct  ........  failed'
  print*;
endif

print*; print *, ' '
print*; print *, 'density and enthalpy, based on the 48-term expression for density'

rho = abs(gsw_rho(sa,ct,p) - 1026.4562376198473d0)
if (rho.lt.rho_ca) then
  print*; print *, 'gsw_rho  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_rho  ........  failed'
  print*;
endif

alpha = abs(gsw_alpha(sa,ct,p) - 2.62460550806784356d-4)
if (alpha.lt.alpha_ca) then
  print*; print *, 'gsw_alpha  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_alpha  ........  failed'
  print*;
endif

beta = abs(gsw_beta(sa,ct,p) - 7.29314455934463365d-4)
if (beta.lt.beta_ca) then
  print*; print *, 'gsw_beta  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_beta  ........  failed'
  print*;
endif

alpha_on_beta = abs(gsw_alpha_on_beta(sa,ct,p) - 0.359872958325632d0)
if (alpha_on_beta.lt.alpha_on_beta_ca) then
  print*; print *, 'gsw_alpha_on_beta  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_alpha_on_beta  ........  failed'
  print*;
endif

call gsw_rho_first_derivatives(sa,ct,p,drho_dsa,drho_dct,drho_dp)
drho_dsa_error = abs(drho_dsa - 0.748609372480258d0)
drho_dct_error = abs(drho_dct + 0.269404269504765d0)
drho_dp_error = abs(drho_dp - 4.287533235942749d-7)
if (drho_dsa_error.lt.drho_dsa_ca.and.drho_dct_error.lt.drho_dct_ca.and.drho_dp_error.lt.drho_dp_ca) then
  print*; print *, 'gsw_rho_first_derivatives  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_rho_first_derivatives  ........  failed'
  print*;
endif

specvol = abs(gsw_specvol(sa,ct,p) - 9.74225654586897711d-4)
if (specvol.lt.specvol_ca) then
  print*; print *, 'gsw_specvol  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_specvol  ........  failed'
  print*;
endif

specvol_anom = abs(gsw_specvol_anom(sa,ct,p) - 2.90948181201264571d-6)
if (specvol_anom.lt.specvol_anom_ca) then
  print*; print *, 'gsw_specvol_anom  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_specvol_anom  ........  failed'
  print*;
endif

sigma0 = abs(gsw_sigma0(sa,ct) - 25.165674636323047d0)
if (sigma0.lt.sigma0_ca) then
  print*; print *, 'gsw_sigma0  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sigma0  ........  failed'
  print*;
endif

sigma1 = abs(gsw_sigma1(sa,ct) - 29.434338510752923d0)
if (sigma1.lt.sigma1_ca) then
  print*; print *, 'gsw_sigma1  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sigma1  ........  failed'
  print*;
endif

sigma2 = abs(gsw_sigma2(sa,ct) - 33.609842926904093d0)
if (sigma2.lt.sigma2_ca) then
  print*; print *, 'gsw_sigma2  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sigma2  ........  failed'
  print*;
endif

sigma3 = abs(gsw_sigma3(sa,ct) - 37.695147569371784d0)
if (sigma3.lt.sigma3_ca) then
  print*; print *, 'gsw_sigma3  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sigma3  ........  failed'
  print*;
endif

sigma4 = abs(gsw_sigma4(sa,ct) - 41.693064726656303d0)
if (sigma4.lt.sigma4_ca) then
  print*; print *, 'gsw_sigma4  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sigma4  ........  failed'
  print*;
endif

sound_speed = abs(gsw_sound_speed(sa,ct,p) - 1527.2011773569989d0)
if (sound_speed.lt.sound_speed_ca) then
  print*; print *, 'gsw_sound_speed  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sound_speed  ........  failed'
  print*;
endif

kappa = abs(gsw_kappa(sa,ct,p) -  4.177024873349404d-10)
if (kappa.lt.kappa_ca) then
  print*; print *, 'gsw_kappa  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_kappa  ........  failed'
  print*;
endif

cabbeling = abs(gsw_cabbeling(sa,ct,p) -  9.463053321129075d-6)
if (cabbeling.lt.cabbeling_ca) then
  print*; print *, 'gsw_cabbeling  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_cabbeling  ........  failed'
endif

thermobaric = abs(gsw_thermobaric(sa,ct,p) -  1.739078662082863d-12)
if (thermobaric.lt.thermobaric_ca) then
  print*; print *, 'gsw_thermobaric  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_thermobaric  ........  failed'
  print*;
endif

internal_energy = abs(gsw_internal_energy(sa,ct,p) - 79740.482561720783d0)
if (internal_energy.lt.internal_energy_ca) then
  print*; print *, 'gsw_internal_energy  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_internal_energy  ........  failed'
  print*;
endif

enthalpy = abs(gsw_enthalpy(sa,ct,p) - 82761.872939932495d0)
if (enthalpy.lt.enthalpy_ca) then
  print*; print *, 'gsw_enthalpy  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_enthalpy  ........  failed'
  print*;
endif

dynamic_enthalpy = abs(gsw_dynamic_enthalpy(sa,ct,p) - 2924.5137975399025d0)
if (dynamic_enthalpy.lt.dynamic_enthalpy_ca) then
  print*; print *, 'gsw_dynamic_enthalpy  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_dynamic_enthalpy  ........  failed'
  print*;
endif

rho = gsw_rho(sa,ct,p)
sa_from_rho = abs(gsw_sa_from_rho(rho,ct,p) - sa) 
if (sa_from_rho.lt.sa_from_rho_ca) then
  print*; print *, 'gsw_sa_from_rho  ........  passed'
else
  gsw_error_flag = 1;
  print*; print *, 'gsw_sa_from_rho  ........  failed'
  print*;
endif

print*; print *, ' '
print*; print *, 'water column properties, based on the 48-term expression for density'

call gsw_nsquared(sa_profile,ct_profile,p_profile,lat_profile,nz,n2,p_mid_n2)

n2_error(1) = abs(n2(1) + 0.070960392693051d-3)
n2_error(2) = abs(n2(2) - 0.175435821615983d-3)
p_mid_n2_error(1) = abs(p_mid_n2(1) - 25d0)
p_mid_n2_error(2) = abs(p_mid_n2(2) - 75d0)
if (n2_error(1).lt.n2_ca.and.n2_error(2).lt.n2_ca.and.p_mid_n2_error(1).lt.p_mid_n2_ca.and.p_mid_n2_error(2).lt.p_mid_n2_ca) then
  print*; print *, 'gsw_nsquared  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_nsquared  ........  failed'
  print*;
endif

call gsw_turner_rsubrho(sa_profile,ct_profile,p_profile,nz,tu,rsubrho,p_mid_tursr)

tu_error(1) = abs(tu(1) + 1.187243981606485d2)
tu_error(2) = abs(tu(2) - 0.494158257088517d2)
rsubrho_error(1) = abs(rsubrho(1) - 3.425146897090065d0)
rsubrho_error(2) = abs(rsubrho(2) - 12.949399443139164d0)
p_mid_tursr_error(1) = abs(p_mid_tursr(1) - 25d0)
p_mid_tursr_error(2) = abs(p_mid_tursr(2) - 75d0)

 if (tu_error(1).lt.tu_ca.and.tu_error(2).lt.tu_ca.and. &
    rsubrho_error(1).lt.rsubrho_ca.and.rsubrho_error(2).lt.rsubrho_ca.and. &
    p_mid_tursr_error(1).lt.p_mid_tursr_ca.and.p_mid_tursr_error(2).lt.p_mid_tursr_ca) then
  print*; print *, 'gsw_turner_rsubrho  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_turner_rsubrho  ........  failed'
  print*;
endif

call gsw_ipv_vs_fnsquared_ratio(sa_profile,ct_profile,p_profile,nz,ipvfn2,p_mid_ipvfn2)

ipvfn2_error(1) = abs(ipvfn2(1) - 0.996783975249010d0)
ipvfn2_error(2) = abs(ipvfn2(2) - 0.992112251478320d0)
p_mid_ipvfn2_error(1) = abs(p_mid_ipvfn2(1) - 25d0)
p_mid_ipvfn2_error(2) = abs(p_mid_ipvfn2(2) - 75d0)
if (ipvfn2_error(1).lt.ipvfn2_ca.and. &
    ipvfn2_error(2).lt.ipvfn2_ca.and. & 
	p_mid_ipvfn2_error(1).lt.p_mid_ipvfn2_ca.and. &
	p_mid_ipvfn2_error(2).lt.p_mid_ipvfn2_ca) then
  print*; print *, 'gsw_ipv_vs_fnsquared_ratio  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_ipv_vs_fnsquared_ratio  ........  failed'
  print*;
endif

print*; print *, ' '
print*; print *, 'freezing temperatures'

ct_freezing = abs(gsw_ct_freezing(sa,p,saturation_fraction) + 2.1801450326174852d0)
if (ct_freezing.lt.ct_freezing_ca) then
  print*; print *, 'gsw_ct_freezing  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_ct_freezing  ........  failed'
  print*;
endif

t_freezing = abs(gsw_t_freezing(sa,p,saturation_fraction) + 2.1765521998023516d0)
if (t_freezing.lt.t_freezing_ca) then
  print*; print *, 'gsw_t_freezing  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_t_freezing  ........  failed'
  print*;
endif

print*; print *, ' '
print*; print *, 'isobaric melting enthalpy and isobaric evaporation enthalpy'
latentheat_melting = abs(gsw_latentheat_melting(sa,p) - 329330.54839618353d0)
if (latentheat_melting.lt.latentheat_melting_ca) then
  print*; print *, 'gsw_latentheat_melting  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_latentheat_melting  ........  failed'
  print*;
endif

latentheat_evap_ct = abs(gsw_latentheat_evap_ct(sa,ct) - 2450871.0228523901d0)
if (latentheat_evap_ct.lt.latentheat_evap_ct_ca) then
  print*; print *, 'gsw_latentheat_evap_ct  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_latentheat_evap_ct  ........  failed'
  print*;
endif

latentheat_evap_t = abs(gsw_latentheat_evap_t(sa,t) - 2462848.2895522709d0)
if (latentheat_evap_t.lt.latentheat_evap_t_ca) then
  print*; print *, 'gsw_latentheat_evap_t  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_latentheat_evap_t  ........  failed'
  print*;
endif

print*; print *, ' '
print*; print *, 'planet Earth properties'
grav = abs(gsw_grav(lat,p) - 9.784910108550843d0)
if (grav.lt.grav_ca) then
  print*; print *, 'gsw_grav  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_grav  ........  failed'
  print*;
endif

print*; print *, ' '
print*; print *, 'basic thermodynamic properties in terms of in-situ t, based on the exact Gibbs function'
rho_t_exact = abs(gsw_rho_t_exact(sa,t,p) - 1027.7128170207150d0)
if (rho_t_exact.lt.rho_t_exact_ca) then
  print*; print *, 'gsw_rho_t_exact  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_rho_t_exact  ........  failed'
  print*;
endif

pot_rho_t_exact = abs(gsw_pot_rho_t_exact(sa,t,p,p_ref) - 1026.8362655887486d0)
if (pot_rho_t_exact.lt.pot_rho_t_exact_ca) then
  print*; print *, 'gsw_pot_rho_t_exact  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_pot_rho_t_exact  ........  failed'
  print*;
endif

alpha_wrt_t_exact = abs(gsw_alpha_wrt_t_exact(sa,t,p) - 2.19066952410728916d-4)
if (alpha_wrt_t_exact.lt.alpha_wrt_t_exact_ca) then
  print*; print *, 'gsw_alpha_wrt_t_exact  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_alpha_wrt_t_exact  ........  failed'
  print*;
endif

beta_const_t_exact = abs(gsw_beta_const_t_exact(sa,t,p) - 7.44744841648729426d-4)
if (beta_const_t_exact.lt.beta_const_t_exact_ca) then
  print*; print *, 'gsw_beta_const_t_exact  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_beta_const_t_exact  ........  failed'
  print*;
endif

specvol_t_exact = abs(gsw_specvol_t_exact(sa,t,p) - 9.73034473676164815d-4)
if (specvol_t_exact.lt.specvol_t_exact_ca) then
  print*; print *, 'gsw_specvol_t_exact  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_specvol_t_exact  ........  failed'
  print*;
endif

sound_speed_t_exact = abs(gsw_sound_speed_t_exact(sa,t,p) - 1512.2053940303056d0)
if (sound_speed_t_exact.lt.sound_speed_t_exact_ca) then
  print*; print *, 'gsw_sound_speed_t_exact  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sound_speed_t_exact  ........  failed'
  print*;
endif

kappa_t_exact = abs(gsw_kappa_t_exact(sa,t,p) - 4.25506953386609075d-010)
if (kappa_t_exact.lt.kappa_t_exact_ca) then
  print*; print *, 'gsw_kappa_t_exact  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_kappa_t_exact  ........  failed'
  print*;
endif

enthalpy_t_exact = abs(gsw_enthalpy_t_exact(sa,t,p) - 62520.680485510929d0)
if (enthalpy_t_exact.lt.enthalpy_t_exact_ca) then
  print*; print *, 'gsw_enthalpy_t_exact  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_enthalpy_t_exact  ........  failed'
  print*;
endif

print*; print *, ' '
print*; print *, 'library functions of the GSW toolbox'

deltasa_atlas = abs(gsw_deltasa_atlas(p,long,lat) - 3.87660373016291727d-3)
if (deltasa_atlas.lt.deltasa_atlas_ca) then
  print*; print *, 'gsw_deltasa_atlas  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_deltasa_atlas  ........  failed'
  print*;
endif

fdelta = abs(gsw_fdelta(p,long,lat) - 1.49916256924158942d-004)
if (fdelta.lt.fdelta_ca) then
  print*; print *, 'gsw_fdelta  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_fdelta  ........  failed'
  print*;
endif

sa_from_sp_baltic = abs(gsw_sa_from_sp_baltic(sp,long_bs,lat_bs) - 35.666154857142850d0)
if (sa_from_sp_baltic.lt.sa_from_sp_ca) then
  print*; print *, 'gsw_sa_from_sp_baltic  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sa_from_sp_baltic  ........  failed'
  print*;
endif

sp_from_sa_baltic = abs(gsw_sp_from_sa_baltic(sa,long_bs,lat_bs) - 35.533769845749660d0)
if (sp_from_sa_baltic.lt.sp_from_sa_ca) then
  print*; print *, 'gsw_sp_from_sa_baltic  ........  passed'
else
  gsw_error_flag = 1;
  print*;
  print*; print *, 'gsw_sp_from_sa_baltic  ........  failed'
  print*;
endif

if (gsw_error_flag.eq.1) then
  print*;
  print*; print *, 'Your installation of the Gibbs SeaWater (GSW) Oceanographic Toolbox has errors !'
else  
  print*;
  print*; print *, 'Well done! The gsw_check_fuctions confirms that the'
  print*; print *, 'Gibbs SeaWater (GSW) Oceanographic Toolbox is installed correctly.'
endif

stop
end

           