module gsw_mod_toolbox

implicit none

public :: gsw_add_barrier
public :: gsw_add_mean
public :: gsw_adiabatic_lapse_rate_from_ct
public :: gsw_adiabatic_lapse_rate_ice
public :: gsw_alpha
public :: gsw_alpha_on_beta
public :: gsw_alpha_wrt_t_exact
public :: gsw_alpha_wrt_t_ice
public :: gsw_beta_const_t_exact
public :: gsw_beta
public :: gsw_brinesa_ct
public :: gsw_brinesa_ct_poly
public :: gsw_brinesa_estimate
public :: gsw_brinesa_t
public :: gsw_brinesa_t_poly
public :: gsw_cabbeling
public :: gsw_c_from_sp
public :: gsw_chem_potential_water_ice
public :: gsw_chem_potential_water_t_exact
public :: gsw_cp_ice
public :: gsw_ct_first_derivatives
public :: gsw_ct_first_derivatives_wrt_t_exact
public :: gsw_ct_freezing_derivative_poly
public :: gsw_ct_freezing_exact
public :: gsw_ct_freezing
public :: gsw_ct_freezing_first_derivatives
public :: gsw_ct_freezing_poly
public :: gsw_ct_from_enthalpy
public :: gsw_ct_from_entropy
public :: gsw_ct_from_pt
public :: gsw_ct_from_rho
public :: gsw_ct_from_t
public :: gsw_ct_maxdensity
public :: gsw_ct_second_derivatives
public :: gsw_deltasa_atlas
public :: gsw_deltasa_from_sp
public :: gsw_dilution_coefficient_t_exact
public :: gsw_dynamic_enthalpy
public :: gsw_enthalpy_diff
public :: gsw_enthalpy
public :: gsw_enthalpy_first_derivatives
public :: gsw_enthalpy_ice
public :: gsw_enthalpy_second_derivatives_ct_exact
public :: gsw_enthalpy_second_derivatives
public :: gsw_enthalpy_sso_0_p
public :: gsw_enthalpy_t_exact
public :: gsw_entropy_first_derivatives
public :: gsw_entropy_from_pt
public :: gsw_entropy_from_t
public :: gsw_entropy_ice
public :: gsw_entropy_part
public :: gsw_entropy_part_zerop
public :: gsw_entropy_second_derivatives
public :: gsw_fdelta
public :: gsw_frazil_ratios
public :: gsw_gibbs
public :: gsw_gibbs_ice
public :: gsw_gibbs_ice_part_t
public :: gsw_gibbs_ice_pt0
public :: gsw_gibbs_ice_pt0_pt0
public :: gsw_gibbs_pt0_pt0
public :: gsw_grav
public :: gsw_helmholtz_energy_ice
public :: gsw_hill_ratio_at_sp2
public :: gsw_ice_fraction_to_freeze_seawater
public :: gsw_internal_energy
public :: gsw_internal_energy_ice
public :: gsw_ipv_vs_fnsquared_ratio
public :: gsw_kappa_const_t_ice
public :: gsw_kappa
public :: gsw_kappa_ice
public :: gsw_kappa_t_exact
public :: gsw_latentheat_evap_ct
public :: gsw_latentheat_evap_t
public :: gsw_latentheat_melting
public :: gsw_melting_ice_equilibrium_sa_ct_ratio
public :: gsw_melting_ice_into_seawater
public :: gsw_melting_ice_sa_ct_ratio
public :: gsw_melting_seaice_equilibrium_sa_ct_ratio
public :: gsw_melting_seaice_into_seawater
public :: gsw_melting_seaice_sa_ct_ratio
public :: gsw_nsquared
public :: gsw_pot_enthalpy_from_pt_ice
public :: gsw_pot_enthalpy_from_pt_ice_poly
public :: gsw_pot_rho_t_exact
public :: gsw_pressure_coefficient_ice
public :: gsw_pressure_freezing_ct
public :: gsw_pt0_cold_ice_poly
public :: gsw_pt0_from_t
public :: gsw_pt0_from_t_ice
public :: gsw_pt_first_derivatives
public :: gsw_pt_from_ct
public :: gsw_pt_from_entropy
public :: gsw_pt_from_pot_enthalpy_ice
public :: gsw_pt_from_pot_enthalpy_ice_poly_dh
public :: gsw_pt_from_pot_enthalpy_ice_poly
public :: gsw_pt_from_t
public :: gsw_pt_from_t_ice
public :: gsw_pt_second_derivatives
public :: gsw_rho_alpha_beta
public :: gsw_rho
public :: gsw_rho_first_derivatives
public :: gsw_rho_ice
public :: gsw_rho_t_exact
public :: gsw_saar
public :: gsw_sa_from_rho
public :: gsw_sa_from_sp_baltic
public :: gsw_sa_from_sp
public :: gsw_sa_from_sstar
public :: gsw_sa_p_inrange
public :: gsw_seaice_fraction_to_freeze_seawater
public :: gsw_sigma0
public :: gsw_sigma1
public :: gsw_sigma2
public :: gsw_sigma3
public :: gsw_sigma4
public :: gsw_sound_speed
public :: gsw_sound_speed_ice
public :: gsw_sound_speed_t_exact
public :: gsw_specvol_anom
public :: gsw_specvol
public :: gsw_specvol_ice
public :: gsw_specvol_sso_0_p
public :: gsw_specvol_t_exact
public :: gsw_sp_from_c
public :: gsw_sp_from_sa_baltic
public :: gsw_sp_from_sa
public :: gsw_sp_from_sk
public :: gsw_sp_from_sr
public :: gsw_sp_from_sstar
public :: gsw_sr_from_sp
public :: gsw_sstar_from_sa
public :: gsw_sstar_from_sp
public :: gsw_t_deriv_chem_potential_water_t_exact
public :: gsw_t_freezing_derivative_poly
public :: gsw_t_freezing_exact
public :: gsw_t_freezing
public :: gsw_t_freezing_first_derivatives
public :: gsw_t_freezing_poly
public :: gsw_t_from_ct
public :: gsw_t_from_pt0_ice
public :: gsw_thermobaric
public :: gsw_turner_rsubrho
public :: gsw_util_indx
public :: gsw_util_xinterp1
public :: gsw_z_from_p

interface

    pure subroutine gsw_add_barrier (input_data, long, lat, long_grid, &
                                  lat_grid, dlong_grid, dlat_grid, output_data)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: long, lat, long_grid, lat_grid, dlong_grid
    real (r14), intent(in) :: dlat_grid
    real (r14), intent(in), dimension(4) :: input_data
    real (r14), intent(out), dimension(4) :: output_data
    end subroutine gsw_add_barrier
    
    pure subroutine gsw_add_mean (data_in, data_out)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in), dimension(4) :: data_in
    real (r14), intent(out), dimension(4) :: data_out
    end subroutine 
    
    elemental function gsw_adiabatic_lapse_rate_from_ct (sa, ct, p) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p 
    real (r14) :: gsw_adiabatic_lapse_rate_from_ct
    end function gsw_adiabatic_lapse_rate_from_ct
    
    elemental function gsw_adiabatic_lapse_rate_ice (t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p
    real (r14) :: gsw_adiabatic_lapse_rate_ice
    end function gsw_adiabatic_lapse_rate_ice
    
    elemental function gsw_alpha (sa, ct, p)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p  
    real (r14) :: gsw_alpha
    end function gsw_alpha
    
    elemental function gsw_alpha_on_beta (sa, ct, p)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p  
    real (r14) :: gsw_alpha_on_beta
    end function gsw_alpha_on_beta
    
    elemental function gsw_alpha_wrt_t_exact (sa, t, p) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t, p 
    real (r14) :: gsw_alpha_wrt_t_exact
    end function gsw_alpha_wrt_t_exact
    
    elemental function gsw_alpha_wrt_t_ice (t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p
    real (r14) :: gsw_alpha_wrt_t_ice
    end function gsw_alpha_wrt_t_ice
    
    elemental function gsw_beta_const_t_exact (sa, t, p)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t, p  
    real (r14) :: gsw_beta_const_t_exact
    end function gsw_beta_const_t_exact
    
    elemental function gsw_beta (sa, ct, p) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p 
    real (r14) :: gsw_beta
    end function gsw_beta
    
    elemental function gsw_brinesa_ct (ct, p, saturation_fraction)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: ct, p, saturation_fraction
    real (r14) :: gsw_brinesa_ct
    end function gsw_brinesa_ct
    
    elemental function gsw_brinesa_ct_poly (ct, p, saturation_fraction)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: ct, p, saturation_fraction
    real (r14) :: gsw_brinesa_ct_poly
    end function gsw_brinesa_ct_poly
    
    elemental function gsw_brinesa_estimate (p, saturation_fraction, ct, t)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: p, saturation_fraction
    real (r14), intent(in), optional :: ct, t
    real (r14) :: gsw_brinesa_estimate
    end function gsw_brinesa_estimate
    
    elemental function gsw_brinesa_t (t, p, saturation_fraction)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p, saturation_fraction
    real (r14) :: gsw_brinesa_t
    end function gsw_brinesa_t
    
    elemental function gsw_brinesa_t_poly (t, p, saturation_fraction)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p, saturation_fraction
    real (r14) :: gsw_brinesa_t_poly
    end function gsw_brinesa_t_poly
    
    elemental function gsw_cabbeling (sa, ct, p)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p  
    real (r14) :: gsw_cabbeling
    end function gsw_cabbeling
    
    elemental function gsw_c_from_sp (sp, t, p)       
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sp, t, p       
    real (r14) :: gsw_c_from_sp
    end function gsw_c_from_sp
    
    elemental function gsw_chem_potential_water_ice (t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p
    real (r14) :: gsw_chem_potential_water_ice
    end function gsw_chem_potential_water_ice
    
    elemental function gsw_chem_potential_water_t_exact (sa, t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t, p
    real (r14) :: gsw_chem_potential_water_t_exact
    end function gsw_chem_potential_water_t_exact
    
    elemental function gsw_cp_ice (t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p
    real (r14) :: gsw_cp_ice
    end function gsw_cp_ice
    
    elemental subroutine gsw_ct_first_derivatives (sa, pt, ct_sa, ct_pt)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, pt
    real (r14), intent(out), optional :: ct_sa, ct_pt
    end subroutine gsw_ct_first_derivatives
    
    elemental subroutine gsw_ct_first_derivatives_wrt_t_exact (sa, t, p, &
    				       ct_sa_wrt_t, ct_t_wrt_t, ct_p_wrt_t)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t, p
    real (r14), intent(out) :: ct_p_wrt_t, ct_sa_wrt_t, ct_t_wrt_t
    end subroutine gsw_ct_first_derivatives_wrt_t_exact
    
    elemental function gsw_ct_freezing_derivative_poly (sa, p, &
                                                        saturation_fraction)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p, saturation_fraction
    real (r14) :: gsw_ct_freezing_derivative_poly
    end function gsw_ct_freezing_derivative_poly
    
    elemental function gsw_ct_freezing_exact (sa, p, saturation_fraction)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p, saturation_fraction
    real (r14) :: gsw_ct_freezing_exact
    end function gsw_ct_freezing_exact
    
    elemental function gsw_ct_freezing (sa, p, saturation_fraction, exact)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p, saturation_fraction
    logical, intent(in), optional :: exact
    real (r14) :: gsw_ct_freezing
    end function gsw_ct_freezing
    
    elemental subroutine gsw_ct_freezing_first_derivatives (sa, p, &
                              saturation_fraction, ctfreezing_sa, ctfreezing_p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p, saturation_fraction
    real (r14), intent(out), optional :: ctfreezing_sa, ctfreezing_p
    end subroutine gsw_ct_freezing_first_derivatives
    
    elemental function gsw_ct_freezing_poly (sa, p, saturation_fraction)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p, saturation_fraction
    real (r14) :: gsw_ct_freezing_poly
    end function gsw_ct_freezing_poly
    
    elemental function gsw_ct_from_enthalpy (sa, h, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, h, p
    real (r14) :: gsw_ct_from_enthalpy
    end function gsw_ct_from_enthalpy
    
    elemental function gsw_ct_from_entropy (sa, entropy)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, entropy
    real (r14) :: gsw_ct_from_entropy
    end function gsw_ct_from_entropy
    
    elemental function gsw_ct_from_pt (sa, pt) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, pt 
    real (r14) :: gsw_ct_from_pt
    end function gsw_ct_from_pt
    
    elemental subroutine gsw_ct_from_rho (rho, sa, p, ct, ct_multiple)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: rho, sa, p
    real (r14), intent(out) :: ct
    real (r14), intent(out), optional :: ct_multiple
    end subroutine gsw_ct_from_rho
    
    elemental function gsw_ct_from_t (sa, t, p) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t, p 
    real (r14) :: gsw_ct_from_t
    end function gsw_ct_from_t
    
    elemental function gsw_ct_maxdensity (sa, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p
    real (r14) :: gsw_ct_maxdensity
    end function gsw_ct_maxdensity
    
    elemental subroutine gsw_ct_second_derivatives (sa, pt, ct_sa_sa, ct_sa_pt, &
                                                    ct_pt_pt)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, pt
    real (r14), intent(out), optional :: ct_sa_sa, ct_sa_pt, ct_pt_pt
    end subroutine gsw_ct_second_derivatives
    
    elemental function gsw_deltasa_atlas (p, long, lat)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: p, long, lat
    real (r14) :: gsw_deltasa_atlas
    end function gsw_deltasa_atlas
    
    elemental function gsw_deltasa_from_sp (sp, p, long, lat) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sp, p, long, lat 
    real (r14) :: gsw_deltasa_from_sp
    end function gsw_deltasa_from_sp
    
    elemental function gsw_dilution_coefficient_t_exact (sa, t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t, p
    real (r14) :: gsw_dilution_coefficient_t_exact
    end function gsw_dilution_coefficient_t_exact
    
    elemental function gsw_dynamic_enthalpy (sa, ct, p) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p 
    real (r14) :: gsw_dynamic_enthalpy
    end function gsw_dynamic_enthalpy
    
    elemental function gsw_enthalpy_diff (sa, ct, p_shallow, p_deep)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p_shallow, p_deep
    real (r14) :: gsw_enthalpy_diff
    end function gsw_enthalpy_diff
    
    elemental function gsw_enthalpy (sa, ct, p)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p  
    real (r14) :: gsw_enthalpy
    end function gsw_enthalpy
    
    elemental subroutine gsw_enthalpy_first_derivatives (sa, ct, p, h_sa, h_ct)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p
    real (r14), intent(out), optional :: h_sa, h_ct
    end subroutine gsw_enthalpy_first_derivatives
    
    elemental function gsw_enthalpy_ice (t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p
    real (r14) :: gsw_enthalpy_ice
    end function gsw_enthalpy_ice
    
    elemental subroutine gsw_enthalpy_second_derivatives_ct_exact (sa, ct, p, &
                                                     h_sa_sa, h_sa_ct, h_ct_ct)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p
    real (r14), intent(out), optional :: h_sa_sa, h_sa_ct, h_ct_ct
    end subroutine gsw_enthalpy_second_derivatives_ct_exact
    
    elemental subroutine gsw_enthalpy_second_derivatives (sa, ct, p, &
                                                     h_sa_sa, h_sa_ct, h_ct_ct)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p
    real (r14), intent(out) :: h_sa_sa, h_sa_ct, h_ct_ct
    end subroutine gsw_enthalpy_second_derivatives
    
    elemental function gsw_enthalpy_sso_0_p (p) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: p 
    real (r14) :: gsw_enthalpy_sso_0_p
    end function gsw_enthalpy_sso_0_p
    
    elemental function gsw_enthalpy_t_exact (sa, t, p) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t, p 
    real (r14) :: gsw_enthalpy_t_exact
    end function gsw_enthalpy_t_exact
    
    elemental subroutine gsw_entropy_first_derivatives (sa, ct, eta_sa, eta_ct)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct
    real (r14), intent(out), optional :: eta_sa, eta_ct
    end subroutine gsw_entropy_first_derivatives
    
    elemental function gsw_entropy_from_pt (sa, pt)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, pt
    real (r14) :: gsw_entropy_from_pt
    end function gsw_entropy_from_pt
    
    elemental function gsw_entropy_from_t (sa, t, p) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t, p 
    real (r14) :: gsw_entropy_from_t
    end function gsw_entropy_from_t
    
    elemental function gsw_entropy_ice (t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p
    real (r14) :: gsw_entropy_ice
    end function gsw_entropy_ice
    
    elemental function gsw_entropy_part (sa, t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t, p
    real (r14) :: gsw_entropy_part
    end function gsw_entropy_part
    
    elemental function gsw_entropy_part_zerop (sa, pt0)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, pt0
    real (r14) :: gsw_entropy_part_zerop
    end function gsw_entropy_part_zerop
    
    elemental subroutine gsw_entropy_second_derivatives (sa, ct, eta_sa_sa, &
                                                         eta_sa_ct, eta_ct_ct)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct
    real (r14), intent(out), optional :: eta_sa_sa, eta_sa_ct, eta_ct_ct
    end subroutine gsw_entropy_second_derivatives
    
    elemental function gsw_fdelta (p, long, lat)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: p, long, lat
    real (r14) :: gsw_fdelta
    end function gsw_fdelta
    
    elemental subroutine gsw_frazil_ratios (sa, p, w_ih, dsa_dct_frazil, &
                                            dsa_dp_frazil, dct_dp_frazil)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p, w_ih
    real (r14), intent(out) :: dsa_dct_frazil, dsa_dp_frazil, dct_dp_frazil
    end subroutine gsw_frazil_ratios
    
    elemental function gsw_gibbs (ns, nt, np, sa, t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    integer, intent(in) :: ns, nt, np
    real (r14), intent(in) :: sa, t, p
    real (r14) :: gsw_gibbs
    end function gsw_gibbs
    
    elemental function gsw_gibbs_ice (nt, np, t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    integer, intent(in) :: nt, np
    real (r14), intent(in) :: t, p
    real (r14) :: gsw_gibbs_ice
    end function gsw_gibbs_ice
    
    elemental function gsw_gibbs_ice_part_t (t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p
    real (r14) :: gsw_gibbs_ice_part_t
    end function gsw_gibbs_ice_part_t
    
    elemental function gsw_gibbs_ice_pt0 (pt0)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: pt0
    real (r14) :: gsw_gibbs_ice_pt0
    end function gsw_gibbs_ice_pt0
    
    elemental function gsw_gibbs_ice_pt0_pt0 (pt0)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: pt0
    real (r14) :: gsw_gibbs_ice_pt0_pt0
    end function gsw_gibbs_ice_pt0_pt0
    
    elemental function gsw_gibbs_pt0_pt0 (sa, pt0)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, pt0
    real (r14) :: gsw_gibbs_pt0_pt0
    end function gsw_gibbs_pt0_pt0
    
    elemental function gsw_grav (lat, p)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: lat, p  
    real (r14) :: gsw_grav
    end function gsw_grav
    
    elemental function gsw_helmholtz_energy_ice (t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p
    real (r14) :: gsw_helmholtz_energy_ice
    end function gsw_helmholtz_energy_ice
    
    elemental function gsw_hill_ratio_at_sp2 (t)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t
    real (r14) :: gsw_hill_ratio_at_sp2
    end function gsw_hill_ratio_at_sp2
    
    elemental subroutine gsw_ice_fraction_to_freeze_seawater (sa, ct, p, &
                         saturation_fraction, t_ih, sa_freeze, ct_freeze, w_ih)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p, saturation_fraction, t_ih
    real (r14), intent(out) :: sa_freeze, ct_freeze, w_ih
    end subroutine gsw_ice_fraction_to_freeze_seawater
    
    elemental function gsw_internal_energy (sa, ct, p)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p  
    real (r14) :: gsw_internal_energy
    end function gsw_internal_energy
    
    elemental function gsw_internal_energy_ice (t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p
    real (r14) :: gsw_internal_energy_ice
    end function gsw_internal_energy_ice
    
    pure subroutine gsw_ipv_vs_fnsquared_ratio (sa, ct, p, p_ref, &
                                                ipv_vs_fnsquared_ratio, p_mid)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa(:), ct(:), p(:), p_ref
    real (r14), intent(out) :: ipv_vs_fnsquared_ratio(:), p_mid(:)
    end subroutine gsw_ipv_vs_fnsquared_ratio
    
    elemental function gsw_kappa_const_t_ice (t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p
    real (r14) :: gsw_kappa_const_t_ice
    end function gsw_kappa_const_t_ice
    
    elemental function gsw_kappa (sa, ct, p)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p  
    real (r14) :: gsw_kappa
    end function gsw_kappa
    
    elemental function gsw_kappa_ice (t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p
    real (r14) :: gsw_kappa_ice
    end function gsw_kappa_ice
    
    elemental function gsw_kappa_t_exact (sa, t, p)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t, p  
    real (r14) :: gsw_kappa_t_exact
    end function gsw_kappa_t_exact
    
    elemental function gsw_latentheat_evap_ct (sa, ct) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct 
    real (r14) :: gsw_latentheat_evap_ct
    end function gsw_latentheat_evap_ct
    
    elemental function gsw_latentheat_evap_t (sa, t)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t  
    real (r14) :: gsw_latentheat_evap_t
    end function gsw_latentheat_evap_t
    
    elemental function gsw_latentheat_melting (sa, p)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p  
    real (r14) :: gsw_latentheat_melting
    end function gsw_latentheat_melting
    
    elemental function gsw_melting_ice_equilibrium_sa_ct_ratio (sa, p, &
                                                           saturation_fraction)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p, saturation_fraction
    real (r14) :: gsw_melting_ice_equilibrium_sa_ct_ratio
    end function gsw_melting_ice_equilibrium_sa_ct_ratio
    
    elemental subroutine gsw_melting_ice_into_seawater (sa, ct, p, &
                           saturation_fraction, w_ih, t_ih, sa_final, ct_final)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p, saturation_fraction, w_ih, t_ih
    real (r14), intent(out) :: sa_final, ct_final
    end subroutine gsw_melting_ice_into_seawater
    
    elemental function gsw_melting_ice_sa_ct_ratio (sa, ct, p, &
                                                    saturation_fraction, t_ih)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p, saturation_fraction, t_ih
    real (r14) :: gsw_melting_ice_sa_ct_ratio
    end function gsw_melting_ice_sa_ct_ratio
    
    elemental function gsw_melting_seaice_equilibrium_sa_ct_ratio (sa, p, &
                                                           saturation_fraction)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p, saturation_fraction
    real (r14) :: gsw_melting_seaice_equilibrium_sa_ct_ratio
    end function gsw_melting_seaice_equilibrium_sa_ct_ratio
    
    elemental subroutine gsw_melting_seaice_into_seawater (sa, ct, p, &
        saturation_fraction, w_seaice, sa_seaice, t_seaice, sa_final, ct_final)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p, saturation_fraction, w_seaice
    real (r14), intent(in) :: sa_seaice, t_seaice
    real (r14), intent(out) :: sa_final, ct_final
    end subroutine gsw_melting_seaice_into_seawater
    
    elemental function gsw_melting_seaice_sa_ct_ratio (sa, ct, p, &
                                      saturation_fraction, sa_seaice, t_seaice)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p, saturation_fraction, sa_seaice
    real (r14), intent(in) :: t_seaice
    real (r14) :: gsw_melting_seaice_sa_ct_ratio
    end function gsw_melting_seaice_sa_ct_ratio
    
    pure subroutine gsw_nsquared (sa, ct, p, lat, n2, p_mid)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa(:), ct(:), p(:), lat(:)
    real (r14), intent(out) :: n2(:), p_mid(:)
    end subroutine 
    
    elemental function gsw_pot_enthalpy_from_pt_ice (pt0_ice)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: pt0_ice
    real (r14) :: gsw_pot_enthalpy_from_pt_ice
    end function gsw_pot_enthalpy_from_pt_ice
    
    elemental function gsw_pot_enthalpy_from_pt_ice_poly (pt0_ice)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: pt0_ice
    real (r14) :: gsw_pot_enthalpy_from_pt_ice_poly
    end function gsw_pot_enthalpy_from_pt_ice_poly
    
    elemental function gsw_pot_rho_t_exact (sa, t, p, p_ref)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t, p, p_ref  
    real (r14) :: gsw_pot_rho_t_exact
    end function gsw_pot_rho_t_exact
    
    elemental function gsw_pressure_coefficient_ice (t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p
    real (r14) :: gsw_pressure_coefficient_ice
    end function gsw_pressure_coefficient_ice
    
    elemental function gsw_pressure_freezing_ct (sa, ct, saturation_fraction)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, saturation_fraction
    real (r14) :: gsw_pressure_freezing_ct
    end function gsw_pressure_freezing_ct
    
    elemental function gsw_pt0_cold_ice_poly (pot_enthalpy_ice)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: pot_enthalpy_ice
    real (r14) :: gsw_pt0_cold_ice_poly
    end function gsw_pt0_cold_ice_poly
    
    elemental function gsw_pt0_from_t (sa, t, p) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t, p 
    real (r14) :: gsw_pt0_from_t
    end function gsw_pt0_from_t
    
    elemental function gsw_pt0_from_t_ice (t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p
    real (r14) :: gsw_pt0_from_t_ice
    end function gsw_pt0_from_t_ice
    
    elemental subroutine gsw_pt_first_derivatives (sa, ct, pt_sa, pt_ct)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct
    real (r14), intent(out), optional :: pt_sa, pt_ct
    end subroutine gsw_pt_first_derivatives
    
    elemental function gsw_pt_from_ct (sa, ct) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct 
    real (r14) :: gsw_pt_from_ct
    end function gsw_pt_from_ct
    
    elemental function gsw_pt_from_entropy (sa, entropy)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, entropy
    real (r14) :: gsw_pt_from_entropy
    end function gsw_pt_from_entropy
    
    elemental function gsw_pt_from_pot_enthalpy_ice (pot_enthalpy_ice)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: pot_enthalpy_ice
    real (r14) :: gsw_pt_from_pot_enthalpy_ice
    end function gsw_pt_from_pot_enthalpy_ice
    
    elemental function gsw_pt_from_pot_enthalpy_ice_poly_dh (pot_enthalpy_ice)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: pot_enthalpy_ice
    real (r14) :: gsw_pt_from_pot_enthalpy_ice_poly_dh
    end function gsw_pt_from_pot_enthalpy_ice_poly_dh
    
    elemental function gsw_pt_from_pot_enthalpy_ice_poly (pot_enthalpy_ice)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: pot_enthalpy_ice
    real (r14) :: gsw_pt_from_pot_enthalpy_ice_poly
    end function gsw_pt_from_pot_enthalpy_ice_poly
    
    elemental function gsw_pt_from_t (sa, t, p, p_ref) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t, p, p_ref 
    real (r14) :: gsw_pt_from_t
    end function gsw_pt_from_t
    
    elemental function gsw_pt_from_t_ice (t, p, p_ref)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p, p_ref
    real (r14) :: gsw_pt_from_t_ice
    end function gsw_pt_from_t_ice
    
    elemental subroutine gsw_pt_second_derivatives (sa, ct, pt_sa_sa, &
                                                    pt_sa_ct, pt_ct_ct)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct
    real (r14), intent(out), optional :: pt_sa_sa, pt_sa_ct, pt_ct_ct
    end subroutine gsw_pt_second_derivatives
    
    elemental subroutine gsw_rho_alpha_beta (sa, ct, p, rho, alpha, beta)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p
    real (r14), intent(out), optional :: rho, alpha, beta
    end subroutine gsw_rho_alpha_beta
    
    elemental function gsw_rho (sa, ct, p) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p 
    real (r14) :: gsw_rho
    end function gsw_rho
    
    elemental subroutine gsw_rho_first_derivatives (sa, ct, p, drho_dsa, &
                                                    drho_dct, drho_dp)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p
    real (r14), intent(out) :: drho_dsa, drho_dct, drho_dp
    end subroutine gsw_rho_first_derivatives
    
    elemental function gsw_rho_ice (t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p
    real (r14) :: gsw_rho_ice
    end function gsw_rho_ice
    
    elemental function gsw_rho_t_exact (sa, t, p) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t, p 
    real (r14) :: gsw_rho_t_exact
    end function gsw_rho_t_exact
    
    elemental function gsw_saar (p, long, lat)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: p, long, lat
    real (r14) :: gsw_saar
    end function gsw_saar
    
    elemental function gsw_sa_from_rho (rho, ct, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: rho, ct, p
    real (r14) :: gsw_sa_from_rho
    end function gsw_sa_from_rho
    
    elemental function gsw_sa_from_sp_baltic (sp, long, lat)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sp, long, lat
    real (r14) :: gsw_sa_from_sp_baltic
    end function gsw_sa_from_sp_baltic
    
    elemental function gsw_sa_from_sp (sp, p, long, lat)       
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sp, p, long, lat       
    real (r14) :: gsw_sa_from_sp
    end function gsw_sa_from_sp
    
    elemental function gsw_sa_from_sstar (sstar, p, long, lat)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sstar, p, long, lat  
    real (r14) :: gsw_sa_from_sstar
    end function gsw_sa_from_sstar
    
    elemental function gsw_sa_p_inrange (sa, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p
    logical :: gsw_sa_p_inrange
    end function gsw_sa_p_inrange
    
    elemental subroutine gsw_seaice_fraction_to_freeze_seawater (sa, ct, p, &
      saturation_fraction, sa_seaice, t_seaice, sa_freeze, ct_freeze, w_seaice)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p, saturation_fraction, sa_seaice
    real (r14), intent(in) :: t_seaice
    real (r14), intent(out) :: sa_freeze, ct_freeze, w_seaice
    end subroutine gsw_seaice_fraction_to_freeze_seawater
    
    elemental function gsw_sigma0 (sa, ct) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct 
    real (r14) :: gsw_sigma0
    end function gsw_sigma0
    
    elemental function gsw_sigma1 (sa, ct) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct 
    real (r14) :: gsw_sigma1
    end function gsw_sigma1
    
    elemental function gsw_sigma2 (sa, ct) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct 
    real (r14) :: gsw_sigma2
    end function gsw_sigma2
    
    elemental function gsw_sigma3 (sa, ct) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct 
    real (r14) :: gsw_sigma3
    end function gsw_sigma3
    
    elemental function gsw_sigma4 (sa, ct) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct 
    real (r14) :: gsw_sigma4
    end function gsw_sigma4
    
    elemental function gsw_sound_speed (sa, ct, p) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p 
    real (r14) :: gsw_sound_speed
    end function gsw_sound_speed
    
    elemental function gsw_sound_speed_ice (t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p
    real (r14) :: gsw_sound_speed_ice
    end function gsw_sound_speed_ice
    
    elemental function gsw_sound_speed_t_exact (sa, t, p)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t, p  
    real (r14) :: gsw_sound_speed_t_exact
    end function gsw_sound_speed_t_exact
    
    elemental function gsw_specvol_anom (sa, ct, p)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p  
    real (r14) :: gsw_specvol_anom
    end function gsw_specvol_anom
    
    elemental function gsw_specvol (sa, ct, p) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p 
    real (r14) :: gsw_specvol
    end function gsw_specvol
    
    elemental function gsw_specvol_ice (t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: t, p
    real (r14) :: gsw_specvol_ice
    end function gsw_specvol_ice
    
    elemental function gsw_specvol_sso_0_p (p) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: p 
    real (r14) :: gsw_specvol_sso_0_p
    end function gsw_specvol_sso_0_p
    
    elemental function gsw_specvol_t_exact (sa, t, p)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t, p  
    real (r14) :: gsw_specvol_t_exact
    end function gsw_specvol_t_exact
    
    elemental function gsw_sp_from_c (c, t, p)       
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: c, t, p       
    real (r14) :: gsw_sp_from_c
    end function gsw_sp_from_c
    
    elemental function gsw_sp_from_sa_baltic (sa, long, lat)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, long, lat
    real (r14) :: gsw_sp_from_sa_baltic
    end function gsw_sp_from_sa_baltic
    
    elemental function gsw_sp_from_sa (sa, p, long, lat) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p, long, lat 
    real (r14) :: gsw_sp_from_sa
    end function gsw_sp_from_sa
    
    elemental function gsw_sp_from_sk (sk)       
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sk       
    real (r14) :: gsw_sp_from_sk
    end function gsw_sp_from_sk
    
    elemental function gsw_sp_from_sr (sr)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sr  
    real (r14) :: gsw_sp_from_sr
    end function gsw_sp_from_sr
    
    elemental function gsw_sp_from_sstar (sstar, p, long, lat)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sstar, p, long, lat  
    real (r14) :: gsw_sp_from_sstar
    end function gsw_sp_from_sstar
    
    elemental function gsw_sr_from_sp (sp) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sp 
    real (r14) :: gsw_sr_from_sp
    end function gsw_sr_from_sp
    
    elemental function gsw_sstar_from_sa (sa, p, long, lat) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p, long, lat 
    real (r14) :: gsw_sstar_from_sa
    end function gsw_sstar_from_sa
    
    elemental function gsw_sstar_from_sp (sp, p, long, lat) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sp, p, long, lat 
    real (r14) :: gsw_sstar_from_sp
    end function gsw_sstar_from_sp
    
    elemental function gsw_t_deriv_chem_potential_water_t_exact (sa, t, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, t, p
    real (r14) :: gsw_t_deriv_chem_potential_water_t_exact
    end function gsw_t_deriv_chem_potential_water_t_exact
    
    elemental function gsw_t_freezing_derivative_poly (sa, p, saturation_fraction)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p, saturation_fraction
    real (r14) :: gsw_t_freezing_derivative_poly
    end function gsw_t_freezing_derivative_poly
    
    elemental function gsw_t_freezing_exact (sa, p, saturation_fraction)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p, saturation_fraction
    real (r14) :: gsw_t_freezing_exact
    end function gsw_t_freezing_exact
    
    elemental function gsw_t_freezing (sa, p, saturation_fraction, exact)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p, saturation_fraction
    logical, intent(in), optional :: exact
    real (r14) :: gsw_t_freezing
    end function gsw_t_freezing
    
    elemental subroutine gsw_t_freezing_first_derivatives (sa, p, &
                                saturation_fraction, tfreezing_sa, tfreezing_p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p, saturation_fraction
    real (r14), intent(out), optional :: tfreezing_sa, tfreezing_p
    end subroutine gsw_t_freezing_first_derivatives
    
    elemental function gsw_t_freezing_poly (sa, p, saturation_fraction)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, p
    real (r14), intent(in), optional :: saturation_fraction
    real (r14) :: gsw_t_freezing_poly
    end function gsw_t_freezing_poly
    
    elemental function gsw_t_from_ct (sa, ct, p)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p  
    real (r14) :: gsw_t_from_ct
    end function gsw_t_from_ct
    
    elemental function gsw_t_from_pt0_ice (pt0_ice, p)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: pt0_ice, p
    real (r14) :: gsw_t_from_pt0_ice
    end function gsw_t_from_pt0_ice
    
    elemental function gsw_thermobaric (sa, ct, p)  
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa, ct, p  
    real (r14) :: gsw_thermobaric
    end function gsw_thermobaric
    
    pure subroutine gsw_turner_rsubrho (sa, ct, p, tu, rsubrho, p_mid)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: sa(:), ct(:), p(:)
    real (r14), intent(out) :: tu(:), rsubrho(:), p_mid(:)
    end subroutine 
    
    pure subroutine gsw_util_indx (x, n, z, k)
    integer, parameter :: r14 = selected_real_kind(14,30)
    integer, intent(in) :: n
    integer, intent(out) :: k
    real (r14), intent(in), dimension(n) :: x
    real (r14), intent(in) :: z
    end subroutine 
    
    pure function gsw_util_xinterp1 (x, y, n, x0)
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    integer, intent(in) :: n
    real (r14), intent(in) :: x0
    real (r14), dimension(n), intent(in) :: x, y
    real (r14) :: gsw_util_xinterp1
    end function gsw_util_xinterp1
    
    elemental function gsw_z_from_p (p, lat) 
    implicit none
    integer, parameter :: r14 = selected_real_kind(14,30)
    real (r14), intent(in) :: p, lat 
    real (r14) :: gsw_z_from_p
    end function gsw_z_from_p
    
end interface

end module gsw_mod_toolbox
