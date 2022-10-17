/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 *********************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <assert.h>
#include <mpi.h>

#include <e3sm_io.h>
#include <e3sm_io_err.h>
#include <e3sm_io_case.hpp>
#include <e3sm_io_driver.hpp>

/*----< get_gattrs() >-------------------------------------------------------*/
/* add global attributes */
static
int get_gattrs(e3sm_io_config &cfg,
               e3sm_io_decom  &decom,
               e3sm_io_driver &driver,
               case_meta      *cmeta,
               int             ncid)
{
    std::string prefix("");
    int err=0, nprocs;

    /* save number of processes as global attributes */
    if (cfg.strategy == blob) {
        if (cfg.api == adios) {
            GET_GATTR_INT("/__pio__/fillmode", 256)
            prefix = "pio_global/";
        }
        else {
            MPI_Comm_size(cfg.io_comm, &nprocs);
            GET_GATTR_INT("global_nprocs", nprocs)
            GET_GATTR_INT("num_decompositions", decom.num_decomp)
            GET_GATTR_INT("num_subfiles", cfg.num_subfiles)
        }
    }

    /* 747 global attributes: */
    GET_GATTR_TXT("initial_file", "/global/cfs/cdirs/e3sm/inputdata/atm/cam/inic/homme/cami_mam3_Linoz_0000-01-ne120np4_L72_c160318.nc")
    GET_GATTR_INT("ne", 120)
    GET_GATTR_INT("np", 9600)
    GET_GATTR_TXT("time_period_freq", "day_5")
    GET_GATTR_TXT("title", "MPAS-Ocean output file information")
    GET_GATTR_TXT("source", "MPAS Ocean")
    GET_GATTR_TXT("source_id", "a79fafdb76")
    GET_GATTR_TXT("product", "model-output")
    GET_GATTR_TXT("realm", "ocean")
    GET_GATTR_TXT("case", "GMPAS-NYF_T62_oRRS18to6v3")
    GET_GATTR_TXT("username", "E3SM")
    GET_GATTR_TXT("hostname", "cori-knl")
    GET_GATTR_TXT("git_version", "a79fafdb76")
    GET_GATTR_TXT("history", "created on 06/04/21 12:46:19")
    GET_GATTR_TXT("Conventions", "CF-1.7")
    GET_GATTR_TXT("institution_id", "E3SM-Project")
    GET_GATTR_TXT ("topography_file", "/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_ne120np4_16xdel2-PFC-consistentSGH.nc")
    GET_GATTR_TXT("institution", "LLNL (Lawrence Livermore National Laboratory, Livermore, CA 94550, USA); ANL (Argonne National Laboratory, Argonne, IL 60439, USA); BNL (Brookhaven National Laboratory, Upton, NY 11973, USA); LANL (Los Alamos National Laboratory, Los Alamos, NM 87545, USA); LBNL (Lawrence Berkeley National Laboratory, Berkeley, CA 94720, USA); ORNL (Oak Ridge National Laboratory, Oak Ridge, TN 37831, USA); PNNL (Pacific Northwest National Laboratory, Richland, WA 99352, USA); SNL (Sandia National Laboratories, Albuquerque,")
    GET_GATTR_TXT("contact", "e3sm-data-support@listserv.llnl.gov")
    GET_GATTR_TXT("on_a_sphere", "YES")
    GET_GATTR_DBL("sphere_radius", 6371229.)
    GET_GATTR_TXT("model_name", "mpas")
    GET_GATTR_TXT("core_name", "ocean")
    GET_GATTR_TXT("parent_id", "rz1tbn0bed\n555vk5hkh9\n0s1lcmezuy\n7z1uysqc5i\nwj661vqvze")
    GET_GATTR_TXT("mesh_spec", "0.0")
    GET_GATTR_TXT("config_ocean_run_mode", "forward")
    GET_GATTR_TXT("config_do_restart", "NO")
    GET_GATTR_TXT("config_restart_timestamp_name", "rpointer.ocn")
    GET_GATTR_TXT("config_start_time", "0001-01-01_00:00:00")
    GET_GATTR_TXT("config_stop_time", "none")
    GET_GATTR_TXT("config_run_duration", "0001-00-00_00:00:00")
    GET_GATTR_TXT("config_calendar_type", "gregorian_noleap")
    GET_GATTR_TXT("config_write_output_on_startup", "NO")
    GET_GATTR_INT("config_pio_num_iotasks", 0)
    GET_GATTR_INT("config_pio_stride", 1)
    GET_GATTR_INT("config_num_halos", 3)
    GET_GATTR_TXT("config_block_decomp_file_prefix", "/global/cfs/cdirs/e3sm/inputdata/ocn/mpas-o/oRRS18to6v3/mpas-o.graph.info.170111.part.")
    GET_GATTR_INT("config_number_of_blocks", 0)
    GET_GATTR_TXT("config_explicit_proc_decomp", "NO")
    GET_GATTR_TXT("config_proc_decomp_file_prefix", "graph.info.part.")
    GET_GATTR_TXT("config_dt", "00:06:00")
    GET_GATTR_TXT("config_time_integrator", "split_explicit")
    GET_GATTR_TXT("config_hmix_scaleWithMesh", "YES")
    GET_GATTR_DBL("config_maxMeshDensity", -1.)
    GET_GATTR_TXT("config_hmix_use_ref_cell_width", "NO")
    GET_GATTR_DBL("config_hmix_ref_cell_width", 30000.)
    GET_GATTR_DBL("config_apvm_scale_factor", 0.)
    GET_GATTR_TXT("config_use_mom_del2", "NO")
    GET_GATTR_DBL("config_mom_del2", 10.)
    GET_GATTR_TXT("config_use_tracer_del2", "NO")
    GET_GATTR_DBL("config_tracer_del2", 10.)
    GET_GATTR_TXT("config_use_mom_del4", "YES")
    GET_GATTR_DBL("config_mom_del4", 3200000000.)
    GET_GATTR_DBL("config_mom_del4_div_factor", 1.)
    GET_GATTR_TXT("config_use_tracer_del4", "NO")
    GET_GATTR_DBL("config_tracer_del4", 0.)
    GET_GATTR_TXT("config_use_Leith_del2", "NO")
    GET_GATTR_DBL("config_Leith_parameter", 1.)
    GET_GATTR_DBL("config_Leith_dx", 15000.)
    GET_GATTR_DBL("config_Leith_visc2_max", 2500.)
    GET_GATTR_TXT("config_eddying_resolution_taper", "ramp")
    GET_GATTR_DBL("config_eddying_resolution_ramp_min", 20000.)
    GET_GATTR_DBL("config_eddying_resolution_ramp_max", 30000.)
    GET_GATTR_TXT("config_use_Redi", "NO")
    GET_GATTR_TXT("config_Redi_closure", "constant")
    GET_GATTR_DBL("config_Redi_constant_kappa", 400.)
    GET_GATTR_DBL("config_Redi_maximum_slope", 0.01)
    GET_GATTR_TXT("config_Redi_use_slope_taper", "YES")
    GET_GATTR_TXT("config_Redi_use_surface_taper", "YES")
    GET_GATTR_TXT("config_Redi_N2_based_taper_enable", "YES")
    GET_GATTR_DBL("config_Redi_N2_based_taper_min", 0.1)
    GET_GATTR_TXT("config_Redi_N2_based_taper_limit_term1", "YES")
    GET_GATTR_TXT("config_use_GM", "NO")
    GET_GATTR_TXT("config_GM_closure", "EdenGreatbatch")
    GET_GATTR_DBL("config_GM_constant_kappa", 900.)
    GET_GATTR_DBL("config_GM_constant_gravWaveSpeed", 0.3)
    GET_GATTR_DBL("config_GM_spatially_variable_min_kappa", 300.)
    GET_GATTR_DBL("config_GM_spatially_variable_max_kappa", 1800.)
    GET_GATTR_DBL("config_GM_spatially_variable_baroclinic_mode", 3.)
    GET_GATTR_DBL("config_GM_Visbeck_alpha", 0.13)
    GET_GATTR_DBL("config_GM_Visbeck_max_depth", 1000.)
    GET_GATTR_DBL("config_GM_EG_riMin", 200.)
    GET_GATTR_DBL("config_GM_EG_kappa_factor", 3.)
    GET_GATTR_DBL("config_GM_EG_Rossby_factor", 2.)
    GET_GATTR_DBL("config_GM_EG_Rhines_factor", 0.3)
    GET_GATTR_TXT("config_Rayleigh_friction", "NO")
    GET_GATTR_DBL("config_Rayleigh_damping_coeff", 0.)
    GET_GATTR_TXT("config_Rayleigh_damping_depth_variable", "NO")
    GET_GATTR_TXT("config_Rayleigh_bottom_friction", "NO")
    GET_GATTR_DBL("config_Rayleigh_bottom_damping_coeff", 0.0001)
    GET_GATTR_TXT("config_use_cvmix", "YES")
    GET_GATTR_DBL("config_cvmix_prandtl_number", 1.)
    GET_GATTR_TXT("config_cvmix_background_scheme", "constant")
    GET_GATTR_DBL("config_cvmix_background_diffusion", 0.)
    GET_GATTR_DBL("config_cvmix_background_viscosity", 0.0001)
    GET_GATTR_DBL("config_cvmix_BryanLewis_bl1", 8.e-05)
    GET_GATTR_DBL("config_cvmix_BryanLewis_bl2", 0.000105)
    GET_GATTR_DBL("config_cvmix_BryanLewis_transitionDepth", 2500.)
    GET_GATTR_DBL("config_cvmix_BryanLewis_transitionWidth", 222.)
    GET_GATTR_TXT("config_use_cvmix_convection", "YES")
    GET_GATTR_DBL("config_cvmix_convective_diffusion", 1.)
    GET_GATTR_DBL("config_cvmix_convective_viscosity", 1.)
    GET_GATTR_TXT("config_cvmix_convective_basedOnBVF", "YES")
    GET_GATTR_DBL("config_cvmix_convective_triggerBVF", 0.)
    GET_GATTR_TXT("config_use_cvmix_shear", "YES")
    GET_GATTR_INT("config_cvmix_num_ri_smooth_loops", 2)
    GET_GATTR_TXT("config_cvmix_use_BLD_smoothing", "YES")
    GET_GATTR_TXT("config_cvmix_shear_mixing_scheme", "KPP")
    GET_GATTR_DBL("config_cvmix_shear_PP_nu_zero", 0.005)
    GET_GATTR_DBL("config_cvmix_shear_PP_alpha", 5.)
    GET_GATTR_DBL("config_cvmix_shear_PP_exp", 2.)
    GET_GATTR_DBL("config_cvmix_shear_KPP_nu_zero", 0.005)
    GET_GATTR_DBL("config_cvmix_shear_KPP_Ri_zero", 0.7)
    GET_GATTR_DBL("config_cvmix_shear_KPP_exp", 3.)
    GET_GATTR_TXT("config_use_cvmix_tidal_mixing", "NO")
    GET_GATTR_TXT("config_use_cvmix_double_diffusion", "NO")
    GET_GATTR_TXT("config_use_cvmix_kpp", "YES")
    GET_GATTR_TXT("config_use_cvmix_fixed_boundary_layer", "NO")
    GET_GATTR_DBL("config_cvmix_kpp_boundary_layer_depth", 30.)
    GET_GATTR_DBL("config_cvmix_kpp_criticalBulkRichardsonNumber", 0.25)
    GET_GATTR_TXT("config_cvmix_kpp_matching", "SimpleShapes")
    GET_GATTR_TXT("config_cvmix_kpp_EkmanOBL", "NO")
    GET_GATTR_TXT("config_cvmix_kpp_MonObOBL", "NO")
    GET_GATTR_TXT("config_cvmix_kpp_interpolationOMLType", "quadratic")
    GET_GATTR_DBL("config_cvmix_kpp_surface_layer_extent", 0.1)
    GET_GATTR_DBL("config_cvmix_kpp_surface_layer_averaging", 5.)
    GET_GATTR_DBL("configure_cvmix_kpp_minimum_OBL_under_sea_ice", 10.)
    GET_GATTR_DBL("config_cvmix_kpp_stop_OBL_search", 100.)
    GET_GATTR_TXT("config_cvmix_kpp_use_enhanced_diff", "YES")
    GET_GATTR_TXT("config_cvmix_kpp_nonlocal_with_implicit_mix", "NO")
    GET_GATTR_TXT("config_cvmix_kpp_use_theory_wave", "NO")
    GET_GATTR_TXT("config_cvmix_kpp_langmuir_mixing_opt", "NONE")
    GET_GATTR_TXT("config_cvmix_kpp_langmuir_entrainment_opt", "NONE")
    GET_GATTR_TXT("config_use_gotm", "NO")
    GET_GATTR_TXT("config_gotm_namelist_file", "gotmturb.nml")
    GET_GATTR_DBL("config_gotm_constant_surface_roughness_length", 0.02)
    GET_GATTR_DBL("config_gotm_constant_bottom_roughness_length", 0.0015)
    GET_GATTR_DBL("config_gotm_constant_bottom_drag_coeff", 0.001)
    GET_GATTR_TXT("config_use_variable_drag", "NO")
    GET_GATTR_TXT("config_use_bulk_wind_stress", "YES")
    GET_GATTR_TXT("config_use_bulk_thickness_flux", "YES")
    GET_GATTR_DBL("config_flux_attenuation_coefficient", 0.001)
    GET_GATTR_DBL("config_flux_attenuation_coefficient_runoff", 10.)
    GET_GATTR_TXT("config_use_time_varying_atmospheric_forcing", "NO")
    GET_GATTR_TXT("config_time_varying_atmospheric_forcing_type", "WINDPRES")
    GET_GATTR_TXT("config_time_varying_atmospheric_forcing_start_time", "0001-01-01_00:00:00")
    GET_GATTR_TXT("config_time_varying_atmospheric_forcing_reference_time", "0001-01-01_00:00:00")
    GET_GATTR_TXT("config_time_varying_atmospheric_forcing_cycle_start", "0001-01-01_00:00:00")
    GET_GATTR_TXT("config_time_varying_atmospheric_forcing_cycle_duration", "2-00-00_00:00:00")
    GET_GATTR_TXT("config_time_varying_atmospheric_forcing_interval", "01:00:00")
    GET_GATTR_DBL("config_time_varying_atmospheric_forcing_ramp", 10.)
    GET_GATTR_DBL("config_time_varying_atmospheric_forcing_ramp_delay", 0.)
    GET_GATTR_TXT("config_use_time_varying_land_ice_forcing", "NO")
    GET_GATTR_TXT("config_time_varying_land_ice_forcing_start_time", "0001-01-01_00:00:00")
    GET_GATTR_TXT("config_time_varying_land_ice_forcing_reference_time", "0001-01-01_00:00:00")
    GET_GATTR_TXT("config_time_varying_land_ice_forcing_cycle_start", "0001-01-01_00:00:00")
    GET_GATTR_TXT("config_time_varying_land_ice_forcing_cycle_duration", "2-00-00_00:00:00")
    GET_GATTR_TXT("config_time_varying_land_ice_forcing_interval", "01:00:00")
    GET_GATTR_DBL("config_ssh_grad_relax_timescale", 0.)
    GET_GATTR_TXT("config_remove_AIS_coupler_runoff", "NO")
    GET_GATTR_TXT("config_sw_absorption_type", "jerlov")
    GET_GATTR_INT("config_jerlov_water_type", 3)
    GET_GATTR_DBL("config_surface_buoyancy_depth", 1.)
    GET_GATTR_TXT("config_enable_shortwave_energy_fixer", "YES")
    GET_GATTR_TXT("config_use_tidal_forcing", "NO")
    GET_GATTR_DBL("config_use_tidal_forcing_tau", 10000.)
    GET_GATTR_TXT("config_tidal_forcing_type", "off")
    GET_GATTR_TXT("config_tidal_forcing_model", "off")
    GET_GATTR_DBL("config_tidal_forcing_monochromatic_amp", 2.)
    GET_GATTR_DBL("config_tidal_forcing_monochromatic_period", 0.5)
    GET_GATTR_DBL("config_tidal_forcing_monochromatic_phaseLag", 0.)
    GET_GATTR_DBL("config_tidal_forcing_monochromatic_baseline", 0.)
    GET_GATTR_TXT("config_use_tidal_potential_forcing", "NO")
    GET_GATTR_TXT("config_tidal_potential_reference_time", "0001-01-01_00:00:00")
    GET_GATTR_TXT("config_use_tidal_potential_forcing_M2", "YES")
    GET_GATTR_TXT("config_use_tidal_potential_forcing_S2", "YES")
    GET_GATTR_TXT("config_use_tidal_potential_forcing_N2", "YES")
    GET_GATTR_TXT("config_use_tidal_potential_forcing_K2", "YES")
    GET_GATTR_TXT("config_use_tidal_potential_forcing_K1", "YES")
    GET_GATTR_TXT("config_use_tidal_potential_forcing_O1", "YES")
    GET_GATTR_TXT("config_use_tidal_potential_forcing_Q1", "YES")
    GET_GATTR_TXT("config_use_tidal_potential_forcing_P1", "YES")
    GET_GATTR_DBL("config_tidal_potential_ramp", 10.)
    GET_GATTR_DBL("config_self_attraction_and_loading_beta", 0.09)
    GET_GATTR_TXT("config_use_vegetation_drag", "NO")
    GET_GATTR_TXT("config_use_vegetation_manning_equation", "NO")
    GET_GATTR_DBL("config_vegetation_drag_coefficient", 1.0900000333786)
    GET_GATTR_TXT("config_use_frazil_ice_formation", "YES")
    GET_GATTR_TXT("config_frazil_in_open_ocean", "YES")
    GET_GATTR_TXT("config_frazil_under_land_ice", "YES")
    GET_GATTR_DBL("config_frazil_heat_of_fusion", 333700.)
    GET_GATTR_DBL("config_frazil_ice_density", 1000.)
    GET_GATTR_DBL("config_frazil_fractional_thickness_limit", 0.1)
    GET_GATTR_DBL("config_specific_heat_sea_water", 3996.)
    GET_GATTR_DBL("config_frazil_maximum_depth", 100.)
    GET_GATTR_DBL("config_frazil_sea_ice_reference_salinity", 4.)
    GET_GATTR_DBL("config_frazil_land_ice_reference_salinity", 0.)
    GET_GATTR_DBL("config_frazil_maximum_freezing_temperature", 0.)
    GET_GATTR_TXT("config_frazil_use_surface_pressure", "NO")
    GET_GATTR_TXT("config_land_ice_flux_mode", "pressure_only")
    GET_GATTR_TXT("config_land_ice_flux_formulation", "Jenkins")
    GET_GATTR_TXT("config_land_ice_flux_useHollandJenkinsAdvDiff", "NO")
    GET_GATTR_DBL("config_land_ice_flux_attenuation_coefficient", 10.)
    GET_GATTR_DBL("config_land_ice_flux_boundaryLayerThickness", 10.)
    GET_GATTR_DBL("config_land_ice_flux_boundaryLayerNeighborWeight", 0.)
    GET_GATTR_DBL("config_land_ice_flux_cp_ice", 2009.)
    GET_GATTR_DBL("config_land_ice_flux_rho_ice", 918.)
    GET_GATTR_DBL("config_land_ice_flux_topDragCoeff", 0.0025)
    GET_GATTR_DBL("config_land_ice_flux_ISOMIP_gammaT", 0.0001)
    GET_GATTR_DBL("config_land_ice_flux_rms_tidal_velocity", 0.05)
    GET_GATTR_DBL("config_land_ice_flux_jenkins_heat_transfer_coefficient", 0.011)
    GET_GATTR_DBL("config_land_ice_flux_jenkins_salt_transfer_coefficient", 0.00031)
    GET_GATTR_TXT("config_vert_tracer_adv", "stencil")
    GET_GATTR_INT("config_vert_tracer_adv_order", 3)
    GET_GATTR_INT("config_horiz_tracer_adv_order", 3)
    GET_GATTR_DBL("config_coef_3rd_order", 0.25)
    GET_GATTR_TXT("config_monotonic", "YES")
    GET_GATTR_TXT("config_use_implicit_bottom_drag", "YES")
    GET_GATTR_DBL("config_implicit_bottom_drag_coeff", 0.001)
    GET_GATTR_TXT("config_use_implicit_bottom_roughness", "NO")
    GET_GATTR_TXT("config_use_implicit_bottom_drag_variable", "NO")
    GET_GATTR_TXT("config_use_implicit_bottom_drag_variable_mannings", "NO")
    GET_GATTR_TXT("config_use_explicit_bottom_drag", "NO")
    GET_GATTR_DBL("config_explicit_bottom_drag_coeff", 0.001)
    GET_GATTR_TXT("config_use_topographic_wave_drag", "NO")
    GET_GATTR_DBL("config_topographic_wave_drag_coeff", 0.0005)
    GET_GATTR_TXT("config_use_wetting_drying", "NO")
    GET_GATTR_TXT("config_prevent_drying", "NO")
    GET_GATTR_DBL("config_drying_min_cell_height", 0.00100000004749745)
    GET_GATTR_TXT("config_zero_drying_velocity", "NO")
    GET_GATTR_TXT("config_verify_not_dry", "NO")
    GET_GATTR_TXT("config_thickness_flux_type", "centered")
    GET_GATTR_DBL("config_drying_safety_height", 0.)
    GET_GATTR_DBL("config_density0", 1026.)
    GET_GATTR_TXT("config_pressure_gradient_type", "Jacobian_from_TS")
    GET_GATTR_DBL("config_common_level_weight", 0.5)
    GET_GATTR_DBL("config_zonal_ssh_grad", 0.)
    GET_GATTR_DBL("config_meridional_ssh_grad", 0.)
    GET_GATTR_TXT("config_eos_type", "jm")
    GET_GATTR_DBL("config_open_ocean_freezing_temperature_coeff_0", 0.)
    GET_GATTR_DBL("config_open_ocean_freezing_temperature_coeff_S", 0.)
    GET_GATTR_DBL("config_open_ocean_freezing_temperature_coeff_p", 0.)
    GET_GATTR_DBL("config_open_ocean_freezing_temperature_coeff_pS", 0.)
    GET_GATTR_DBL("config_open_ocean_freezing_temperature_coeff_mushy_az1_liq", -18.48)
    GET_GATTR_DBL("config_land_ice_cavity_freezing_temperature_coeff_0", 0.0622)
    GET_GATTR_DBL("config_land_ice_cavity_freezing_temperature_coeff_S", -0.0563)
    GET_GATTR_DBL("config_land_ice_cavity_freezing_temperature_coeff_p", -7.43e-08)
    GET_GATTR_DBL("config_land_ice_cavity_freezing_temperature_coeff_pS", -1.74e-10)
    GET_GATTR_DBL("config_eos_linear_alpha", 0.2)
    GET_GATTR_DBL("config_eos_linear_beta", 0.8)
    GET_GATTR_DBL("config_eos_linear_Tref", 5.)
    GET_GATTR_DBL("config_eos_linear_Sref", 35.)
    GET_GATTR_DBL("config_eos_linear_densityref", 1000.)
    GET_GATTR_DBL("config_eos_wright_ref_pressure", 0.)
    GET_GATTR_INT("config_n_ts_iter", 2)
    GET_GATTR_INT("config_n_bcl_iter_beg", 1)
    GET_GATTR_INT("config_n_bcl_iter_mid", 2)
    GET_GATTR_INT("config_n_bcl_iter_end", 2)
    GET_GATTR_TXT("config_btr_dt", "0000_00:00:12")
    GET_GATTR_INT("config_n_btr_cor_iter", 2)
    GET_GATTR_TXT("config_vel_correction", "YES")
    GET_GATTR_INT("config_btr_subcycle_loop_factor", 2)
    GET_GATTR_DBL("config_btr_gam1_velWt1", 0.5)
    GET_GATTR_DBL("config_btr_gam2_SSHWt1", 1.)
    GET_GATTR_DBL("config_btr_gam3_velWt2", 1.)
    GET_GATTR_TXT("config_btr_solve_SSH2", "NO")
    GET_GATTR_TXT("config_btr_si_preconditioner", "ras")
    GET_GATTR_DBL("config_btr_si_tolerance", 1.e-09)
    GET_GATTR_INT("config_n_btr_si_outer_iter", 2)
    GET_GATTR_TXT("config_btr_si_partition_match_mode", "NO")
    GET_GATTR_TXT("config_vert_coord_movement", "uniform_stretching")
    GET_GATTR_TXT("config_ALE_thickness_proportionality", "restingThickness_times_weights")
    GET_GATTR_DBL("config_vert_taper_weight_depth_1", 250.)
    GET_GATTR_DBL("config_vert_taper_weight_depth_2", 500.)
    GET_GATTR_TXT("config_use_min_max_thickness", "NO")
    GET_GATTR_DBL("config_min_thickness", 1.)
    GET_GATTR_DBL("config_max_thickness_factor", 6.)
    GET_GATTR_TXT("config_dzdk_positive", "NO")
    GET_GATTR_TXT("config_use_freq_filtered_thickness", "NO")
    GET_GATTR_DBL("config_thickness_filter_timescale", 5.)
    GET_GATTR_TXT("config_use_highFreqThick_restore", "NO")
    GET_GATTR_DBL("config_highFreqThick_restore_time", 30.)
    GET_GATTR_TXT("config_use_highFreqThick_del2", "NO")
    GET_GATTR_DBL("config_highFreqThick_del2", 100.)
    GET_GATTR_TXT("config_check_zlevel_consistency", "NO")
    GET_GATTR_TXT("config_check_ssh_consistency", "YES")
    GET_GATTR_TXT("config_filter_btr_mode", "NO")
    GET_GATTR_TXT("config_prescribe_velocity", "NO")
    GET_GATTR_TXT("config_prescribe_thickness", "NO")
    GET_GATTR_TXT("config_include_KE_vertex", "YES")
    GET_GATTR_TXT("config_check_tracer_monotonicity", "NO")
    GET_GATTR_TXT("config_compute_active_tracer_budgets", "YES")
    GET_GATTR_TXT("config_disable_thick_all_tend", "NO")
    GET_GATTR_TXT("config_disable_thick_hadv", "NO")
    GET_GATTR_TXT("config_disable_thick_vadv", "NO")
    GET_GATTR_TXT("config_disable_thick_sflux", "NO")
    GET_GATTR_TXT("config_disable_vel_all_tend", "NO")
    GET_GATTR_TXT("config_disable_vel_coriolis", "NO")
    GET_GATTR_TXT("config_disable_vel_pgrad", "NO")
    GET_GATTR_TXT("config_disable_vel_hmix", "NO")
    GET_GATTR_TXT("config_disable_vel_surface_stress", "NO")
    GET_GATTR_TXT("config_disable_vel_topographic_wave_drag", "NO")
    GET_GATTR_TXT("config_disable_vel_explicit_bottom_drag", "NO")
    GET_GATTR_TXT("config_disable_vel_vmix", "NO")
    GET_GATTR_TXT("config_disable_vel_vadv", "NO")
    GET_GATTR_TXT("config_disable_tr_all_tend", "NO")
    GET_GATTR_TXT("config_disable_tr_adv", "NO")
    GET_GATTR_TXT("config_disable_tr_hmix", "NO")
    GET_GATTR_TXT("config_disable_tr_vmix", "NO")
    GET_GATTR_TXT("config_disable_tr_sflux", "NO")
    GET_GATTR_TXT("config_disable_tr_nonlocalflux", "NO")
    GET_GATTR_TXT("config_disable_redi_k33", "NO")
    GET_GATTR_TXT("config_read_nearest_restart", "NO")
    GET_GATTR_TXT("config_conduct_tests", "NO")
    GET_GATTR_TXT("config_test_tensors", "NO")
    GET_GATTR_TXT("config_tensor_test_function", "sph_uCosCos")
    GET_GATTR_INT("config_vert_levels", -1)
    GET_GATTR_TXT("config_use_activeTracers", "YES")
    GET_GATTR_TXT("config_use_activeTracers_surface_bulk_forcing", "YES")
    GET_GATTR_TXT("config_use_activeTracers_surface_restoring", "YES")
    GET_GATTR_TXT("config_use_activeTracers_interior_restoring", "NO")
    GET_GATTR_TXT("config_use_activeTracers_exponential_decay", "NO")
    GET_GATTR_TXT("config_use_activeTracers_idealAge_forcing", "NO")
    GET_GATTR_TXT("config_use_activeTracers_ttd_forcing", "NO")
    GET_GATTR_TXT("config_use_surface_salinity_monthly_restoring", "YES")
    GET_GATTR_TXT("config_surface_salinity_monthly_restoring_compute_interval", "0000-00-01_00:00:00")
    GET_GATTR_DBL("config_salinity_restoring_constant_piston_velocity", 1.585e-06)
    GET_GATTR_DBL("config_salinity_restoring_max_difference", 100.)
    GET_GATTR_TXT("config_salinity_restoring_under_sea_ice", "NO")
    GET_GATTR_TXT("config_use_debugTracers", "NO")
    GET_GATTR_TXT("config_reset_debugTracers_near_surface", "NO")
    GET_GATTR_INT("config_reset_debugTracers_top_nLayers", 20)
    GET_GATTR_TXT("config_use_debugTracers_surface_bulk_forcing", "NO")
    GET_GATTR_TXT("config_use_debugTracers_surface_restoring", "NO")
    GET_GATTR_TXT("config_use_debugTracers_interior_restoring", "NO")
    GET_GATTR_TXT("config_use_debugTracers_exponential_decay", "NO")
    GET_GATTR_TXT("config_use_debugTracers_idealAge_forcing", "NO")
    GET_GATTR_TXT("config_use_debugTracers_ttd_forcing", "NO")
    GET_GATTR_TXT("config_use_ecosysTracers", "NO")
    GET_GATTR_TXT("config_ecosys_atm_co2_option", "constant")
    GET_GATTR_TXT("config_ecosys_atm_alt_co2_option", "constant")
    GET_GATTR_TXT("config_ecosys_atm_alt_co2_use_eco", "NO")
    GET_GATTR_DBL("config_ecosys_atm_co2_constant_value", 379.)
    GET_GATTR_TXT("config_use_ecosysTracers_surface_bulk_forcing", "NO")
    GET_GATTR_TXT("config_use_ecosysTracers_surface_restoring", "NO")
    GET_GATTR_TXT("config_use_ecosysTracers_interior_restoring", "NO")
    GET_GATTR_TXT("config_use_ecosysTracers_exponential_decay", "NO")
    GET_GATTR_TXT("config_use_ecosysTracers_idealAge_forcing", "NO")
    GET_GATTR_TXT("config_use_ecosysTracers_ttd_forcing", "NO")
    GET_GATTR_TXT("config_use_ecosysTracers_surface_value", "NO")
    GET_GATTR_TXT("config_use_ecosysTracers_sea_ice_coupling", "NO")
    GET_GATTR_TXT("config_ecosysTracers_diagnostic_fields_level1", "NO")
    GET_GATTR_TXT("config_ecosysTracers_diagnostic_fields_level2", "NO")
    GET_GATTR_TXT("config_ecosysTracers_diagnostic_fields_level3", "NO")
    GET_GATTR_TXT("config_ecosysTracers_diagnostic_fields_level4", "NO")
    GET_GATTR_TXT("config_ecosysTracers_diagnostic_fields_level5", "NO")
    GET_GATTR_TXT("config_use_DMSTracers", "NO")
    GET_GATTR_TXT("config_use_DMSTracers_surface_bulk_forcing", "NO")
    GET_GATTR_TXT("config_use_DMSTracers_surface_restoring", "NO")
    GET_GATTR_TXT("config_use_DMSTracers_interior_restoring", "NO")
    GET_GATTR_TXT("config_use_DMSTracers_exponential_decay", "NO")
    GET_GATTR_TXT("config_use_DMSTracers_idealAge_forcing", "NO")
    GET_GATTR_TXT("config_use_DMSTracers_ttd_forcing", "NO")
    GET_GATTR_TXT("config_use_DMSTracers_surface_value", "NO")
    GET_GATTR_TXT("config_use_DMSTracers_sea_ice_coupling", "NO")
    GET_GATTR_TXT("config_use_MacroMoleculesTracers", "NO")
    GET_GATTR_TXT("config_use_MacroMoleculesTracers_surface_bulk_forcing", "NO")
    GET_GATTR_TXT("config_use_MacroMoleculesTracers_surface_restoring", "NO")
    GET_GATTR_TXT("config_use_MacroMoleculesTracers_interior_restoring", "NO")
    GET_GATTR_TXT("config_use_MacroMoleculesTracers_exponential_decay", "NO")
    GET_GATTR_TXT("config_use_MacroMoleculesTracers_idealAge_forcing", "NO")
    GET_GATTR_TXT("config_use_MacroMoleculesTracers_ttd_forcing", "NO")
    GET_GATTR_TXT("config_use_MacroMoleculesTracers_surface_value", "NO")
    GET_GATTR_TXT("config_use_MacroMoleculesTracers_sea_ice_coupling", "NO")
    GET_GATTR_TXT("config_AM_globalStats_enable", "YES")
    GET_GATTR_TXT("config_AM_globalStats_compute_interval", "output_interval")
    GET_GATTR_TXT("config_AM_globalStats_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_globalStats_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_globalStats_text_file", "NO")
    GET_GATTR_TXT("config_AM_globalStats_directory", "analysis_members")
    GET_GATTR_TXT("config_AM_globalStats_output_stream", "globalStatsOutput")
    GET_GATTR_TXT("config_AM_surfaceAreaWeightedAverages_enable", "YES")
    GET_GATTR_TXT("config_AM_surfaceAreaWeightedAverages_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_surfaceAreaWeightedAverages_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_surfaceAreaWeightedAverages_compute_interval", "0000-00-00_01:00:00")
    GET_GATTR_TXT("config_AM_surfaceAreaWeightedAverages_output_stream", "surfaceAreaWeightedAveragesOutput")
    GET_GATTR_TXT("config_AM_waterMassCensus_enable", "NO")
    GET_GATTR_TXT("config_AM_waterMassCensus_compute_interval", "0000-00-00_01:00:00")
    GET_GATTR_TXT("config_AM_waterMassCensus_output_stream", "waterMassCensusOutput")
    GET_GATTR_TXT("config_AM_waterMassCensus_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_waterMassCensus_write_on_startup", "NO")
    GET_GATTR_DBL("config_AM_waterMassCensus_minTemperature", -2.)
    GET_GATTR_DBL("config_AM_waterMassCensus_maxTemperature", 30.)
    GET_GATTR_DBL("config_AM_waterMassCensus_minSalinity", 32.)
    GET_GATTR_DBL("config_AM_waterMassCensus_maxSalinity", 37.)
    GET_GATTR_TXT("config_AM_waterMassCensus_compute_predefined_regions", "YES")
    GET_GATTR_TXT("config_AM_waterMassCensus_region_group", "")
    GET_GATTR_TXT("config_AM_layerVolumeWeightedAverage_enable", "YES")
    GET_GATTR_TXT("config_AM_layerVolumeWeightedAverage_compute_interval", "0000-00-00_01:00:00")
    GET_GATTR_TXT("config_AM_layerVolumeWeightedAverage_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_layerVolumeWeightedAverage_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_layerVolumeWeightedAverage_output_stream", "layerVolumeWeightedAverageOutput")
    GET_GATTR_TXT("config_AM_zonalMean_enable", "NO")
    GET_GATTR_TXT("config_AM_zonalMean_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_zonalMean_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_zonalMean_compute_interval", "0000-00-00_01:00:00")
    GET_GATTR_TXT("config_AM_zonalMean_output_stream", "zonalMeanOutput")
    GET_GATTR_INT("config_AM_zonalMean_num_bins", 180)
    GET_GATTR_DBL("config_AM_zonalMean_min_bin", -1.e+34)
    GET_GATTR_DBL("config_AM_zonalMean_max_bin", -1.e+34)
    GET_GATTR_TXT("config_AM_okuboWeiss_enable", "NO")
    GET_GATTR_TXT("config_AM_okuboWeiss_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_okuboWeiss_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_okuboWeiss_compute_interval", "0000-00-00_01:00:00")
    GET_GATTR_TXT("config_AM_okuboWeiss_output_stream", "okuboWeissOutput")
    GET_GATTR_TXT("config_AM_okuboWeiss_directory", "analysis_members")
    GET_GATTR_DBL("config_AM_okuboWeiss_threshold_value", -0.2)
    GET_GATTR_DBL("config_AM_okuboWeiss_normalization", 1.e-10)
    GET_GATTR_DBL("config_AM_okuboWeiss_lambda2_normalization", 1.e-10)
    GET_GATTR_TXT("config_AM_okuboWeiss_use_lat_lon_coords", "YES")
    GET_GATTR_TXT("config_AM_okuboWeiss_compute_eddy_census", "YES")
    GET_GATTR_INT("config_AM_okuboWeiss_eddy_min_cells", 20)
    GET_GATTR_TXT("config_AM_meridionalHeatTransport_enable", "YES")
    GET_GATTR_TXT("config_AM_meridionalHeatTransport_compute_interval", "0000-00-00_01:00:00")
    GET_GATTR_TXT("config_AM_meridionalHeatTransport_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_meridionalHeatTransport_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_meridionalHeatTransport_output_stream", "meridionalHeatTransportOutput")
    GET_GATTR_INT("config_AM_meridionalHeatTransport_num_bins", 180)
    GET_GATTR_DBL("config_AM_meridionalHeatTransport_min_bin", -1.e+34)
    GET_GATTR_DBL("config_AM_meridionalHeatTransport_max_bin", -1.e+34)
    GET_GATTR_TXT("config_AM_meridionalHeatTransport_region_group", "")
    GET_GATTR_TXT("config_AM_testComputeInterval_enable", "NO")
    GET_GATTR_TXT("config_AM_testComputeInterval_compute_interval", "00-00-01_00:00:00")
    GET_GATTR_TXT("config_AM_testComputeInterval_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_testComputeInterval_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_testComputeInterval_output_stream", "testComputeIntervalOutput")
    GET_GATTR_TXT("config_AM_highFrequencyOutput_enable", "YES")
    GET_GATTR_TXT("config_AM_highFrequencyOutput_compute_interval", "output_interval")
    GET_GATTR_TXT("config_AM_highFrequencyOutput_output_stream", "highFrequencyOutput")
    GET_GATTR_TXT("config_AM_highFrequencyOutput_compute_on_startup", "NO")
    GET_GATTR_TXT("config_AM_highFrequencyOutput_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_timeFilters_enable", "NO")
    GET_GATTR_TXT("config_AM_timeFilters_compute_interval", "dt")
    GET_GATTR_TXT("config_AM_timeFilters_output_stream", "timeFiltersOutput")
    GET_GATTR_TXT("config_AM_timeFilters_restart_stream", "timeFiltersRestart")
    GET_GATTR_TXT("config_AM_timeFilters_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_timeFilters_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_timeFilters_initialize_filters", "YES")
    GET_GATTR_TXT("config_AM_timeFilters_tau", "90_00:00:00")
    GET_GATTR_TXT("config_AM_timeFilters_compute_cell_centered_values", "YES")
    GET_GATTR_TXT("config_AM_lagrPartTrack_enable", "NO")
    GET_GATTR_TXT("config_AM_lagrPartTrack_compute_interval", "dt")
    GET_GATTR_TXT("config_AM_lagrPartTrack_compute_on_startup", "NO")
    GET_GATTR_TXT("config_AM_lagrPartTrack_output_stream", "lagrPartTrackOutput")
    GET_GATTR_TXT("config_AM_lagrPartTrack_restart_stream", "lagrPartTrackRestart")
    GET_GATTR_TXT("config_AM_lagrPartTrack_input_stream", "lagrPartTrackInput")
    GET_GATTR_TXT("config_AM_lagrPartTrack_write_on_startup", "NO")
    GET_GATTR_INT("config_AM_lagrPartTrack_filter_number", 0)
    GET_GATTR_INT("config_AM_lagrPartTrack_timeIntegration", 2)
    GET_GATTR_TXT("config_AM_lagrPartTrack_reset_criteria", "none")
    GET_GATTR_TXT("config_AM_lagrPartTrack_reset_global_timestamp", "0000_00:00:00")
    GET_GATTR_TXT("config_AM_lagrPartTrack_region_stream", "lagrPartTrackRegions")
    GET_GATTR_TXT("config_AM_lagrPartTrack_reset_if_outside_region", "NO")
    GET_GATTR_TXT("config_AM_lagrPartTrack_reset_if_inside_region", "NO")
    GET_GATTR_TXT("config_AM_lagrPartTrack_sample_horizontal_interp", "YES")
    GET_GATTR_TXT("config_AM_lagrPartTrack_sample_temperature", "YES")
    GET_GATTR_TXT("config_AM_lagrPartTrack_sample_salinity", "YES")
    GET_GATTR_TXT("config_AM_lagrPartTrack_sample_DIC", "NO")
    GET_GATTR_TXT("config_AM_lagrPartTrack_sample_ALK", "NO")
    GET_GATTR_TXT("config_AM_lagrPartTrack_sample_PO4", "NO")
    GET_GATTR_TXT("config_AM_lagrPartTrack_sample_NO3", "NO")
    GET_GATTR_TXT("config_AM_lagrPartTrack_sample_SiO3", "NO")
    GET_GATTR_TXT("config_AM_lagrPartTrack_sample_NH4", "NO")
    GET_GATTR_TXT("config_AM_lagrPartTrack_sample_Fe", "NO")
    GET_GATTR_TXT("config_AM_lagrPartTrack_sample_O2", "NO")
    GET_GATTR_TXT("config_AM_eliassenPalm_enable", "NO")
    GET_GATTR_TXT("config_AM_eliassenPalm_compute_interval", "output_interval")
    GET_GATTR_TXT("config_AM_eliassenPalm_output_stream", "eliassenPalmOutput")
    GET_GATTR_TXT("config_AM_eliassenPalm_restart_stream", "eliassenPalmRestart")
    GET_GATTR_TXT("config_AM_eliassenPalm_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_eliassenPalm_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_eliassenPalm_debug", "NO")
    GET_GATTR_INT("config_AM_eliassenPalm_nBuoyancyLayers", 45)
    GET_GATTR_DBL("config_AM_eliassenPalm_rhomin_buoycoor", 900.)
    GET_GATTR_DBL("config_AM_eliassenPalm_rhomax_buoycoor", 1080.)
    GET_GATTR_TXT("config_AM_mixedLayerDepths_enable", "YES")
    GET_GATTR_TXT("config_AM_mixedLayerDepths_compute_interval", "dt")
    GET_GATTR_TXT("config_AM_mixedLayerDepths_output_stream", "mixedLayerDepthsOutput")
    GET_GATTR_TXT("config_AM_mixedLayerDepths_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_mixedLayerDepths_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_mixedLayerDepths_Tthreshold", "YES")
    GET_GATTR_TXT("config_AM_mixedLayerDepths_Dthreshold", "YES")
    GET_GATTR_DBL("config_AM_mixedLayerDepths_crit_temp_threshold", 0.2)
    GET_GATTR_DBL("config_AM_mixedLayerDepths_crit_dens_threshold", 0.03)
    GET_GATTR_DBL("config_AM_mixedLayerDepths_reference_pressure", 100000.)
    GET_GATTR_TXT("config_AM_mixedLayerDepths_Tgradient", "NO")
    GET_GATTR_TXT("config_AM_mixedLayerDepths_Dgradient", "NO")
    GET_GATTR_DBL("config_AM_mixedLayerDepths_temp_gradient_threshold", 5.e-07)
    GET_GATTR_DBL("config_AM_mixedLayerDepths_den_gradient_threshold", 5.e-08)
    GET_GATTR_INT("config_AM_mixedLayerDepths_interp_method", 1)
    GET_GATTR_TXT("config_AM_regionalStatsDaily_enable", "NO")
    GET_GATTR_TXT("config_AM_regionalStatsDaily_compute_on_startup", "NO")
    GET_GATTR_TXT("config_AM_regionalStatsDaily_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_regionalStatsDaily_compute_interval", "output_interval")
    GET_GATTR_TXT("config_AM_regionalStatsDaily_output_stream", "regionalStatsDailyOutput")
    GET_GATTR_TXT("config_AM_regionalStatsDaily_restart_stream", "regionalMasksInput")
    GET_GATTR_TXT("config_AM_regionalStatsDaily_input_stream", "regionalMasksInput")
    GET_GATTR_TXT("config_AM_regionalStatsDaily_operation", "avg")
    GET_GATTR_TXT("config_AM_regionalStatsDaily_region_type", "cell")
    GET_GATTR_TXT("config_AM_regionalStatsDaily_region_group", "all")
    GET_GATTR_TXT("config_AM_regionalStatsDaily_1d_weighting_function", "mul")
    GET_GATTR_TXT("config_AM_regionalStatsDaily_2d_weighting_function", "mul")
    GET_GATTR_TXT("config_AM_regionalStatsDaily_1d_weighting_field", "areaCell")
    GET_GATTR_TXT("config_AM_regionalStatsDaily_2d_weighting_field", "volumeCell")
    GET_GATTR_TXT("config_AM_regionalStatsDaily_vertical_mask", "cellMask")
    GET_GATTR_TXT("config_AM_regionalStatsDaily_vertical_dimension", "nVertLevels")
    GET_GATTR_TXT("config_AM_regionalStatsWeekly_enable", "NO")
    GET_GATTR_TXT("config_AM_regionalStatsWeekly_compute_on_startup", "NO")
    GET_GATTR_TXT("config_AM_regionalStatsWeekly_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_regionalStatsWeekly_compute_interval", "output_interval")
    GET_GATTR_TXT("config_AM_regionalStatsWeekly_output_stream", "regionalStatsWeeklyOutput")
    GET_GATTR_TXT("config_AM_regionalStatsWeekly_restart_stream", "regionalMasksInput")
    GET_GATTR_TXT("config_AM_regionalStatsWeekly_input_stream", "regionalMasksInput")
    GET_GATTR_TXT("config_AM_regionalStatsWeekly_operation", "avg")
    GET_GATTR_TXT("config_AM_regionalStatsWeekly_region_type", "cell")
    GET_GATTR_TXT("config_AM_regionalStatsWeekly_region_group", "all")
    GET_GATTR_TXT("config_AM_regionalStatsWeekly_1d_weighting_function", "mul")
    GET_GATTR_TXT("config_AM_regionalStatsWeekly_2d_weighting_function", "mul")
    GET_GATTR_TXT("config_AM_regionalStatsWeekly_1d_weighting_field", "areaCell")
    GET_GATTR_TXT("config_AM_regionalStatsWeekly_2d_weighting_field", "volumeCell")
    GET_GATTR_TXT("config_AM_regionalStatsWeekly_vertical_mask", "cellMask")
    GET_GATTR_TXT("config_AM_regionalStatsWeekly_vertical_dimension", "nVertLevels")
    GET_GATTR_TXT("config_AM_regionalStatsMonthly_enable", "NO")
    GET_GATTR_TXT("config_AM_regionalStatsMonthly_compute_on_startup", "NO")
    GET_GATTR_TXT("config_AM_regionalStatsMonthly_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_regionalStatsMonthly_compute_interval", "output_interval")
    GET_GATTR_TXT("config_AM_regionalStatsMonthly_output_stream", "regionalStatsMonthlyOutput")
    GET_GATTR_TXT("config_AM_regionalStatsMonthly_restart_stream", "regionalMasksInput")
    GET_GATTR_TXT("config_AM_regionalStatsMonthly_input_stream", "regionalMasksInput")
    GET_GATTR_TXT("config_AM_regionalStatsMonthly_operation", "avg")
    GET_GATTR_TXT("config_AM_regionalStatsMonthly_region_type", "cell")
    GET_GATTR_TXT("config_AM_regionalStatsMonthly_region_group", "all")
    GET_GATTR_TXT("config_AM_regionalStatsMonthly_1d_weighting_function", "mul")
    GET_GATTR_TXT("config_AM_regionalStatsMonthly_2d_weighting_function", "mul")
    GET_GATTR_TXT("config_AM_regionalStatsMonthly_1d_weighting_field", "areaCell")
    GET_GATTR_TXT("config_AM_regionalStatsMonthly_2d_weighting_field", "volumeCell")
    GET_GATTR_TXT("config_AM_regionalStatsMonthly_vertical_mask", "cellMask")
    GET_GATTR_TXT("config_AM_regionalStatsMonthly_vertical_dimension", "nVertLevels")
    GET_GATTR_TXT("config_AM_regionalStatsCustom_enable", "NO")
    GET_GATTR_TXT("config_AM_regionalStatsCustom_compute_on_startup", "NO")
    GET_GATTR_TXT("config_AM_regionalStatsCustom_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_regionalStatsCustom_compute_interval", "output_interval")
    GET_GATTR_TXT("config_AM_regionalStatsCustom_output_stream", "regionalStatsCustomOutput")
    GET_GATTR_TXT("config_AM_regionalStatsCustom_restart_stream", "regionalMasksInput")
    GET_GATTR_TXT("config_AM_regionalStatsCustom_input_stream", "regionalMasksInput")
    GET_GATTR_TXT("config_AM_regionalStatsCustom_operation", "avg")
    GET_GATTR_TXT("config_AM_regionalStatsCustom_region_type", "cell")
    GET_GATTR_TXT("config_AM_regionalStatsCustom_region_group", "all")
    GET_GATTR_TXT("config_AM_regionalStatsCustom_1d_weighting_function", "mul")
    GET_GATTR_TXT("config_AM_regionalStatsCustom_2d_weighting_function", "mul")
    GET_GATTR_TXT("config_AM_regionalStatsCustom_1d_weighting_field", "areaCell")
    GET_GATTR_TXT("config_AM_regionalStatsCustom_2d_weighting_field", "volumeCell")
    GET_GATTR_TXT("config_AM_regionalStatsCustom_vertical_mask", "cellMask")
    GET_GATTR_TXT("config_AM_regionalStatsCustom_vertical_dimension", "nVertLevels")
    GET_GATTR_TXT("config_AM_timeSeriesStatsDaily_enable", "NO")
    GET_GATTR_TXT("config_AM_timeSeriesStatsDaily_compute_on_startup", "NO")
    GET_GATTR_TXT("config_AM_timeSeriesStatsDaily_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_timeSeriesStatsDaily_compute_interval", "00-00-00_01:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsDaily_output_stream", "timeSeriesStatsDailyOutput")
    GET_GATTR_TXT("config_AM_timeSeriesStatsDaily_restart_stream", "none")
    GET_GATTR_TXT("config_AM_timeSeriesStatsDaily_operation", "avg")
    GET_GATTR_TXT("config_AM_timeSeriesStatsDaily_reference_times", "initial_time")
    GET_GATTR_TXT("config_AM_timeSeriesStatsDaily_duration_intervals", "repeat_interval")
    GET_GATTR_TXT("config_AM_timeSeriesStatsDaily_repeat_intervals", "reset_interval")
    GET_GATTR_TXT("config_AM_timeSeriesStatsDaily_reset_intervals", "00-00-01_00:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsDaily_backward_output_offset", "00-00-01_00:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthly_enable", "YES")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthly_compute_on_startup", "NO")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthly_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthly_compute_interval", "dt")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthly_output_stream", "timeSeriesStatsMonthlyOutput")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthly_restart_stream", "none")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthly_operation", "avg")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthly_reference_times", "initial_time")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthly_duration_intervals", "repeat_interval")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthly_repeat_intervals", "reset_interval")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthly_reset_intervals", "00-01-00_00:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthly_backward_output_offset", "00-01-00_00:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsClimatology_enable", "NO")
    GET_GATTR_TXT("config_AM_timeSeriesStatsClimatology_compute_on_startup", "NO")
    GET_GATTR_TXT("config_AM_timeSeriesStatsClimatology_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_timeSeriesStatsClimatology_compute_interval", "00-00-00_01:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsClimatology_output_stream", "timeSeriesStatsClimatologyOutput")
    GET_GATTR_TXT("config_AM_timeSeriesStatsClimatology_restart_stream", "timeSeriesStatsClimatologyRestart")
    GET_GATTR_TXT("config_AM_timeSeriesStatsClimatology_operation", "avg")
    GET_GATTR_TXT("config_AM_timeSeriesStatsClimatology_reference_times", "00-03-01_00:00:00;00-06-01_00:00:00;00-09-01_00:00:00;00-12-01_00:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsClimatology_duration_intervals", "00-03-00_00:00:00;00-03-00_00:00:00;00-03-00_00:00:00;00-03-00_00:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsClimatology_repeat_intervals", "01-00-00_00:00:00;01-00-00_00:00:00;01-00-00_00:00:00;01-00-00_00:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsClimatology_reset_intervals", "1000-00-00_00:00:00;1000-00-00_00:00:00;1000-00-00_00:00:00;1000-00-00_00:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsClimatology_backward_output_offset", "00-03-00_00:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMax_enable", "YES")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMax_compute_on_startup", "NO")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMax_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMax_compute_interval", "00-00-00_01:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMax_output_stream", "timeSeriesStatsMonthlyMaxOutput")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMax_restart_stream", "none")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMax_operation", "max")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMax_reference_times", "initial_time")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMax_duration_intervals", "repeat_interval")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMax_repeat_intervals", "reset_interval")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMax_reset_intervals", "00-01-00_00:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMax_backward_output_offset", "00-01-00_00:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMin_enable", "YES")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMin_compute_on_startup", "NO")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMin_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMin_compute_interval", "00-00-00_01:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMin_output_stream", "timeSeriesStatsMonthlyMinOutput")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMin_restart_stream", "none")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMin_operation", "min")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMin_reference_times", "initial_time")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMin_duration_intervals", "repeat_interval")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMin_repeat_intervals", "reset_interval")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMin_reset_intervals", "00-01-00_00:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsMonthlyMin_backward_output_offset", "00-01-00_00:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsCustom_enable", "NO")
    GET_GATTR_TXT("config_AM_timeSeriesStatsCustom_compute_on_startup", "NO")
    GET_GATTR_TXT("config_AM_timeSeriesStatsCustom_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_timeSeriesStatsCustom_compute_interval", "00-00-00_01:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsCustom_output_stream", "timeSeriesStatsCustomOutput")
    GET_GATTR_TXT("config_AM_timeSeriesStatsCustom_restart_stream", "timeSeriesStatsCustomRestart")
    GET_GATTR_TXT("config_AM_timeSeriesStatsCustom_operation", "avg")
    GET_GATTR_TXT("config_AM_timeSeriesStatsCustom_reference_times", "initial_time")
    GET_GATTR_TXT("config_AM_timeSeriesStatsCustom_duration_intervals", "repeat_interval")
    GET_GATTR_TXT("config_AM_timeSeriesStatsCustom_repeat_intervals", "reset_interval")
    GET_GATTR_TXT("config_AM_timeSeriesStatsCustom_reset_intervals", "00-00-07_00:00:00")
    GET_GATTR_TXT("config_AM_timeSeriesStatsCustom_backward_output_offset", "00-00-01_00:00:00")
    GET_GATTR_TXT("config_AM_pointwiseStats_enable", "NO")
    GET_GATTR_TXT("config_AM_pointwiseStats_compute_interval", "output_interval")
    GET_GATTR_TXT("config_AM_pointwiseStats_output_stream", "pointwiseStatsOutput")
    GET_GATTR_TXT("config_AM_pointwiseStats_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_pointwiseStats_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_debugDiagnostics_enable", "NO")
    GET_GATTR_TXT("config_AM_debugDiagnostics_compute_interval", "dt")
    GET_GATTR_TXT("config_AM_debugDiagnostics_output_stream", "debugDiagnosticsOutput")
    GET_GATTR_TXT("config_AM_debugDiagnostics_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_debugDiagnostics_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_debugDiagnostics_check_state", "YES")
    GET_GATTR_TXT("config_AM_rpnCalculator_enable", "NO")
    GET_GATTR_TXT("config_AM_rpnCalculator_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_rpnCalculator_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_rpnCalculator_compute_interval", "0010-00-00_00:00:00")
    GET_GATTR_TXT("config_AM_rpnCalculator_output_stream", "none")
    GET_GATTR_TXT("config_AM_rpnCalculator_variable_a", "layerThickness")
    GET_GATTR_TXT("config_AM_rpnCalculator_variable_b", "areaCell")
    GET_GATTR_TXT("config_AM_rpnCalculator_variable_c", "none")
    GET_GATTR_TXT("config_AM_rpnCalculator_variable_d", "none")
    GET_GATTR_TXT("config_AM_rpnCalculator_variable_e", "none")
    GET_GATTR_TXT("config_AM_rpnCalculator_variable_f", "none")
    GET_GATTR_TXT("config_AM_rpnCalculator_variable_g", "none")
    GET_GATTR_TXT("config_AM_rpnCalculator_variable_h", "none")
    GET_GATTR_TXT("config_AM_rpnCalculator_expression_1", "a b *")
    GET_GATTR_TXT("config_AM_rpnCalculator_expression_2", "none")
    GET_GATTR_TXT("config_AM_rpnCalculator_expression_3", "none")
    GET_GATTR_TXT("config_AM_rpnCalculator_expression_4", "none")
    GET_GATTR_TXT("config_AM_rpnCalculator_output_name_1", "volumeCell")
    GET_GATTR_TXT("config_AM_rpnCalculator_output_name_2", "none")
    GET_GATTR_TXT("config_AM_rpnCalculator_output_name_3", "none")
    GET_GATTR_TXT("config_AM_rpnCalculator_output_name_4", "none")
    GET_GATTR_TXT("config_AM_transectTransport_enable", "NO")
    GET_GATTR_TXT("config_AM_transectTransport_compute_interval", "output_interval")
    GET_GATTR_TXT("config_AM_transectTransport_output_stream", "transectTransportOutput")
    GET_GATTR_TXT("config_AM_transectTransport_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_transectTransport_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_transectTransport_transect_group", "all")
    GET_GATTR_TXT("config_AM_eddyProductVariables_enable", "YES")
    GET_GATTR_TXT("config_AM_eddyProductVariables_compute_interval", "0000-00-00_01:00:00")
    GET_GATTR_TXT("config_AM_eddyProductVariables_output_stream", "eddyProductVariablesOutput")
    GET_GATTR_TXT("config_AM_eddyProductVariables_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_eddyProductVariables_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_mocStreamfunction_enable", "YES")
    GET_GATTR_TXT("config_AM_mocStreamfunction_compute_interval", "0000-00-00_01:00:00")
    GET_GATTR_TXT("config_AM_mocStreamfunction_output_stream", "mocStreamfunctionOutput")
    GET_GATTR_TXT("config_AM_mocStreamfunction_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_mocStreamfunction_write_on_startup", "NO")
    GET_GATTR_DBL("config_AM_mocStreamfunction_min_bin", -1.e+34)
    GET_GATTR_DBL("config_AM_mocStreamfunction_max_bin", -1.e+34)
    GET_GATTR_INT("config_AM_mocStreamfunction_num_bins", 180)
    GET_GATTR_TXT("config_AM_mocStreamfunction_region_group", "all")
    GET_GATTR_TXT("config_AM_mocStreamfunction_transect_group", "all")
    GET_GATTR_TXT("config_AM_oceanHeatContent_enable", "YES")
    GET_GATTR_TXT("config_AM_oceanHeatContent_compute_interval", "0000-00-00_01:00:00")
    GET_GATTR_TXT("config_AM_oceanHeatContent_output_stream", "oceanHeatContentOutput")
    GET_GATTR_TXT("config_AM_oceanHeatContent_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_oceanHeatContent_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_mixedLayerHeatBudget_enable", "YES")
    GET_GATTR_TXT("config_AM_mixedLayerHeatBudget_compute_interval", "0000-00-00_01:00:00")
    GET_GATTR_TXT("config_AM_mixedLayerHeatBudget_output_stream", "mixedLayerHeatBudgetOutput")
    GET_GATTR_TXT("config_AM_mixedLayerHeatBudget_compute_on_startup", "NO")
    GET_GATTR_TXT("config_AM_mixedLayerHeatBudget_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_sedimentFluxIndex_enable", "NO")
    GET_GATTR_TXT("config_AM_sedimentFluxIndex_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_sedimentFluxIndex_write_on_startup", "YES")
    GET_GATTR_TXT("config_AM_sedimentFluxIndex_compute_interval", "output_interval")
    GET_GATTR_TXT("config_AM_sedimentFluxIndex_output_stream", "sedimentFluxIndexOutput")
    GET_GATTR_TXT("config_AM_sedimentFluxIndex_directory", "analysis_members")
    GET_GATTR_TXT("config_AM_sedimentFluxIndex_use_lat_lon_coords", "YES")
    GET_GATTR_TXT("config_AM_sedimentTransport_enable", "NO")
    GET_GATTR_TXT("config_AM_sedimentTransport_compute_on_startup", "YES")
    GET_GATTR_TXT("config_AM_sedimentTransport_write_on_startup", "YES")
    GET_GATTR_TXT("config_AM_sedimentTransport_compute_interval", "output_interval")
    GET_GATTR_TXT("config_AM_sedimentTransport_output_stream", "sedimentTransportOutput")
    GET_GATTR_TXT("config_AM_sedimentTransport_directory", "analysis_members")
    GET_GATTR_DBL("config_AM_sedimentTransport_grain_size", 0.00025)
    GET_GATTR_TXT("config_AM_sedimentTransport_ws_formula", "VanRijn1993")
    GET_GATTR_TXT("config_AM_sedimentTransport_bedld_formula", "Soulsby-Damgaard")
    GET_GATTR_TXT("config_AM_sedimentTransport_SSC_ref_formula", "Lee2004")
    GET_GATTR_DBL("config_AM_sedimentTransport_drag_coefficient", 0.0025)
    GET_GATTR_DBL("config_AM_sedimentTransport_erate", 0.0005)
    GET_GATTR_DBL("config_AM_sedimentTransport_tau_ce", 0.1)
    GET_GATTR_DBL("config_AM_sedimentTransport_tau_cd", 0.1)
    GET_GATTR_DBL("config_AM_sedimentTransport_Manning_coef", 0.022)
    GET_GATTR_DBL("config_AM_sedimentTransport_grain_porosity", 0.5)
    GET_GATTR_DBL("config_AM_sedimentTransport_water_density", 1020.)
    GET_GATTR_DBL("config_AM_sedimentTransport_grain_density", 2650.)
    GET_GATTR_DBL("config_AM_sedimentTransport_alpha", 0.01)
    GET_GATTR_DBL("config_AM_sedimentTransport_kinematic_viscosity", 1.e-06)
    GET_GATTR_DBL("config_AM_sedimentTransport_vertical_diffusion_coefficient", 0.01)
    GET_GATTR_TXT("config_AM_sedimentTransport_bedload", "YES")
    GET_GATTR_TXT("config_AM_sedimentTransport_suspended", "YES")
    GET_GATTR_TXT("config_AM_sedimentTransport_use_lat_lon_coords", "YES")
    GET_GATTR_TXT("config_AM_harmonicAnalysis_enable", "NO")
    GET_GATTR_TXT("config_AM_harmonicAnalysis_compute_interval", "output_interval")
    GET_GATTR_TXT("config_AM_harmonicAnalysis_start", "0001-01-01_00:00:00")
    GET_GATTR_TXT("config_AM_harmonicAnalysis_end", "0001-01-01_00:00:00")
    GET_GATTR_TXT("config_AM_harmonicAnalysis_output_stream", "harmonicAnalysisOutput")
    GET_GATTR_TXT("config_AM_harmonicAnalysis_restart_stream", "harmonicAnalysisRestart")
    GET_GATTR_TXT("config_AM_harmonicAnalysis_compute_on_startup", "NO")
    GET_GATTR_TXT("config_AM_harmonicAnalysis_write_on_startup", "NO")
    GET_GATTR_TXT("config_AM_harmonicAnalysis_use_M2", "YES")
    GET_GATTR_TXT("config_AM_harmonicAnalysis_use_S2", "YES")
    GET_GATTR_TXT("config_AM_harmonicAnalysis_use_N2", "YES")
    GET_GATTR_TXT("config_AM_harmonicAnalysis_use_K2", "YES")
    GET_GATTR_TXT("config_AM_harmonicAnalysis_use_K1", "YES")
    GET_GATTR_TXT("config_AM_harmonicAnalysis_use_O1", "YES")
    GET_GATTR_TXT("config_AM_harmonicAnalysis_use_Q1", "YES")
    GET_GATTR_TXT("config_AM_harmonicAnalysis_use_P1", "YES")
    GET_GATTR_TXT("file_id", "6270j0pl9q")

err_out:
    return err;
}

/*----< inq_G_case() >-------------------------------------------------------*/
int e3sm_io_case::inq_G_case(e3sm_io_config   &cfg,
                             e3sm_io_decom    &decom,
                             e3sm_io_driver   &driver,
                             int               ncid)    /* file ID */
{
    /* Total 52 variables */
    int i, err, nprocs, dimids[4], nvars_decomp;
    int dim_time, dim_nCells, dim_nEdges, dim_nVertices, dim_StrLen;
    int dim_nVertLevels, dim_nVertLevelsP1, dim_nblobs;
    int dim_nelems[MAX_NUM_DECOMP], dim_max_nreqs[MAX_NUM_DECOMP];
    int g_dimids[MAX_NUM_DECOMP][4];
    std::map<int, std::string> dnames;
    var_meta *varp;
    case_meta *cmeta;

    if (cfg.run_case == F) {
        if (cfg.hist == h0) cmeta = &cfg.F_case_h0;
        else                cmeta = &cfg.F_case_h1;
    } else if (cfg.run_case == I) {
        if (cfg.hist == h0) cmeta = &cfg.I_case_h0;
        else                cmeta = &cfg.I_case_h1;
    } else if (cfg.run_case == G)
        cmeta = &cfg.G_case;


    /* add global attributes */
    err = get_gattrs(cfg, decom, driver, cmeta, ncid);
    CHECK_ERR

    /* define dimensions */

    /* GMPAS-NYF_T62_oRRS18to6v3.mpaso.hist.0001-01-01_00000
      nCells        =  3693225 ;
      nEdges        = 11135652 ;
      nVertices     =  7441216 ;
      nVertLevelsP1 =       81 ;
      nVertLevels   =       80 ;
      StrLen        =       64 ;

      num_decomp = 6 ;
      decomp_nprocs = 9600 ;
      D1.total_nreqs = 3692863 ;     3693225         nCells
      D2.total_nreqs = 9804239 ;    11135652         nEdges
      D3.total_nreqs = 3693225 ;     3693225 x 80    nCells    x nVertLevels
      D4.total_nreqs = 11135652 ;   11135652 x 80    nEdges    x nVertLevels
      D5.total_nreqs = 7441216 ;     7441216 x 80    nVertices x nVertLevels
      D6.total_nreqs = 3693225 ;     3693225 x 81    nCells    x nVertLevelsP1
    */

    INQ_DIM("Time",          NC_UNLIMITED,     &dim_time)
    INQ_DIM("nCells",        decom.dims[0][0], &dim_nCells)
    INQ_DIM("nEdges",        decom.dims[1][0], &dim_nEdges)
    INQ_DIM("nVertices",     decom.dims[4][0], &dim_nVertices)
    INQ_DIM("nVertLevelsP1", decom.dims[5][1], &dim_nVertLevelsP1)
    INQ_DIM("nVertLevels",   decom.dims[2][1], &dim_nVertLevels)
    INQ_DIM("StrLen",        64    ,           &dim_StrLen)

    if (cfg.strategy == blob && cfg.api != adios) {
        char name[64];

        g_dimids[0][0] = dim_nCells;
        g_dimids[1][0] = dim_nEdges;
        g_dimids[2][0] = dim_nCells;    g_dimids[2][1] = dim_nVertLevels;
        g_dimids[3][0] = dim_nEdges;    g_dimids[3][1] = dim_nVertLevels;
        g_dimids[4][0] = dim_nVertices; g_dimids[4][1] = dim_nVertLevels;
        g_dimids[5][0] = dim_nCells;    g_dimids[5][1] = dim_nVertLevelsP1;

        /* additional dimensions to be used by decomposition variables */
        MPI_Comm_size(cfg.sub_comm, &nprocs);
        INQ_DIM("nblobs", (MPI_Offset)nprocs, &dim_nblobs)

        for (i=0; i<decom.num_decomp; i++) {
            sprintf(name, "D%d.nelems", i+1);
            INQ_DIM(name, decom.nelems[i], &dim_nelems[i])
        }

        for (i=0; i<decom.num_decomp; i++) {
            sprintf(name, "D%d.max_nreqs", i+1);
            INQ_DIM(name, (MPI_Offset)decom.max_nreqs[i], &dim_max_nreqs[i])
        }

        for (i=0; i<decom.num_decomp; i++) {
            fix_dimids[i]    = dim_nelems[i];
            rec_dimids[i][0] = dim_time;
            rec_dimids[i][1] = dim_nelems[i];
        }
    }

    /* define variables related to decompositions */
    nvars_decomp = 0;
    if (cfg.strategy == blob) {
        if (cfg.api == adios)
            nvars_decomp = decom.num_decomp + 1;
        else
            nvars_decomp = NVARS_DECOMP * decom.num_decomp;

        err = inq_var_decomp(cfg, decom, driver, cmeta, ncid, dim_time,
                             dim_nblobs, dim_max_nreqs, g_dimids);
        CHECK_ERR
    }

    /* There are 52 climate variables:
     *    0 scalar variables     + 52 array variables
     *   11 fixed-size variables + 41 record variables
     *   11 not partitioned      + 41 partitioned
     */

    varp = vars + nvars_decomp - 1;

    /* double salinitySurfaceRestoringTendency(Time, nCells) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    INQ_VAR("salinitySurfaceRestoringTendency", NC_DOUBLE, 2, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("units", "m PSU/s")
    GET_ATTR_TXT("long_name", "salinity tendency due to surface restoring")

    /* double vertTransportVelocityTop(Time, nCells, nVertLevelsP1) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevelsP1;
    INQ_VAR("vertTransportVelocityTop", NC_DOUBLE, 3, dimids, REC_ITYPE, 5)
    GET_ATTR_TXT("units", "m s^{-1}")
    GET_ATTR_TXT("long_name", "vertical tracer-transport velocity defined at center (horizontally) and top (vertically) of cell.  This is not the vertical ALE transport, but is Eulerian (fixed-frame) in the vertical, and computed from the continuity equation from the horizontal total tracer-transport velocity.")

    /* double vertGMBolusVelocityTop(Time, nCells, nVertLevelsP1) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevelsP1;
    INQ_VAR("vertGMBolusVelocityTop", NC_DOUBLE, 3, dimids, REC_ITYPE, 5)
    GET_ATTR_TXT("units", "m s^{-1}")
    GET_ATTR_TXT("long_name", "vertical tracer-transport velocity defined at center (horizontally) and top (vertically) of cell.  This is not the vertical ALE transport, but is Eulerian (fixed-frame) in the vertical, and computed from the continuity equation from the horizontal GM Bolus velocity.")

    /* double vertAleTransportTop(Time, nCells, nVertLevelsP1) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevelsP1;
    INQ_VAR("vertAleTransportTop", NC_DOUBLE, 3, dimids, REC_ITYPE, 5)
    GET_ATTR_TXT("units", "m s^{-1}")
    GET_ATTR_TXT("long_name", "vertical transport through the layer interface at the top of the cell")

    /* double tendSSH(Time, nCells) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    INQ_VAR("tendSSH", NC_DOUBLE, 2, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("units", "m s^{-1}")
    GET_ATTR_TXT("long_name", "time tendency of sea-surface height")

    /* double layerThickness(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("layerThickness", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("long_name", "layer thickness")

    /* double normalVelocity(Time, nEdges, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nEdges;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("normalVelocity", NC_DOUBLE, 3, dimids, REC_ITYPE, 3)
    GET_ATTR_TXT("units", "m s^{-1}")
    GET_ATTR_TXT("long_name", "horizonal velocity, normal component to an edge")

    /* double ssh(Time, nCells) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    INQ_VAR("ssh", NC_DOUBLE, 2, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("long_name", "sea surface height")

    /* int maxLevelEdgeTop(nEdges) */
    INQ_VAR("maxLevelEdgeTop", NC_INT, 1, &dim_nEdges, MPI_INT, 1)
    GET_ATTR_TXT("units", "unitless")
    GET_ATTR_TXT("long_name", "Index to the last edge in a column with active ocean cells on both sides of it.")

    /* double vertCoordMovementWeights(nVertLevels) */
    INQ_VAR("vertCoordMovementWeights", NC_DOUBLE, 1, &dim_nVertLevels, REC_ITYPE, -1)
    GET_ATTR_TXT("units", "unitless")
    GET_ATTR_TXT("long_name", "Weights used for distribution of sea surface height perturbations through multiple vertical levels.")

    /* int edgeMask(nEdges, nVertLevels) */
    dimids[0] = dim_nEdges;
    dimids[1] = dim_nVertLevels;
    INQ_VAR("edgeMask", NC_INT, 2, dimids, MPI_INT, 3)
    GET_ATTR_TXT("units", "unitless")
    GET_ATTR_TXT("long_name", "Mask on edges that determines if computations should be done on edges.")

    /* int cellMask(nCells, nVertLevels) */
    dimids[0] = dim_nCells;
    dimids[1] = dim_nVertLevels;
    INQ_VAR("cellMask", NC_INT, 2, dimids, MPI_INT, 2)
    GET_ATTR_TXT("units", "unitless")
    GET_ATTR_TXT("long_name", "Mask on cells that determines if computations should be done on cells.")

    /* int vertexMask(nVertices, nVertLevels) */
    dimids[0] = dim_nVertices;
    dimids[1] = dim_nVertLevels;
    INQ_VAR("vertexMask", NC_INT, 2, dimids, MPI_INT, 4)
    GET_ATTR_TXT("units", "unitless")
    GET_ATTR_TXT("long_name", "Mask on vertices that determines if computations should be done on vertices.")

    /* double refZMid(nVertLevels) */
    INQ_VAR("refZMid", NC_DOUBLE, 1, &dim_nVertLevels, REC_ITYPE, -1)
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("long_name", "Reference mid z-coordinate of ocean for each vertical level. This has a negative value.")

    /* double refLayerThickness(nVertLevels) */
    INQ_VAR("refLayerThickness", NC_DOUBLE, 1, &dim_nVertLevels, REC_ITYPE, -1)
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("long_name", "Reference layerThickness of ocean for each vertical level.")

    /* char xtime(Time, StrLen) */
    dimids[0] = dim_time;
    dimids[1] = dim_StrLen;
    INQ_VAR("xtime", NC_CHAR, 2, dimids, MPI_CHAR, -1)
    GET_ATTR_TXT("units", "unitless")
    GET_ATTR_TXT("long_name", "model time, with format \'YYYY-MM-DD_HH:MM:SS\'")

    /* double kineticEnergyCell(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("kineticEnergyCell", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("units", "m^2 s^{-2}")
    GET_ATTR_TXT("long_name", "kinetic energy of horizontal velocity on cells")

    /* double relativeVorticityCell(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("relativeVorticityCell", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("units", "s^{-1}")
    GET_ATTR_TXT("long_name", "curl of horizontal velocity, averaged from vertices to cell centers")

    /* double relativeVorticity(Time, nVertices, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nVertices;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("relativeVorticity", NC_DOUBLE, 3, dimids, REC_ITYPE, 4)
    GET_ATTR_TXT("units", "s^{-1}")
    GET_ATTR_TXT("long_name", "curl of horizontal velocity, defined at vertices")

    /* double divergence(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("divergence", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("units", "s^{-1}")
    GET_ATTR_TXT("long_name", "divergence of horizontal velocity")

    /* double areaCellGlobal(Time) */
    INQ_VAR("areaCellGlobal", NC_DOUBLE, 1, &dim_time, REC_ITYPE, -1)
    GET_ATTR_TXT("units", "m^2")
    GET_ATTR_TXT("long_name", "sum of the areaCell variable over the full domain, used to normalize global statistics")

    /* double areaEdgeGlobal(Time) */
    INQ_VAR("areaEdgeGlobal", NC_DOUBLE, 1, &dim_time, REC_ITYPE, -1)
    GET_ATTR_TXT("units", "m^2")
    GET_ATTR_TXT("long_name", "sum of the areaEdge variable over the full domain, used to normalize global statistics")

    /* double areaTriangleGlobal(Time) */
    INQ_VAR("areaTriangleGlobal", NC_DOUBLE, 1, &dim_time, REC_ITYPE, -1)
    GET_ATTR_TXT("units", "m^2")
    GET_ATTR_TXT("long_name", "sum of the areaTriangle variable over the full domain, used to normalize global statistics")

    /* double volumeCellGlobal(Time) */
    INQ_VAR("volumeCellGlobal", NC_DOUBLE, 1, &dim_time, REC_ITYPE, -1)
    GET_ATTR_TXT("units", "m^3")
    GET_ATTR_TXT("long_name", "sum of the volumeCell variable over the full domain, used to normalize global statistics")

    /* double volumeEdgeGlobal(Time) */
    INQ_VAR("volumeEdgeGlobal", NC_DOUBLE, 1, &dim_time, REC_ITYPE, -1)
    GET_ATTR_TXT("units", "m^3")
    GET_ATTR_TXT("long_name", "sum of the volumeEdge variable over the full domain, used to normalize global statistics")

    /* double CFLNumberGlobal(Time) */
    INQ_VAR("CFLNumberGlobal", NC_DOUBLE, 1, &dim_time, REC_ITYPE, -1)
    GET_ATTR_TXT("units", "unitless")
    GET_ATTR_TXT("long_name", "maximum CFL number over the full domain")

    /* double BruntVaisalaFreqTop(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("BruntVaisalaFreqTop", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("units", "s^{-2}")
    GET_ATTR_TXT("long_name", "Brunt Vaisala frequency defined at the center (horizontally) and top (vertically) of cell")

    /* double vertVelocityTop(Time, nCells, nVertLevelsP1) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevelsP1;
    INQ_VAR("vertVelocityTop", NC_DOUBLE, 3, dimids, REC_ITYPE, 5)
    GET_ATTR_TXT("units", "m s^{-1}")
    GET_ATTR_TXT("long_name", "vertical velocity defined at center (horizontally) and top (vertically) of cell")

    /* double velocityZonal(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("velocityZonal", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("units", "m s^{-1}")
    GET_ATTR_TXT("long_name", "component of horizontal velocity in the eastward direction")

    /* double velocityMeridional(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("velocityMeridional", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("units", "m s^{-1}")
    GET_ATTR_TXT("long_name", "component of horizontal velocity in the northward direction")

    /* double displacedDensity(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("displacedDensity", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("units", "kg m^{-3}")
    GET_ATTR_TXT("long_name", "Density displaced adiabatically to the mid-depth one layer deeper.  That is, layer k has been displaced to the depth of layer k+1.")

    /* double potentialDensity(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("potentialDensity", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("units", "kg m^{-3}")
    GET_ATTR_TXT("long_name", "potential density: density displaced adiabatically to the mid-depth of top layer")

    /* double pressure(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("pressure", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("units", "N m^{-2}")
    GET_ATTR_TXT("long_name", "pressure used in the momentum equation")

    /* double refBottomDepth(nVertLevels) */
    INQ_VAR("refBottomDepth", NC_DOUBLE, 1, &dim_nVertLevels, REC_ITYPE, -1)
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("long_name", "Reference depth of ocean for each vertical level. Used in \'z-level\' type runs.")

    /* double zMid(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("zMid", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("long_name", "z-coordinate of the mid-depth of the layer")

    /* double bottomDepth(nCells) */
    INQ_VAR("bottomDepth", NC_DOUBLE, 1, &dim_nCells, REC_ITYPE, 0)
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("long_name", "Depth of the bottom of the ocean. Given as a positive distance from sea level.")

    /* int maxLevelCell(nCells) */
    INQ_VAR("maxLevelCell", NC_INT, 1, &dim_nCells, MPI_INT, 0)
    GET_ATTR_TXT("units", "unitless")
    GET_ATTR_TXT("long_name", "Index to the last active ocean cell in each column.")

    /* int maxLevelEdgeBot(nEdges) */
    INQ_VAR("maxLevelEdgeBot", NC_INT, 1, &dim_nEdges, MPI_INT, 1)
    GET_ATTR_TXT("units", "unitless")
    GET_ATTR_TXT("long_name", "Index to the last edge in a column with at least one active ocean cell on either side of it.")

    /* double columnIntegratedSpeed(Time, nCells) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    INQ_VAR("columnIntegratedSpeed", NC_DOUBLE, 2, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("units", "m^2 s^{-1}")
    GET_ATTR_TXT("long_name", "speed = sum(h*sqrt(2*ke)), where ke is kineticEnergyCell and the sum is over the full column at cell centers.")

    /* double temperatureHorizontalAdvectionTendency(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("temperatureHorizontalAdvectionTendency", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("long_name", "potential temperature tendency due to horizontal advection")
    GET_ATTR_TXT("units", "degrees Celsius per second")

    /* double salinityHorizontalAdvectionTendency(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("salinityHorizontalAdvectionTendency", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("long_name", "salinity tendency due to horizontal advection")
    GET_ATTR_TXT("units", "PSU per second")

    /* double temperatureVerticalAdvectionTendency(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("temperatureVerticalAdvectionTendency", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("long_name", "potential temperature tendency due to vertical advection")
    GET_ATTR_TXT("units", "degrees Celsius per second")

    /* double salinityVerticalAdvectionTendency(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("salinityVerticalAdvectionTendency", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("long_name", "salinity tendency due to vertical advection")
    GET_ATTR_TXT("units", "PSU per second")

    /* double temperatureVertMixTendency(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("temperatureVertMixTendency", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("long_name", "potential temperature tendency due to vertical mixing")
    GET_ATTR_TXT("units", "degrees Celsius per second")

    /* double salinityVertMixTendency(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("salinityVertMixTendency", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("long_name", "salinity tendency due to vertical mixing")
    GET_ATTR_TXT("units", "PSU per second")

    /* double temperatureSurfaceFluxTendency(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("temperatureSurfaceFluxTendency", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("long_name", "potential temperature tendency due to surface fluxes")
    GET_ATTR_TXT("units", "degrees Celsius per second")

    /* double salinitySurfaceFluxTendency(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("salinitySurfaceFluxTendency", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("long_name", "salinity tendency due to surface fluxes")
    GET_ATTR_TXT("units", "PSU per second")

    /* double temperatureShortWaveTendency(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("temperatureShortWaveTendency", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("units", "degrees Celsius per second")
    GET_ATTR_TXT("long_name", "potential temperature tendency due to penetrating shortwave")

    /* double temperatureNonLocalTendency(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("temperatureNonLocalTendency", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("long_name", "potential temperature tendency due to kpp non-local flux")
    GET_ATTR_TXT("units", "degrees Celsius per second")

    /* double salinityNonLocalTendency(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("salinityNonLocalTendency", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("long_name", "salinity tendency due to kpp non-local flux")
    GET_ATTR_TXT("units", "PSU per second")

    /* double temperature(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("temperature", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("long_name", "potential temperature")
    GET_ATTR_TXT("units", "degrees Celsius")

    /* double salinity(Time, nCells, nVertLevels) */
    dimids[0] = dim_time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;
    INQ_VAR("salinity", NC_DOUBLE, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("long_name", "salinity")
    GET_ATTR_TXT("units", "grams salt per kilogram seawater")

    assert(varp - vars + 1 == cfg.nvars + nvars_decomp);

err_out:
    return err;
}

