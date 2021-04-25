/*********************************************************************
 *
 * Copyright (C) 2018, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 *********************************************************************/

#include <assert.h>
#include <e3sm_io.h>

#define INQ_VID(A, B, D, E, F, C) ncmpi_inq_varid (A, B, C)
#define NOP(A, B, D, E, C)        NC_NOERR
#define NOP2(A, B, D, E, F, C)    NC_NOERR

static int define_global_attributes (int ncid) {
    int err, nerrs = 0, iattr;
    double dattr;

    /* 963 global attributes: */
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "on_a_sphere", 3, "YES");
    ERR
    dattr = 6371229.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "sphere_radius", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "model_name", 4, "mpas");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "core_name", 5, "ocean");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "history", 28, "mpirun -n 9600 ./ocean_model");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "source", 4, "MPAS");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "Conventions", 4, "MPAS");
    ERR
    err = ncmpi_put_att_text (
        ncid, NC_GLOBAL, "parent_id", 65,
        "mt6vdoeok9\nrz1tbn0bed\n555vk5hkh9\n0s1lcmezuy\n7z1uysqc5i\nwj661vqvze");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "mesh_spec", 3, "0.0");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "git_version", 16, "MPAS_GIT_VERSION");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_ocean_run_mode", 7, "forward");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_do_restart", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_restart_timestamp_name", 12, "rpointer.ocn");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_start_time", 19, "0001-01-01_00:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_stop_time", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_run_duration", 19, "0001-00-00_00:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_calendar_type", 16, "gregorian_noleap");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_write_output_on_startup", 2, "NO");
    ERR
    iattr = 0;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_pio_num_iotasks", NC_INT, 1, &iattr);
    ERR
    iattr = 1;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_pio_stride", NC_INT, 1, &iattr);
    ERR
    iattr = 3;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_num_halos", NC_INT, 1, &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_block_decomp_file_prefix", 89,
                              "/project/projectdirs/acme/inputdata/ocn/mpas-o/oRRS18to6v3/"
                              "mpas-o.graph.info.170111.part.");
    ERR
    iattr = 0;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_number_of_blocks", NC_INT, 1, &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_explicit_proc_decomp", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_proc_decomp_file_prefix", 16,
                              "graph.info.part.");
    ERR
    iattr = -1;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_init_configuration", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_expand_sphere", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_realistic_coriolis_parameter", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_write_cull_cell_mask", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_vertical_grid", 7, "uniform");
    ERR
    dattr = 1.077;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_1dCVTgenerator_stretch1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.0275;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_1dCVTgenerator_stretch2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.2;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_1dCVTgenerator_dzSeed", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_iterative_init_variable", 15,
                              "landIcePressure");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_dt", 8, "00:06:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_time_integrator", 14, "split_explicit");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_vert_coord_movement", 18,
                              "uniform_stretching");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_min_max_thickness", 2, "NO");
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_min_thickness", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 6.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_max_thickness_factor", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_dzdk_positive", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_freq_filtered_thickness", 2, "NO");
    ERR
    dattr = 5.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_thickness_filter_timescale", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_highFreqThick_restore", 2, "NO");
    ERR
    dattr = 30.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_highFreqThick_restore_time", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_highFreqThick_del2", 2, "NO");
    ERR
    dattr = 100.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_highFreqThick_del2", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_alter_ICs_for_pbcs", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_pbc_alteration_type", 9, "full_cell");
    ERR
    dattr = 0.1;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_min_pbc_fraction", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_hmix_scaleWithMesh", 3, "YES");
    ERR
    dattr = -1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_maxMeshDensity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_apvm_scale_factor", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_mom_del2", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_tracer_del2", 2, "NO");
    ERR
    dattr = 10.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_mom_del2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 10.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_tracer_del2", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_mom_del4", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_tracer_del4", 2, "NO");
    ERR
    dattr = 3200000000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_mom_del4", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_mom_del4_div_factor", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_tracer_del4", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_Leith_del2", 2, "NO");
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_Leith_parameter", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 15000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_Leith_dx", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2500.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_Leith_visc2_max", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_standardGM", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_Redi_surface_layer_tapering", 2, "NO");
    ERR
    dattr = 0.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_Redi_surface_layer_tapering_extent", NC_DOUBLE, 1,
                         &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_Redi_bottom_layer_tapering", 2, "NO");
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_Redi_bottom_layer_tapering_depth", NC_DOUBLE, 1,
                         &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_GM_Bolus_kappa_function", 8, "constant");
    ERR
    dattr = 1800.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_standardGM_tracer_kappa", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_GM_Bolus_kappa_min", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 600.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_GM_Bolus_kappa_max", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 20000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_GM_Bolus_cell_size_min", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 30000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_GM_Bolus_cell_size_max", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_Redi_kappa", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.3;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_gravWaveSpeed_trunc", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.01;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_max_relative_slope", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_mom_del2_tensor", 2, "NO");
    ERR
    dattr = 10.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_mom_del2_tensor", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_mom_del4_tensor", 2, "NO");
    ERR
    dattr = 50000000000000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_mom_del4_tensor", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_Rayleigh_friction", 2, "NO");
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_Rayleigh_damping_coeff", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_Rayleigh_bottom_friction", 2, "NO");
    ERR
    dattr = 0.0001;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_Rayleigh_bottom_damping_coeff", NC_DOUBLE, 1,
                         &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_cvmix", 3, "YES");
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_prandtl_number", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_cvmix_background", 3, "YES");
    ERR
    dattr = 0.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_background_diffusion", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.0001;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_background_viscosity", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_cvmix_convection", 3, "YES");
    ERR
    dattr = 1.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_convective_diffusion", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_convective_viscosity", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_cvmix_convective_basedOnBVF", 3, "YES");
    ERR
    dattr = 0.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_convective_triggerBVF", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_cvmix_shear", 3, "YES");
    ERR
    iattr = 2;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_num_ri_smooth_loops", NC_INT, 1, &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_cvmix_use_BLD_smoothing", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_cvmix_shear_mixing_scheme", 3, "KPP");
    ERR
    dattr = 0.005;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_shear_PP_nu_zero", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 5.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_shear_PP_alpha", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_shear_PP_exp", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.005;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_shear_KPP_nu_zero", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.7;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_shear_KPP_Ri_zero", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 3.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_shear_KPP_exp", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_cvmix_tidal_mixing", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_cvmix_double_diffusion", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_cvmix_kpp", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_cvmix_fixed_boundary_layer", 2, "NO");
    ERR
    dattr = 30.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_kpp_boundary_layer_depth", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 0.25;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_kpp_criticalBulkRichardsonNumber",
                         NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_cvmix_kpp_matching", 12, "SimpleShapes");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_cvmix_kpp_EkmanOBL", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_cvmix_kpp_MonObOBL", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_cvmix_kpp_interpolationOMLType", 9,
                              "quadratic");
    ERR
    dattr = 0.1;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_kpp_surface_layer_extent", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 5.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_kpp_surface_layer_averaging", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 10.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "configure_cvmix_kpp_minimum_OBL_under_sea_ice",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 100.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_kpp_stop_OBL_search", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_cvmix_kpp_use_enhanced_diff", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_cvmix_kpp_nonlocal_with_implicit_mix", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_const_visc", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_const_diff", 2, "NO");
    ERR
    dattr = 0.0001;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_vert_visc", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.e-05;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_vert_diff", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_rich_visc", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_rich_diff", 2, "NO");
    ERR
    dattr = 0.0001;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_bkrd_vert_visc", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.e-05;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_bkrd_vert_diff", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.005;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_rich_mix", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_convective_visc", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_convective_diff", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_tanh_visc", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_tanh_diff", 2, "NO");
    ERR
    dattr = 0.25;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_max_visc_tanh", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.0001;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_min_visc_tanh", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.025;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_max_diff_tanh", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.e-05;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_min_diff_tanh", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -100.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_zMid_tanh", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 100.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_zWidth_tanh", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_bulk_wind_stress", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_bulk_thickness_flux", 3, "YES");
    ERR
    dattr = 0.001;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_flux_attenuation_coefficient", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 10.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_flux_attenuation_coefficient_runoff", NC_DOUBLE,
                         1, &dattr);
    ERR
    dattr = 86400.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_ssh_grad_relax_timescale", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_remove_AIS_coupler_runoff", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_sw_absorption_type", 6, "jerlov");
    ERR
    iattr = 3;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_jerlov_water_type", NC_INT, 1, &iattr);
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_surface_buoyancy_depth", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_frazil_ice_formation", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_frazil_in_open_ocean", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_frazil_under_land_ice", 3, "YES");
    ERR
    dattr = 333700.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_frazil_heat_of_fusion", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_frazil_ice_density", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.1;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_frazil_fractional_thickness_limit", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 3996.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_specific_heat_sea_water", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 100.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_frazil_maximum_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 4.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_frazil_sea_ice_reference_salinity", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 0.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_frazil_land_ice_reference_salinity", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_frazil_maximum_freezing_temperature", NC_DOUBLE,
                         1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_frazil_use_surface_pressure", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_land_ice_flux_mode", 13, "pressure_only");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_land_ice_flux_formulation", 7, "Jenkins");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_land_ice_flux_useHollandJenkinsAdvDiff", 2,
                              "NO");
    ERR
    dattr = 10.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_land_ice_flux_attenuation_coefficient", NC_DOUBLE,
                         1, &dattr);
    ERR
    dattr = 10.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_land_ice_flux_boundaryLayerThickness", NC_DOUBLE,
                         1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_land_ice_flux_boundaryLayerNeighborWeight",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2009.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_land_ice_flux_cp_ice", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 918.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_land_ice_flux_rho_ice", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.0025;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_land_ice_flux_topDragCoeff", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.0001;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_land_ice_flux_ISOMIP_gammaT", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.05;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_land_ice_flux_rms_tidal_velocity", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 0.011;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_land_ice_flux_jenkins_heat_transfer_coefficient",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.00031;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_land_ice_flux_jenkins_salt_transfer_coefficient",
                         NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_vert_tracer_adv", 7, "stencil");
    ERR
    iattr = 3;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_vert_tracer_adv_order", NC_INT, 1, &iattr);
    ERR
    iattr = 3;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_horiz_tracer_adv_order", NC_INT, 1, &iattr);
    ERR
    dattr = 0.25;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_coef_3rd_order", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_monotonic", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_implicit_bottom_drag", 3, "YES");
    ERR
    dattr = 0.001;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_implicit_bottom_drag_coeff", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_explicit_bottom_drag", 2, "NO");
    ERR
    dattr = 0.001;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_explicit_bottom_drag_coeff", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1026.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_density0", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_pressure_gradient_type", 16,
                              "Jacobian_from_TS");
    ERR
    dattr = 0.5;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_common_level_weight", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_eos_type", 2, "jm");
    ERR
    dattr = -1.8;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_open_ocean_freezing_temperature_coeff_0",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_open_ocean_freezing_temperature_coeff_S",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_open_ocean_freezing_temperature_coeff_p",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_open_ocean_freezing_temperature_coeff_pS",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_open_ocean_freezing_temperature_reference_pressure",
                       NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.0622;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_land_ice_cavity_freezing_temperature_coeff_0",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -0.0563;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_land_ice_cavity_freezing_temperature_coeff_S",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -7.43e-08;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_land_ice_cavity_freezing_temperature_coeff_p",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -1.74e-10;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_land_ice_cavity_freezing_temperature_coeff_pS",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL,
                         "config_land_ice_cavity_freezing_temperature_reference_pressure",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.2;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_eos_linear_alpha", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.8;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_eos_linear_beta", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 5.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_eos_linear_Tref", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 35.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_eos_linear_Sref", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_eos_linear_densityref", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 2;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_n_ts_iter", NC_INT, 1, &iattr);
    ERR
    iattr = 1;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_n_bcl_iter_beg", NC_INT, 1, &iattr);
    ERR
    iattr = 2;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_n_bcl_iter_mid", NC_INT, 1, &iattr);
    ERR
    iattr = 2;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_n_bcl_iter_end", NC_INT, 1, &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_btr_dt", 13, "0000_00:00:12");
    ERR
    iattr = 2;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_n_btr_cor_iter", NC_INT, 1, &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_vel_correction", 3, "YES");
    ERR
    iattr = 2;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_btr_subcycle_loop_factor", NC_INT, 1, &iattr);
    ERR
    dattr = 0.5;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_btr_gam1_velWt1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_btr_gam2_SSHWt1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_btr_gam3_velWt2", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_btr_solve_SSH2", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_conduct_tests", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_test_tensors", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_tensor_test_function", 11, "sph_uCosCos");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_redi_k33", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_redi_horizontal_term1", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_redi_horizontal_term2", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_redi_horizontal_term3", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_check_zlevel_consistency", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_check_ssh_consistency", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_filter_btr_mode", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_prescribe_velocity", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_prescribe_thickness", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_include_KE_vertex", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_check_tracer_monotonicity", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_compute_active_tracer_budgets", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_thick_all_tend", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_thick_hadv", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_thick_vadv", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_thick_sflux", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_vel_all_tend", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_vel_coriolis", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_vel_pgrad", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_vel_hmix", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_vel_surface_stress", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_vel_explicit_bottom_drag", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_vel_vmix", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_vel_vadv", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_tr_all_tend", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_tr_adv", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_tr_hmix", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_tr_vmix", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_tr_sflux", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_disable_tr_nonlocalflux", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_read_nearest_restart", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_rx1_constraint", 2, "NO");
    ERR
    iattr = 20;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_rx1_outer_iter_count", NC_INT, 1, &iattr);
    ERR
    iattr = 10;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_rx1_inner_iter_count", NC_INT, 1, &iattr);
    ERR
    dattr = 0.1;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_rx1_init_inner_weight", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 5.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_rx1_max", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_rx1_horiz_smooth_weight", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_rx1_vert_smooth_weight", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.1;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_rx1_slope_weight", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_rx1_zstar_weight", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 20;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_rx1_horiz_smooth_open_ocean_cells", NC_INT, 1,
                         &iattr);
    ERR
    iattr = 3;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_rx1_min_levels", NC_INT, 1, &iattr);
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_rx1_min_layer_thickness", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 20;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_baroclinic_channel_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_baroclinic_channel_use_distances", 2, "NO");
    ERR
    dattr = 13.1000003814697;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_baroclinic_channel_surface_temperature",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 10.1000003814697;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_baroclinic_channel_bottom_temperature", NC_DOUBLE,
                         1, &dattr);
    ERR
    dattr = 1.20000004768372;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_baroclinic_channel_temperature_difference",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.0799999982118607;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_baroclinic_channel_gradient_width_frac",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 40000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_baroclinic_channel_gradient_width_dist",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_baroclinic_channel_bottom_depth", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 35.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_baroclinic_channel_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -0.000119999996968545;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_baroclinic_channel_coriolis_parameter", NC_DOUBLE,
                         1, &dattr);
    ERR
    iattr = 20;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_lock_exchange_vert_levels", NC_INT, 1, &iattr);
    ERR
    dattr = 20.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_lock_exchange_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 5.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_lock_exchange_cold_temperature", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 30.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_lock_exchange_warm_temperature", NC_DOUBLE, 1,
                         &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_lock_exchange_direction", 1, "y");
    ERR
    dattr = 35.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_lock_exchange_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_lock_exchange_layer_type", 7, "z-level");
    ERR
    dattr = 0.00999999977648258;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_lock_exchange_isopycnal_min_thickness", NC_DOUBLE,
                         1, &dattr);
    ERR
    iattr = 20;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_internal_waves_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_internal_waves_use_distances", 2, "NO");
    ERR
    dattr = 20.1000003814697;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_internal_waves_surface_temperature", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 10.1000003814697;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_internal_waves_bottom_temperature", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 2.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_internal_waves_temperature_difference", NC_DOUBLE,
                         1, &dattr);
    ERR
    dattr = 0.330000013113022;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_internal_waves_amplitude_width_frac", NC_DOUBLE,
                         1, &dattr);
    ERR
    dattr = 50000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_internal_waves_amplitude_width_dist", NC_DOUBLE,
                         1, &dattr);
    ERR
    dattr = 500.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_internal_waves_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 35.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_internal_waves_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_internal_waves_layer_type", 7, "z-level");
    ERR
    dattr = 125.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_internal_waves_isopycnal_displacement", NC_DOUBLE,
                         1, &dattr);
    ERR
    iattr = 100;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_overflow_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_overflow_use_distances", 2, "NO");
    ERR
    dattr = 2000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_overflow_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 500.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_overflow_ridge_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 10.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_overflow_plug_temperature", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 20.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_overflow_domain_temperature", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 35.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_overflow_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.100000001490116;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_overflow_plug_width_frac", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.200000002980232;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_overflow_slope_center_frac", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.0500000007450581;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_overflow_slope_width_frac", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 20000.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_overflow_plug_width_dist", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 40000.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_overflow_slope_center_dist", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 7000.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_overflow_slope_width_dist", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_overflow_layer_type", 7, "z-level");
    ERR
    dattr = 0.00999999977648258;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_overflow_isopycnal_min_thickness", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 15.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_global_ocean_minimum_depth", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_depth_file", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_depth_dimname", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_depth_varname", 4, "none");
    ERR
    dattr = 1.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_global_ocean_depth_conversion_factor", NC_DOUBLE,
                         1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_temperature_file", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_salinity_file", 4, "none");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_tracer_nlat_dimname", 4, "none");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_tracer_nlon_dimname", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_tracer_ndepth_dimname", 4,
                              "none");
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_global_ocean_tracer_depth_conversion_factor",
                         NC_DOUBLE, 1, &dattr);
    ERR
    iattr = -1;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_global_ocean_tracer_vert_levels", NC_INT, 1,
                         &iattr);
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_temperature_varname", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_salinity_varname", 4, "none");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_tracer_latlon_degrees", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_tracer_lat_varname", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_tracer_lon_varname", 4, "none");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_tracer_depth_varname", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_tracer_method", 22,
                              "bilinear_interpolation");
    ERR
    iattr = 0;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_global_ocean_smooth_TS_iterations", NC_INT, 1,
                         &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_swData_file", 4, "none");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_swData_nlat_dimname", 4, "none");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_swData_nlon_dimname", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_swData_lat_varname", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_swData_lon_varname", 4, "none");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_swData_latlon_degrees", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_swData_method", 22,
                              "bilinear_interpolation");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_chlorophyll_varname", 4, "none");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_zenithAngle_varname", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_clearSky_varname", 4, "none");
    ERR
    dattr = 4.99999987368938e-05;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_global_ocean_piston_velocity", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 0.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_global_ocean_interior_restore_rate", NC_DOUBLE, 1,
                         &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_topography_file", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_topography_nlat_dimname", 4,
                              "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_topography_nlon_dimname", 4,
                              "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_topography_latlon_degrees", 3,
                              "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_topography_lat_varname", 4,
                              "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_topography_lon_varname", 4,
                              "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_topography_varname", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_topography_has_ocean_frac", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_topography_ocean_frac_varname",
                              4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_topography_method", 22,
                              "bilinear_interpolation");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_smooth_topography", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_deepen_critical_passages", 3,
                              "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_depress_by_land_ice", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_file", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_nlat_dimname", 4,
                              "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_nlon_dimname", 4,
                              "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_latlon_degrees",
                              3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_lat_varname", 4,
                              "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_lon_varname", 4,
                              "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_global_ocean_land_ice_topo_thickness_varname", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_draft_varname", 4,
                              "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_ice_frac_varname",
                              4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_global_ocean_land_ice_topo_grounded_frac_varname", 4, "none");
    ERR
    err = ncmpi_put_att_text (
        ncid, NC_GLOBAL, "config_global_ocean_use_constant_land_ice_cavity_temperature", 2, "NO");
    ERR
    dattr = -1.79999995231628;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_global_ocean_constant_land_ice_cavity_temperature",
                       NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_cull_inland_seas", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_windstress_file", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_windstress_nlat_dimname", 4,
                              "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_windstress_nlon_dimname", 4,
                              "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_windstress_latlon_degrees", 3,
                              "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_windstress_lat_varname", 4,
                              "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_windstress_lon_varname", 4,
                              "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_windstress_zonal_varname", 4,
                              "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_windstress_meridional_varname",
                              4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_windstress_method", 22,
                              "bilinear_interpolation");
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_global_ocean_windstress_conversion_factor",
                         NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_ecosys_file", 7, "unknown");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_ecosys_forcing_file", 7,
                              "unknown");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_ecosys_nlat_dimname", 4, "none");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_ecosys_nlon_dimname", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_ecosys_ndepth_dimname", 4,
                              "none");
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_global_ocean_ecosys_depth_conversion_factor",
                         NC_DOUBLE, 1, &dattr);
    ERR
    iattr = -1;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_global_ocean_ecosys_vert_levels", NC_INT, 1,
                         &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_ecosys_lat_varname", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_ecosys_lon_varname", 4, "none");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_ecosys_depth_varname", 4, "none");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_ecosys_latlon_degrees", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_ecosys_method", 22,
                              "bilinear_interpolation");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_global_ocean_ecosys_forcing_time_dimname", 4,
                              "none");
    ERR
    iattr = 0;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_global_ocean_smooth_ecosys_iterations", NC_INT, 1,
                         &iattr);
    ERR
    iattr = 100;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_vert_levels", NC_INT, 1, &iattr);
    ERR
    dattr = 15.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_surface_temperature", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 35.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_surface_salinity", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 15.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_surface_restoring_temperature",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 35.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_surface_restoring_salinity",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 3.99999998990097e-06;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_temperature_piston_velocity",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 3.99999998990097e-06;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_salinity_piston_velocity", NC_DOUBLE,
                         1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_sensible_heat_flux", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_latent_heat_flux", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_shortwave_heat_flux", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_rain_flux", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_evaporation_flux", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 9.99999997475243e-07;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_interior_temperature_restoring_rate",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 9.99999997475243e-07;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_interior_salinity_restoring_rate",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.00999999977648258;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_temperature_gradient", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_salinity_gradient", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_temperature_gradient_mixed_layer",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_salinity_gradient_mixed_layer",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_mixed_layer_depth_temperature",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_mixed_layer_depth_salinity",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_mixed_layer_temperature_change",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_mixed_layer_salinity_change",
                         NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_vertical_grid", 7, "uniform");
    ERR
    dattr = 400.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.100000001490116;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_max_windstress", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 9.99999974737875e-05;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_coriolis_parameter", NC_DOUBLE, 1,
                         &dattr);
    ERR
    iattr = 100;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_vert_levels", NC_INT, 1, &iattr);
    ERR
    dattr = 4000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_main_channel_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -50.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_north_wall_lat", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -70.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_south_wall_lat", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_iso_ridge_flag", 3, "YES");
    ERR
    dattr = 180.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_ridge_center_lon", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_ridge_height", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2000000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_ridge_width", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_iso_plateau_flag", 3, "YES");
    ERR
    dattr = 300.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_plateau_center_lon", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -58.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_plateau_center_lat", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_plateau_height", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 200000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_plateau_radius", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1000000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_plateau_slope_width", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_iso_shelf_flag", 3, "YES");
    ERR
    dattr = 500.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_shelf_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 120000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_shelf_width", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_iso_cont_slope_flag", 3, "YES");
    ERR
    dattr = 0.00999999977648258;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_max_cont_slope", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_iso_embayment_flag", 3, "YES");
    ERR
    dattr = 60.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_embayment_center_lon", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -71.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_embayment_center_lat", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 500000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_embayment_radius", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_embayment_depth", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_iso_depression_flag", 3, "YES");
    ERR
    dattr = 60.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_depression_center_lon", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -72.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_depression_south_lat", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -65.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_depression_north_lat", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 480000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_depression_width", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 800.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_depression_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 35.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.00999999977648258;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_wind_stress_max", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.200000002980232;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_acc_wind", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -0.0500000007450581;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_asf_wind", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -65.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_wind_trans", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -5.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_heat_flux_south", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 10.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_heat_flux_middle", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -5.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_heat_flux_north", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -70.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_heat_flux_lat_ss", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -65.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_heat_flux_lat_sm", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -53.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_heat_flux_lat_mn", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 60.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_region1_center_lon", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -75.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_region1_center_lat", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 150.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_region2_center_lon", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -71.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_region2_center_lat", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 240.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_region3_center_lon", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -71.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_region3_center_lat", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 330.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_region4_center_lon", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -71.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_region4_center_lat", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_iso_heat_flux_region1_flag", 2, "NO");
    ERR
    dattr = -5.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_heat_flux_region1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 300000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_heat_flux_region1_radius", NC_DOUBLE, 1,
                         &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_iso_heat_flux_region2_flag", 2, "NO");
    ERR
    dattr = -5.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_heat_flux_region2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 240000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_heat_flux_region2_radius", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 5.80000014451798e-05;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_surface_temperature_piston_velocity",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 3.5;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_initial_temp_t1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 4.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_initial_temp_t2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1200.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_initial_temp_h0", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 500.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_initial_temp_h1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 7.50000035623088e-05;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_initial_temp_mt", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -75.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_initial_temp_latS", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -50.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_initial_temp_latN", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 10.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_temperature_sponge_t1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1000.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_temperature_sponge_h1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 120000.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_temperature_sponge_l1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 10.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_temperature_sponge_tau1", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_iso_temperature_restore_region1_flag", 3,
                              "YES");
    ERR
    dattr = -1.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_temperature_restore_t1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 600000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_temperature_restore_lcx1", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 600000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_temperature_restore_lcy1", NC_DOUBLE, 1,
                         &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_iso_temperature_restore_region2_flag", 3,
                              "YES");
    ERR
    dattr = -1.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_temperature_restore_t2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 600000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_temperature_restore_lcx2", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 250000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_temperature_restore_lcy2", NC_DOUBLE, 1,
                         &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_iso_temperature_restore_region3_flag", 3,
                              "YES");
    ERR
    dattr = -1.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_temperature_restore_t3", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 600000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_temperature_restore_lcx3", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 250000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_temperature_restore_lcy3", NC_DOUBLE, 1,
                         &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_iso_temperature_restore_region4_flag", 3,
                              "YES");
    ERR
    dattr = -1.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_temperature_restore_t4", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 600000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_temperature_restore_lcx4", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 250000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_iso_temperature_restore_lcy4", NC_DOUBLE, 1,
                         &dattr);
    ERR
    iattr = 100;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_soma_vert_levels", NC_INT, 1, &iattr);
    ERR
    dattr = 1250000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_soma_domain_width", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 35.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_soma_center_latitude", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_soma_center_longitude", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.100000001490116;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_soma_phi", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2500.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_soma_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -0.400000005960464;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_soma_shelf_width", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 100.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_soma_shelf_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_soma_ref_density", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 4.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_soma_density_difference", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 300.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_soma_thermocline_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.0500000007450581;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_soma_density_difference_linear", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 20.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_soma_surface_temperature", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 33.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_soma_surface_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_soma_use_surface_temp_restoring", 2, "NO");
    ERR
    dattr = 7.5;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_soma_surface_temp_restoring_at_center_latitude",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.5;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_soma_surface_temp_restoring_latitude_gradient",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 9.99999974737875e-06;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_soma_restoring_temp_piston_vel", NC_DOUBLE, 1,
                         &dattr);
    ERR
    iattr = 100;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_ziso_add_easterly_wind_stress_ASF", 2, "NO");
    ERR
    dattr = 800000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_wind_transition_position", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 600000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_antarctic_shelf_front_width", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = -0.0500000007450581;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_wind_stress_shelf_front_max", NC_DOUBLE, 1,
                         &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_ziso_use_slopping_bathymetry", 2, "NO");
    ERR
    dattr = 2000000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_meridional_extent", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1000000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_zonal_extent", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2500.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 500.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_shelf_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 100000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_slope_half_width", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 500000.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_slope_center_position", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -9.99999974737875e-05;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_reference_coriolis", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_coriolis_gradient", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.200000002980232;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_wind_stress_max", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 3.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_mean_restoring_temp", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_restoring_temp_dev_ta", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_restoring_temp_dev_tb", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 30.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_restoring_temp_tau", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.89999991562217e-05;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_restoring_temp_piston_vel", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 1250.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_restoring_temp_ze", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 80000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_restoring_sponge_l", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 6.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_initial_temp_t1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 3.59999990463257;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_initial_temp_t2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 300.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_initial_temp_h1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 7.50000035623088e-05;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_initial_temp_mt", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_ziso_frazil_enable", 2, "NO");
    ERR
    dattr = -3.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ziso_frazil_temperature_anomaly", NC_DOUBLE, 1,
                         &dattr);
    ERR
    iattr = 20;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_vert_levels", NC_INT, 1, &iattr);
    ERR
    dattr = 2000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_bottom_depth", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 25.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_cavity_thickness", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 500.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_slope_height", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 15000.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_edge_width", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 30000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_y1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 60000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_y2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_temperature", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 34.5;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_surface_salinity", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 34.7000007629395;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_bottom_salinity", NC_DOUBLE, 1,
                         &dattr);
    ERR
    iattr = 100;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_periodic_planar_vert_levels", NC_INT, 1, &iattr);
    ERR
    dattr = 2500.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_periodic_planar_bottom_depth", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 1.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_periodic_planar_velocity_strength", NC_DOUBLE, 1,
                         &dattr);
    ERR
    iattr = 100;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ecosys_column_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_ecosys_column_vertical_grid", 14,
                              "100layerACMEv1");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_ecosys_column_TS_filename", 7, "unknown");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_ecosys_column_ecosys_filename", 7, "unknown");
    ERR
    dattr = 6000.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_ecosys_column_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 10;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sea_mount_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_sea_mount_layer_type", 5, "sigma");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_sea_mount_stratification_type", 11,
                              "exponential");
    ERR
    dattr = 1024.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sea_mount_density_coef_linear", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 1028.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_sea_mount_density_coef_exp", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.100000001490116;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_sea_mount_density_gradient_linear", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 3.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sea_mount_density_gradient_exp", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 4500.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sea_mount_density_depth_linear", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 500.;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_sea_mount_density_depth_exp", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1028.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sea_mount_density_ref", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 5.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sea_mount_density_Tref", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.200000002980232;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sea_mount_density_alpha", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 5000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sea_mount_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 4500.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sea_mount_height", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 10000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sea_mount_radius", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 40000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sea_mount_width", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 35.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sea_mount_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -9.99999974737875e-05;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_sea_mount_coriolis_parameter", NC_DOUBLE, 1,
                         &dattr);
    ERR
    iattr = 30;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_isomip_vertical_level_distribution", 8,
                              "constant");
    ERR
    dattr = -900.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -1.89999997615814;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_temperature", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 34.4000015258789;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -1.89999997615814;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_restoring_temperature", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 1.20000004244503e-05;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_temperature_piston_velocity", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 34.4000015258789;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_restoring_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.20000004244503e-05;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_salinity_piston_velocity", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = -0.00014000000373926;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_coriolis_parameter", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_southern_boundary", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1000000.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_northern_boundary", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_western_boundary", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 500000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_eastern_boundary", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_y1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -700.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_z1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_ice_fraction1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 400000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_y2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -200.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_z2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_ice_fraction2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1000000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_y3", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -200.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_z3", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_ice_fraction3", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 36;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_isomip_plus_vertical_level_distribution", 8,
                              "constant");
    ERR
    dattr = -720.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_max_bottom_depth", NC_DOUBLE, 1,
                         &dattr);
    ERR
    iattr = 3;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_minimum_levels", NC_INT, 1, &iattr);
    ERR
    dattr = 10.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_min_column_thickness", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 0.5;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_min_ocean_fraction", NC_DOUBLE, 1,
                         &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_isomip_plus_topography_file", 27,
                              "input_geometry_processed.nc");
    ERR
    dattr = -1.89999997615814;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_init_top_temp", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -1.89999997615814;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_init_bot_temp", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 33.7999992370605;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_init_top_sal", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 34.5;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_init_bot_sal", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -1.89999997615814;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_restore_top_temp", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 1.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_restore_bot_temp", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 33.7999992370605;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_restore_top_sal", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 34.7000007629395;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_restore_bot_sal", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 10.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_restore_rate", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 200.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_restore_evap_rate", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 790000.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_restore_xMin", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 800000.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_restore_xMax", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -0.000140999996801838;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_coriolis_parameter", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 1026.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_isomip_plus_effective_density", NC_DOUBLE, 1,
                         &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_activeTracers", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_activeTracers_surface_bulk_forcing", 3,
                              "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_activeTracers_surface_restoring", 3,
                              "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_activeTracers_interior_restoring", 2,
                              "NO");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_activeTracers_exponential_decay", 2, "NO");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_activeTracers_idealAge_forcing", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_activeTracers_ttd_forcing", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_surface_salinity_monthly_restoring", 3,
                              "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_surface_salinity_monthly_restoring_compute_interval", 19,
                              "0000-00-01_00:00:00");
    ERR
    dattr = 1.585e-06;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_salinity_restoring_constant_piston_velocity",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 100;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_salinity_restoring_max_difference", NC_DOUBLE, 1,
                         &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_salinity_restoring_under_sea_ice", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_debugTracers", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_debugTracers_surface_bulk_forcing", 2,
                              "NO");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_debugTracers_surface_restoring", 2, "NO");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_debugTracers_interior_restoring", 2, "NO");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_debugTracers_exponential_decay", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_debugTracers_idealAge_forcing", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_debugTracers_ttd_forcing", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_ecosysTracers", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_ecosys_atm_co2_option", 8, "constant");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_ecosys_atm_alt_co2_option", 8, "constant");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_ecosys_atm_alt_co2_use_eco", 2, "NO");
    ERR
    dattr = 379.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_ecosys_atm_co2_constant_value", NC_DOUBLE, 1,
                         &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_ecosysTracers_surface_bulk_forcing", 2,
                              "NO");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_ecosysTracers_surface_restoring", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_ecosysTracers_interior_restoring", 2,
                              "NO");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_ecosysTracers_exponential_decay", 2, "NO");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_ecosysTracers_idealAge_forcing", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_ecosysTracers_ttd_forcing", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_ecosysTracers_surface_value", 2, "NO");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_ecosysTracers_sea_ice_coupling", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_ecosysTracers_diagnostic_fields_level1", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_ecosysTracers_diagnostic_fields_level2", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_ecosysTracers_diagnostic_fields_level3", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_ecosysTracers_diagnostic_fields_level4", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_ecosysTracers_diagnostic_fields_level5", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_DMSTracers", 2, "NO");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_DMSTracers_surface_bulk_forcing", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_DMSTracers_surface_restoring", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_DMSTracers_interior_restoring", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_DMSTracers_exponential_decay", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_DMSTracers_idealAge_forcing", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_DMSTracers_ttd_forcing", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_DMSTracers_surface_value", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_DMSTracers_sea_ice_coupling", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_MacroMoleculesTracers", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_use_MacroMoleculesTracers_surface_bulk_forcing", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_MacroMoleculesTracers_surface_restoring",
                              2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_use_MacroMoleculesTracers_interior_restoring", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_MacroMoleculesTracers_exponential_decay",
                              2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_MacroMoleculesTracers_idealAge_forcing",
                              2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_MacroMoleculesTracers_ttd_forcing", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_MacroMoleculesTracers_surface_value", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_use_MacroMoleculesTracers_sea_ice_coupling",
                              2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_globalStats_enable", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_globalStats_compute_interval", 15,
                              "output_interval");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_globalStats_compute_on_startup", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_globalStats_write_on_startup", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_globalStats_text_file", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_globalStats_directory", 16,
                              "analysis_members");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_globalStats_output_stream", 17,
                              "globalStatsOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_surfaceAreaWeightedAverages_enable", 3,
                              "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_surfaceAreaWeightedAverages_compute_on_startup", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_surfaceAreaWeightedAverages_write_on_startup", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_surfaceAreaWeightedAverages_compute_interval", 19,
                              "0000-00-00_01:00:00");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_surfaceAreaWeightedAverages_output_stream",
                            33, "surfaceAreaWeightedAveragesOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_waterMassCensus_enable", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_waterMassCensus_compute_interval", 19,
                              "0000-00-00_01:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_waterMassCensus_output_stream", 21,
                              "waterMassCensusOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_waterMassCensus_compute_on_startup", 3,
                              "YES");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_waterMassCensus_write_on_startup", 2, "NO");
    ERR
    dattr = -2.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_waterMassCensus_minTemperature", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 30.;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_waterMassCensus_maxTemperature", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 32.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_waterMassCensus_minSalinity", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 37.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_waterMassCensus_maxSalinity", NC_DOUBLE, 1,
                         &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_waterMassCensus_compute_predefined_regions", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_waterMassCensus_region_group", 0, "");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_layerVolumeWeightedAverage_enable", 3,
                              "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_layerVolumeWeightedAverage_compute_interval", 19,
                              "0000-00-00_01:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_layerVolumeWeightedAverage_compute_on_startup", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_layerVolumeWeightedAverage_write_on_startup", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_layerVolumeWeightedAverage_output_stream",
                              32, "layerVolumeWeightedAverageOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_zonalMean_enable", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_zonalMean_compute_on_startup", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_zonalMean_write_on_startup", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_zonalMean_compute_interval", 19,
                              "0000-00-00_01:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_zonalMean_output_stream", 15,
                              "zonalMeanOutput");
    ERR
    iattr = 180;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_zonalMean_num_bins", NC_INT, 1, &iattr);
    ERR
    dattr = -1.e+34;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_zonalMean_min_bin", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -1.e+34;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_zonalMean_max_bin", NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_okuboWeiss_enable", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_okuboWeiss_compute_on_startup", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_okuboWeiss_write_on_startup", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_okuboWeiss_compute_interval", 19,
                              "0000-00-00_01:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_okuboWeiss_output_stream", 16,
                              "okuboWeissOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_okuboWeiss_directory", 16,
                              "analysis_members");
    ERR
    dattr = -0.2;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_okuboWeiss_threshold_value", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 1.e-10;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_okuboWeiss_normalization", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.e-10;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_okuboWeiss_lambda2_normalization", NC_DOUBLE,
                         1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_okuboWeiss_use_lat_lon_coords", 3, "YES");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_okuboWeiss_compute_eddy_census", 3, "YES");
    ERR
    iattr = 20;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_okuboWeiss_eddy_min_cells", NC_INT, 1, &iattr);
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_meridionalHeatTransport_enable", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_meridionalHeatTransport_compute_interval",
                              19, "0000-00-00_01:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_meridionalHeatTransport_compute_on_startup", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_meridionalHeatTransport_write_on_startup",
                              2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_meridionalHeatTransport_output_stream",
                              29, "meridionalHeatTransportOutput");
    ERR
    iattr = 180;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_meridionalHeatTransport_num_bins", NC_INT, 1,
                         &iattr);
    ERR
    dattr = -1.e+34;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_meridionalHeatTransport_min_bin", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = -1.e+34;
    err = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_meridionalHeatTransport_max_bin", NC_DOUBLE, 1,
                         &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_meridionalHeatTransport_region_group", 0,
                              "");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_testComputeInterval_enable", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_testComputeInterval_compute_interval", 17,
                              "00-00-01_00:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_testComputeInterval_compute_on_startup",
                              3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_testComputeInterval_write_on_startup", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_testComputeInterval_output_stream", 25,
                              "testComputeIntervalOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_highFrequencyOutput_enable", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_highFrequencyOutput_compute_interval", 15,
                              "output_interval");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_highFrequencyOutput_output_stream", 19,
                              "highFrequencyOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_highFrequencyOutput_compute_on_startup",
                              2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_highFrequencyOutput_write_on_startup", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeFilters_enable", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeFilters_compute_interval", 2, "dt");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeFilters_output_stream", 17,
                              "timeFiltersOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeFilters_restart_stream", 18,
                              "timeFiltersRestart");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeFilters_compute_on_startup", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeFilters_write_on_startup", 2, "NO");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeFilters_initialize_filters", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeFilters_tau", 11, "90_00:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeFilters_compute_cell_centered_values",
                              3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_enable", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_compute_interval", 2, "dt");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_compute_on_startup", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_output_stream", 19,
                              "lagrPartTrackOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_restart_stream", 20,
                              "lagrPartTrackRestart");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_input_stream", 18,
                              "lagrPartTrackInput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_write_on_startup", 2, "NO");
    ERR
    iattr = 0;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_filter_number", NC_INT, 1, &iattr);
    ERR
    iattr = 2;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_timeIntegration", NC_INT, 1,
                         &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_reset_criteria", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_reset_global_timestamp", 13,
                              "0000_00:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_region_stream", 20,
                              "lagrPartTrackRegions");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_reset_if_outside_region", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_reset_if_inside_region", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_eliassenPalm_enable", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_eliassenPalm_compute_interval", 15,
                              "output_interval");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_eliassenPalm_output_stream", 18,
                              "eliassenPalmOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_eliassenPalm_restart_stream", 19,
                              "eliassenPalmRestart");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_eliassenPalm_compute_on_startup", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_eliassenPalm_write_on_startup", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_eliassenPalm_debug", 2, "NO");
    ERR
    iattr = 45;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_eliassenPalm_nBuoyancyLayers", NC_INT, 1,
                         &iattr);
    ERR
    dattr = 900.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_eliassenPalm_rhomin_buoycoor", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = 1080.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_eliassenPalm_rhomax_buoycoor", NC_DOUBLE, 1,
                         &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_enable", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_compute_interval", 19,
                              "0000-00-00_01:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_output_stream", 22,
                              "mixedLayerDepthsOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_write_on_startup", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_compute_on_startup", 3,
                              "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_Tthreshold", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_Dthreshold", 3, "YES");
    ERR
    dattr = 0.2;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_crit_temp_threshold",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.03;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_crit_dens_threshold",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 100000.;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_reference_pressure",
                         NC_DOUBLE, 1, &dattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_Tgradient", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_Dgradient", 2, "NO");
    ERR
    dattr = 5.e-07;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_temp_gradient_threshold",
                         NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 5.e-08;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_den_gradient_threshold",
                         NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 1;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_interp_method", NC_INT, 1,
                         &iattr);
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_enable", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_compute_on_startup", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_write_on_startup", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_compute_interval", 15,
                              "output_interval");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_output_stream", 24,
                              "regionalStatsDailyOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_restart_stream", 18,
                              "regionalMasksInput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_input_stream", 18,
                              "regionalMasksInput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_operation", 3, "avg");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_region_type", 4, "cell");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_region_group", 3, "all");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_1d_weighting_function",
                              3, "mul");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_2d_weighting_function",
                              3, "mul");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_1d_weighting_field", 8,
                              "areaCell");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_2d_weighting_field",
                              10, "volumeCell");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_vertical_mask", 8,
                              "cellMask");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_vertical_dimension",
                              11, "nVertLevels");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_enable", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_compute_on_startup",
                              2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_write_on_startup", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_compute_interval", 15,
                              "output_interval");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_output_stream", 25,
                              "regionalStatsWeeklyOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_restart_stream", 18,
                              "regionalMasksInput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_input_stream", 18,
                              "regionalMasksInput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_operation", 3, "avg");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_region_type", 4,
                              "cell");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_region_group", 3,
                              "all");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_regionalStatsWeekly_1d_weighting_function", 3, "mul");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_regionalStatsWeekly_2d_weighting_function", 3, "mul");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_1d_weighting_field",
                              8, "areaCell");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_2d_weighting_field",
                              10, "volumeCell");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_vertical_mask", 8,
                              "cellMask");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_vertical_dimension",
                              11, "nVertLevels");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_enable", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_compute_on_startup",
                              2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_write_on_startup", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_compute_interval",
                              15, "output_interval");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_output_stream", 26,
                              "regionalStatsMonthlyOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_restart_stream", 18,
                              "regionalMasksInput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_input_stream", 18,
                              "regionalMasksInput");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_operation", 3, "avg");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_region_type", 4,
                              "cell");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_region_group", 3,
                              "all");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_regionalStatsMonthly_1d_weighting_function", 3, "mul");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_regionalStatsMonthly_2d_weighting_function", 3, "mul");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_1d_weighting_field",
                              8, "areaCell");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_2d_weighting_field",
                              10, "volumeCell");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_vertical_mask", 8,
                              "cellMask");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_vertical_dimension",
                              11, "nVertLevels");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_enable", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_compute_on_startup",
                              2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_write_on_startup", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_compute_interval", 15,
                              "output_interval");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_output_stream", 25,
                              "regionalStatsCustomOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_restart_stream", 18,
                              "regionalMasksInput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_input_stream", 18,
                              "regionalMasksInput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_operation", 3, "avg");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_region_type", 4,
                              "cell");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_region_group", 3,
                              "all");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_regionalStatsCustom_1d_weighting_function", 3, "mul");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_regionalStatsCustom_2d_weighting_function", 3, "mul");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_1d_weighting_field",
                              8, "areaCell");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_2d_weighting_field",
                              10, "volumeCell");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_vertical_mask", 8,
                              "cellMask");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_vertical_dimension",
                              11, "nVertLevels");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_enable", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_compute_on_startup",
                              2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_write_on_startup", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_compute_interval",
                              17, "00-00-00_01:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_output_stream", 26,
                              "timeSeriesStatsDailyOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_restart_stream", 4,
                              "none");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_operation", 3, "avg");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_reference_times", 12,
                              "initial_time");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_duration_intervals",
                              15, "repeat_interval");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_repeat_intervals",
                              14, "reset_interval");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_reset_intervals", 17,
                              "00-00-01_00:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_timeSeriesStatsDaily_backward_output_offset", 17,
                              "00-00-01_00:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_enable", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_timeSeriesStatsMonthly_compute_on_startup", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_write_on_startup",
                              2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_compute_interval",
                              17, "00-00-00_01:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_output_stream", 28,
                              "timeSeriesStatsMonthlyOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_restart_stream", 4,
                              "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_operation", 3,
                              "avg");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_reference_times",
                              12, "initial_time");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_duration_intervals",
                            15, "repeat_interval");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_repeat_intervals",
                              14, "reset_interval");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_reset_intervals",
                              17, "00-01-00_00:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_timeSeriesStatsMonthly_backward_output_offset", 17,
                              "00-01-00_00:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_enable", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_timeSeriesStatsClimatology_compute_on_startup", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_timeSeriesStatsClimatology_write_on_startup", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_timeSeriesStatsClimatology_compute_interval", 17,
                              "00-00-00_01:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_output_stream",
                              32, "timeSeriesStatsClimatologyOutput");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_restart_stream",
                            33, "timeSeriesStatsClimatologyRestart");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_operation", 3,
                              "avg");
    ERR
    err = ncmpi_put_att_text (
        ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_reference_times", 71,
        "00-03-01_00:00:00;00-06-01_00:00:00;00-09-01_00:00:00;00-12-01_00:00:00");
    ERR
    err = ncmpi_put_att_text (
        ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_duration_intervals", 71,
        "00-03-00_00:00:00;00-03-00_00:00:00;00-03-00_00:00:00;00-03-00_00:00:00");
    ERR
    err = ncmpi_put_att_text (
        ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_repeat_intervals", 71,
        "01-00-00_00:00:00;01-00-00_00:00:00;01-00-00_00:00:00;01-00-00_00:00:00");
    ERR
    err = ncmpi_put_att_text (
        ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_reset_intervals", 79,
        "1000-00-00_00:00:00;1000-00-00_00:00:00;1000-00-00_00:00:00;1000-00-00_00:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_timeSeriesStatsClimatology_backward_output_offset", 17,
                              "00-03-00_00:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_enable", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_compute_on_startup",
                              2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_write_on_startup",
                              2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_compute_interval",
                              17, "00-00-00_01:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_output_stream", 27,
                              "timeSeriesStatsCustomOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_restart_stream", 28,
                              "timeSeriesStatsCustomRestart");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_operation", 3, "avg");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_reference_times",
                              12, "initial_time");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_duration_intervals",
                              15, "repeat_interval");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_repeat_intervals",
                              14, "reset_interval");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_reset_intervals",
                              17, "00-00-07_00:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL,
                              "config_AM_timeSeriesStatsCustom_backward_output_offset", 17,
                              "00-00-01_00:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_pointwiseStats_enable", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_pointwiseStats_compute_interval", 15,
                              "output_interval");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_pointwiseStats_output_stream", 20,
                              "pointwiseStatsOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_pointwiseStats_compute_on_startup", 3,
                              "YES");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_pointwiseStats_write_on_startup", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_debugDiagnostics_enable", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_debugDiagnostics_compute_interval", 2,
                              "dt");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_debugDiagnostics_output_stream", 22,
                              "debugDiagnosticsOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_debugDiagnostics_compute_on_startup", 3,
                              "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_debugDiagnostics_write_on_startup", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_debugDiagnostics_check_state", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_enable", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_compute_on_startup", 3,
                              "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_write_on_startup", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_compute_interval", 19,
                              "0010-00-00_00:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_output_stream", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_variable_a", 14,
                              "layerThickness");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_variable_b", 8, "areaCell");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_variable_c", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_variable_d", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_variable_e", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_variable_f", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_variable_g", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_variable_h", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_expression_1", 5, "a b *");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_expression_2", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_expression_3", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_expression_4", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_output_name_1", 10,
                              "volumeCell");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_output_name_2", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_output_name_3", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_rpnCalculator_output_name_4", 4, "none");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_transectTransport_enable", 2, "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_transectTransport_compute_interval", 15,
                              "output_interval");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_transectTransport_output_stream", 23,
                              "transectTransportOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_transectTransport_compute_on_startup", 3,
                              "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_transectTransport_write_on_startup", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_transectTransport_transect_group", 3,
                              "all");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_eddyProductVariables_enable", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_eddyProductVariables_compute_interval",
                              19, "0000-00-00_01:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_eddyProductVariables_output_stream", 26,
                              "eddyProductVariablesOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_eddyProductVariables_compute_on_startup",
                              3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_eddyProductVariables_write_on_startup", 2,
                              "NO");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_enable", 3, "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_compute_interval", 19,
                              "0000-00-00_01:00:00");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_output_stream", 23,
                              "mocStreamfunctionOutput");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_compute_on_startup", 3,
                              "YES");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_write_on_startup", 2,
                              "NO");
    ERR
    dattr = -1.e+34;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_min_bin", NC_DOUBLE, 1,
                         &dattr);
    ERR
    dattr = -1.e+34;
    err   = ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_max_bin", NC_DOUBLE, 1,
                         &dattr);
    ERR
    iattr = 180;
    err =
        ncmpi_put_att (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_num_bins", NC_INT, 1, &iattr);
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_vertical_velocity_value",
                            15, "vertVelocityTop");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_normal_velocity_value",
                              14, "normalVelocity");
    ERR
    err =
        ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_region_group", 3, "all");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_transect_group", 3,
                              "all");
    ERR
    err = ncmpi_put_att_text (ncid, NC_GLOBAL, "file_id", 10, "3q7jdjett8");
    ERR

fn_exit:
    return nerrs;
}

/*----< def_G_case_h0() >----------------------------------------------------*/
int def_G_case_h0 (int ncid,                    /* file ID */
                   const MPI_Offset dims_D1[1], /* dimension sizes of decomposition 1 */
                   const MPI_Offset dims_D2[1], /* dimension sizes of decomposition 2 */
                   const MPI_Offset dims_D3[2], /* dimension sizes of decomposition 3 */
                   const MPI_Offset dims_D4[2], /* dimension sizes of decomposition 4 */
                   const MPI_Offset dims_D5[2], /* dimension sizes of decomposition 5 */
                   const MPI_Offset dims_D6[2], /* dimension sizes of decomposition 6 */
                   int nvars,                   /* number of variables */
                   int *varids)                 /* variable IDs */
{
    /* Total 52 variables */
    int salinitySurfaceRestoringTendency, vertTransportVelocityTop, vertGMBolusVelocityTop,
        vertAleTransportTop, tendSSH, layerThickness, normalVelocity, ssh, maxLevelEdgeTop,
        vertCoordMovementWeights, edgeMask, cellMask, vertexMask, refZMid, refLayerThickness, xtime,
        kineticEnergyCell, relativeVorticityCell, relativeVorticity, divergence, areaCellGlobal,
        areaEdgeGlobal, areaTriangleGlobal, volumeCellGlobal, volumeEdgeGlobal, CFLNumberGlobal,
        BruntVaisalaFreqTop, vertVelocityTop, velocityZonal, velocityMeridional, displacedDensity,
        potentialDensity, pressure, refBottomDepth, zMid, bottomDepth, maxLevelCell,
        maxLevelEdgeBot, columnIntegratedSpeed, temperatureHorizontalAdvectionTendency,
        salinityHorizontalAdvectionTendency, temperatureVerticalAdvectionTendency,
        salinityVerticalAdvectionTendency, temperatureVertMixTendency, salinityVertMixTendency,
        temperatureSurfaceFluxTendency, salinitySurfaceFluxTendency, temperatureShortWaveTendency,
        temperatureNonLocalTendency, salinityNonLocalTendency, temperature, salinity;

    int i, err, nerrs = 0, dimids[3];
    int dim_nVertLevelsP1, dim_nCells, dim_Time, dim_nVertLevels, dim_nEdges, dim_nVertices,
        dim_StrLen;

    err = define_global_attributes (ncid);
    ERR

    /* define dimensions */
    err = ncmpi_def_dim (ncid, "nCells", dims_D1[0], &dim_nCells);
    ERR
    err = ncmpi_def_dim (ncid, "Time", NC_UNLIMITED, &dim_Time);
    ERR
    err = ncmpi_def_dim (ncid, "nVertLevelsP1", dims_D6[1], &dim_nVertLevelsP1);
    ERR
    err = ncmpi_def_dim (ncid, "nVertLevels", dims_D3[1], &dim_nVertLevels);
    ERR
    err = ncmpi_def_dim (ncid, "nEdges", dims_D2[0], &dim_nEdges);
    ERR
    err = ncmpi_def_dim (ncid, "nVertices", dims_D5[0], &dim_nVertices);
    ERR
    err = ncmpi_def_dim (ncid, "StrLen", 64, &dim_StrLen);
    ERR

    i = 0;

    /* define variables */
    /* 1 double (Time, nCells) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;

    err = ncmpi_def_var (ncid, "salinitySurfaceRestoringTendency", NC_DOUBLE, 2, dimids,
                         &salinitySurfaceRestoringTendency);
    ERR
    err = ncmpi_put_att_text (ncid, salinitySurfaceRestoringTendency, "units", 7, "m PSU/s");
    ERR
    err = ncmpi_put_att_text (ncid, salinitySurfaceRestoringTendency, "long_name", 42,
                              "salinity tendency due to surface restoring");
    ERR
    varids[i++] = salinitySurfaceRestoringTendency;

    /* 3 double (Time, nCells, nVertLevelsP1) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevelsP1;

    err = ncmpi_def_var (ncid, "vertTransportVelocityTop", NC_DOUBLE, 3, dimids,
                         &vertTransportVelocityTop);
    ERR
    err = ncmpi_put_att_text (ncid, vertTransportVelocityTop, "units", 8, "m s^{-1}");
    ERR
    err = ncmpi_put_att_text (
        ncid, vertTransportVelocityTop, "long_name", 280,
        "vertical tracer-transport velocity defined at center (horizontally) "
        "and top (vertically) of cell.  This is not the vertical ALE transport, "
        "but is Eulerian (fixed-frame) in the vertical, and computed from the "
        "continuity equation from the horizontal total tracer-transport velocity.");
    ERR
    varids[i++] = vertTransportVelocityTop;

    err = ncmpi_def_var (ncid, "vertGMBolusVelocityTop", NC_DOUBLE, 3, dimids,
                         &vertGMBolusVelocityTop);
    ERR
    err = ncmpi_put_att_text (ncid, vertGMBolusVelocityTop, "units", 8, "m s^{-1}");
    ERR
    err = ncmpi_put_att_text (
        ncid, vertGMBolusVelocityTop, "long_name", 266,
        "vertical tracer-transport velocity defined at center (horizontally) "
        "and top (vertically) of cell.  This is not the vertical ALE transport, "
        "but is Eulerian (fixed-frame) in the vertical, and computed from the "
        "continuity equation from the horizontal GM Bolus velocity.");
    ERR
    varids[i++] = vertGMBolusVelocityTop;

    err = ncmpi_def_var (ncid, "vertAleTransportTop", NC_DOUBLE, 3, dimids, &vertAleTransportTop);
    ERR
    err = ncmpi_put_att_text (ncid, vertAleTransportTop, "units", 8, "m s^{-1}");
    ERR
    err = ncmpi_put_att_text (ncid, vertAleTransportTop, "long_name", 69,
                              "vertical transport through "
                              "the layer interface at the top of the cell");
    ERR
    varids[i++] = vertAleTransportTop;

    /* 1 double (Time, nCells) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;

    err = ncmpi_def_var (ncid, "tendSSH", NC_DOUBLE, 2, dimids, &tendSSH);
    ERR
    err = ncmpi_put_att_text (ncid, tendSSH, "units", 8, "m s^{-1}");
    ERR
    err =
        ncmpi_put_att_text (ncid, tendSSH, "long_name", 35, "time tendency of sea-surface height");
    ERR
    varids[i++] = tendSSH;

    /* 1 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = ncmpi_def_var (ncid, "layerThickness", NC_DOUBLE, 3, dimids, &layerThickness);
    ERR
    err = ncmpi_put_att_text (ncid, layerThickness, "units", 1, "m");
    ERR
    err = ncmpi_put_att_text (ncid, layerThickness, "long_name", 15, "layer thickness");
    ERR
    varids[i++] = layerThickness;

    /* 1 double (Time, nEdges, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nEdges;
    dimids[2] = dim_nVertLevels;

    err = ncmpi_def_var (ncid, "normalVelocity", NC_DOUBLE, 3, dimids, &normalVelocity);
    ERR
    err = ncmpi_put_att_text (ncid, normalVelocity, "units", 8, "m s^{-1}");
    ERR
    err = ncmpi_put_att_text (ncid, normalVelocity, "long_name", 47,
                              "horizonal velocity, "
                              "normal component to an edge");
    ERR
    varids[i++] = normalVelocity;

    /* 1 double (Time, nCells) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;

    err = ncmpi_def_var (ncid, "ssh", NC_DOUBLE, 2, dimids, &ssh);
    ERR
    err = ncmpi_put_att_text (ncid, ssh, "units", 1, "m");
    ERR
    err = ncmpi_put_att_text (ncid, ssh, "long_name", 18, "sea surface height");
    ERR
    varids[i++] = ssh;

    /* 1 int (nEdges) */
    dimids[0] = dim_nEdges;

    err = ncmpi_def_var (ncid, "maxLevelEdgeTop", NC_INT, 1, dimids, &maxLevelEdgeTop);
    ERR
    err = ncmpi_put_att_text (ncid, maxLevelEdgeTop, "units", 8, "unitless");
    ERR
    err = ncmpi_put_att_text (ncid, maxLevelEdgeTop, "long_name", 79,
                              "Index to the last edge "
                              "in a column with active ocean cells on both sides of it.");
    ERR
    varids[i++] = maxLevelEdgeTop;

    /* 1 double (nVertLevels) */
    dimids[0] = dim_nVertLevels;

    err = ncmpi_def_var (ncid, "vertCoordMovementWeights", NC_DOUBLE, 1, dimids,
                         &vertCoordMovementWeights);
    ERR
    err = ncmpi_put_att_text (ncid, vertCoordMovementWeights, "units", 8, "unitless");
    ERR
    err = ncmpi_put_att_text (ncid, vertCoordMovementWeights, "long_name", 98,
                              "Weights used "
                              "for distribution of sea surface heigh purturbations through "
                              "multiple vertical levels.");
    ERR
    varids[i++] = vertCoordMovementWeights;

    /* 1 int (nEdges, nVertLevels) */
    dimids[0] = dim_nEdges;
    dimids[1] = dim_nVertLevels;

    err = ncmpi_def_var (ncid, "edgeMask", NC_INT, 2, dimids, &edgeMask);
    ERR
    err = ncmpi_put_att_text (ncid, edgeMask, "units", 8, "unitless");
    ERR
    err = ncmpi_put_att_text (ncid, edgeMask, "long_name", 69,
                              "Mask on edges that determines "
                              "if computations should be done on edge.");
    ERR
    varids[i++] = edgeMask;

    /* 1 int (nCells, nVertLevels) */
    dimids[0] = dim_nCells;
    dimids[1] = dim_nVertLevels;

    err = ncmpi_def_var (ncid, "cellMask", NC_INT, 2, dimids, &cellMask);
    ERR
    err = ncmpi_put_att_text (ncid, cellMask, "units", 8, "unitless");
    ERR
    err = ncmpi_put_att_text (ncid, cellMask, "long_name", 69,
                              "Mask on cells that determines "
                              "if computations should be done on cell.");
    ERR
    varids[i++] = cellMask;

    /* 1 int (nVertices, nVertLevels) */
    dimids[0] = dim_nVertices;
    dimids[1] = dim_nVertLevels;

    err = ncmpi_def_var (ncid, "vertexMask", NC_INT, 2, dimids, &vertexMask);
    ERR
    err = ncmpi_put_att_text (ncid, vertexMask, "units", 8, "unitless");
    ERR
    err = ncmpi_put_att_text (ncid, vertexMask, "long_name", 75,
                              "Mask on vertices that determines "
                              "if computations should be done on vertice.");
    ERR
    varids[i++] = vertexMask;

    /* 2 double (nVertLevels) */
    dimids[0] = dim_nVertLevels;

    err = ncmpi_def_var (ncid, "refZMid", NC_DOUBLE, 1, dimids, &refZMid);
    ERR
    err = ncmpi_put_att_text (ncid, refZMid, "units", 1, "m");
    ERR
    err = ncmpi_put_att_text (ncid, refZMid, "long_name", 87,
                              "Reference mid z-coordinate of ocean "
                              "for each vertical level. This has a negative value.");
    ERR
    varids[i++] = refZMid;

    err = ncmpi_def_var (ncid, "refLayerThickness", NC_DOUBLE, 1, dimids, &refLayerThickness);
    ERR
    err = ncmpi_put_att_text (ncid, refLayerThickness, "units", 1, "m");
    ERR
    err = ncmpi_put_att_text (ncid, refLayerThickness, "long_name", 58,
                              "Reference layerThickness "
                              "of ocean for each vertical level.");
    ERR
    varids[i++] = refLayerThickness;

    /* 1 char (Time, StrLen) */
    dimids[0] = dim_Time;
    dimids[1] = dim_StrLen;

    err = ncmpi_def_var (ncid, "xtime", NC_CHAR, 2, dimids, &xtime);
    ERR
    err = ncmpi_put_att_text (ncid, xtime, "units", 8, "unitless");
    ERR
    err = ncmpi_put_att_text (ncid, xtime, "long_name", 45,
                              "model time, with format \'YYYY-MM-DD_HH:MM:SS\'");
    ERR
    varids[i++] = xtime;

    /* 2 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = ncmpi_def_var (ncid, "kineticEnergyCell", NC_DOUBLE, 3, dimids, &kineticEnergyCell);
    ERR
    err = ncmpi_put_att_text (ncid, kineticEnergyCell, "units", 10, "m^2 s^{-2}");
    ERR
    err = ncmpi_put_att_text (ncid, kineticEnergyCell, "long_name", 45,
                              "kinetic energy of horizonal "
                              "velocity on cells");
    ERR
    varids[i++] = kineticEnergyCell;

    err =
        ncmpi_def_var (ncid, "relativeVorticityCell", NC_DOUBLE, 3, dimids, &relativeVorticityCell);
    ERR
    err = ncmpi_put_att_text (ncid, relativeVorticityCell, "units", 6, "s^{-1}");
    ERR
    err = ncmpi_put_att_text (ncid, relativeVorticityCell, "long_name", 67,
                              "curl of horizontal velocity, "
                              "averaged from vertices to cell centers");
    ERR
    varids[i++] = relativeVorticityCell;

    /* 1 double (Time, nVertices, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nVertices;
    dimids[2] = dim_nVertLevels;

    err = ncmpi_def_var (ncid, "relativeVorticity", NC_DOUBLE, 3, dimids, &relativeVorticity);
    ERR
    err = ncmpi_put_att_text (ncid, relativeVorticity, "units", 6, "s^{-1}");
    ERR
    err = ncmpi_put_att_text (ncid, relativeVorticity, "long_name", 48,
                              "curl of horizontal velocity, "
                              "defined at vertices");
    ERR
    varids[i++] = relativeVorticity;

    /* 1 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = ncmpi_def_var (ncid, "divergence", NC_DOUBLE, 3, dimids, &divergence);
    ERR
    err = ncmpi_put_att_text (ncid, divergence, "units", 6, "s^{-1}");
    ERR
    err =
        ncmpi_put_att_text (ncid, divergence, "long_name", 32, "divergence of horizonal velocity");
    ERR
    varids[i++] = divergence;

    /* 6 double (Time) */
    dimids[0] = dim_Time;

    err = ncmpi_def_var (ncid, "areaCellGlobal", NC_DOUBLE, 1, dimids, &areaCellGlobal);
    ERR
    err = ncmpi_put_att_text (ncid, areaCellGlobal, "units", 3, "m^2");
    ERR
    err = ncmpi_put_att_text (ncid, areaCellGlobal, "long_name", 86,
                              "sum of the areaCell variable over "
                              "the full domain, used to normalize global statistics");
    ERR
    varids[i++] = areaCellGlobal;

    err = ncmpi_def_var (ncid, "areaEdgeGlobal", NC_DOUBLE, 1, dimids, &areaEdgeGlobal);
    ERR
    err = ncmpi_put_att_text (ncid, areaEdgeGlobal, "units", 3, "m^2");
    ERR
    err = ncmpi_put_att_text (ncid, areaEdgeGlobal, "long_name", 86,
                              "sum of the areaEdge variable over "
                              "the full domain, used to normalize global statistics");
    ERR
    varids[i++] = areaEdgeGlobal;

    err = ncmpi_def_var (ncid, "areaTriangleGlobal", NC_DOUBLE, 1, dimids, &areaTriangleGlobal);
    ERR
    err = ncmpi_put_att_text (ncid, areaTriangleGlobal, "units", 3, "m^2");
    ERR
    err = ncmpi_put_att_text (ncid, areaTriangleGlobal, "long_name", 90,
                              "sum of the areaTriangle variable "
                              "over the full domain, used to normalize global statistics");
    ERR
    varids[i++] = areaTriangleGlobal;

    err = ncmpi_def_var (ncid, "volumeCellGlobal", NC_DOUBLE, 1, dimids, &volumeCellGlobal);
    ERR
    err = ncmpi_put_att_text (ncid, volumeCellGlobal, "units", 3, "m^3");
    ERR
    err = ncmpi_put_att_text (ncid, volumeCellGlobal, "long_name", 88,
                              "sum of the volumeCell variable over "
                              "the full domain, used to normalize global statistics");
    ERR
    varids[i++] = volumeCellGlobal;

    err = ncmpi_def_var (ncid, "volumeEdgeGlobal", NC_DOUBLE, 1, dimids, &volumeEdgeGlobal);
    ERR
    err = ncmpi_put_att_text (ncid, volumeEdgeGlobal, "units", 3, "m^3");
    ERR
    err = ncmpi_put_att_text (ncid, volumeEdgeGlobal, "long_name", 88,
                              "sum of the volumeEdge variable over "
                              "the full domain, used to normalize global statistics");
    ERR
    varids[i++] = volumeEdgeGlobal;

    err = ncmpi_def_var (ncid, "CFLNumberGlobal", NC_DOUBLE, 1, dimids, &CFLNumberGlobal);
    ERR
    err = ncmpi_put_att_text (ncid, CFLNumberGlobal, "units", 8, "unitless");
    ERR
    err = ncmpi_put_att_text (ncid, CFLNumberGlobal, "long_name", 39,
                              "maximum CFL number over the full domain");
    ERR
    varids[i++] = CFLNumberGlobal;

    /* 1 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = ncmpi_def_var (ncid, "BruntVaisalaFreqTop", NC_DOUBLE, 3, dimids, &BruntVaisalaFreqTop);
    ERR
    err = ncmpi_put_att_text (ncid, BruntVaisalaFreqTop, "units", 6, "s^{-2}");
    ERR
    err = ncmpi_put_att_text (ncid, BruntVaisalaFreqTop, "long_name", 89,
                              "Brunt Vaisala frequency defined at "
                              "the center (horizontally) and top (vertically) of cell");
    ERR
    varids[i++] = BruntVaisalaFreqTop;

    /* 1 double (Time, nCells, nVertLevelsP1) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevelsP1;

    err = ncmpi_def_var (ncid, "vertVelocityTop", NC_DOUBLE, 3, dimids, &vertVelocityTop);
    ERR
    err = ncmpi_put_att_text (ncid, vertVelocityTop, "units", 8, "m s^{-1}");
    ERR
    err = ncmpi_put_att_text (ncid, vertVelocityTop, "long_name", 79,
                              "vertical velocity defined at center "
                              "(horizontally) and top (vertically) of cell");
    ERR
    varids[i++] = vertVelocityTop;

    /* 5 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = ncmpi_def_var (ncid, "velocityZonal", NC_DOUBLE, 3, dimids, &velocityZonal);
    ERR
    err = ncmpi_put_att_text (ncid, velocityZonal, "units", 8, "m s^{-1}");
    ERR
    err = ncmpi_put_att_text (ncid, velocityZonal, "long_name", 58,
                              "component of horizontal velocity in "
                              "the eastward direction");
    ERR
    varids[i++] = velocityZonal;

    err = ncmpi_def_var (ncid, "velocityMeridional", NC_DOUBLE, 3, dimids, &velocityMeridional);
    ERR
    err = ncmpi_put_att_text (ncid, velocityMeridional, "units", 8, "m s^{-1}");
    ERR
    err = ncmpi_put_att_text (ncid, velocityMeridional, "long_name", 59,
                              "component of horizontal velocity in "
                              "the northward direction");
    ERR
    varids[i++] = velocityMeridional;

    err = ncmpi_def_var (ncid, "displacedDensity", NC_DOUBLE, 3, dimids, &displacedDensity);
    ERR
    err = ncmpi_put_att_text (ncid, displacedDensity, "units", 9, "kg m^{-3}");
    ERR
    err = ncmpi_put_att_text (
        ncid, displacedDensity, "long_name", 130,
        "Density displaced adiabatically to "
        "the mid-depth one layer deeper.  That is, layer k has been displaced to the "
        "depth of layer k+1.");
    ERR
    varids[i++] = displacedDensity;

    err = ncmpi_def_var (ncid, "potentialDensity", NC_DOUBLE, 3, dimids, &potentialDensity);
    ERR
    err = ncmpi_put_att_text (ncid, potentialDensity, "units", 9, "kg m^{-3}");
    ERR
    err = ncmpi_put_att_text (ncid, potentialDensity, "long_name", 80,
                              "potential density: density displaced "
                              "adiabatically to the mid-depth of top layer");
    ERR
    varids[i++] = potentialDensity;

    err = ncmpi_def_var (ncid, "pressure", NC_DOUBLE, 3, dimids, &pressure);
    ERR
    err = ncmpi_put_att_text (ncid, pressure, "units", 8, "N m^{-2}");
    ERR
    err = ncmpi_put_att_text (ncid, pressure, "long_name", 38,
                              "pressure used in the momentum equation");
    ERR
    varids[i++] = pressure;

    /* 1 double (nVertLevels) */
    dimids[0] = dim_nVertLevels;

    err = ncmpi_def_var (ncid, "refBottomDepth", NC_DOUBLE, 1, dimids, &refBottomDepth);
    ERR
    err = ncmpi_put_att_text (ncid, refBottomDepth, "units", 1, "m");
    ERR
    err = ncmpi_put_att_text (ncid, refBottomDepth, "long_name", 78,
                              "Reference depth of ocean for each "
                              "vertical level. Used in \'z-level\' type runs.");
    ERR
    varids[i++] = refBottomDepth;

    /* 1 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = ncmpi_def_var (ncid, "zMid", NC_DOUBLE, 3, dimids, &zMid);
    ERR
    err = ncmpi_put_att_text (ncid, zMid, "units", 1, "m");
    ERR
    err = ncmpi_put_att_text (ncid, zMid, "long_name", 42,
                              "z-coordinate of the mid-depth of the layer");
    ERR
    varids[i++] = zMid;

    /* 1 double (nCells) */
    dimids[0] = dim_nCells;

    err = ncmpi_def_var (ncid, "bottomDepth", NC_DOUBLE, 1, dimids, &bottomDepth);
    ERR
    err = ncmpi_put_att_text (ncid, bottomDepth, "units", 1, "m");
    ERR
    err = ncmpi_put_att_text (ncid, bottomDepth, "long_name", 78,
                              "Depth of the bottom of the ocean. Given "
                              "as a positive distance from sea level.");
    ERR
    varids[i++] = bottomDepth;

    /* 1 int (nCells) */
    dimids[0] = dim_nCells;

    err = ncmpi_def_var (ncid, "maxLevelCell", NC_INT, 1, dimids, &maxLevelCell);
    ERR
    err = ncmpi_put_att_text (ncid, maxLevelCell, "units", 8, "unitless");
    ERR
    err = ncmpi_put_att_text (ncid, maxLevelCell, "long_name", 51,
                              "Index to the last active ocean cell in each column.");
    ERR
    varids[i++] = maxLevelCell;

    /* 1 int (nEdges) */
    dimids[0] = dim_nEdges;

    err = ncmpi_def_var (ncid, "maxLevelEdgeBot", NC_INT, 1, dimids, &maxLevelEdgeBot);
    ERR
    err = ncmpi_put_att_text (ncid, maxLevelEdgeBot, "units", 8, "unitless");
    ERR
    err = ncmpi_put_att_text (ncid, maxLevelEdgeBot, "long_name", 92,
                              "Index to the last edge in a column with at "
                              "least one active ocean cell on either side of it.");
    ERR
    varids[i++] = maxLevelEdgeBot;

    /* 1 double (Time, nCells) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;

    err =
        ncmpi_def_var (ncid, "columnIntegratedSpeed", NC_DOUBLE, 2, dimids, &columnIntegratedSpeed);
    ERR
    err = ncmpi_put_att_text (ncid, columnIntegratedSpeed, "units", 10, "m^2 s^{-1}");
    ERR
    err = ncmpi_put_att_text (
        ncid, columnIntegratedSpeed, "long_name", 109,
        "speed = sum(h*sqrt(2*ke)), where ke "
        "is kineticEnergyCell and the sum is over the full column at cell centers.");
    ERR
    varids[i++] = columnIntegratedSpeed;

    /* 13 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = ncmpi_def_var (ncid, "temperatureHorizontalAdvectionTendency", NC_DOUBLE, 3, dimids,
                         &temperatureHorizontalAdvectionTendency);
    ERR
    err = ncmpi_put_att_text (ncid, temperatureHorizontalAdvectionTendency, "long_name", 58,
                              "potential temperature "
                              "tendency due to horizontal advection");
    ERR
    err = ncmpi_put_att_text (ncid, temperatureHorizontalAdvectionTendency, "units", 26,
                              "degrees Celsius per second");
    ERR
    varids[i++] = temperatureHorizontalAdvectionTendency;

    err = ncmpi_def_var (ncid, "salinityHorizontalAdvectionTendency", NC_DOUBLE, 3, dimids,
                         &salinityHorizontalAdvectionTendency);
    ERR
    err = ncmpi_put_att_text (ncid, salinityHorizontalAdvectionTendency, "long_name", 45,
                              "salinity tendency due "
                              "to horizontal advection");
    ERR
    err = ncmpi_put_att_text (ncid, salinityHorizontalAdvectionTendency, "units", 14,
                              "PSU per second");
    ERR
    varids[i++] = salinityHorizontalAdvectionTendency;

    err = ncmpi_def_var (ncid, "temperatureVerticalAdvectionTendency", NC_DOUBLE, 3, dimids,
                         &temperatureVerticalAdvectionTendency);
    ERR
    err = ncmpi_put_att_text (ncid, temperatureVerticalAdvectionTendency, "long_name", 56,
                              "potential temperature "
                              "tendency due to vertical advection");
    ERR
    err = ncmpi_put_att_text (ncid, temperatureVerticalAdvectionTendency, "units", 26,
                              "degrees Celsius per second");
    ERR
    varids[i++] = temperatureVerticalAdvectionTendency;

    err = ncmpi_def_var (ncid, "salinityVerticalAdvectionTendency", NC_DOUBLE, 3, dimids,
                         &salinityVerticalAdvectionTendency);
    ERR
    err = ncmpi_put_att_text (ncid, salinityVerticalAdvectionTendency, "long_name", 43,
                              "salinity tendency due "
                              "to vertical advection");
    ERR
    err =
        ncmpi_put_att_text (ncid, salinityVerticalAdvectionTendency, "units", 14, "PSU per second");
    ERR
    varids[i++] = salinityVerticalAdvectionTendency;

    err = ncmpi_def_var (ncid, "temperatureVertMixTendency", NC_DOUBLE, 3, dimids,
                         &temperatureVertMixTendency);
    ERR
    err = ncmpi_put_att_text (ncid, temperatureVertMixTendency, "long_name", 53,
                              "potential temperature tendency "
                              "due to vertical mixing");
    ERR
    err = ncmpi_put_att_text (ncid, temperatureVertMixTendency, "units", 26,
                              "degrees Celsius per second");
    ERR
    varids[i++] = temperatureVertMixTendency;

    err = ncmpi_def_var (ncid, "salinityVertMixTendency", NC_DOUBLE, 3, dimids,
                         &salinityVertMixTendency);
    ERR
    err = ncmpi_put_att_text (ncid, salinityVertMixTendency, "long_name", 40,
                              "salinity tendency due to vertical mixing");
    ERR
    err = ncmpi_put_att_text (ncid, salinityVertMixTendency, "units", 14, "PSU per second");
    ERR
    varids[i++] = salinityVertMixTendency;

    err = ncmpi_def_var (ncid, "temperatureSurfaceFluxTendency", NC_DOUBLE, 3, dimids,
                         &temperatureSurfaceFluxTendency);
    ERR
    err = ncmpi_put_att_text (ncid, temperatureSurfaceFluxTendency, "long_name", 52,
                              "potential temperature tendency "
                              "due to surface fluxes");
    ERR
    err = ncmpi_put_att_text (ncid, temperatureSurfaceFluxTendency, "units", 26,
                              "degrees Celsius per second");
    ERR
    varids[i++] = temperatureSurfaceFluxTendency;

    err = ncmpi_def_var (ncid, "salinitySurfaceFluxTendency", NC_DOUBLE, 3, dimids,
                         &salinitySurfaceFluxTendency);
    ERR
    err = ncmpi_put_att_text (ncid, salinitySurfaceFluxTendency, "long_name", 39,
                              "salinity tendency due to surface fluxes");
    ERR
    err = ncmpi_put_att_text (ncid, salinitySurfaceFluxTendency, "units", 14, "PSU per second");
    ERR
    varids[i++] = salinitySurfaceFluxTendency;

    err = ncmpi_def_var (ncid, "temperatureShortWaveTendency", NC_DOUBLE, 3, dimids,
                         &temperatureShortWaveTendency);
    ERR
    err = ncmpi_put_att_text (ncid, temperatureShortWaveTendency, "units", 26,
                              "degrees Celsius per second");
    ERR
    err = ncmpi_put_att_text (ncid, temperatureShortWaveTendency, "long_name", 59,
                              "potential temperature tendency due "
                              "to penetrating shortwave");
    ERR
    varids[i++] = temperatureShortWaveTendency;

    err = ncmpi_def_var (ncid, "temperatureNonLocalTendency", NC_DOUBLE, 3, dimids,
                         &temperatureNonLocalTendency);
    ERR
    err = ncmpi_put_att_text (ncid, temperatureNonLocalTendency, "long_name", 56,
                              "potential temperature tendency due "
                              "to kpp non-local flux");
    ERR
    err = ncmpi_put_att_text (ncid, temperatureNonLocalTendency, "units", 26,
                              "degrees Celsius per second");
    ERR
    varids[i++] = temperatureNonLocalTendency;

    err = ncmpi_def_var (ncid, "salinityNonLocalTendency", NC_DOUBLE, 3, dimids,
                         &salinityNonLocalTendency);
    ERR
    err = ncmpi_put_att_text (ncid, salinityNonLocalTendency, "long_name", 43,
                              "salinity tendency due to kpp non-local flux");
    ERR
    err = ncmpi_put_att_text (ncid, salinityNonLocalTendency, "units", 14, "PSU per second");
    ERR
    varids[i++] = salinityNonLocalTendency;

    err = ncmpi_def_var (ncid, "temperature", NC_DOUBLE, 3, dimids, &temperature);
    ERR
    err = ncmpi_put_att_text (ncid, temperature, "long_name", 21, "potential temperature");
    ERR
    err = ncmpi_put_att_text (ncid, temperature, "units", 15, "degrees Celsius");
    ERR
    varids[i++] = temperature;

    err = ncmpi_def_var (ncid, "salinity", NC_DOUBLE, 3, dimids, &salinity);
    ERR
    err = ncmpi_put_att_text (ncid, salinity, "long_name", 8, "salinity");
    ERR
    err = ncmpi_put_att_text (ncid, salinity, "units", 32, "grams salt per kilogram seawater");
    ERR
    varids[i++] = salinity;

    assert (i == nvars);

fn_exit:
    return nerrs;
}

#define GET_ATT_TEXT(A, B, C, D, E) ncmpi_get_att_text (A, B, C, txtbuf)
#define GET_ATT(A, B, C, D, E, F)   ncmpi_get_att (A, B, C, F)

static int inq_global_attributes (int ncid) {
    int err, nerrs = 0, iattr;
    double dattr;
    char txtbuf[1024];

    /* 963 global attributes: */
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "on_a_sphere", 3, "YES");
    ERR
    dattr = 6371229.;
    err   = GET_ATT (ncid, NC_GLOBAL, "sphere_radius", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "model_name", 4, "mpas");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "core_name", 5, "ocean");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "history", 28, "mpirun -n 9600 ./ocean_model");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "source", 4, "MPAS");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "Conventions", 4, "MPAS");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "parent_id", 65,
                        "mt6vdoeok9\nrz1tbn0bed\n555vk5hkh9\n0s1lcmezuy\n7z1uysqc5i\nwj661vqvze");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "mesh_spec", 3, "0.0");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "git_version", 16, "MPAS_GIT_VERSION");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_ocean_run_mode", 7, "forward");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_do_restart", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_restart_timestamp_name", 12, "rpointer.ocn");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_start_time", 19, "0001-01-01_00:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_stop_time", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_run_duration", 19, "0001-00-00_00:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_calendar_type", 16, "gregorian_noleap");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_write_output_on_startup", 2, "NO");
    ERR
    iattr = 0;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_pio_num_iotasks", NC_INT, 1, &iattr);
    ERR
    iattr = 1;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_pio_stride", NC_INT, 1, &iattr);
    ERR
    iattr = 3;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_num_halos", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_block_decomp_file_prefix", 89,
                        "/project/projectdirs/acme/inputdata/ocn/mpas-o/oRRS18to6v3/"
                        "mpas-o.graph.info.170111.part.");
    ERR
    iattr = 0;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_number_of_blocks", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_explicit_proc_decomp", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_proc_decomp_file_prefix", 16, "graph.info.part.");
    ERR
    iattr = -1;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_init_configuration", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_expand_sphere", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_realistic_coriolis_parameter", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_write_cull_cell_mask", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_vertical_grid", 7, "uniform");
    ERR
    dattr = 1.077;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_1dCVTgenerator_stretch1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.0275;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_1dCVTgenerator_stretch2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.2;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_1dCVTgenerator_dzSeed", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_iterative_init_variable", 15, "landIcePressure");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_dt", 8, "00:06:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_time_integrator", 14, "split_explicit");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_vert_coord_movement", 18, "uniform_stretching");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_min_max_thickness", 2, "NO");
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_min_thickness", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 6.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_max_thickness_factor", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_dzdk_positive", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_freq_filtered_thickness", 2, "NO");
    ERR
    dattr = 5.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_thickness_filter_timescale", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_highFreqThick_restore", 2, "NO");
    ERR
    dattr = 30.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_highFreqThick_restore_time", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_highFreqThick_del2", 2, "NO");
    ERR
    dattr = 100.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_highFreqThick_del2", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_alter_ICs_for_pbcs", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_pbc_alteration_type", 9, "full_cell");
    ERR
    dattr = 0.1;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_min_pbc_fraction", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_hmix_scaleWithMesh", 3, "YES");
    ERR
    dattr = -1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_maxMeshDensity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_apvm_scale_factor", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_mom_del2", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_tracer_del2", 2, "NO");
    ERR
    dattr = 10.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_mom_del2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 10.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_tracer_del2", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_mom_del4", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_tracer_del4", 2, "NO");
    ERR
    dattr = 3200000000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_mom_del4", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_mom_del4_div_factor", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_tracer_del4", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_Leith_del2", 2, "NO");
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_Leith_parameter", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 15000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_Leith_dx", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2500.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_Leith_visc2_max", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_standardGM", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_Redi_surface_layer_tapering", 2, "NO");
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_Redi_surface_layer_tapering_extent", NC_DOUBLE, 1,
                   &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_Redi_bottom_layer_tapering", 2, "NO");
    ERR
    dattr = 0.;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_Redi_bottom_layer_tapering_depth", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_GM_Bolus_kappa_function", 8, "constant");
    ERR
    dattr = 1800.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_standardGM_tracer_kappa", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_GM_Bolus_kappa_min", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 600.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_GM_Bolus_kappa_max", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 20000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_GM_Bolus_cell_size_min", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 30000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_GM_Bolus_cell_size_max", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_Redi_kappa", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.3;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_gravWaveSpeed_trunc", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.01;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_max_relative_slope", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_mom_del2_tensor", 2, "NO");
    ERR
    dattr = 10.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_mom_del2_tensor", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_mom_del4_tensor", 2, "NO");
    ERR
    dattr = 50000000000000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_mom_del4_tensor", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_Rayleigh_friction", 2, "NO");
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_Rayleigh_damping_coeff", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_Rayleigh_bottom_friction", 2, "NO");
    ERR
    dattr = 0.0001;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_Rayleigh_bottom_damping_coeff", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_cvmix", 3, "YES");
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_prandtl_number", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_cvmix_background", 3, "YES");
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_background_diffusion", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.0001;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_background_viscosity", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_cvmix_convection", 3, "YES");
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_convective_diffusion", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_convective_viscosity", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_cvmix_convective_basedOnBVF", 3, "YES");
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_convective_triggerBVF", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_cvmix_shear", 3, "YES");
    ERR
    iattr = 2;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_num_ri_smooth_loops", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_cvmix_use_BLD_smoothing", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_cvmix_shear_mixing_scheme", 3, "KPP");
    ERR
    dattr = 0.005;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_shear_PP_nu_zero", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 5.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_shear_PP_alpha", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_shear_PP_exp", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.005;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_shear_KPP_nu_zero", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.7;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_shear_KPP_Ri_zero", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 3.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_shear_KPP_exp", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_cvmix_tidal_mixing", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_cvmix_double_diffusion", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_cvmix_kpp", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_cvmix_fixed_boundary_layer", 2, "NO");
    ERR
    dattr = 30.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_kpp_boundary_layer_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.25;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_kpp_criticalBulkRichardsonNumber", NC_DOUBLE, 1,
                   &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_cvmix_kpp_matching", 12, "SimpleShapes");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_cvmix_kpp_EkmanOBL", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_cvmix_kpp_MonObOBL", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_cvmix_kpp_interpolationOMLType", 9, "quadratic");
    ERR
    dattr = 0.1;
    err = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_kpp_surface_layer_extent", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 5.;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_cvmix_kpp_surface_layer_averaging", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 10.;
    err   = GET_ATT (ncid, NC_GLOBAL, "configure_cvmix_kpp_minimum_OBL_under_sea_ice", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 100.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_kpp_stop_OBL_search", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_cvmix_kpp_use_enhanced_diff", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_cvmix_kpp_nonlocal_with_implicit_mix", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_const_visc", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_const_diff", 2, "NO");
    ERR
    dattr = 0.0001;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_vert_visc", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.e-05;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_vert_diff", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_rich_visc", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_rich_diff", 2, "NO");
    ERR
    dattr = 0.0001;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_bkrd_vert_visc", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.e-05;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_bkrd_vert_diff", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.005;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_rich_mix", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_convective_visc", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_convective_diff", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_tanh_visc", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_tanh_diff", 2, "NO");
    ERR
    dattr = 0.25;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_max_visc_tanh", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.0001;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_min_visc_tanh", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.025;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_max_diff_tanh", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.e-05;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_min_diff_tanh", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -100.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_zMid_tanh", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 100.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_zWidth_tanh", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_bulk_wind_stress", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_bulk_thickness_flux", 3, "YES");
    ERR
    dattr = 0.001;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_flux_attenuation_coefficient", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 10.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_flux_attenuation_coefficient_runoff", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 86400.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ssh_grad_relax_timescale", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_remove_AIS_coupler_runoff", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_sw_absorption_type", 6, "jerlov");
    ERR
    iattr = 3;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_jerlov_water_type", NC_INT, 1, &iattr);
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_surface_buoyancy_depth", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_frazil_ice_formation", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_frazil_in_open_ocean", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_frazil_under_land_ice", 3, "YES");
    ERR
    dattr = 333700.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_frazil_heat_of_fusion", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_frazil_ice_density", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.1;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_frazil_fractional_thickness_limit", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 3996.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_specific_heat_sea_water", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 100.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_frazil_maximum_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 4.;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_frazil_sea_ice_reference_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_frazil_land_ice_reference_salinity", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_frazil_maximum_freezing_temperature", NC_DOUBLE, 1,
                   &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_frazil_use_surface_pressure", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_land_ice_flux_mode", 13, "pressure_only");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_land_ice_flux_formulation", 7, "Jenkins");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_land_ice_flux_useHollandJenkinsAdvDiff", 2, "NO");
    ERR
    dattr = 10.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_land_ice_flux_attenuation_coefficient", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 10.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_land_ice_flux_boundaryLayerThickness", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_land_ice_flux_boundaryLayerNeighborWeight", NC_DOUBLE,
                   1, &dattr);
    ERR
    dattr = 2009.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_land_ice_flux_cp_ice", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 918.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_land_ice_flux_rho_ice", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.0025;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_land_ice_flux_topDragCoeff", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.0001;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_land_ice_flux_ISOMIP_gammaT", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.05;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_land_ice_flux_rms_tidal_velocity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.011;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_land_ice_flux_jenkins_heat_transfer_coefficient",
                   NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.00031;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_land_ice_flux_jenkins_salt_transfer_coefficient",
                   NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_vert_tracer_adv", 7, "stencil");
    ERR
    iattr = 3;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_vert_tracer_adv_order", NC_INT, 1, &iattr);
    ERR
    iattr = 3;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_horiz_tracer_adv_order", NC_INT, 1, &iattr);
    ERR
    dattr = 0.25;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_coef_3rd_order", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_monotonic", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_implicit_bottom_drag", 3, "YES");
    ERR
    dattr = 0.001;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_implicit_bottom_drag_coeff", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_explicit_bottom_drag", 2, "NO");
    ERR
    dattr = 0.001;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_explicit_bottom_drag_coeff", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1026.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_density0", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_pressure_gradient_type", 16, "Jacobian_from_TS");
    ERR
    dattr = 0.5;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_common_level_weight", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_eos_type", 2, "jm");
    ERR
    dattr = -1.8;
    err = GET_ATT (ncid, NC_GLOBAL, "config_open_ocean_freezing_temperature_coeff_0", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 0.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_open_ocean_freezing_temperature_coeff_S", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 0.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_open_ocean_freezing_temperature_coeff_p", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 0.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_open_ocean_freezing_temperature_coeff_pS", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_open_ocean_freezing_temperature_reference_pressure",
                   NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.0622;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_land_ice_cavity_freezing_temperature_coeff_0",
                   NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -0.0563;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_land_ice_cavity_freezing_temperature_coeff_S",
                   NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -7.43e-08;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_land_ice_cavity_freezing_temperature_coeff_p",
                   NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -1.74e-10;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_land_ice_cavity_freezing_temperature_coeff_pS",
                   NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_land_ice_cavity_freezing_temperature_reference_pressure",
                 NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.2;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_eos_linear_alpha", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.8;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_eos_linear_beta", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 5.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_eos_linear_Tref", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 35.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_eos_linear_Sref", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_eos_linear_densityref", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 2;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_n_ts_iter", NC_INT, 1, &iattr);
    ERR
    iattr = 1;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_n_bcl_iter_beg", NC_INT, 1, &iattr);
    ERR
    iattr = 2;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_n_bcl_iter_mid", NC_INT, 1, &iattr);
    ERR
    iattr = 2;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_n_bcl_iter_end", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_btr_dt", 13, "0000_00:00:12");
    ERR
    iattr = 2;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_n_btr_cor_iter", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_vel_correction", 3, "YES");
    ERR
    iattr = 2;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_btr_subcycle_loop_factor", NC_INT, 1, &iattr);
    ERR
    dattr = 0.5;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_btr_gam1_velWt1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_btr_gam2_SSHWt1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_btr_gam3_velWt2", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_btr_solve_SSH2", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_conduct_tests", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_test_tensors", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_tensor_test_function", 11, "sph_uCosCos");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_redi_k33", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_redi_horizontal_term1", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_redi_horizontal_term2", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_redi_horizontal_term3", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_check_zlevel_consistency", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_check_ssh_consistency", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_filter_btr_mode", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_prescribe_velocity", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_prescribe_thickness", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_include_KE_vertex", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_check_tracer_monotonicity", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_compute_active_tracer_budgets", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_thick_all_tend", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_thick_hadv", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_thick_vadv", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_thick_sflux", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_vel_all_tend", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_vel_coriolis", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_vel_pgrad", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_vel_hmix", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_vel_surface_stress", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_vel_explicit_bottom_drag", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_vel_vmix", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_vel_vadv", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_tr_all_tend", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_tr_adv", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_tr_hmix", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_tr_vmix", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_tr_sflux", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_disable_tr_nonlocalflux", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_read_nearest_restart", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_rx1_constraint", 2, "NO");
    ERR
    iattr = 20;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_rx1_outer_iter_count", NC_INT, 1, &iattr);
    ERR
    iattr = 10;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_rx1_inner_iter_count", NC_INT, 1, &iattr);
    ERR
    dattr = 0.1;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_rx1_init_inner_weight", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 5.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_rx1_max", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_rx1_horiz_smooth_weight", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_rx1_vert_smooth_weight", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.1;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_rx1_slope_weight", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_rx1_zstar_weight", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 20;
    err = GET_ATT (ncid, NC_GLOBAL, "config_rx1_horiz_smooth_open_ocean_cells", NC_INT, 1, &iattr);
    ERR
    iattr = 3;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_rx1_min_levels", NC_INT, 1, &iattr);
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_rx1_min_layer_thickness", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 20;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_baroclinic_channel_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_baroclinic_channel_use_distances", 2, "NO");
    ERR
    dattr = 13.1000003814697;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_baroclinic_channel_surface_temperature", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 10.1000003814697;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_baroclinic_channel_bottom_temperature", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 1.20000004768372;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_baroclinic_channel_temperature_difference", NC_DOUBLE,
                   1, &dattr);
    ERR
    dattr = 0.0799999982118607;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_baroclinic_channel_gradient_width_frac", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 40000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_baroclinic_channel_gradient_width_dist", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 1000.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_baroclinic_channel_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 35.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_baroclinic_channel_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -0.000119999996968545;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_baroclinic_channel_coriolis_parameter", NC_DOUBLE, 1,
                   &dattr);
    ERR
    iattr = 20;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_lock_exchange_vert_levels", NC_INT, 1, &iattr);
    ERR
    dattr = 20.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_lock_exchange_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 5.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_lock_exchange_cold_temperature", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 30.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_lock_exchange_warm_temperature", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_lock_exchange_direction", 1, "y");
    ERR
    dattr = 35.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_lock_exchange_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_lock_exchange_layer_type", 7, "z-level");
    ERR
    dattr = 0.00999999977648258;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_lock_exchange_isopycnal_min_thickness", NC_DOUBLE, 1,
                   &dattr);
    ERR
    iattr = 20;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_internal_waves_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_internal_waves_use_distances", 2, "NO");
    ERR
    dattr = 20.1000003814697;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_internal_waves_surface_temperature", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 10.1000003814697;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_internal_waves_bottom_temperature", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_internal_waves_temperature_difference", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 0.330000013113022;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_internal_waves_amplitude_width_frac", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 50000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_internal_waves_amplitude_width_dist", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 500.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_internal_waves_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 35.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_internal_waves_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_internal_waves_layer_type", 7, "z-level");
    ERR
    dattr = 125.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_internal_waves_isopycnal_displacement", NC_DOUBLE, 1,
                   &dattr);
    ERR
    iattr = 100;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_overflow_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_overflow_use_distances", 2, "NO");
    ERR
    dattr = 2000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_overflow_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 500.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_overflow_ridge_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 10.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_overflow_plug_temperature", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 20.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_overflow_domain_temperature", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 35.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_overflow_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.100000001490116;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_overflow_plug_width_frac", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.200000002980232;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_overflow_slope_center_frac", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.0500000007450581;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_overflow_slope_width_frac", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 20000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_overflow_plug_width_dist", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 40000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_overflow_slope_center_dist", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 7000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_overflow_slope_width_dist", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_overflow_layer_type", 7, "z-level");
    ERR
    dattr = 0.00999999977648258;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_overflow_isopycnal_min_thickness", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 15.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_global_ocean_minimum_depth", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_depth_file", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_depth_dimname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_depth_varname", 4, "none");
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_global_ocean_depth_conversion_factor", NC_DOUBLE, 1,
                   &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_temperature_file", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_salinity_file", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_tracer_nlat_dimname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_tracer_nlon_dimname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_tracer_ndepth_dimname", 4, "none");
    ERR
    dattr = 1.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_global_ocean_tracer_depth_conversion_factor", NC_DOUBLE,
                   1, &dattr);
    ERR
    iattr = -1;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_global_ocean_tracer_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_temperature_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_salinity_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_tracer_latlon_degrees", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_tracer_lat_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_tracer_lon_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_tracer_depth_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_tracer_method", 22,
                        "bilinear_interpolation");
    ERR
    iattr = 0;
    err = GET_ATT (ncid, NC_GLOBAL, "config_global_ocean_smooth_TS_iterations", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_swData_file", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_swData_nlat_dimname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_swData_nlon_dimname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_swData_lat_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_swData_lon_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_swData_latlon_degrees", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_swData_method", 22,
                        "bilinear_interpolation");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_chlorophyll_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_zenithAngle_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_clearSky_varname", 4, "none");
    ERR
    dattr = 4.99999987368938e-05;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_global_ocean_piston_velocity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_global_ocean_interior_restore_rate", NC_DOUBLE, 1,
                   &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_topography_file", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_topography_nlat_dimname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_topography_nlon_dimname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_topography_latlon_degrees", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_topography_lat_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_topography_lon_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_topography_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_topography_has_ocean_frac", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_topography_ocean_frac_varname", 4,
                        "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_topography_method", 22,
                        "bilinear_interpolation");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_smooth_topography", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_deepen_critical_passages", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_depress_by_land_ice", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_file", 4, "none");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_nlat_dimname", 4, "none");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_nlon_dimname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_latlon_degrees", 3,
                        "YES");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_lat_varname", 4, "none");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_lon_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_thickness_varname", 4,
                        "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_draft_varname", 4,
                        "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_ice_frac_varname", 4,
                        "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_land_ice_topo_grounded_frac_varname",
                        4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL,
                        "config_global_ocean_use_constant_land_ice_cavity_temperature", 2, "NO");
    ERR
    dattr = -1.79999995231628;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_global_ocean_constant_land_ice_cavity_temperature",
                   NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_cull_inland_seas", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_windstress_file", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_windstress_nlat_dimname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_windstress_nlon_dimname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_windstress_latlon_degrees", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_windstress_lat_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_windstress_lon_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_windstress_zonal_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_windstress_meridional_varname", 4,
                        "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_windstress_method", 22,
                        "bilinear_interpolation");
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_global_ocean_windstress_conversion_factor", NC_DOUBLE,
                   1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_ecosys_file", 7, "unknown");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_ecosys_forcing_file", 7, "unknown");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_ecosys_nlat_dimname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_ecosys_nlon_dimname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_ecosys_ndepth_dimname", 4, "none");
    ERR
    dattr = 1.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_global_ocean_ecosys_depth_conversion_factor", NC_DOUBLE,
                   1, &dattr);
    ERR
    iattr = -1;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_global_ocean_ecosys_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_ecosys_lat_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_ecosys_lon_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_ecosys_depth_varname", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_ecosys_latlon_degrees", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_ecosys_method", 22,
                        "bilinear_interpolation");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_global_ocean_ecosys_forcing_time_dimname", 4,
                        "none");
    ERR
    iattr = 0;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_global_ocean_smooth_ecosys_iterations", NC_INT, 1,
                   &iattr);
    ERR
    iattr = 100;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_vert_levels", NC_INT, 1, &iattr);
    ERR
    dattr = 15.;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_surface_temperature", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 35.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_surface_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 15.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_surface_restoring_temperature", NC_DOUBLE,
                   1, &dattr);
    ERR
    dattr = 35.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_surface_restoring_salinity", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 3.99999998990097e-06;
    err = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_temperature_piston_velocity", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 3.99999998990097e-06;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_salinity_piston_velocity", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 0.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_sensible_heat_flux", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_latent_heat_flux", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_shortwave_heat_flux", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_rain_flux", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_evaporation_flux", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 9.99999997475243e-07;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_interior_temperature_restoring_rate",
                   NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 9.99999997475243e-07;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_interior_salinity_restoring_rate",
                   NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.00999999977648258;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_temperature_gradient", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_salinity_gradient", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_temperature_gradient_mixed_layer",
                   NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_salinity_gradient_mixed_layer", NC_DOUBLE,
                   1, &dattr);
    ERR
    dattr = 0.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_mixed_layer_depth_temperature", NC_DOUBLE,
                   1, &dattr);
    ERR
    dattr = 0.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_mixed_layer_depth_salinity", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 0.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_mixed_layer_temperature_change", NC_DOUBLE,
                   1, &dattr);
    ERR
    dattr = 0.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_mixed_layer_salinity_change", NC_DOUBLE, 1,
                   &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_vertical_grid", 7, "uniform");
    ERR
    dattr = 400.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.100000001490116;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_max_windstress", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 9.99999974737875e-05;
    err = GET_ATT (ncid, NC_GLOBAL, "config_cvmix_WSwSBF_coriolis_parameter", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 100;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_vert_levels", NC_INT, 1, &iattr);
    ERR
    dattr = 4000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_main_channel_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -50.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_north_wall_lat", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -70.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_south_wall_lat", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_iso_ridge_flag", 3, "YES");
    ERR
    dattr = 180.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_ridge_center_lon", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_ridge_height", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2000000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_ridge_width", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_iso_plateau_flag", 3, "YES");
    ERR
    dattr = 300.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_plateau_center_lon", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -58.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_plateau_center_lat", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_plateau_height", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 200000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_plateau_radius", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1000000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_plateau_slope_width", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_iso_shelf_flag", 3, "YES");
    ERR
    dattr = 500.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_shelf_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 120000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_shelf_width", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_iso_cont_slope_flag", 3, "YES");
    ERR
    dattr = 0.00999999977648258;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_max_cont_slope", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_iso_embayment_flag", 3, "YES");
    ERR
    dattr = 60.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_embayment_center_lon", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -71.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_embayment_center_lat", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 500000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_embayment_radius", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_embayment_depth", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_iso_depression_flag", 3, "YES");
    ERR
    dattr = 60.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_depression_center_lon", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -72.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_depression_south_lat", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -65.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_depression_north_lat", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 480000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_depression_width", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 800.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_depression_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 35.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.00999999977648258;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_wind_stress_max", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.200000002980232;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_acc_wind", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -0.0500000007450581;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_asf_wind", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -65.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_wind_trans", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -5.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_heat_flux_south", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 10.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_heat_flux_middle", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -5.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_heat_flux_north", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -70.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_heat_flux_lat_ss", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -65.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_heat_flux_lat_sm", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -53.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_heat_flux_lat_mn", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 60.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_region1_center_lon", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -75.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_region1_center_lat", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 150.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_region2_center_lon", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -71.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_region2_center_lat", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 240.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_region3_center_lon", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -71.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_region3_center_lat", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 330.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_region4_center_lon", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -71.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_region4_center_lat", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_iso_heat_flux_region1_flag", 2, "NO");
    ERR
    dattr = -5.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_heat_flux_region1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 300000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_heat_flux_region1_radius", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_iso_heat_flux_region2_flag", 2, "NO");
    ERR
    dattr = -5.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_heat_flux_region2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 240000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_heat_flux_region2_radius", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 5.80000014451798e-05;
    err = GET_ATT (ncid, NC_GLOBAL, "config_iso_surface_temperature_piston_velocity", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 3.5;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_initial_temp_t1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 4.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_initial_temp_t2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1200.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_initial_temp_h0", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 500.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_initial_temp_h1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 7.50000035623088e-05;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_initial_temp_mt", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -75.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_initial_temp_latS", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -50.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_initial_temp_latN", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 10.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_temperature_sponge_t1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_temperature_sponge_h1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 120000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_temperature_sponge_l1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 10.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_temperature_sponge_tau1", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_iso_temperature_restore_region1_flag", 3, "YES");
    ERR
    dattr = -1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_temperature_restore_t1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 600000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_temperature_restore_lcx1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 600000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_temperature_restore_lcy1", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_iso_temperature_restore_region2_flag", 3, "YES");
    ERR
    dattr = -1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_temperature_restore_t2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 600000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_temperature_restore_lcx2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 250000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_temperature_restore_lcy2", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_iso_temperature_restore_region3_flag", 3, "YES");
    ERR
    dattr = -1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_temperature_restore_t3", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 600000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_temperature_restore_lcx3", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 250000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_temperature_restore_lcy3", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_iso_temperature_restore_region4_flag", 3, "YES");
    ERR
    dattr = -1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_temperature_restore_t4", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 600000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_temperature_restore_lcx4", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 250000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_iso_temperature_restore_lcy4", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 100;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_soma_vert_levels", NC_INT, 1, &iattr);
    ERR
    dattr = 1250000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_soma_domain_width", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 35.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_soma_center_latitude", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_soma_center_longitude", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.100000001490116;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_soma_phi", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2500.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_soma_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -0.400000005960464;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_soma_shelf_width", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 100.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_soma_shelf_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_soma_ref_density", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 4.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_soma_density_difference", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 300.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_soma_thermocline_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.0500000007450581;
    err = GET_ATT (ncid, NC_GLOBAL, "config_soma_density_difference_linear", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 20.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_soma_surface_temperature", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 33.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_soma_surface_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_soma_use_surface_temp_restoring", 2, "NO");
    ERR
    dattr = 7.5;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_soma_surface_temp_restoring_at_center_latitude",
                   NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.5;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_soma_surface_temp_restoring_latitude_gradient",
                   NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 9.99999974737875e-06;
    err = GET_ATT (ncid, NC_GLOBAL, "config_soma_restoring_temp_piston_vel", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 100;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_ziso_add_easterly_wind_stress_ASF", 2, "NO");
    ERR
    dattr = 800000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_wind_transition_position", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 600000.;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_ziso_antarctic_shelf_front_width", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -0.0500000007450581;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_ziso_wind_stress_shelf_front_max", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_ziso_use_slopping_bathymetry", 2, "NO");
    ERR
    dattr = 2000000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_meridional_extent", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1000000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_zonal_extent", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2500.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 500.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_shelf_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 100000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_slope_half_width", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 500000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_slope_center_position", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -9.99999974737875e-05;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_reference_coriolis", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_coriolis_gradient", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.200000002980232;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_wind_stress_max", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 3.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_mean_restoring_temp", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_restoring_temp_dev_ta", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 2.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_restoring_temp_dev_tb", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 30.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_restoring_temp_tau", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.89999991562217e-05;
    err = GET_ATT (ncid, NC_GLOBAL, "config_ziso_restoring_temp_piston_vel", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1250.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_restoring_temp_ze", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 80000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_restoring_sponge_l", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 6.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_initial_temp_t1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 3.59999990463257;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_initial_temp_t2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 300.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_initial_temp_h1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 7.50000035623088e-05;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ziso_initial_temp_mt", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_ziso_frazil_enable", 2, "NO");
    ERR
    dattr = -3.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_ziso_frazil_temperature_anomaly", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 20;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_vert_levels", NC_INT, 1, &iattr);
    ERR
    dattr = 2000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 25.;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_cavity_thickness", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 500.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_slope_height", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 15000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_edge_width", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 30000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_y1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 60000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_y2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_temperature", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 34.5;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_surface_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 34.7000007629395;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_sub_ice_shelf_2D_bottom_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 100;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_periodic_planar_vert_levels", NC_INT, 1, &iattr);
    ERR
    dattr = 2500.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_periodic_planar_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_periodic_planar_velocity_strength", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 100;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ecosys_column_vert_levels", NC_INT, 1, &iattr);
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_ecosys_column_vertical_grid", 14, "100layerACMEv1");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_ecosys_column_TS_filename", 7, "unknown");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_ecosys_column_ecosys_filename", 7, "unknown");
    ERR
    dattr = 6000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ecosys_column_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 10;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sea_mount_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_sea_mount_layer_type", 5, "sigma");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_sea_mount_stratification_type", 11, "exponential");
    ERR
    dattr = 1024.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sea_mount_density_coef_linear", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1028.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sea_mount_density_coef_exp", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.100000001490116;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_sea_mount_density_gradient_linear", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 3.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_sea_mount_density_gradient_exp", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 4500.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_sea_mount_density_depth_linear", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 500.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sea_mount_density_depth_exp", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1028.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sea_mount_density_ref", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 5.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sea_mount_density_Tref", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.200000002980232;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sea_mount_density_alpha", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 5000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sea_mount_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 4500.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sea_mount_height", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 10000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sea_mount_radius", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 40000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sea_mount_width", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 35.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sea_mount_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -9.99999974737875e-05;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_sea_mount_coriolis_parameter", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 30;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_vert_levels", NC_INT, 1, &iattr);
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_isomip_vertical_level_distribution", 8, "constant");
    ERR
    dattr = -900.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -1.89999997615814;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_temperature", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 34.4000015258789;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -1.89999997615814;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_restoring_temperature", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.20000004244503e-05;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_temperature_piston_velocity", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 34.4000015258789;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_restoring_salinity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.20000004244503e-05;
    err = GET_ATT (ncid, NC_GLOBAL, "config_isomip_salinity_piston_velocity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -0.00014000000373926;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_coriolis_parameter", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_southern_boundary", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1000000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_northern_boundary", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_western_boundary", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 500000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_eastern_boundary", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_y1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -700.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_z1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_ice_fraction1", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 400000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_y2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -200.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_z2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_ice_fraction2", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1000000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_y3", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -200.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_z3", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_ice_fraction3", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 36;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_vert_levels", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_isomip_plus_vertical_level_distribution", 8,
                        "constant");
    ERR
    dattr = -720.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_max_bottom_depth", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 3;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_minimum_levels", NC_INT, 1, &iattr);
    ERR
    dattr = 10.;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_min_column_thickness", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 0.5;
    err = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_min_ocean_fraction", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_isomip_plus_topography_file", 27,
                        "input_geometry_processed.nc");
    ERR
    dattr = -1.89999997615814;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_init_top_temp", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -1.89999997615814;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_init_bot_temp", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 33.7999992370605;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_init_top_sal", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 34.5;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_init_bot_sal", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -1.89999997615814;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_restore_top_temp", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_restore_bot_temp", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 33.7999992370605;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_restore_top_sal", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 34.7000007629395;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_restore_bot_sal", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 10.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_restore_rate", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 200.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_restore_evap_rate", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 790000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_restore_xMin", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 800000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_restore_xMax", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -0.000140999996801838;
    err = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_coriolis_parameter", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1026.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_isomip_plus_effective_density", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_activeTracers", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_activeTracers_surface_bulk_forcing", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_activeTracers_surface_restoring", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_activeTracers_interior_restoring", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_activeTracers_exponential_decay", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_activeTracers_idealAge_forcing", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_activeTracers_ttd_forcing", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_surface_salinity_monthly_restoring", 3, "YES");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_surface_salinity_monthly_restoring_compute_interval",
                      19, "0000-00-01_00:00:00");
    ERR
    dattr = 1.585e-06;
    err = GET_ATT (ncid, NC_GLOBAL, "config_salinity_restoring_constant_piston_velocity", NC_DOUBLE,
                   1, &dattr);
    ERR
    dattr = 100;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_salinity_restoring_max_difference", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_salinity_restoring_under_sea_ice", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_debugTracers", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_debugTracers_surface_bulk_forcing", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_debugTracers_surface_restoring", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_debugTracers_interior_restoring", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_debugTracers_exponential_decay", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_debugTracers_idealAge_forcing", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_debugTracers_ttd_forcing", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_ecosysTracers", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_ecosys_atm_co2_option", 8, "constant");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_ecosys_atm_alt_co2_option", 8, "constant");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_ecosys_atm_alt_co2_use_eco", 2, "NO");
    ERR
    dattr = 379.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_ecosys_atm_co2_constant_value", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_ecosysTracers_surface_bulk_forcing", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_ecosysTracers_surface_restoring", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_ecosysTracers_interior_restoring", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_ecosysTracers_exponential_decay", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_ecosysTracers_idealAge_forcing", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_ecosysTracers_ttd_forcing", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_ecosysTracers_surface_value", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_ecosysTracers_sea_ice_coupling", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_ecosysTracers_diagnostic_fields_level1", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_ecosysTracers_diagnostic_fields_level2", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_ecosysTracers_diagnostic_fields_level3", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_ecosysTracers_diagnostic_fields_level4", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_ecosysTracers_diagnostic_fields_level5", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_DMSTracers", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_DMSTracers_surface_bulk_forcing", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_DMSTracers_surface_restoring", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_DMSTracers_interior_restoring", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_DMSTracers_exponential_decay", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_DMSTracers_idealAge_forcing", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_DMSTracers_ttd_forcing", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_DMSTracers_surface_value", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_DMSTracers_sea_ice_coupling", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_MacroMoleculesTracers", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_MacroMoleculesTracers_surface_bulk_forcing", 2,
                        "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_MacroMoleculesTracers_surface_restoring", 2,
                        "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_MacroMoleculesTracers_interior_restoring", 2,
                        "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_MacroMoleculesTracers_exponential_decay", 2,
                        "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_MacroMoleculesTracers_idealAge_forcing", 2,
                        "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_MacroMoleculesTracers_ttd_forcing", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_MacroMoleculesTracers_surface_value", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_use_MacroMoleculesTracers_sea_ice_coupling", 2,
                        "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_globalStats_enable", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_globalStats_compute_interval", 15,
                        "output_interval");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_globalStats_compute_on_startup", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_globalStats_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_globalStats_text_file", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_globalStats_directory", 16, "analysis_members");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_globalStats_output_stream", 17,
                        "globalStatsOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_surfaceAreaWeightedAverages_enable", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_surfaceAreaWeightedAverages_compute_on_startup",
                        3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_surfaceAreaWeightedAverages_write_on_startup",
                        2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_surfaceAreaWeightedAverages_compute_interval",
                        19, "0000-00-00_01:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_surfaceAreaWeightedAverages_output_stream", 33,
                        "surfaceAreaWeightedAveragesOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_waterMassCensus_enable", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_waterMassCensus_compute_interval", 19,
                        "0000-00-00_01:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_waterMassCensus_output_stream", 21,
                        "waterMassCensusOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_waterMassCensus_compute_on_startup", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_waterMassCensus_write_on_startup", 2, "NO");
    ERR
    dattr = -2.;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_AM_waterMassCensus_minTemperature", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 30.;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_AM_waterMassCensus_maxTemperature", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 32.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_AM_waterMassCensus_minSalinity", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 37.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_AM_waterMassCensus_maxSalinity", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_waterMassCensus_compute_predefined_regions", 3,
                        "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_waterMassCensus_region_group", 0, "");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_layerVolumeWeightedAverage_enable", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_layerVolumeWeightedAverage_compute_interval",
                        19, "0000-00-00_01:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_layerVolumeWeightedAverage_compute_on_startup",
                        3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_layerVolumeWeightedAverage_write_on_startup", 2,
                        "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_layerVolumeWeightedAverage_output_stream", 32,
                        "layerVolumeWeightedAverageOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_zonalMean_enable", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_zonalMean_compute_on_startup", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_zonalMean_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_zonalMean_compute_interval", 19,
                        "0000-00-00_01:00:00");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_zonalMean_output_stream", 15, "zonalMeanOutput");
    ERR
    iattr = 180;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_AM_zonalMean_num_bins", NC_INT, 1, &iattr);
    ERR
    dattr = -1.e+34;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_AM_zonalMean_min_bin", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -1.e+34;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_AM_zonalMean_max_bin", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_okuboWeiss_enable", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_okuboWeiss_compute_on_startup", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_okuboWeiss_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_okuboWeiss_compute_interval", 19,
                        "0000-00-00_01:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_okuboWeiss_output_stream", 16,
                        "okuboWeissOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_okuboWeiss_directory", 16, "analysis_members");
    ERR
    dattr = -0.2;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_AM_okuboWeiss_threshold_value", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.e-10;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_AM_okuboWeiss_normalization", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1.e-10;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_AM_okuboWeiss_lambda2_normalization", NC_DOUBLE, 1,
                   &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_okuboWeiss_use_lat_lon_coords", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_okuboWeiss_compute_eddy_census", 3, "YES");
    ERR
    iattr = 20;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_AM_okuboWeiss_eddy_min_cells", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_meridionalHeatTransport_enable", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_meridionalHeatTransport_compute_interval", 19,
                        "0000-00-00_01:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_meridionalHeatTransport_compute_on_startup", 3,
                        "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_meridionalHeatTransport_write_on_startup", 2,
                        "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_meridionalHeatTransport_output_stream", 29,
                        "meridionalHeatTransportOutput");
    ERR
    iattr = 180;
    err =
        GET_ATT (ncid, NC_GLOBAL, "config_AM_meridionalHeatTransport_num_bins", NC_INT, 1, &iattr);
    ERR
    dattr = -1.e+34;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_AM_meridionalHeatTransport_min_bin", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = -1.e+34;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_AM_meridionalHeatTransport_max_bin", NC_DOUBLE, 1,
                   &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_meridionalHeatTransport_region_group", 0, "");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_testComputeInterval_enable", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_testComputeInterval_compute_interval", 17,
                        "00-00-01_00:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_testComputeInterval_compute_on_startup", 3,
                        "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_testComputeInterval_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_testComputeInterval_output_stream", 25,
                        "testComputeIntervalOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_highFrequencyOutput_enable", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_highFrequencyOutput_compute_interval", 15,
                        "output_interval");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_highFrequencyOutput_output_stream", 19,
                        "highFrequencyOutput");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_highFrequencyOutput_compute_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_highFrequencyOutput_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeFilters_enable", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeFilters_compute_interval", 2, "dt");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeFilters_output_stream", 17,
                        "timeFiltersOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeFilters_restart_stream", 18,
                        "timeFiltersRestart");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeFilters_compute_on_startup", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeFilters_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeFilters_initialize_filters", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeFilters_tau", 11, "90_00:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeFilters_compute_cell_centered_values", 3,
                        "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_enable", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_compute_interval", 2, "dt");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_compute_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_output_stream", 19,
                        "lagrPartTrackOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_restart_stream", 20,
                        "lagrPartTrackRestart");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_input_stream", 18,
                        "lagrPartTrackInput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_write_on_startup", 2, "NO");
    ERR
    iattr = 0;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_filter_number", NC_INT, 1, &iattr);
    ERR
    iattr = 2;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_timeIntegration", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_reset_criteria", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_reset_global_timestamp", 13,
                        "0000_00:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_region_stream", 20,
                        "lagrPartTrackRegions");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_reset_if_outside_region", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_lagrPartTrack_reset_if_inside_region", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_eliassenPalm_enable", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_eliassenPalm_compute_interval", 15,
                        "output_interval");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_eliassenPalm_output_stream", 18,
                        "eliassenPalmOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_eliassenPalm_restart_stream", 19,
                        "eliassenPalmRestart");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_eliassenPalm_compute_on_startup", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_eliassenPalm_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_eliassenPalm_debug", 2, "NO");
    ERR
    iattr = 45;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_AM_eliassenPalm_nBuoyancyLayers", NC_INT, 1, &iattr);
    ERR
    dattr = 900.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_AM_eliassenPalm_rhomin_buoycoor", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = 1080.;
    err = GET_ATT (ncid, NC_GLOBAL, "config_AM_eliassenPalm_rhomax_buoycoor", NC_DOUBLE, 1, &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_enable", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_compute_interval", 19,
                        "0000-00-00_01:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_output_stream", 22,
                        "mixedLayerDepthsOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_compute_on_startup", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_Tthreshold", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_Dthreshold", 3, "YES");
    ERR
    dattr = 0.2;
    err = GET_ATT (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_crit_temp_threshold", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 0.03;
    err = GET_ATT (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_crit_dens_threshold", NC_DOUBLE, 1,
                   &dattr);
    ERR
    dattr = 100000.;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_reference_pressure", NC_DOUBLE, 1,
                   &dattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_Tgradient", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_Dgradient", 2, "NO");
    ERR
    dattr = 5.e-07;
    err = GET_ATT (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_temp_gradient_threshold", NC_DOUBLE,
                   1, &dattr);
    ERR
    dattr = 5.e-08;
    err = GET_ATT (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_den_gradient_threshold", NC_DOUBLE,
                   1, &dattr);
    ERR
    iattr = 1;
    err = GET_ATT (ncid, NC_GLOBAL, "config_AM_mixedLayerDepths_interp_method", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_enable", 2, "NO");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_compute_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_compute_interval", 15,
                        "output_interval");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_output_stream", 24,
                        "regionalStatsDailyOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_restart_stream", 18,
                        "regionalMasksInput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_input_stream", 18,
                        "regionalMasksInput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_operation", 3, "avg");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_region_type", 4, "cell");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_region_group", 3, "all");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_1d_weighting_function", 3,
                        "mul");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_2d_weighting_function", 3,
                        "mul");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_1d_weighting_field", 8,
                        "areaCell");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_2d_weighting_field", 10,
                        "volumeCell");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_vertical_mask", 8, "cellMask");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsDaily_vertical_dimension", 11,
                        "nVertLevels");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_enable", 2, "NO");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_compute_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_compute_interval", 15,
                        "output_interval");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_output_stream", 25,
                        "regionalStatsWeeklyOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_restart_stream", 18,
                        "regionalMasksInput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_input_stream", 18,
                        "regionalMasksInput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_operation", 3, "avg");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_region_type", 4, "cell");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_region_group", 3, "all");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_1d_weighting_function", 3,
                        "mul");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_2d_weighting_function", 3,
                        "mul");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_1d_weighting_field", 8,
                        "areaCell");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_2d_weighting_field", 10,
                        "volumeCell");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_vertical_mask", 8,
                        "cellMask");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsWeekly_vertical_dimension", 11,
                        "nVertLevels");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_enable", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_compute_on_startup", 2,
                        "NO");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_compute_interval", 15,
                        "output_interval");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_output_stream", 26,
                        "regionalStatsMonthlyOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_restart_stream", 18,
                        "regionalMasksInput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_input_stream", 18,
                        "regionalMasksInput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_operation", 3, "avg");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_region_type", 4, "cell");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_region_group", 3, "all");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_1d_weighting_function", 3,
                        "mul");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_2d_weighting_function", 3,
                        "mul");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_1d_weighting_field", 8,
                        "areaCell");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_2d_weighting_field", 10,
                        "volumeCell");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_vertical_mask", 8,
                        "cellMask");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsMonthly_vertical_dimension", 11,
                        "nVertLevels");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_enable", 2, "NO");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_compute_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_compute_interval", 15,
                        "output_interval");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_output_stream", 25,
                        "regionalStatsCustomOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_restart_stream", 18,
                        "regionalMasksInput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_input_stream", 18,
                        "regionalMasksInput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_operation", 3, "avg");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_region_type", 4, "cell");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_region_group", 3, "all");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_1d_weighting_function", 3,
                        "mul");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_2d_weighting_function", 3,
                        "mul");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_1d_weighting_field", 8,
                        "areaCell");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_2d_weighting_field", 10,
                        "volumeCell");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_vertical_mask", 8,
                        "cellMask");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_regionalStatsCustom_vertical_dimension", 11,
                        "nVertLevels");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_enable", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_compute_on_startup", 2,
                        "NO");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_compute_interval", 17,
                        "00-00-00_01:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_output_stream", 26,
                        "timeSeriesStatsDailyOutput");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_restart_stream", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_operation", 3, "avg");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_reference_times", 12,
                        "initial_time");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_duration_intervals", 15,
                        "repeat_interval");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_repeat_intervals", 14,
                        "reset_interval");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_reset_intervals", 17,
                        "00-00-01_00:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsDaily_backward_output_offset",
                        17, "00-00-01_00:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_enable", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_compute_on_startup", 2,
                        "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_write_on_startup", 2,
                        "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_compute_interval", 17,
                        "00-00-00_01:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_output_stream", 28,
                        "timeSeriesStatsMonthlyOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_restart_stream", 4,
                        "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_operation", 3, "avg");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_reference_times", 12,
                        "initial_time");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_duration_intervals", 15,
                        "repeat_interval");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_repeat_intervals", 14,
                        "reset_interval");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_reset_intervals", 17,
                        "00-01-00_00:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsMonthly_backward_output_offset",
                        17, "00-01-00_00:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_enable", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_compute_on_startup",
                        2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_write_on_startup", 2,
                        "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_compute_interval",
                        17, "00-00-00_01:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_output_stream", 32,
                        "timeSeriesStatsClimatologyOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_restart_stream", 33,
                        "timeSeriesStatsClimatologyRestart");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_operation", 3, "avg");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_reference_times", 71,
                        "00-03-01_00:00:00;00-06-01_00:00:00;00-09-01_00:00:00;00-12-01_00:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_duration_intervals",
                        71,
                        "00-03-00_00:00:00;00-03-00_00:00:00;00-03-00_00:00:00;00-03-00_00:00:00");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_repeat_intervals", 71,
                      "01-00-00_00:00:00;01-00-00_00:00:00;01-00-00_00:00:00;01-00-00_00:00:00");
    ERR
    err = GET_ATT_TEXT (
        ncid, NC_GLOBAL, "config_AM_timeSeriesStatsClimatology_reset_intervals", 79,
        "1000-00-00_00:00:00;1000-00-00_00:00:00;1000-00-00_00:00:00;1000-00-00_00:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL,
                        "config_AM_timeSeriesStatsClimatology_backward_output_offset", 17,
                        "00-03-00_00:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_enable", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_compute_on_startup", 2,
                        "NO");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_compute_interval", 17,
                        "00-00-00_01:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_output_stream", 27,
                        "timeSeriesStatsCustomOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_restart_stream", 28,
                        "timeSeriesStatsCustomRestart");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_operation", 3, "avg");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_reference_times", 12,
                        "initial_time");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_duration_intervals", 15,
                        "repeat_interval");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_repeat_intervals", 14,
                        "reset_interval");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_reset_intervals", 17,
                        "00-00-07_00:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_timeSeriesStatsCustom_backward_output_offset",
                        17, "00-00-01_00:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_pointwiseStats_enable", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_pointwiseStats_compute_interval", 15,
                        "output_interval");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_pointwiseStats_output_stream", 20,
                        "pointwiseStatsOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_pointwiseStats_compute_on_startup", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_pointwiseStats_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_debugDiagnostics_enable", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_debugDiagnostics_compute_interval", 2, "dt");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_debugDiagnostics_output_stream", 22,
                        "debugDiagnosticsOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_debugDiagnostics_compute_on_startup", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_debugDiagnostics_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_debugDiagnostics_check_state", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_enable", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_compute_on_startup", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_compute_interval", 19,
                        "0010-00-00_00:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_output_stream", 4, "none");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_variable_a", 14, "layerThickness");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_variable_b", 8, "areaCell");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_variable_c", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_variable_d", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_variable_e", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_variable_f", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_variable_g", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_variable_h", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_expression_1", 5, "a b *");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_expression_2", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_expression_3", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_expression_4", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_output_name_1", 10, "volumeCell");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_output_name_2", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_output_name_3", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_rpnCalculator_output_name_4", 4, "none");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_transectTransport_enable", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_transectTransport_compute_interval", 15,
                        "output_interval");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_transectTransport_output_stream", 23,
                        "transectTransportOutput");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_transectTransport_compute_on_startup", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_transectTransport_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_transectTransport_transect_group", 3, "all");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_eddyProductVariables_enable", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_eddyProductVariables_compute_interval", 19,
                        "0000-00-00_01:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_eddyProductVariables_output_stream", 26,
                        "eddyProductVariablesOutput");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_eddyProductVariables_compute_on_startup", 3,
                        "YES");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_eddyProductVariables_write_on_startup", 2, "NO");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_enable", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_compute_interval", 19,
                        "0000-00-00_01:00:00");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_output_stream", 23,
                        "mocStreamfunctionOutput");
    ERR
    err =
        GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_compute_on_startup", 3, "YES");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_write_on_startup", 2, "NO");
    ERR
    dattr = -1.e+34;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_min_bin", NC_DOUBLE, 1, &dattr);
    ERR
    dattr = -1.e+34;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_max_bin", NC_DOUBLE, 1, &dattr);
    ERR
    iattr = 180;
    err   = GET_ATT (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_num_bins", NC_INT, 1, &iattr);
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_vertical_velocity_value", 15,
                        "vertVelocityTop");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_normal_velocity_value", 14,
                        "normalVelocity");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_region_group", 3, "all");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "config_AM_mocStreamfunction_transect_group", 3, "all");
    ERR
    err = GET_ATT_TEXT (ncid, NC_GLOBAL, "file_id", 10, "3q7jdjett8");
    ERR

fn_exit:
    return nerrs;
}

/*----< def_G_case_h0() >----------------------------------------------------*/
int inq_G_case_h0 (int ncid,              /* file ID */
                   MPI_Offset dims_D1[1], /* dimension sizes of decomposition 1 */
                   MPI_Offset dims_D2[1], /* dimension sizes of decomposition 2 */
                   MPI_Offset dims_D3[2], /* dimension sizes of decomposition 3 */
                   MPI_Offset dims_D4[2], /* dimension sizes of decomposition 4 */
                   MPI_Offset dims_D5[2], /* dimension sizes of decomposition 5 */
                   MPI_Offset dims_D6[2], /* dimension sizes of decomposition 6 */
                   int nvars,             /* number of variables */
                   int *varids)           /* variable IDs */
{
    /* Total 52 variables */
    int salinitySurfaceRestoringTendency, vertTransportVelocityTop, vertGMBolusVelocityTop,
        vertAleTransportTop, tendSSH, layerThickness, normalVelocity, ssh, maxLevelEdgeTop,
        vertCoordMovementWeights, edgeMask, cellMask, vertexMask, refZMid, refLayerThickness, xtime,
        kineticEnergyCell, relativeVorticityCell, relativeVorticity, divergence, areaCellGlobal,
        areaEdgeGlobal, areaTriangleGlobal, volumeCellGlobal, volumeEdgeGlobal, CFLNumberGlobal,
        BruntVaisalaFreqTop, vertVelocityTop, velocityZonal, velocityMeridional, displacedDensity,
        potentialDensity, pressure, refBottomDepth, zMid, bottomDepth, maxLevelCell,
        maxLevelEdgeBot, columnIntegratedSpeed, temperatureHorizontalAdvectionTendency,
        salinityHorizontalAdvectionTendency, temperatureVerticalAdvectionTendency,
        salinityVerticalAdvectionTendency, temperatureVertMixTendency, salinityVertMixTendency,
        temperatureSurfaceFluxTendency, salinitySurfaceFluxTendency, temperatureShortWaveTendency,
        temperatureNonLocalTendency, salinityNonLocalTendency, temperature, salinity;

    int i, err, nerrs = 0, dimids[3];
    int dim_nVertLevelsP1, dim_nCells, dim_Time, dim_nVertLevels, dim_nEdges, dim_nVertices,
        dim_StrLen;

    // err = define_global_attributes(ncid); ERR

    /* define dimensions */
    /*
    err = ncmpi_def_dim(ncid, "nCells", dims_D1[0], &dim_nCells); ERR
    err = ncmpi_def_dim(ncid, "Time", NC_UNLIMITED, &dim_Time); ERR
    err = ncmpi_def_dim(ncid, "nVertLevelsP1", dims_D6[1], &dim_nVertLevelsP1); ERR
    err = ncmpi_def_dim(ncid, "nVertLevels", dims_D3[1], &dim_nVertLevels); ERR
    err = ncmpi_def_dim(ncid, "nEdges", dims_D2[0], &dim_nEdges); ERR
    err = ncmpi_def_dim(ncid, "nVertices", dims_D5[0], &dim_nVertices); ERR
    err = ncmpi_def_dim(ncid, "StrLen", 64, &dim_StrLen); ERR
    */

    err = ncmpi_inq_dimid (ncid, "nCells", &dim_nCells);
    ERR
    err = ncmpi_inq_dimid (ncid, "nVertLevelsP1", &dim_nVertLevelsP1);
    ERR
    err = ncmpi_inq_dimid (ncid, "nVertLevels", &dim_nVertLevels);
    ERR
    err = ncmpi_inq_dimid (ncid, "nEdges", &dim_nEdges);
    ERR
    err = ncmpi_inq_dimid (ncid, "nVertices", &dim_nVertices);
    ERR

    err = ncmpi_inq_dimlen (ncid, dim_nCells, &(((MPI_Offset *)dims_D1)[0]));
    ERR
    err = ncmpi_inq_dimlen (ncid, dim_nVertLevelsP1, &(((MPI_Offset *)dims_D6)[1]));
    ERR
    err = ncmpi_inq_dimlen (ncid, dim_nVertLevels, &(((MPI_Offset *)dims_D3)[1]));
    ERR
    err = ncmpi_inq_dimlen (ncid, dim_nEdges, &(((MPI_Offset *)dims_D2)[0]));
    ERR
    err = ncmpi_inq_dimlen (ncid, dim_nVertices, &(((MPI_Offset *)dims_D5)[0]));
    ERR

    i = 0;

    /* define variables */
    /* 1 double (Time, nCells) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;

    err = INQ_VID (ncid, "salinitySurfaceRestoringTendency", NC_DOUBLE, 2, dimids,
                   &salinitySurfaceRestoringTendency);
    ERR
    err = NOP (ncid, salinitySurfaceRestoringTendency, "units", 7, "m PSU/s");
    ERR
    err = NOP (ncid, salinitySurfaceRestoringTendency, "long_name", 42,
               "salinity tendency due to surface restoring");
    ERR
    varids[i++] = salinitySurfaceRestoringTendency;

    /* 3 double (Time, nCells, nVertLevelsP1) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevelsP1;

    err =
        INQ_VID (ncid, "vertTransportVelocityTop", NC_DOUBLE, 3, dimids, &vertTransportVelocityTop);
    ERR
    err = NOP (ncid, vertTransportVelocityTop, "units", 8, "m s^{-1}");
    ERR
    err = NOP (ncid, vertTransportVelocityTop, "long_name", 280,
               "vertical tracer-transport velocity defined at center (horizontally) "
               "and top (vertically) of cell.  This is not the vertical ALE transport, "
               "but is Eulerian (fixed-frame) in the vertical, and computed from the "
               "continuity equation from the horizontal total tracer-transport velocity.");
    ERR
    varids[i++] = vertTransportVelocityTop;

    err = INQ_VID (ncid, "vertGMBolusVelocityTop", NC_DOUBLE, 3, dimids, &vertGMBolusVelocityTop);
    ERR
    err = NOP (ncid, vertGMBolusVelocityTop, "units", 8, "m s^{-1}");
    ERR
    err = NOP (ncid, vertGMBolusVelocityTop, "long_name", 266,
               "vertical tracer-transport velocity defined at center (horizontally) "
               "and top (vertically) of cell.  This is not the vertical ALE transport, "
               "but is Eulerian (fixed-frame) in the vertical, and computed from the "
               "continuity equation from the horizontal GM Bolus velocity.");
    ERR
    varids[i++] = vertGMBolusVelocityTop;

    err = INQ_VID (ncid, "vertAleTransportTop", NC_DOUBLE, 3, dimids, &vertAleTransportTop);
    ERR
    err = NOP (ncid, vertAleTransportTop, "units", 8, "m s^{-1}");
    ERR
    err = NOP (ncid, vertAleTransportTop, "long_name", 69,
               "vertical transport through "
               "the layer interface at the top of the cell");
    ERR
    varids[i++] = vertAleTransportTop;

    /* 1 double (Time, nCells) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;

    err = INQ_VID (ncid, "tendSSH", NC_DOUBLE, 2, dimids, &tendSSH);
    ERR
    err = NOP (ncid, tendSSH, "units", 8, "m s^{-1}");
    ERR
    err = NOP (ncid, tendSSH, "long_name", 35, "time tendency of sea-surface height");
    ERR
    varids[i++] = tendSSH;

    /* 1 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "layerThickness", NC_DOUBLE, 3, dimids, &layerThickness);
    ERR
    err = NOP (ncid, layerThickness, "units", 1, "m");
    ERR
    err = NOP (ncid, layerThickness, "long_name", 15, "layer thickness");
    ERR
    varids[i++] = layerThickness;

    /* 1 double (Time, nEdges, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nEdges;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "normalVelocity", NC_DOUBLE, 3, dimids, &normalVelocity);
    ERR
    err = NOP (ncid, normalVelocity, "units", 8, "m s^{-1}");
    ERR
    err = NOP (ncid, normalVelocity, "long_name", 47,
               "horizonal velocity, "
               "normal component to an edge");
    ERR
    varids[i++] = normalVelocity;

    /* 1 double (Time, nCells) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;

    err = INQ_VID (ncid, "ssh", NC_DOUBLE, 2, dimids, &ssh);
    ERR
    err = NOP (ncid, ssh, "units", 1, "m");
    ERR
    err = NOP (ncid, ssh, "long_name", 18, "sea surface height");
    ERR
    varids[i++] = ssh;

    /* 1 int (nEdges) */
    dimids[0] = dim_nEdges;

    err = INQ_VID (ncid, "maxLevelEdgeTop", NC_INT, 1, dimids, &maxLevelEdgeTop);
    ERR
    err = NOP (ncid, maxLevelEdgeTop, "units", 8, "unitless");
    ERR
    err = NOP (ncid, maxLevelEdgeTop, "long_name", 79,
               "Index to the last edge "
               "in a column with active ocean cells on both sides of it.");
    ERR
    varids[i++] = maxLevelEdgeTop;

    /* 1 double (nVertLevels) */
    dimids[0] = dim_nVertLevels;

    err =
        INQ_VID (ncid, "vertCoordMovementWeights", NC_DOUBLE, 1, dimids, &vertCoordMovementWeights);
    ERR
    err = NOP (ncid, vertCoordMovementWeights, "units", 8, "unitless");
    ERR
    err = NOP (ncid, vertCoordMovementWeights, "long_name", 98,
               "Weights used "
               "for distribution of sea surface heigh purturbations through "
               "multiple vertical levels.");
    ERR
    varids[i++] = vertCoordMovementWeights;

    /* 1 int (nEdges, nVertLevels) */
    dimids[0] = dim_nEdges;
    dimids[1] = dim_nVertLevels;

    err = INQ_VID (ncid, "edgeMask", NC_INT, 2, dimids, &edgeMask);
    ERR
    err = NOP (ncid, edgeMask, "units", 8, "unitless");
    ERR
    err = NOP (ncid, edgeMask, "long_name", 69,
               "Mask on edges that determines "
               "if computations should be done on edge.");
    ERR
    varids[i++] = edgeMask;

    /* 1 int (nCells, nVertLevels) */
    dimids[0] = dim_nCells;
    dimids[1] = dim_nVertLevels;

    err = INQ_VID (ncid, "cellMask", NC_INT, 2, dimids, &cellMask);
    ERR
    err = NOP (ncid, cellMask, "units", 8, "unitless");
    ERR
    err = NOP (ncid, cellMask, "long_name", 69,
               "Mask on cells that determines "
               "if computations should be done on cell.");
    ERR
    varids[i++] = cellMask;

    /* 1 int (nVertices, nVertLevels) */
    dimids[0] = dim_nVertices;
    dimids[1] = dim_nVertLevels;

    err = INQ_VID (ncid, "vertexMask", NC_INT, 2, dimids, &vertexMask);
    ERR
    err = NOP (ncid, vertexMask, "units", 8, "unitless");
    ERR
    err = NOP (ncid, vertexMask, "long_name", 75,
               "Mask on vertices that determines "
               "if computations should be done on vertice.");
    ERR
    varids[i++] = vertexMask;

    /* 2 double (nVertLevels) */
    dimids[0] = dim_nVertLevels;

    err = INQ_VID (ncid, "refZMid", NC_DOUBLE, 1, dimids, &refZMid);
    ERR
    err = NOP (ncid, refZMid, "units", 1, "m");
    ERR
    err = NOP (ncid, refZMid, "long_name", 87,
               "Reference mid z-coordinate of ocean "
               "for each vertical level. This has a negative value.");
    ERR
    varids[i++] = refZMid;

    err = INQ_VID (ncid, "refLayerThickness", NC_DOUBLE, 1, dimids, &refLayerThickness);
    ERR
    err = NOP (ncid, refLayerThickness, "units", 1, "m");
    ERR
    err = NOP (ncid, refLayerThickness, "long_name", 58,
               "Reference layerThickness "
               "of ocean for each vertical level.");
    ERR
    varids[i++] = refLayerThickness;

    /* 1 char (Time, StrLen) */
    dimids[0] = dim_Time;
    dimids[1] = dim_StrLen;

    err = INQ_VID (ncid, "xtime", NC_CHAR, 2, dimids, &xtime);
    ERR
    err = NOP (ncid, xtime, "units", 8, "unitless");
    ERR
    err = NOP (ncid, xtime, "long_name", 45, "model time, with format \'YYYY-MM-DD_HH:MM:SS\'");
    ERR
    varids[i++] = xtime;

    /* 2 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "kineticEnergyCell", NC_DOUBLE, 3, dimids, &kineticEnergyCell);
    ERR
    err = NOP (ncid, kineticEnergyCell, "units", 10, "m^2 s^{-2}");
    ERR
    err = NOP (ncid, kineticEnergyCell, "long_name", 45,
               "kinetic energy of horizonal "
               "velocity on cells");
    ERR
    varids[i++] = kineticEnergyCell;

    err = INQ_VID (ncid, "relativeVorticityCell", NC_DOUBLE, 3, dimids, &relativeVorticityCell);
    ERR
    err = NOP (ncid, relativeVorticityCell, "units", 6, "s^{-1}");
    ERR
    err = NOP (ncid, relativeVorticityCell, "long_name", 67,
               "curl of horizontal velocity, "
               "averaged from vertices to cell centers");
    ERR
    varids[i++] = relativeVorticityCell;

    /* 1 double (Time, nVertices, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nVertices;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "relativeVorticity", NC_DOUBLE, 3, dimids, &relativeVorticity);
    ERR
    err = NOP (ncid, relativeVorticity, "units", 6, "s^{-1}");
    ERR
    err = NOP (ncid, relativeVorticity, "long_name", 48,
               "curl of horizontal velocity, "
               "defined at vertices");
    ERR
    varids[i++] = relativeVorticity;

    /* 1 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "divergence", NC_DOUBLE, 3, dimids, &divergence);
    ERR
    err = NOP (ncid, divergence, "units", 6, "s^{-1}");
    ERR
    err = NOP (ncid, divergence, "long_name", 32, "divergence of horizonal velocity");
    ERR
    varids[i++] = divergence;

    /* 6 double (Time) */
    dimids[0] = dim_Time;

    err = INQ_VID (ncid, "areaCellGlobal", NC_DOUBLE, 1, dimids, &areaCellGlobal);
    ERR
    err = NOP (ncid, areaCellGlobal, "units", 3, "m^2");
    ERR
    err = NOP (ncid, areaCellGlobal, "long_name", 86,
               "sum of the areaCell variable over "
               "the full domain, used to normalize global statistics");
    ERR
    varids[i++] = areaCellGlobal;

    err = INQ_VID (ncid, "areaEdgeGlobal", NC_DOUBLE, 1, dimids, &areaEdgeGlobal);
    ERR
    err = NOP (ncid, areaEdgeGlobal, "units", 3, "m^2");
    ERR
    err = NOP (ncid, areaEdgeGlobal, "long_name", 86,
               "sum of the areaEdge variable over "
               "the full domain, used to normalize global statistics");
    ERR
    varids[i++] = areaEdgeGlobal;

    err = INQ_VID (ncid, "areaTriangleGlobal", NC_DOUBLE, 1, dimids, &areaTriangleGlobal);
    ERR
    err = NOP (ncid, areaTriangleGlobal, "units", 3, "m^2");
    ERR
    err = NOP (ncid, areaTriangleGlobal, "long_name", 90,
               "sum of the areaTriangle variable "
               "over the full domain, used to normalize global statistics");
    ERR
    varids[i++] = areaTriangleGlobal;

    err = INQ_VID (ncid, "volumeCellGlobal", NC_DOUBLE, 1, dimids, &volumeCellGlobal);
    ERR
    err = NOP (ncid, volumeCellGlobal, "units", 3, "m^3");
    ERR
    err = NOP (ncid, volumeCellGlobal, "long_name", 88,
               "sum of the volumeCell variable over "
               "the full domain, used to normalize global statistics");
    ERR
    varids[i++] = volumeCellGlobal;

    err = INQ_VID (ncid, "volumeEdgeGlobal", NC_DOUBLE, 1, dimids, &volumeEdgeGlobal);
    ERR
    err = NOP (ncid, volumeEdgeGlobal, "units", 3, "m^3");
    ERR
    err = NOP (ncid, volumeEdgeGlobal, "long_name", 88,
               "sum of the volumeEdge variable over "
               "the full domain, used to normalize global statistics");
    ERR
    varids[i++] = volumeEdgeGlobal;

    err = INQ_VID (ncid, "CFLNumberGlobal", NC_DOUBLE, 1, dimids, &CFLNumberGlobal);
    ERR
    err = NOP (ncid, CFLNumberGlobal, "units", 8, "unitless");
    ERR
    err = NOP (ncid, CFLNumberGlobal, "long_name", 39, "maximum CFL number over the full domain");
    ERR
    varids[i++] = CFLNumberGlobal;

    /* 1 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "BruntVaisalaFreqTop", NC_DOUBLE, 3, dimids, &BruntVaisalaFreqTop);
    ERR
    err = NOP (ncid, BruntVaisalaFreqTop, "units", 6, "s^{-2}");
    ERR
    err = NOP (ncid, BruntVaisalaFreqTop, "long_name", 89,
               "Brunt Vaisala frequency defined at "
               "the center (horizontally) and top (vertically) of cell");
    ERR
    varids[i++] = BruntVaisalaFreqTop;

    /* 1 double (Time, nCells, nVertLevelsP1) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevelsP1;

    err = INQ_VID (ncid, "vertVelocityTop", NC_DOUBLE, 3, dimids, &vertVelocityTop);
    ERR
    err = NOP (ncid, vertVelocityTop, "units", 8, "m s^{-1}");
    ERR
    err = NOP (ncid, vertVelocityTop, "long_name", 79,
               "vertical velocity defined at center "
               "(horizontally) and top (vertically) of cell");
    ERR
    varids[i++] = vertVelocityTop;

    /* 5 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "velocityZonal", NC_DOUBLE, 3, dimids, &velocityZonal);
    ERR
    err = NOP (ncid, velocityZonal, "units", 8, "m s^{-1}");
    ERR
    err = NOP (ncid, velocityZonal, "long_name", 58,
               "component of horizontal velocity in "
               "the eastward direction");
    ERR
    varids[i++] = velocityZonal;

    err = INQ_VID (ncid, "velocityMeridional", NC_DOUBLE, 3, dimids, &velocityMeridional);
    ERR
    err = NOP (ncid, velocityMeridional, "units", 8, "m s^{-1}");
    ERR
    err = NOP (ncid, velocityMeridional, "long_name", 59,
               "component of horizontal velocity in "
               "the northward direction");
    ERR
    varids[i++] = velocityMeridional;

    err = INQ_VID (ncid, "displacedDensity", NC_DOUBLE, 3, dimids, &displacedDensity);
    ERR
    err = NOP (ncid, displacedDensity, "units", 9, "kg m^{-3}");
    ERR
    err = NOP (ncid, displacedDensity, "long_name", 130,
               "Density displaced adiabatically to "
               "the mid-depth one layer deeper.  That is, layer k has been displaced to the "
               "depth of layer k+1.");
    ERR
    varids[i++] = displacedDensity;

    err = INQ_VID (ncid, "potentialDensity", NC_DOUBLE, 3, dimids, &potentialDensity);
    ERR
    err = NOP (ncid, potentialDensity, "units", 9, "kg m^{-3}");
    ERR
    err = NOP (ncid, potentialDensity, "long_name", 80,
               "potential density: density displaced "
               "adiabatically to the mid-depth of top layer");
    ERR
    varids[i++] = potentialDensity;

    err = INQ_VID (ncid, "pressure", NC_DOUBLE, 3, dimids, &pressure);
    ERR
    err = NOP (ncid, pressure, "units", 8, "N m^{-2}");
    ERR
    err = NOP (ncid, pressure, "long_name", 38, "pressure used in the momentum equation");
    ERR
    varids[i++] = pressure;

    /* 1 double (nVertLevels) */
    dimids[0] = dim_nVertLevels;

    err = INQ_VID (ncid, "refBottomDepth", NC_DOUBLE, 1, dimids, &refBottomDepth);
    ERR
    err = NOP (ncid, refBottomDepth, "units", 1, "m");
    ERR
    err = NOP (ncid, refBottomDepth, "long_name", 78,
               "Reference depth of ocean for each "
               "vertical level. Used in \'z-level\' type runs.");
    ERR
    varids[i++] = refBottomDepth;

    /* 1 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "zMid", NC_DOUBLE, 3, dimids, &zMid);
    ERR
    err = NOP (ncid, zMid, "units", 1, "m");
    ERR
    err = NOP (ncid, zMid, "long_name", 42, "z-coordinate of the mid-depth of the layer");
    ERR
    varids[i++] = zMid;

    /* 1 double (nCells) */
    dimids[0] = dim_nCells;

    err = INQ_VID (ncid, "bottomDepth", NC_DOUBLE, 1, dimids, &bottomDepth);
    ERR
    err = NOP (ncid, bottomDepth, "units", 1, "m");
    ERR
    err = NOP (ncid, bottomDepth, "long_name", 78,
               "Depth of the bottom of the ocean. Given "
               "as a positive distance from sea level.");
    ERR
    varids[i++] = bottomDepth;

    /* 1 int (nCells) */
    dimids[0] = dim_nCells;

    err = INQ_VID (ncid, "maxLevelCell", NC_INT, 1, dimids, &maxLevelCell);
    ERR
    err = NOP (ncid, maxLevelCell, "units", 8, "unitless");
    ERR
    err = NOP (ncid, maxLevelCell, "long_name", 51,
               "Index to the last active ocean cell in each column.");
    ERR
    varids[i++] = maxLevelCell;

    /* 1 int (nEdges) */
    dimids[0] = dim_nEdges;

    err = INQ_VID (ncid, "maxLevelEdgeBot", NC_INT, 1, dimids, &maxLevelEdgeBot);
    ERR
    err = NOP (ncid, maxLevelEdgeBot, "units", 8, "unitless");
    ERR
    err = NOP (ncid, maxLevelEdgeBot, "long_name", 92,
               "Index to the last edge in a column with at "
               "least one active ocean cell on either side of it.");
    ERR
    varids[i++] = maxLevelEdgeBot;

    /* 1 double (Time, nCells) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;

    err = INQ_VID (ncid, "columnIntegratedSpeed", NC_DOUBLE, 2, dimids, &columnIntegratedSpeed);
    ERR
    err = NOP (ncid, columnIntegratedSpeed, "units", 10, "m^2 s^{-1}");
    ERR
    err = NOP (ncid, columnIntegratedSpeed, "long_name", 109,
               "speed = sum(h*sqrt(2*ke)), where ke "
               "is kineticEnergyCell and the sum is over the full column at cell centers.");
    ERR
    varids[i++] = columnIntegratedSpeed;

    /* 13 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "temperatureHorizontalAdvectionTendency", NC_DOUBLE, 3, dimids,
                   &temperatureHorizontalAdvectionTendency);
    ERR
    err = NOP (ncid, temperatureHorizontalAdvectionTendency, "long_name", 58,
               "potential temperature "
               "tendency due to horizontal advection");
    ERR
    err = NOP (ncid, temperatureHorizontalAdvectionTendency, "units", 26,
               "degrees Celsius per second");
    ERR
    varids[i++] = temperatureHorizontalAdvectionTendency;

    err = INQ_VID (ncid, "salinityHorizontalAdvectionTendency", NC_DOUBLE, 3, dimids,
                   &salinityHorizontalAdvectionTendency);
    ERR
    err = NOP (ncid, salinityHorizontalAdvectionTendency, "long_name", 45,
               "salinity tendency due "
               "to horizontal advection");
    ERR
    err = NOP (ncid, salinityHorizontalAdvectionTendency, "units", 14, "PSU per second");
    ERR
    varids[i++] = salinityHorizontalAdvectionTendency;

    err = INQ_VID (ncid, "temperatureVerticalAdvectionTendency", NC_DOUBLE, 3, dimids,
                   &temperatureVerticalAdvectionTendency);
    ERR
    err = NOP (ncid, temperatureVerticalAdvectionTendency, "long_name", 56,
               "potential temperature "
               "tendency due to vertical advection");
    ERR
    err =
        NOP (ncid, temperatureVerticalAdvectionTendency, "units", 26, "degrees Celsius per second");
    ERR
    varids[i++] = temperatureVerticalAdvectionTendency;

    err = INQ_VID (ncid, "salinityVerticalAdvectionTendency", NC_DOUBLE, 3, dimids,
                   &salinityVerticalAdvectionTendency);
    ERR
    err = NOP (ncid, salinityVerticalAdvectionTendency, "long_name", 43,
               "salinity tendency due "
               "to vertical advection");
    ERR
    err = NOP (ncid, salinityVerticalAdvectionTendency, "units", 14, "PSU per second");
    ERR
    varids[i++] = salinityVerticalAdvectionTendency;

    err = INQ_VID (ncid, "temperatureVertMixTendency", NC_DOUBLE, 3, dimids,
                   &temperatureVertMixTendency);
    ERR
    err = NOP (ncid, temperatureVertMixTendency, "long_name", 53,
               "potential temperature tendency "
               "due to vertical mixing");
    ERR
    err = NOP (ncid, temperatureVertMixTendency, "units", 26, "degrees Celsius per second");
    ERR
    varids[i++] = temperatureVertMixTendency;

    err = INQ_VID (ncid, "salinityVertMixTendency", NC_DOUBLE, 3, dimids, &salinityVertMixTendency);
    ERR
    err = NOP (ncid, salinityVertMixTendency, "long_name", 40,
               "salinity tendency due to vertical mixing");
    ERR
    err = NOP (ncid, salinityVertMixTendency, "units", 14, "PSU per second");
    ERR
    varids[i++] = salinityVertMixTendency;

    err = INQ_VID (ncid, "temperatureSurfaceFluxTendency", NC_DOUBLE, 3, dimids,
                   &temperatureSurfaceFluxTendency);
    ERR
    err = NOP (ncid, temperatureSurfaceFluxTendency, "long_name", 52,
               "potential temperature tendency "
               "due to surface fluxes");
    ERR
    err = NOP (ncid, temperatureSurfaceFluxTendency, "units", 26, "degrees Celsius per second");
    ERR
    varids[i++] = temperatureSurfaceFluxTendency;

    err = INQ_VID (ncid, "salinitySurfaceFluxTendency", NC_DOUBLE, 3, dimids,
                   &salinitySurfaceFluxTendency);
    ERR
    err = NOP (ncid, salinitySurfaceFluxTendency, "long_name", 39,
               "salinity tendency due to surface fluxes");
    ERR
    err = NOP (ncid, salinitySurfaceFluxTendency, "units", 14, "PSU per second");
    ERR
    varids[i++] = salinitySurfaceFluxTendency;

    err = INQ_VID (ncid, "temperatureShortWaveTendency", NC_DOUBLE, 3, dimids,
                   &temperatureShortWaveTendency);
    ERR
    err = NOP (ncid, temperatureShortWaveTendency, "units", 26, "degrees Celsius per second");
    ERR
    err = NOP (ncid, temperatureShortWaveTendency, "long_name", 59,
               "potential temperature tendency due "
               "to penetrating shortwave");
    ERR
    varids[i++] = temperatureShortWaveTendency;

    err = INQ_VID (ncid, "temperatureNonLocalTendency", NC_DOUBLE, 3, dimids,
                   &temperatureNonLocalTendency);
    ERR
    err = NOP (ncid, temperatureNonLocalTendency, "long_name", 56,
               "potential temperature tendency due "
               "to kpp non-local flux");
    ERR
    err = NOP (ncid, temperatureNonLocalTendency, "units", 26, "degrees Celsius per second");
    ERR
    varids[i++] = temperatureNonLocalTendency;

    err =
        INQ_VID (ncid, "salinityNonLocalTendency", NC_DOUBLE, 3, dimids, &salinityNonLocalTendency);
    ERR
    err = NOP (ncid, salinityNonLocalTendency, "long_name", 43,
               "salinity tendency due to kpp non-local flux");
    ERR
    err = NOP (ncid, salinityNonLocalTendency, "units", 14, "PSU per second");
    ERR
    varids[i++] = salinityNonLocalTendency;

    err = INQ_VID (ncid, "temperature", NC_DOUBLE, 3, dimids, &temperature);
    ERR
    err = NOP (ncid, temperature, "long_name", 21, "potential temperature");
    ERR
    err = NOP (ncid, temperature, "units", 15, "degrees Celsius");
    ERR
    varids[i++] = temperature;

    err = INQ_VID (ncid, "salinity", NC_DOUBLE, 3, dimids, &salinity);
    ERR
    err = NOP (ncid, salinity, "long_name", 8, "salinity");
    ERR
    err = NOP (ncid, salinity, "units", 32, "grams salt per kilogram seawater");
    ERR
    varids[i++] = salinity;

    assert (i == nvars);

fn_exit:
    return nerrs;
}
