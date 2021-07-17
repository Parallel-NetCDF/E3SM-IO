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
//
#include <assert.h>
#include <string.h>
//
#include <mpi.h>
//
#include <e3sm_io.h>
#include <e3sm_io_err.h>

#include <e3sm_io_case_scorpio.hpp>
#include <e3sm_io_driver.hpp>

#define CHECK_VAR_ERR(varid)                                                                  \
    {                                                                                         \
        if (err != 0) {                                                                       \
            char var_name[64];                                                                \
            driver.inq_var_name (ncid, (varid).data, var_name);                               \
            printf ("Error in %s:%d: %s() var %s\n", __FILE__, __LINE__, __func__, var_name); \
            goto err_out;                                                                     \
        }                                                                                     \
    }
#define DEF_DIM(name, num, dimid)                                                  \
    {                                                                              \
        err = e3sm_io_scorpio_define_dim (driver, ncid, name, num, dnames, dimid); \
        CHECK_ERR                                                                  \
    }
#define DEF_VAR(name, type, ndims, dimids)                                                      \
    {                                                                                           \
        varid++;                                                                                \
        err = e3sm_io_scorpio_define_var (driver, cfg, dnames, decom, decomids[varid - varids], \
                                          ncid, name, type, ndims, dimids, varid);              \
        if (err != 0) {                                                                         \
            printf ("Error in %s line %d: def_var %s\n", __FILE__, __LINE__, name);             \
            goto err_out;                                                                       \
        }                                                                                       \
    }
#define PUT_GATTR_TXT(name, buf)                                                          \
    {                                                                                     \
        err = e3sm_io_scorpio_put_att (driver, ncid, E3SM_IO_GLOBAL_ATTR, name, MPI_CHAR, \
                                       strlen (buf), (void *)(buf));                      \
        CHECK_ERR                                                                         \
    }
#define PUT_GATTR_INT(name, val)                                                                   \
    {                                                                                              \
        int buf = val;                                                                             \
        err = e3sm_io_scorpio_put_att (driver, ncid, E3SM_IO_GLOBAL_ATTR, name, MPI_INT, 1, &buf); \
        CHECK_ERR                                                                                  \
    }
#define PUT_GATTR_DBL(name, val)                                                               \
    {                                                                                          \
        double buf = val;                                                                      \
        err = e3sm_io_scorpio_put_att (driver, ncid, E3SM_IO_GLOBAL_ATTR, name, MPI_DOUBLE, 1, \
                                       &buf);                                                  \
        CHECK_ERR                                                                              \
    }
#define PUT_ATTR_TXT(name, buf)                                                            \
    {                                                                                      \
        err = e3sm_io_scorpio_put_att (driver, ncid, *varid, name, MPI_CHAR, strlen (buf), \
                                       (void *)(buf));                                     \
        CHECK_VAR_ERR (*varid)                                                             \
    }
#define PUT_ATTR_INT(name, num, buf)                                                   \
    {                                                                                  \
        err = e3sm_io_scorpio_put_att (driver, ncid, *varid, name, MPI_INT, num, buf); \
        CHECK_VAR_ERR (*varid)                                                         \
    }
#define PUT_ATTR_FLOAT(name, num, buf)                                                   \
    {                                                                                    \
        err = e3sm_io_scorpio_put_att (driver, ncid, *varid, name, MPI_FLOAT, num, buf); \
        CHECK_VAR_ERR (*varid)                                                           \
    }
#define PUT_ATTR_INT64(name, num, buf)                                                       \
    {                                                                                        \
        err = e3sm_io_scorpio_put_att (driver, ncid, *varid, name, MPI_LONG_LONG, num, buf); \
        CHECK_VAR_ERR (*varid)                                                               \
    }
#define PUT_ATTR_DECOMP(D, ndims, dimids)                                                        \
    {                                                                                            \
        if (cfg.strategy == blob) {                                                              \
            err = e3sm_io_scorpio_put_att (driver, ncid, *varid, "decomposition_ID", MPI_INT, 1, \
                                           &D);                                                  \
            CHECK_VAR_ERR (*varid)                                                               \
            err = e3sm_io_scorpio_put_att (driver, ncid, *varid, "global_dimids", MPI_INT,       \
                                           (MPI_Offset)ndims, (void *)(dimids));                 \
            CHECK_VAR_ERR (*varid)                                                               \
        }                                                                                        \
    }

/*----< add_gattrs() >-------------------------------------------------------*/
static int add_gattrs_scorpio (e3sm_io_config &cfg,
                               e3sm_io_decom &decom,
                               e3sm_io_driver &driver,
                               int ncid) {
    int err = 0, nprocs;
    int k;

    // PIO attributes
    k   = 256;
    err = driver.put_att (ncid, E3SM_IO_GLOBAL_ATTR, "/__pio__/fillmode", MPI_INT, 1, &k);
    CHECK_ERR

    /* 742 global attributes: */
    PUT_GATTR_TXT ("title", "MPAS-Ocean output file information")
    PUT_GATTR_TXT ("source", "MPAS Ocean")
    PUT_GATTR_TXT ("source_id", "a79fafdb76")
    PUT_GATTR_TXT ("product", "model-output")
    PUT_GATTR_TXT ("realm", "ocean")
    PUT_GATTR_TXT ("case", "GMPAS-NYF_T62_oRRS18to6v3")
    PUT_GATTR_TXT ("username", "dqwu")
    PUT_GATTR_TXT ("hostname", "cori-knl")
    PUT_GATTR_TXT ("git_version", "a79fafdb76")
    PUT_GATTR_TXT ("history", "created on 06/04/21 12:46:19")
    PUT_GATTR_TXT ("Conventions", "CF-1.7")
    PUT_GATTR_TXT ("institution_id", "E3SM-Project")
    PUT_GATTR_TXT (
        "institution",
        "LLNL (Lawrence Livermore National Laboratory, Livermore, CA 94550, USA); ANL (Argonne "
        "National Laboratory, Argonne, IL 60439, USA); BNL (Brookhaven National Laboratory, Upton, "
        "NY 11973, USA); LANL (Los Alamos National Laboratory, Los Alamos, NM 87545, USA); LBNL "
        "(Lawrence Berkeley National Laboratory, Berkeley, CA 94720, USA); ORNL (Oak Ridge "
        "National Laboratory, Oak Ridge, TN 37831, USA); PNNL (Pacific Northwest National "
        "Laboratory, Richland, WA 99352, USA); SNL (Sandia National Laboratories, Albuquerque,")
    PUT_GATTR_TXT ("contact", "e3sm-data-support@listserv.llnl.gov")
    PUT_GATTR_TXT ("on_a_sphere", "YES")
    PUT_GATTR_DBL ("sphere_radius", 6371229.)
    PUT_GATTR_TXT ("model_name", "mpas")
    PUT_GATTR_TXT ("core_name", "ocean")
    PUT_GATTR_TXT ("parent_id", "rz1tbn0bed\n555vk5hkh9\n0s1lcmezuy\n7z1uysqc5i\nwj661vqvze")
    PUT_GATTR_TXT ("mesh_spec", "0.0")
    PUT_GATTR_TXT ("config_ocean_run_mode", "forward")
    PUT_GATTR_TXT ("config_do_restart", "NO")
    PUT_GATTR_TXT ("config_restart_timestamp_name", "rpointer.ocn")
    PUT_GATTR_TXT ("config_start_time", "0001-01-01_00:00:00")
    PUT_GATTR_TXT ("config_stop_time", "none")
    PUT_GATTR_TXT ("config_run_duration", "0001-00-00_00:00:00")
    PUT_GATTR_TXT ("config_calendar_type", "gregorian_noleap")
    PUT_GATTR_TXT ("config_write_output_on_startup", "NO")
    PUT_GATTR_INT ("config_pio_num_iotasks", 0)
    PUT_GATTR_INT ("config_pio_stride", 1)
    PUT_GATTR_INT ("config_num_halos", 3)
    PUT_GATTR_TXT (
        "config_block_decomp_file_prefix",
        "/global/cfs/cdirs/e3sm/inputdata/ocn/mpas-o/oRRS18to6v3/mpas-o.graph.info.170111.part.")
    PUT_GATTR_INT ("config_number_of_blocks", 0)
    PUT_GATTR_TXT ("config_explicit_proc_decomp", "NO")
    PUT_GATTR_TXT ("config_proc_decomp_file_prefix", "graph.info.part.")
    PUT_GATTR_TXT ("config_dt", "00:06:00")
    PUT_GATTR_TXT ("config_time_integrator", "split_explicit")
    PUT_GATTR_TXT ("config_hmix_scaleWithMesh", "YES")
    PUT_GATTR_DBL ("config_maxMeshDensity", -1.)
    PUT_GATTR_TXT ("config_hmix_use_ref_cell_width", "NO")
    PUT_GATTR_DBL ("config_hmix_ref_cell_width", 30000.)
    PUT_GATTR_DBL ("config_apvm_scale_factor", 0.)
    PUT_GATTR_TXT ("config_use_mom_del2", "NO")
    PUT_GATTR_DBL ("config_mom_del2", 10.)
    PUT_GATTR_TXT ("config_use_tracer_del2", "NO")
    PUT_GATTR_DBL ("config_tracer_del2", 10.)
    PUT_GATTR_TXT ("config_use_mom_del4", "YES")
    PUT_GATTR_DBL ("config_mom_del4", 3200000000.)
    PUT_GATTR_DBL ("config_mom_del4_div_factor", 1.)
    PUT_GATTR_TXT ("config_use_tracer_del4", "NO")
    PUT_GATTR_DBL ("config_tracer_del4", 0.)
    PUT_GATTR_TXT ("config_use_Leith_del2", "NO")
    PUT_GATTR_DBL ("config_Leith_parameter", 1.)
    PUT_GATTR_DBL ("config_Leith_dx", 15000.)
    PUT_GATTR_DBL ("config_Leith_visc2_max", 2500.)
    PUT_GATTR_TXT ("config_eddying_resolution_taper", "ramp")
    PUT_GATTR_DBL ("config_eddying_resolution_ramp_min", 20000.)
    PUT_GATTR_DBL ("config_eddying_resolution_ramp_max", 30000.)
    PUT_GATTR_TXT ("config_use_Redi", "NO")
    PUT_GATTR_TXT ("config_Redi_closure", "constant")
    PUT_GATTR_DBL ("config_Redi_constant_kappa", 400.)
    PUT_GATTR_DBL ("config_Redi_maximum_slope", 0.01)
    PUT_GATTR_TXT ("config_Redi_use_slope_taper", "YES")
    PUT_GATTR_TXT ("config_Redi_use_surface_taper", "YES")
    PUT_GATTR_TXT ("config_Redi_N2_based_taper_enable", "YES")
    PUT_GATTR_DBL ("config_Redi_N2_based_taper_min", 0.1)
    PUT_GATTR_TXT ("config_Redi_N2_based_taper_limit_term1", "YES")
    PUT_GATTR_TXT ("config_use_GM", "NO")
    PUT_GATTR_TXT ("config_GM_closure", "EdenGreatbatch")
    PUT_GATTR_DBL ("config_GM_constant_kappa", 900.)
    PUT_GATTR_DBL ("config_GM_constant_gravWaveSpeed", 0.3)
    PUT_GATTR_DBL ("config_GM_spatially_variable_min_kappa", 300.)
    PUT_GATTR_DBL ("config_GM_spatially_variable_max_kappa", 1800.)
    PUT_GATTR_DBL ("config_GM_spatially_variable_baroclinic_mode", 3.)
    PUT_GATTR_DBL ("config_GM_Visbeck_alpha", 0.13)
    PUT_GATTR_DBL ("config_GM_Visbeck_max_depth", 1000.)
    PUT_GATTR_DBL ("config_GM_EG_riMin", 200.)
    PUT_GATTR_DBL ("config_GM_EG_kappa_factor", 3.)
    PUT_GATTR_DBL ("config_GM_EG_Rossby_factor", 2.)
    PUT_GATTR_DBL ("config_GM_EG_Rhines_factor", 0.3)
    PUT_GATTR_TXT ("config_Rayleigh_friction", "NO")
    PUT_GATTR_DBL ("config_Rayleigh_damping_coeff", 0.)
    PUT_GATTR_TXT ("config_Rayleigh_damping_depth_variable", "NO")
    PUT_GATTR_TXT ("config_Rayleigh_bottom_friction", "NO")
    PUT_GATTR_DBL ("config_Rayleigh_bottom_damping_coeff", 0.0001)
    PUT_GATTR_TXT ("config_use_cvmix", "YES")
    PUT_GATTR_DBL ("config_cvmix_prandtl_number", 1.)
    PUT_GATTR_TXT ("config_cvmix_background_scheme", "constant")
    PUT_GATTR_DBL ("config_cvmix_background_diffusion", 0.)
    PUT_GATTR_DBL ("config_cvmix_background_viscosity", 0.0001)
    PUT_GATTR_DBL ("config_cvmix_BryanLewis_bl1", 8.e-05)
    PUT_GATTR_DBL ("config_cvmix_BryanLewis_bl2", 0.000105)
    PUT_GATTR_DBL ("config_cvmix_BryanLewis_transitionDepth", 2500.)
    PUT_GATTR_DBL ("config_cvmix_BryanLewis_transitionWidth", 222.)
    PUT_GATTR_TXT ("config_use_cvmix_convection", "YES")
    PUT_GATTR_DBL ("config_cvmix_convective_diffusion", 1.)
    PUT_GATTR_DBL ("config_cvmix_convective_viscosity", 1.)
    PUT_GATTR_TXT ("config_cvmix_convective_basedOnBVF", "YES")
    PUT_GATTR_DBL ("config_cvmix_convective_triggerBVF", 0.)
    PUT_GATTR_TXT ("config_use_cvmix_shear", "YES")
    PUT_GATTR_INT ("config_cvmix_num_ri_smooth_loops", 2)
    PUT_GATTR_TXT ("config_cvmix_use_BLD_smoothing", "YES")
    PUT_GATTR_TXT ("config_cvmix_shear_mixing_scheme", "KPP")
    PUT_GATTR_DBL ("config_cvmix_shear_PP_nu_zero", 0.005)
    PUT_GATTR_DBL ("config_cvmix_shear_PP_alpha", 5.)
    PUT_GATTR_DBL ("config_cvmix_shear_PP_exp", 2.)
    PUT_GATTR_DBL ("config_cvmix_shear_KPP_nu_zero", 0.005)
    PUT_GATTR_DBL ("config_cvmix_shear_KPP_Ri_zero", 0.7)
    PUT_GATTR_DBL ("config_cvmix_shear_KPP_exp", 3.)
    PUT_GATTR_TXT ("config_use_cvmix_tidal_mixing", "NO")
    PUT_GATTR_TXT ("config_use_cvmix_double_diffusion", "NO")
    PUT_GATTR_TXT ("config_use_cvmix_kpp", "YES")
    PUT_GATTR_TXT ("config_use_cvmix_fixed_boundary_layer", "NO")
    PUT_GATTR_DBL ("config_cvmix_kpp_boundary_layer_depth", 30.)
    PUT_GATTR_DBL ("config_cvmix_kpp_criticalBulkRichardsonNumber", 0.25)
    PUT_GATTR_TXT ("config_cvmix_kpp_matching", "SimpleShapes")
    PUT_GATTR_TXT ("config_cvmix_kpp_EkmanOBL", "NO")
    PUT_GATTR_TXT ("config_cvmix_kpp_MonObOBL", "NO")
    PUT_GATTR_TXT ("config_cvmix_kpp_interpolationOMLType", "quadratic")
    PUT_GATTR_DBL ("config_cvmix_kpp_surface_layer_extent", 0.1)
    PUT_GATTR_DBL ("config_cvmix_kpp_surface_layer_averaging", 5.)
    PUT_GATTR_DBL ("configure_cvmix_kpp_minimum_OBL_under_sea_ice", 10.)
    PUT_GATTR_DBL ("config_cvmix_kpp_stop_OBL_search", 100.)
    PUT_GATTR_TXT ("config_cvmix_kpp_use_enhanced_diff", "YES")
    PUT_GATTR_TXT ("config_cvmix_kpp_nonlocal_with_implicit_mix", "NO")
    PUT_GATTR_TXT ("config_cvmix_kpp_use_theory_wave", "NO")
    PUT_GATTR_TXT ("config_cvmix_kpp_langmuir_mixing_opt", "NONE")
    PUT_GATTR_TXT ("config_cvmix_kpp_langmuir_entrainment_opt", "NONE")
    PUT_GATTR_TXT ("config_use_gotm", "NO")
    PUT_GATTR_TXT ("config_gotm_namelist_file", "gotmturb.nml")
    PUT_GATTR_DBL ("config_gotm_constant_surface_roughness_length", 0.02)
    PUT_GATTR_DBL ("config_gotm_constant_bottom_roughness_length", 0.0015)
    PUT_GATTR_DBL ("config_gotm_constant_bottom_drag_coeff", 0.001)
    PUT_GATTR_TXT ("config_use_variable_drag", "NO")
    PUT_GATTR_TXT ("config_use_bulk_wind_stress", "YES")
    PUT_GATTR_TXT ("config_use_bulk_thickness_flux", "YES")
    PUT_GATTR_DBL ("config_flux_attenuation_coefficient", 0.001)
    PUT_GATTR_DBL ("config_flux_attenuation_coefficient_runoff", 10.)
    PUT_GATTR_TXT ("config_use_time_varying_atmospheric_forcing", "NO")
    PUT_GATTR_TXT ("config_time_varying_atmospheric_forcing_type", "WINDPRES")
    PUT_GATTR_TXT ("config_time_varying_atmospheric_forcing_start_time", "0001-01-01_00:00:00")
    PUT_GATTR_TXT ("config_time_varying_atmospheric_forcing_reference_time", "0001-01-01_00:00:00")
    PUT_GATTR_TXT ("config_time_varying_atmospheric_forcing_cycle_start", "0001-01-01_00:00:00")
    PUT_GATTR_TXT ("config_time_varying_atmospheric_forcing_cycle_duration", "2-00-00_00:00:00")
    PUT_GATTR_TXT ("config_time_varying_atmospheric_forcing_interval", "01:00:00")
    PUT_GATTR_DBL ("config_time_varying_atmospheric_forcing_ramp", 10.)
    PUT_GATTR_DBL ("config_time_varying_atmospheric_forcing_ramp_delay", 0.)
    PUT_GATTR_TXT ("config_use_time_varying_land_ice_forcing", "NO")
    PUT_GATTR_TXT ("config_time_varying_land_ice_forcing_start_time", "0001-01-01_00:00:00")
    PUT_GATTR_TXT ("config_time_varying_land_ice_forcing_reference_time", "0001-01-01_00:00:00")
    PUT_GATTR_TXT ("config_time_varying_land_ice_forcing_cycle_start", "0001-01-01_00:00:00")
    PUT_GATTR_TXT ("config_time_varying_land_ice_forcing_cycle_duration", "2-00-00_00:00:00")
    PUT_GATTR_TXT ("config_time_varying_land_ice_forcing_interval", "01:00:00")
    PUT_GATTR_DBL ("config_ssh_grad_relax_timescale", 0.)
    PUT_GATTR_TXT ("config_remove_AIS_coupler_runoff", "NO")
    PUT_GATTR_TXT ("config_sw_absorption_type", "jerlov")
    PUT_GATTR_INT ("config_jerlov_water_type", 3)
    PUT_GATTR_DBL ("config_surface_buoyancy_depth", 1.)
    PUT_GATTR_TXT ("config_enable_shortwave_energy_fixer", "YES")
    PUT_GATTR_TXT ("config_use_tidal_forcing", "NO")
    PUT_GATTR_DBL ("config_use_tidal_forcing_tau", 10000.)
    PUT_GATTR_TXT ("config_tidal_forcing_type", "off")
    PUT_GATTR_TXT ("config_tidal_forcing_model", "off")
    PUT_GATTR_DBL ("config_tidal_forcing_monochromatic_amp", 2.)
    PUT_GATTR_DBL ("config_tidal_forcing_monochromatic_period", 0.5)
    PUT_GATTR_DBL ("config_tidal_forcing_monochromatic_phaseLag", 0.)
    PUT_GATTR_DBL ("config_tidal_forcing_monochromatic_baseline", 0.)
    PUT_GATTR_TXT ("config_use_tidal_potential_forcing", "NO")
    PUT_GATTR_TXT ("config_tidal_potential_reference_time", "0001-01-01_00:00:00")
    PUT_GATTR_TXT ("config_use_tidal_potential_forcing_M2", "YES")
    PUT_GATTR_TXT ("config_use_tidal_potential_forcing_S2", "YES")
    PUT_GATTR_TXT ("config_use_tidal_potential_forcing_N2", "YES")
    PUT_GATTR_TXT ("config_use_tidal_potential_forcing_K2", "YES")
    PUT_GATTR_TXT ("config_use_tidal_potential_forcing_K1", "YES")
    PUT_GATTR_TXT ("config_use_tidal_potential_forcing_O1", "YES")
    PUT_GATTR_TXT ("config_use_tidal_potential_forcing_Q1", "YES")
    PUT_GATTR_TXT ("config_use_tidal_potential_forcing_P1", "YES")
    PUT_GATTR_DBL ("config_tidal_potential_ramp", 10.)
    PUT_GATTR_DBL ("config_self_attraction_and_loading_beta", 0.09)
    PUT_GATTR_TXT ("config_use_vegetation_drag", "NO")
    PUT_GATTR_TXT ("config_use_vegetation_manning_equation", "NO")
    PUT_GATTR_DBL ("config_vegetation_drag_coefficient", 1.0900000333786)
    PUT_GATTR_TXT ("config_use_frazil_ice_formation", "YES")
    PUT_GATTR_TXT ("config_frazil_in_open_ocean", "YES")
    PUT_GATTR_TXT ("config_frazil_under_land_ice", "YES")
    PUT_GATTR_DBL ("config_frazil_heat_of_fusion", 333700.)
    PUT_GATTR_DBL ("config_frazil_ice_density", 1000.)
    PUT_GATTR_DBL ("config_frazil_fractional_thickness_limit", 0.1)
    PUT_GATTR_DBL ("config_specific_heat_sea_water", 3996.)
    PUT_GATTR_DBL ("config_frazil_maximum_depth", 100.)
    PUT_GATTR_DBL ("config_frazil_sea_ice_reference_salinity", 4.)
    PUT_GATTR_DBL ("config_frazil_land_ice_reference_salinity", 0.)
    PUT_GATTR_DBL ("config_frazil_maximum_freezing_temperature", 0.)
    PUT_GATTR_TXT ("config_frazil_use_surface_pressure", "NO")
    PUT_GATTR_TXT ("config_land_ice_flux_mode", "pressure_only")
    PUT_GATTR_TXT ("config_land_ice_flux_formulation", "Jenkins")
    PUT_GATTR_TXT ("config_land_ice_flux_useHollandJenkinsAdvDiff", "NO")
    PUT_GATTR_DBL ("config_land_ice_flux_attenuation_coefficient", 10.)
    PUT_GATTR_DBL ("config_land_ice_flux_boundaryLayerThickness", 10.)
    PUT_GATTR_DBL ("config_land_ice_flux_boundaryLayerNeighborWeight", 0.)
    PUT_GATTR_DBL ("config_land_ice_flux_cp_ice", 2009.)
    PUT_GATTR_DBL ("config_land_ice_flux_rho_ice", 918.)
    PUT_GATTR_DBL ("config_land_ice_flux_topDragCoeff", 0.0025)
    PUT_GATTR_DBL ("config_land_ice_flux_ISOMIP_gammaT", 0.0001)
    PUT_GATTR_DBL ("config_land_ice_flux_rms_tidal_velocity", 0.05)
    PUT_GATTR_DBL ("config_land_ice_flux_jenkins_heat_transfer_coefficient", 0.011)
    PUT_GATTR_DBL ("config_land_ice_flux_jenkins_salt_transfer_coefficient", 0.00031)
    PUT_GATTR_TXT ("config_vert_tracer_adv", "stencil")
    PUT_GATTR_INT ("config_vert_tracer_adv_order", 3)
    PUT_GATTR_INT ("config_horiz_tracer_adv_order", 3)
    PUT_GATTR_DBL ("config_coef_3rd_order", 0.25)
    PUT_GATTR_TXT ("config_monotonic", "YES")
    PUT_GATTR_TXT ("config_use_implicit_bottom_drag", "YES")
    PUT_GATTR_DBL ("config_implicit_bottom_drag_coeff", 0.001)
    PUT_GATTR_TXT ("config_use_implicit_bottom_roughness", "NO")
    PUT_GATTR_TXT ("config_use_implicit_bottom_drag_variable", "NO")
    PUT_GATTR_TXT ("config_use_implicit_bottom_drag_variable_mannings", "NO")
    PUT_GATTR_TXT ("config_use_explicit_bottom_drag", "NO")
    PUT_GATTR_DBL ("config_explicit_bottom_drag_coeff", 0.001)
    PUT_GATTR_TXT ("config_use_topographic_wave_drag", "NO")
    PUT_GATTR_DBL ("config_topographic_wave_drag_coeff", 0.0005)
    PUT_GATTR_TXT ("config_use_wetting_drying", "NO")
    PUT_GATTR_TXT ("config_prevent_drying", "NO")
    PUT_GATTR_DBL ("config_drying_min_cell_height", 0.00100000004749745)
    PUT_GATTR_TXT ("config_zero_drying_velocity", "NO")
    PUT_GATTR_TXT ("config_verify_not_dry", "NO")
    PUT_GATTR_TXT ("config_thickness_flux_type", "centered")
    PUT_GATTR_DBL ("config_drying_safety_height", 0.)
    PUT_GATTR_DBL ("config_density0", 1026.)
    PUT_GATTR_TXT ("config_pressure_gradient_type", "Jacobian_from_TS")
    PUT_GATTR_DBL ("config_common_level_weight", 0.5)
    PUT_GATTR_DBL ("config_zonal_ssh_grad", 0.)
    PUT_GATTR_DBL ("config_meridional_ssh_grad", 0.)
    PUT_GATTR_TXT ("config_eos_type", "jm")
    PUT_GATTR_DBL ("config_open_ocean_freezing_temperature_coeff_0", 0.)
    PUT_GATTR_DBL ("config_open_ocean_freezing_temperature_coeff_S", 0.)
    PUT_GATTR_DBL ("config_open_ocean_freezing_temperature_coeff_p", 0.)
    PUT_GATTR_DBL ("config_open_ocean_freezing_temperature_coeff_pS", 0.)
    PUT_GATTR_DBL ("config_open_ocean_freezing_temperature_coeff_mushy_az1_liq", -18.48)
    PUT_GATTR_DBL ("config_land_ice_cavity_freezing_temperature_coeff_0", 0.0622)
    PUT_GATTR_DBL ("config_land_ice_cavity_freezing_temperature_coeff_S", -0.0563)
    PUT_GATTR_DBL ("config_land_ice_cavity_freezing_temperature_coeff_p", -7.43e-08)
    PUT_GATTR_DBL ("config_land_ice_cavity_freezing_temperature_coeff_pS", -1.74e-10)
    PUT_GATTR_DBL ("config_eos_linear_alpha", 0.2)
    PUT_GATTR_DBL ("config_eos_linear_beta", 0.8)
    PUT_GATTR_DBL ("config_eos_linear_Tref", 5.)
    PUT_GATTR_DBL ("config_eos_linear_Sref", 35.)
    PUT_GATTR_DBL ("config_eos_linear_densityref", 1000.)
    PUT_GATTR_DBL ("config_eos_wright_ref_pressure", 0.)
    PUT_GATTR_INT ("config_n_ts_iter", 2)
    PUT_GATTR_INT ("config_n_bcl_iter_beg", 1)
    PUT_GATTR_INT ("config_n_bcl_iter_mid", 2)
    PUT_GATTR_INT ("config_n_bcl_iter_end", 2)
    PUT_GATTR_TXT ("config_btr_dt", "0000_00:00:12")
    PUT_GATTR_INT ("config_n_btr_cor_iter", 2)
    PUT_GATTR_TXT ("config_vel_correction", "YES")
    PUT_GATTR_INT ("config_btr_subcycle_loop_factor", 2)
    PUT_GATTR_DBL ("config_btr_gam1_velWt1", 0.5)
    PUT_GATTR_DBL ("config_btr_gam2_SSHWt1", 1.)
    PUT_GATTR_DBL ("config_btr_gam3_velWt2", 1.)
    PUT_GATTR_TXT ("config_btr_solve_SSH2", "NO")
    PUT_GATTR_TXT ("config_btr_si_preconditioner", "ras")
    PUT_GATTR_DBL ("config_btr_si_tolerance", 1.e-09)
    PUT_GATTR_INT ("config_n_btr_si_outer_iter", 2)
    PUT_GATTR_TXT ("config_btr_si_partition_match_mode", "NO")
    PUT_GATTR_TXT ("config_vert_coord_movement", "uniform_stretching")
    PUT_GATTR_TXT ("config_ALE_thickness_proportionality", "restingThickness_times_weights")
    PUT_GATTR_DBL ("config_vert_taper_weight_depth_1", 250.)
    PUT_GATTR_DBL ("config_vert_taper_weight_depth_2", 500.)
    PUT_GATTR_TXT ("config_use_min_max_thickness", "NO")
    PUT_GATTR_DBL ("config_min_thickness", 1.)
    PUT_GATTR_DBL ("config_max_thickness_factor", 6.)
    PUT_GATTR_TXT ("config_dzdk_positive", "NO")
    PUT_GATTR_TXT ("config_use_freq_filtered_thickness", "NO")
    PUT_GATTR_DBL ("config_thickness_filter_timescale", 5.)
    PUT_GATTR_TXT ("config_use_highFreqThick_restore", "NO")
    PUT_GATTR_DBL ("config_highFreqThick_restore_time", 30.)
    PUT_GATTR_TXT ("config_use_highFreqThick_del2", "NO")
    PUT_GATTR_DBL ("config_highFreqThick_del2", 100.)
    PUT_GATTR_TXT ("config_check_zlevel_consistency", "NO")
    PUT_GATTR_TXT ("config_check_ssh_consistency", "YES")
    PUT_GATTR_TXT ("config_filter_btr_mode", "NO")
    PUT_GATTR_TXT ("config_prescribe_velocity", "NO")
    PUT_GATTR_TXT ("config_prescribe_thickness", "NO")
    PUT_GATTR_TXT ("config_include_KE_vertex", "YES")
    PUT_GATTR_TXT ("config_check_tracer_monotonicity", "NO")
    PUT_GATTR_TXT ("config_compute_active_tracer_budgets", "YES")
    PUT_GATTR_TXT ("config_disable_thick_all_tend", "NO")
    PUT_GATTR_TXT ("config_disable_thick_hadv", "NO")
    PUT_GATTR_TXT ("config_disable_thick_vadv", "NO")
    PUT_GATTR_TXT ("config_disable_thick_sflux", "NO")
    PUT_GATTR_TXT ("config_disable_vel_all_tend", "NO")
    PUT_GATTR_TXT ("config_disable_vel_coriolis", "NO")
    PUT_GATTR_TXT ("config_disable_vel_pgrad", "NO")
    PUT_GATTR_TXT ("config_disable_vel_hmix", "NO")
    PUT_GATTR_TXT ("config_disable_vel_surface_stress", "NO")
    PUT_GATTR_TXT ("config_disable_vel_topographic_wave_drag", "NO")
    PUT_GATTR_TXT ("config_disable_vel_explicit_bottom_drag", "NO")
    PUT_GATTR_TXT ("config_disable_vel_vmix", "NO")
    PUT_GATTR_TXT ("config_disable_vel_vadv", "NO")
    PUT_GATTR_TXT ("config_disable_tr_all_tend", "NO")
    PUT_GATTR_TXT ("config_disable_tr_adv", "NO")
    PUT_GATTR_TXT ("config_disable_tr_hmix", "NO")
    PUT_GATTR_TXT ("config_disable_tr_vmix", "NO")
    PUT_GATTR_TXT ("config_disable_tr_sflux", "NO")
    PUT_GATTR_TXT ("config_disable_tr_nonlocalflux", "NO")
    PUT_GATTR_TXT ("config_disable_redi_k33", "NO")
    PUT_GATTR_TXT ("config_read_nearest_restart", "NO")
    PUT_GATTR_TXT ("config_conduct_tests", "NO")
    PUT_GATTR_TXT ("config_test_tensors", "NO")
    PUT_GATTR_TXT ("config_tensor_test_function", "sph_uCosCos")
    PUT_GATTR_INT ("config_vert_levels", -1)
    PUT_GATTR_TXT ("config_use_activeTracers", "YES")
    PUT_GATTR_TXT ("config_use_activeTracers_surface_bulk_forcing", "YES")
    PUT_GATTR_TXT ("config_use_activeTracers_surface_restoring", "YES")
    PUT_GATTR_TXT ("config_use_activeTracers_interior_restoring", "NO")
    PUT_GATTR_TXT ("config_use_activeTracers_exponential_decay", "NO")
    PUT_GATTR_TXT ("config_use_activeTracers_idealAge_forcing", "NO")
    PUT_GATTR_TXT ("config_use_activeTracers_ttd_forcing", "NO")
    PUT_GATTR_TXT ("config_use_surface_salinity_monthly_restoring", "YES")
    PUT_GATTR_TXT ("config_surface_salinity_monthly_restoring_compute_interval",
                   "0000-00-01_00:00:00")
    PUT_GATTR_DBL ("config_salinity_restoring_constant_piston_velocity", 1.585e-06)
    PUT_GATTR_DBL ("config_salinity_restoring_max_difference", 100.)
    PUT_GATTR_TXT ("config_salinity_restoring_under_sea_ice", "NO")
    PUT_GATTR_TXT ("config_use_debugTracers", "NO")
    PUT_GATTR_TXT ("config_reset_debugTracers_near_surface", "NO")
    PUT_GATTR_INT ("config_reset_debugTracers_top_nLayers", 20)
    PUT_GATTR_TXT ("config_use_debugTracers_surface_bulk_forcing", "NO")
    PUT_GATTR_TXT ("config_use_debugTracers_surface_restoring", "NO")
    PUT_GATTR_TXT ("config_use_debugTracers_interior_restoring", "NO")
    PUT_GATTR_TXT ("config_use_debugTracers_exponential_decay", "NO")
    PUT_GATTR_TXT ("config_use_debugTracers_idealAge_forcing", "NO")
    PUT_GATTR_TXT ("config_use_debugTracers_ttd_forcing", "NO")
    PUT_GATTR_TXT ("config_use_ecosysTracers", "NO")
    PUT_GATTR_TXT ("config_ecosys_atm_co2_option", "constant")
    PUT_GATTR_TXT ("config_ecosys_atm_alt_co2_option", "constant")
    PUT_GATTR_TXT ("config_ecosys_atm_alt_co2_use_eco", "NO")
    PUT_GATTR_DBL ("config_ecosys_atm_co2_constant_value", 379.)
    PUT_GATTR_TXT ("config_use_ecosysTracers_surface_bulk_forcing", "NO")
    PUT_GATTR_TXT ("config_use_ecosysTracers_surface_restoring", "NO")
    PUT_GATTR_TXT ("config_use_ecosysTracers_interior_restoring", "NO")
    PUT_GATTR_TXT ("config_use_ecosysTracers_exponential_decay", "NO")
    PUT_GATTR_TXT ("config_use_ecosysTracers_idealAge_forcing", "NO")
    PUT_GATTR_TXT ("config_use_ecosysTracers_ttd_forcing", "NO")
    PUT_GATTR_TXT ("config_use_ecosysTracers_surface_value", "NO")
    PUT_GATTR_TXT ("config_use_ecosysTracers_sea_ice_coupling", "NO")
    PUT_GATTR_TXT ("config_ecosysTracers_diagnostic_fields_level1", "NO")
    PUT_GATTR_TXT ("config_ecosysTracers_diagnostic_fields_level2", "NO")
    PUT_GATTR_TXT ("config_ecosysTracers_diagnostic_fields_level3", "NO")
    PUT_GATTR_TXT ("config_ecosysTracers_diagnostic_fields_level4", "NO")
    PUT_GATTR_TXT ("config_ecosysTracers_diagnostic_fields_level5", "NO")
    PUT_GATTR_TXT ("config_use_DMSTracers", "NO")
    PUT_GATTR_TXT ("config_use_DMSTracers_surface_bulk_forcing", "NO")
    PUT_GATTR_TXT ("config_use_DMSTracers_surface_restoring", "NO")
    PUT_GATTR_TXT ("config_use_DMSTracers_interior_restoring", "NO")
    PUT_GATTR_TXT ("config_use_DMSTracers_exponential_decay", "NO")
    PUT_GATTR_TXT ("config_use_DMSTracers_idealAge_forcing", "NO")
    PUT_GATTR_TXT ("config_use_DMSTracers_ttd_forcing", "NO")
    PUT_GATTR_TXT ("config_use_DMSTracers_surface_value", "NO")
    PUT_GATTR_TXT ("config_use_DMSTracers_sea_ice_coupling", "NO")
    PUT_GATTR_TXT ("config_use_MacroMoleculesTracers", "NO")
    PUT_GATTR_TXT ("config_use_MacroMoleculesTracers_surface_bulk_forcing", "NO")
    PUT_GATTR_TXT ("config_use_MacroMoleculesTracers_surface_restoring", "NO")
    PUT_GATTR_TXT ("config_use_MacroMoleculesTracers_interior_restoring", "NO")
    PUT_GATTR_TXT ("config_use_MacroMoleculesTracers_exponential_decay", "NO")
    PUT_GATTR_TXT ("config_use_MacroMoleculesTracers_idealAge_forcing", "NO")
    PUT_GATTR_TXT ("config_use_MacroMoleculesTracers_ttd_forcing", "NO")
    PUT_GATTR_TXT ("config_use_MacroMoleculesTracers_surface_value", "NO")
    PUT_GATTR_TXT ("config_use_MacroMoleculesTracers_sea_ice_coupling", "NO")
    PUT_GATTR_TXT ("config_AM_globalStats_enable", "YES")
    PUT_GATTR_TXT ("config_AM_globalStats_compute_interval", "output_interval")
    PUT_GATTR_TXT ("config_AM_globalStats_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_globalStats_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_globalStats_text_file", "NO")
    PUT_GATTR_TXT ("config_AM_globalStats_directory", "analysis_members")
    PUT_GATTR_TXT ("config_AM_globalStats_output_stream", "globalStatsOutput")
    PUT_GATTR_TXT ("config_AM_surfaceAreaWeightedAverages_enable", "YES")
    PUT_GATTR_TXT ("config_AM_surfaceAreaWeightedAverages_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_surfaceAreaWeightedAverages_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_surfaceAreaWeightedAverages_compute_interval", "0000-00-00_01:00:00")
    PUT_GATTR_TXT ("config_AM_surfaceAreaWeightedAverages_output_stream",
                   "surfaceAreaWeightedAveragesOutput")
    PUT_GATTR_TXT ("config_AM_waterMassCensus_enable", "NO")
    PUT_GATTR_TXT ("config_AM_waterMassCensus_compute_interval", "0000-00-00_01:00:00")
    PUT_GATTR_TXT ("config_AM_waterMassCensus_output_stream", "waterMassCensusOutput")
    PUT_GATTR_TXT ("config_AM_waterMassCensus_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_waterMassCensus_write_on_startup", "NO")
    PUT_GATTR_DBL ("config_AM_waterMassCensus_minTemperature", -2.)
    PUT_GATTR_DBL ("config_AM_waterMassCensus_maxTemperature", 30.)
    PUT_GATTR_DBL ("config_AM_waterMassCensus_minSalinity", 32.)
    PUT_GATTR_DBL ("config_AM_waterMassCensus_maxSalinity", 37.)
    PUT_GATTR_TXT ("config_AM_waterMassCensus_compute_predefined_regions", "YES")
    PUT_GATTR_TXT ("config_AM_waterMassCensus_region_group", "")
    PUT_GATTR_TXT ("config_AM_layerVolumeWeightedAverage_enable", "YES")
    PUT_GATTR_TXT ("config_AM_layerVolumeWeightedAverage_compute_interval", "0000-00-00_01:00:00")
    PUT_GATTR_TXT ("config_AM_layerVolumeWeightedAverage_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_layerVolumeWeightedAverage_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_layerVolumeWeightedAverage_output_stream",
                   "layerVolumeWeightedAverageOutput")
    PUT_GATTR_TXT ("config_AM_zonalMean_enable", "NO")
    PUT_GATTR_TXT ("config_AM_zonalMean_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_zonalMean_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_zonalMean_compute_interval", "0000-00-00_01:00:00")
    PUT_GATTR_TXT ("config_AM_zonalMean_output_stream", "zonalMeanOutput")
    PUT_GATTR_INT ("config_AM_zonalMean_num_bins", 180)
    PUT_GATTR_DBL ("config_AM_zonalMean_min_bin", -1.e+34)
    PUT_GATTR_DBL ("config_AM_zonalMean_max_bin", -1.e+34)
    PUT_GATTR_TXT ("config_AM_okuboWeiss_enable", "NO")
    PUT_GATTR_TXT ("config_AM_okuboWeiss_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_okuboWeiss_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_okuboWeiss_compute_interval", "0000-00-00_01:00:00")
    PUT_GATTR_TXT ("config_AM_okuboWeiss_output_stream", "okuboWeissOutput")
    PUT_GATTR_TXT ("config_AM_okuboWeiss_directory", "analysis_members")
    PUT_GATTR_DBL ("config_AM_okuboWeiss_threshold_value", -0.2)
    PUT_GATTR_DBL ("config_AM_okuboWeiss_normalization", 1.e-10)
    PUT_GATTR_DBL ("config_AM_okuboWeiss_lambda2_normalization", 1.e-10)
    PUT_GATTR_TXT ("config_AM_okuboWeiss_use_lat_lon_coords", "YES")
    PUT_GATTR_TXT ("config_AM_okuboWeiss_compute_eddy_census", "YES")
    PUT_GATTR_INT ("config_AM_okuboWeiss_eddy_min_cells", 20)
    PUT_GATTR_TXT ("config_AM_meridionalHeatTransport_enable", "YES")
    PUT_GATTR_TXT ("config_AM_meridionalHeatTransport_compute_interval", "0000-00-00_01:00:00")
    PUT_GATTR_TXT ("config_AM_meridionalHeatTransport_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_meridionalHeatTransport_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_meridionalHeatTransport_output_stream",
                   "meridionalHeatTransportOutput")
    PUT_GATTR_INT ("config_AM_meridionalHeatTransport_num_bins", 180)
    PUT_GATTR_DBL ("config_AM_meridionalHeatTransport_min_bin", -1.e+34)
    PUT_GATTR_DBL ("config_AM_meridionalHeatTransport_max_bin", -1.e+34)
    PUT_GATTR_TXT ("config_AM_meridionalHeatTransport_region_group", "")
    PUT_GATTR_TXT ("config_AM_testComputeInterval_enable", "NO")
    PUT_GATTR_TXT ("config_AM_testComputeInterval_compute_interval", "00-00-01_00:00:00")
    PUT_GATTR_TXT ("config_AM_testComputeInterval_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_testComputeInterval_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_testComputeInterval_output_stream", "testComputeIntervalOutput")
    PUT_GATTR_TXT ("config_AM_highFrequencyOutput_enable", "YES")
    PUT_GATTR_TXT ("config_AM_highFrequencyOutput_compute_interval", "output_interval")
    PUT_GATTR_TXT ("config_AM_highFrequencyOutput_output_stream", "highFrequencyOutput")
    PUT_GATTR_TXT ("config_AM_highFrequencyOutput_compute_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_highFrequencyOutput_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_timeFilters_enable", "NO")
    PUT_GATTR_TXT ("config_AM_timeFilters_compute_interval", "dt")
    PUT_GATTR_TXT ("config_AM_timeFilters_output_stream", "timeFiltersOutput")
    PUT_GATTR_TXT ("config_AM_timeFilters_restart_stream", "timeFiltersRestart")
    PUT_GATTR_TXT ("config_AM_timeFilters_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_timeFilters_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_timeFilters_initialize_filters", "YES")
    PUT_GATTR_TXT ("config_AM_timeFilters_tau", "90_00:00:00")
    PUT_GATTR_TXT ("config_AM_timeFilters_compute_cell_centered_values", "YES")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_enable", "NO")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_compute_interval", "dt")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_compute_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_output_stream", "lagrPartTrackOutput")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_restart_stream", "lagrPartTrackRestart")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_input_stream", "lagrPartTrackInput")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_write_on_startup", "NO")
    PUT_GATTR_INT ("config_AM_lagrPartTrack_filter_number", 0)
    PUT_GATTR_INT ("config_AM_lagrPartTrack_timeIntegration", 2)
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_reset_criteria", "none")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_reset_global_timestamp", "0000_00:00:00")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_region_stream", "lagrPartTrackRegions")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_reset_if_outside_region", "NO")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_reset_if_inside_region", "NO")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_sample_horizontal_interp", "YES")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_sample_temperature", "YES")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_sample_salinity", "YES")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_sample_DIC", "NO")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_sample_ALK", "NO")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_sample_PO4", "NO")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_sample_NO3", "NO")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_sample_SiO3", "NO")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_sample_NH4", "NO")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_sample_Fe", "NO")
    PUT_GATTR_TXT ("config_AM_lagrPartTrack_sample_O2", "NO")
    PUT_GATTR_TXT ("config_AM_eliassenPalm_enable", "NO")
    PUT_GATTR_TXT ("config_AM_eliassenPalm_compute_interval", "output_interval")
    PUT_GATTR_TXT ("config_AM_eliassenPalm_output_stream", "eliassenPalmOutput")
    PUT_GATTR_TXT ("config_AM_eliassenPalm_restart_stream", "eliassenPalmRestart")
    PUT_GATTR_TXT ("config_AM_eliassenPalm_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_eliassenPalm_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_eliassenPalm_debug", "NO")
    PUT_GATTR_INT ("config_AM_eliassenPalm_nBuoyancyLayers", 45)
    PUT_GATTR_DBL ("config_AM_eliassenPalm_rhomin_buoycoor", 900.)
    PUT_GATTR_DBL ("config_AM_eliassenPalm_rhomax_buoycoor", 1080.)
    PUT_GATTR_TXT ("config_AM_mixedLayerDepths_enable", "YES")
    PUT_GATTR_TXT ("config_AM_mixedLayerDepths_compute_interval", "dt")
    PUT_GATTR_TXT ("config_AM_mixedLayerDepths_output_stream", "mixedLayerDepthsOutput")
    PUT_GATTR_TXT ("config_AM_mixedLayerDepths_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_mixedLayerDepths_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_mixedLayerDepths_Tthreshold", "YES")
    PUT_GATTR_TXT ("config_AM_mixedLayerDepths_Dthreshold", "YES")
    PUT_GATTR_DBL ("config_AM_mixedLayerDepths_crit_temp_threshold", 0.2)
    PUT_GATTR_DBL ("config_AM_mixedLayerDepths_crit_dens_threshold", 0.03)
    PUT_GATTR_DBL ("config_AM_mixedLayerDepths_reference_pressure", 100000.)
    PUT_GATTR_TXT ("config_AM_mixedLayerDepths_Tgradient", "NO")
    PUT_GATTR_TXT ("config_AM_mixedLayerDepths_Dgradient", "NO")
    PUT_GATTR_DBL ("config_AM_mixedLayerDepths_temp_gradient_threshold", 5.e-07)
    PUT_GATTR_DBL ("config_AM_mixedLayerDepths_den_gradient_threshold", 5.e-08)
    PUT_GATTR_INT ("config_AM_mixedLayerDepths_interp_method", 1)
    PUT_GATTR_TXT ("config_AM_regionalStatsDaily_enable", "NO")
    PUT_GATTR_TXT ("config_AM_regionalStatsDaily_compute_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_regionalStatsDaily_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_regionalStatsDaily_compute_interval", "output_interval")
    PUT_GATTR_TXT ("config_AM_regionalStatsDaily_output_stream", "regionalStatsDailyOutput")
    PUT_GATTR_TXT ("config_AM_regionalStatsDaily_restart_stream", "regionalMasksInput")
    PUT_GATTR_TXT ("config_AM_regionalStatsDaily_input_stream", "regionalMasksInput")
    PUT_GATTR_TXT ("config_AM_regionalStatsDaily_operation", "avg")
    PUT_GATTR_TXT ("config_AM_regionalStatsDaily_region_type", "cell")
    PUT_GATTR_TXT ("config_AM_regionalStatsDaily_region_group", "all")
    PUT_GATTR_TXT ("config_AM_regionalStatsDaily_1d_weighting_function", "mul")
    PUT_GATTR_TXT ("config_AM_regionalStatsDaily_2d_weighting_function", "mul")
    PUT_GATTR_TXT ("config_AM_regionalStatsDaily_1d_weighting_field", "areaCell")
    PUT_GATTR_TXT ("config_AM_regionalStatsDaily_2d_weighting_field", "volumeCell")
    PUT_GATTR_TXT ("config_AM_regionalStatsDaily_vertical_mask", "cellMask")
    PUT_GATTR_TXT ("config_AM_regionalStatsDaily_vertical_dimension", "nVertLevels")
    PUT_GATTR_TXT ("config_AM_regionalStatsWeekly_enable", "NO")
    PUT_GATTR_TXT ("config_AM_regionalStatsWeekly_compute_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_regionalStatsWeekly_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_regionalStatsWeekly_compute_interval", "output_interval")
    PUT_GATTR_TXT ("config_AM_regionalStatsWeekly_output_stream", "regionalStatsWeeklyOutput")
    PUT_GATTR_TXT ("config_AM_regionalStatsWeekly_restart_stream", "regionalMasksInput")
    PUT_GATTR_TXT ("config_AM_regionalStatsWeekly_input_stream", "regionalMasksInput")
    PUT_GATTR_TXT ("config_AM_regionalStatsWeekly_operation", "avg")
    PUT_GATTR_TXT ("config_AM_regionalStatsWeekly_region_type", "cell")
    PUT_GATTR_TXT ("config_AM_regionalStatsWeekly_region_group", "all")
    PUT_GATTR_TXT ("config_AM_regionalStatsWeekly_1d_weighting_function", "mul")
    PUT_GATTR_TXT ("config_AM_regionalStatsWeekly_2d_weighting_function", "mul")
    PUT_GATTR_TXT ("config_AM_regionalStatsWeekly_1d_weighting_field", "areaCell")
    PUT_GATTR_TXT ("config_AM_regionalStatsWeekly_2d_weighting_field", "volumeCell")
    PUT_GATTR_TXT ("config_AM_regionalStatsWeekly_vertical_mask", "cellMask")
    PUT_GATTR_TXT ("config_AM_regionalStatsWeekly_vertical_dimension", "nVertLevels")
    PUT_GATTR_TXT ("config_AM_regionalStatsMonthly_enable", "NO")
    PUT_GATTR_TXT ("config_AM_regionalStatsMonthly_compute_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_regionalStatsMonthly_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_regionalStatsMonthly_compute_interval", "output_interval")
    PUT_GATTR_TXT ("config_AM_regionalStatsMonthly_output_stream", "regionalStatsMonthlyOutput")
    PUT_GATTR_TXT ("config_AM_regionalStatsMonthly_restart_stream", "regionalMasksInput")
    PUT_GATTR_TXT ("config_AM_regionalStatsMonthly_input_stream", "regionalMasksInput")
    PUT_GATTR_TXT ("config_AM_regionalStatsMonthly_operation", "avg")
    PUT_GATTR_TXT ("config_AM_regionalStatsMonthly_region_type", "cell")
    PUT_GATTR_TXT ("config_AM_regionalStatsMonthly_region_group", "all")
    PUT_GATTR_TXT ("config_AM_regionalStatsMonthly_1d_weighting_function", "mul")
    PUT_GATTR_TXT ("config_AM_regionalStatsMonthly_2d_weighting_function", "mul")
    PUT_GATTR_TXT ("config_AM_regionalStatsMonthly_1d_weighting_field", "areaCell")
    PUT_GATTR_TXT ("config_AM_regionalStatsMonthly_2d_weighting_field", "volumeCell")
    PUT_GATTR_TXT ("config_AM_regionalStatsMonthly_vertical_mask", "cellMask")
    PUT_GATTR_TXT ("config_AM_regionalStatsMonthly_vertical_dimension", "nVertLevels")
    PUT_GATTR_TXT ("config_AM_regionalStatsCustom_enable", "NO")
    PUT_GATTR_TXT ("config_AM_regionalStatsCustom_compute_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_regionalStatsCustom_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_regionalStatsCustom_compute_interval", "output_interval")
    PUT_GATTR_TXT ("config_AM_regionalStatsCustom_output_stream", "regionalStatsCustomOutput")
    PUT_GATTR_TXT ("config_AM_regionalStatsCustom_restart_stream", "regionalMasksInput")
    PUT_GATTR_TXT ("config_AM_regionalStatsCustom_input_stream", "regionalMasksInput")
    PUT_GATTR_TXT ("config_AM_regionalStatsCustom_operation", "avg")
    PUT_GATTR_TXT ("config_AM_regionalStatsCustom_region_type", "cell")
    PUT_GATTR_TXT ("config_AM_regionalStatsCustom_region_group", "all")
    PUT_GATTR_TXT ("config_AM_regionalStatsCustom_1d_weighting_function", "mul")
    PUT_GATTR_TXT ("config_AM_regionalStatsCustom_2d_weighting_function", "mul")
    PUT_GATTR_TXT ("config_AM_regionalStatsCustom_1d_weighting_field", "areaCell")
    PUT_GATTR_TXT ("config_AM_regionalStatsCustom_2d_weighting_field", "volumeCell")
    PUT_GATTR_TXT ("config_AM_regionalStatsCustom_vertical_mask", "cellMask")
    PUT_GATTR_TXT ("config_AM_regionalStatsCustom_vertical_dimension", "nVertLevels")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsDaily_enable", "NO")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsDaily_compute_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsDaily_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsDaily_compute_interval", "00-00-00_01:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsDaily_output_stream", "timeSeriesStatsDailyOutput")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsDaily_restart_stream", "none")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsDaily_operation", "avg")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsDaily_reference_times", "initial_time")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsDaily_duration_intervals", "repeat_interval")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsDaily_repeat_intervals", "reset_interval")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsDaily_reset_intervals", "00-00-01_00:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsDaily_backward_output_offset", "00-00-01_00:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthly_enable", "YES")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthly_compute_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthly_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthly_compute_interval", "dt")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthly_output_stream", "timeSeriesStatsMonthlyOutput")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthly_restart_stream", "none")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthly_operation", "avg")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthly_reference_times", "initial_time")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthly_duration_intervals", "repeat_interval")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthly_repeat_intervals", "reset_interval")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthly_reset_intervals", "00-01-00_00:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthly_backward_output_offset", "00-01-00_00:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsClimatology_enable", "NO")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsClimatology_compute_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsClimatology_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsClimatology_compute_interval", "00-00-00_01:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsClimatology_output_stream",
                   "timeSeriesStatsClimatologyOutput")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsClimatology_restart_stream",
                   "timeSeriesStatsClimatologyRestart")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsClimatology_operation", "avg")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsClimatology_reference_times",
                   "00-03-01_00:00:00;00-06-01_00:00:00;00-09-01_00:00:00;00-12-01_00:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsClimatology_duration_intervals",
                   "00-03-00_00:00:00;00-03-00_00:00:00;00-03-00_00:00:00;00-03-00_00:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsClimatology_repeat_intervals",
                   "01-00-00_00:00:00;01-00-00_00:00:00;01-00-00_00:00:00;01-00-00_00:00:00")
    PUT_GATTR_TXT (
        "config_AM_timeSeriesStatsClimatology_reset_intervals",
        "1000-00-00_00:00:00;1000-00-00_00:00:00;1000-00-00_00:00:00;1000-00-00_00:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsClimatology_backward_output_offset",
                   "00-03-00_00:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMax_enable", "YES")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMax_compute_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMax_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMax_compute_interval", "00-00-00_01:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMax_output_stream",
                   "timeSeriesStatsMonthlyMaxOutput")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMax_restart_stream", "none")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMax_operation", "max")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMax_reference_times", "initial_time")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMax_duration_intervals", "repeat_interval")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMax_repeat_intervals", "reset_interval")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMax_reset_intervals", "00-01-00_00:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMax_backward_output_offset",
                   "00-01-00_00:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMin_enable", "YES")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMin_compute_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMin_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMin_compute_interval", "00-00-00_01:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMin_output_stream",
                   "timeSeriesStatsMonthlyMinOutput")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMin_restart_stream", "none")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMin_operation", "min")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMin_reference_times", "initial_time")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMin_duration_intervals", "repeat_interval")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMin_repeat_intervals", "reset_interval")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMin_reset_intervals", "00-01-00_00:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsMonthlyMin_backward_output_offset",
                   "00-01-00_00:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsCustom_enable", "NO")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsCustom_compute_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsCustom_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsCustom_compute_interval", "00-00-00_01:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsCustom_output_stream", "timeSeriesStatsCustomOutput")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsCustom_restart_stream", "timeSeriesStatsCustomRestart")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsCustom_operation", "avg")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsCustom_reference_times", "initial_time")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsCustom_duration_intervals", "repeat_interval")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsCustom_repeat_intervals", "reset_interval")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsCustom_reset_intervals", "00-00-07_00:00:00")
    PUT_GATTR_TXT ("config_AM_timeSeriesStatsCustom_backward_output_offset", "00-00-01_00:00:00")
    PUT_GATTR_TXT ("config_AM_pointwiseStats_enable", "NO")
    PUT_GATTR_TXT ("config_AM_pointwiseStats_compute_interval", "output_interval")
    PUT_GATTR_TXT ("config_AM_pointwiseStats_output_stream", "pointwiseStatsOutput")
    PUT_GATTR_TXT ("config_AM_pointwiseStats_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_pointwiseStats_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_debugDiagnostics_enable", "NO")
    PUT_GATTR_TXT ("config_AM_debugDiagnostics_compute_interval", "dt")
    PUT_GATTR_TXT ("config_AM_debugDiagnostics_output_stream", "debugDiagnosticsOutput")
    PUT_GATTR_TXT ("config_AM_debugDiagnostics_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_debugDiagnostics_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_debugDiagnostics_check_state", "YES")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_enable", "NO")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_compute_interval", "0010-00-00_00:00:00")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_output_stream", "none")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_variable_a", "layerThickness")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_variable_b", "areaCell")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_variable_c", "none")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_variable_d", "none")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_variable_e", "none")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_variable_f", "none")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_variable_g", "none")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_variable_h", "none")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_expression_1", "a b *")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_expression_2", "none")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_expression_3", "none")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_expression_4", "none")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_output_name_1", "volumeCell")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_output_name_2", "none")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_output_name_3", "none")
    PUT_GATTR_TXT ("config_AM_rpnCalculator_output_name_4", "none")
    PUT_GATTR_TXT ("config_AM_transectTransport_enable", "NO")
    PUT_GATTR_TXT ("config_AM_transectTransport_compute_interval", "output_interval")
    PUT_GATTR_TXT ("config_AM_transectTransport_output_stream", "transectTransportOutput")
    PUT_GATTR_TXT ("config_AM_transectTransport_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_transectTransport_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_transectTransport_transect_group", "all")
    PUT_GATTR_TXT ("config_AM_eddyProductVariables_enable", "YES")
    PUT_GATTR_TXT ("config_AM_eddyProductVariables_compute_interval", "0000-00-00_01:00:00")
    PUT_GATTR_TXT ("config_AM_eddyProductVariables_output_stream", "eddyProductVariablesOutput")
    PUT_GATTR_TXT ("config_AM_eddyProductVariables_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_eddyProductVariables_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_mocStreamfunction_enable", "YES")
    PUT_GATTR_TXT ("config_AM_mocStreamfunction_compute_interval", "0000-00-00_01:00:00")
    PUT_GATTR_TXT ("config_AM_mocStreamfunction_output_stream", "mocStreamfunctionOutput")
    PUT_GATTR_TXT ("config_AM_mocStreamfunction_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_mocStreamfunction_write_on_startup", "NO")
    PUT_GATTR_DBL ("config_AM_mocStreamfunction_min_bin", -1.e+34)
    PUT_GATTR_DBL ("config_AM_mocStreamfunction_max_bin", -1.e+34)
    PUT_GATTR_INT ("config_AM_mocStreamfunction_num_bins", 180)
    PUT_GATTR_TXT ("config_AM_mocStreamfunction_region_group", "all")
    PUT_GATTR_TXT ("config_AM_mocStreamfunction_transect_group", "all")
    PUT_GATTR_TXT ("config_AM_oceanHeatContent_enable", "YES")
    PUT_GATTR_TXT ("config_AM_oceanHeatContent_compute_interval", "0000-00-00_01:00:00")
    PUT_GATTR_TXT ("config_AM_oceanHeatContent_output_stream", "oceanHeatContentOutput")
    PUT_GATTR_TXT ("config_AM_oceanHeatContent_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_oceanHeatContent_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_mixedLayerHeatBudget_enable", "YES")
    PUT_GATTR_TXT ("config_AM_mixedLayerHeatBudget_compute_interval", "0000-00-00_01:00:00")
    PUT_GATTR_TXT ("config_AM_mixedLayerHeatBudget_output_stream", "mixedLayerHeatBudgetOutput")
    PUT_GATTR_TXT ("config_AM_mixedLayerHeatBudget_compute_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_mixedLayerHeatBudget_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_sedimentFluxIndex_enable", "NO")
    PUT_GATTR_TXT ("config_AM_sedimentFluxIndex_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_sedimentFluxIndex_write_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_sedimentFluxIndex_compute_interval", "output_interval")
    PUT_GATTR_TXT ("config_AM_sedimentFluxIndex_output_stream", "sedimentFluxIndexOutput")
    PUT_GATTR_TXT ("config_AM_sedimentFluxIndex_directory", "analysis_members")
    PUT_GATTR_TXT ("config_AM_sedimentFluxIndex_use_lat_lon_coords", "YES")
    PUT_GATTR_TXT ("config_AM_sedimentTransport_enable", "NO")
    PUT_GATTR_TXT ("config_AM_sedimentTransport_compute_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_sedimentTransport_write_on_startup", "YES")
    PUT_GATTR_TXT ("config_AM_sedimentTransport_compute_interval", "output_interval")
    PUT_GATTR_TXT ("config_AM_sedimentTransport_output_stream", "sedimentTransportOutput")
    PUT_GATTR_TXT ("config_AM_sedimentTransport_directory", "analysis_members")
    PUT_GATTR_DBL ("config_AM_sedimentTransport_grain_size", 0.00025)
    PUT_GATTR_TXT ("config_AM_sedimentTransport_ws_formula", "VanRijn1993")
    PUT_GATTR_TXT ("config_AM_sedimentTransport_bedld_formula", "Soulsby-Damgaard")
    PUT_GATTR_TXT ("config_AM_sedimentTransport_SSC_ref_formula", "Lee2004")
    PUT_GATTR_DBL ("config_AM_sedimentTransport_drag_coefficient", 0.0025)
    PUT_GATTR_DBL ("config_AM_sedimentTransport_erate", 0.0005)
    PUT_GATTR_DBL ("config_AM_sedimentTransport_tau_ce", 0.1)
    PUT_GATTR_DBL ("config_AM_sedimentTransport_tau_cd", 0.1)
    PUT_GATTR_DBL ("config_AM_sedimentTransport_Manning_coef", 0.022)
    PUT_GATTR_DBL ("config_AM_sedimentTransport_grain_porosity", 0.5)
    PUT_GATTR_DBL ("config_AM_sedimentTransport_water_density", 1020.)
    PUT_GATTR_DBL ("config_AM_sedimentTransport_grain_density", 2650.)
    PUT_GATTR_DBL ("config_AM_sedimentTransport_alpha", 0.01)
    PUT_GATTR_DBL ("config_AM_sedimentTransport_kinematic_viscosity", 1.e-06)
    PUT_GATTR_DBL ("config_AM_sedimentTransport_vertical_diffusion_coefficient", 0.01)
    PUT_GATTR_TXT ("config_AM_sedimentTransport_bedload", "YES")
    PUT_GATTR_TXT ("config_AM_sedimentTransport_suspended", "YES")
    PUT_GATTR_TXT ("config_AM_sedimentTransport_use_lat_lon_coords", "YES")
    PUT_GATTR_TXT ("config_AM_harmonicAnalysis_enable", "NO")
    PUT_GATTR_TXT ("config_AM_harmonicAnalysis_compute_interval", "output_interval")
    PUT_GATTR_TXT ("config_AM_harmonicAnalysis_start", "0001-01-01_00:00:00")
    PUT_GATTR_TXT ("config_AM_harmonicAnalysis_end", "0001-01-01_00:00:00")
    PUT_GATTR_TXT ("config_AM_harmonicAnalysis_output_stream", "harmonicAnalysisOutput")
    PUT_GATTR_TXT ("config_AM_harmonicAnalysis_restart_stream", "harmonicAnalysisRestart")
    PUT_GATTR_TXT ("config_AM_harmonicAnalysis_compute_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_harmonicAnalysis_write_on_startup", "NO")
    PUT_GATTR_TXT ("config_AM_harmonicAnalysis_use_M2", "YES")
    PUT_GATTR_TXT ("config_AM_harmonicAnalysis_use_S2", "YES")
    PUT_GATTR_TXT ("config_AM_harmonicAnalysis_use_N2", "YES")
    PUT_GATTR_TXT ("config_AM_harmonicAnalysis_use_K2", "YES")
    PUT_GATTR_TXT ("config_AM_harmonicAnalysis_use_K1", "YES")
    PUT_GATTR_TXT ("config_AM_harmonicAnalysis_use_O1", "YES")
    PUT_GATTR_TXT ("config_AM_harmonicAnalysis_use_Q1", "YES")
    PUT_GATTR_TXT ("config_AM_harmonicAnalysis_use_P1", "YES")
    PUT_GATTR_TXT ("file_id", "6270j0pl9q")

err_out:
    return err;
}

/*----< def_G_case() >-------------------------------------------------------*/
int def_G_case_scorpio (e3sm_io_config &cfg,
                        e3sm_io_decom &decom,
                        e3sm_io_driver &driver,
                        int ncid, /* file ID */
                        std::vector<int> &decomids,
                        e3sm_io_scorpio_var *varids, /* variable IDs */
                        int *scorpiovars) {
    /* Total 52 variables */
    int one = 1, two = 2, three = 3, four = 4, five = 5, six = 6;
    int i, j, k, err, nprocs, nprocs_ID, ndims, dimids[2];
    e3sm_io_scorpio_var *varid;
    int nelems_D[MAX_NUM_DECOMP], max_nreqs_dimid[MAX_NUM_DECOMP];
    int dim_Time, dim_nCells, dim_nEdges, dim_nVertices, dim_StrLen;
    int dim_nVertLevels, dim_nVertLevelsP1, dimids_D[MAX_NUM_DECOMP][3];
    int fix_D[MAX_NUM_DECOMP][3], rec_D[MAX_NUM_DECOMP][3];
    std::map<int, std::string> dnames;

    char name[64];

    err = add_gattrs_scorpio (cfg, decom, driver, ncid);
    CHECK_ERR

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

    /* define dimensions */
    DEF_DIM ("Time", NC_UNLIMITED, &dim_Time)
    DEF_DIM ("nCells", decom.dims[0][0], &dim_nCells)
    DEF_DIM ("nEdges", decom.dims[1][0], &dim_nEdges)
    DEF_DIM ("nVertices", decom.dims[4][0], &dim_nVertices)
    DEF_DIM ("nVertLevelsP1", decom.dims[5][1], &dim_nVertLevelsP1)
    DEF_DIM ("nVertLevels", decom.dims[2][1], &dim_nVertLevels)
    DEF_DIM ("StrLen", 64, &dim_StrLen)

    dimids_D[0][0] = dim_Time;
    dimids_D[1][0] = dim_Time;
    dimids_D[2][0] = dim_Time;
    dimids_D[3][0] = dim_Time;
    dimids_D[4][0] = dim_Time;
    dimids_D[5][0] = dim_Time;

    dimids_D[0][1] = dim_nCells;
    dimids_D[1][1] = dim_nEdges;
    dimids_D[2][1] = dim_nCells;
    dimids_D[2][2] = dim_nVertLevels;
    dimids_D[3][1] = dim_nEdges;
    dimids_D[3][2] = dim_nVertLevels;
    dimids_D[4][1] = dim_nVertices;
    dimids_D[4][2] = dim_nVertLevels;
    dimids_D[5][1] = dim_nCells;
    dimids_D[5][2] = dim_nVertLevelsP1;

    fix_D[0][0] = dim_nCells;
    fix_D[1][0] = dim_nEdges;
    fix_D[2][0] = dim_nCells;
    fix_D[2][1] = dim_nVertLevels;
    fix_D[3][0] = dim_nEdges;
    fix_D[3][1] = dim_nVertLevels;
    fix_D[4][0] = dim_nVertices;
    fix_D[4][1] = dim_nVertLevels;
    fix_D[5][0] = dim_nCells;
    fix_D[5][1] = dim_nVertLevelsP1;

    rec_D[0][0] = dim_Time;
    rec_D[1][0] = dim_Time;
    rec_D[2][0] = dim_Time;
    rec_D[3][0] = dim_Time;
    rec_D[4][0] = dim_Time;
    rec_D[5][0] = dim_Time;
    rec_D[0][1] = dim_nCells;
    rec_D[1][1] = dim_nEdges;
    rec_D[2][1] = dim_nCells;
    rec_D[2][2] = dim_nVertLevels;
    rec_D[3][1] = dim_nEdges;
    rec_D[3][2] = dim_nVertLevels;
    rec_D[4][1] = dim_nVertices;
    rec_D[4][2] = dim_nVertLevels;
    rec_D[5][1] = dim_nCells;
    rec_D[5][2] = dim_nVertLevelsP1;

    varid = varids - 1;

    /* define scorpio decom map variables */
    for (j = 0; j < 6; j++) {
        int piodecomid[] = {0, 1, 2, 3, 4, 5};
        int piodims[5];

        sprintf (name, "/__pio__/decomp/%d", (j + 512));
        err = driver.def_local_var (ncid, name, MPI_LONG_LONG, 1, decom.raw_nreqs + piodecomid[j],
                                    scorpiovars + j);
        CHECK_ERR

        for (i = 0; i < decom.ndims[piodecomid[j]]; i++){
            piodims[i] = (int)decom.dims[piodecomid[j]][i];
        }
        err = driver.put_att (ncid, scorpiovars[j], "dimlen", MPI_INT, decom.ndims[piodecomid[j]],
                              piodims);
        CHECK_ERR

        err =
            driver.put_att (ncid, scorpiovars[j], "ndims", MPI_INT, 1, decom.ndims + piodecomid[j]);
        CHECK_ERR

        k   = 6;
        err = driver.put_att (ncid, scorpiovars[j], "piotype", MPI_INT, 1, &k);
        CHECK_ERR
    }

    // TODO: only the first subfile contain nproc
    err = driver.def_local_var (ncid, "/__pio__/info/nproc", MPI_INT, 0, NULL, scorpiovars + (j++));
    CHECK_ERR

    /* define climate variables */

    /* double salinitySurfaceRestoringTendency(Time, nCells) */
    DEF_VAR ("salinitySurfaceRestoringTendency", MPI_DOUBLE, 2, rec_D[0])
    PUT_ATTR_TXT ("units", "m PSU/s")
    PUT_ATTR_TXT ("long_name", "salinity tendency due to surface restoring")
    PUT_ATTR_DECOMP (one, 2, dimids_D[0])

    /* double vertTransportVelocityTop(Time, nCells, nVertLevelsP1) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("vertTransportVelocityTop", MPI_DOUBLE, ndims, rec_D[5])
    PUT_ATTR_TXT ("units", "m s^{-1}")
    PUT_ATTR_TXT ("long_name",
                  "vertical tracer-transport velocity defined at center (horizontally) and top "
                  "(vertically) of cell.  This is not the vertical ALE transport, but is "
                  "Eulerian (fixed-frame) in the vertical, and computed from the continuity "
                  "equation from the horizontal total tracer-transport velocity.")
    PUT_ATTR_DECOMP (six, 3, dimids_D[5])

    /* double vertGMBolusVelocityTop(Time, nCells, nVertLevelsP1) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("vertGMBolusVelocityTop", MPI_DOUBLE, ndims, rec_D[5])
    PUT_ATTR_TXT ("units", "m s^{-1}")
    PUT_ATTR_TXT ("long_name",
                  "vertical tracer-transport velocity defined at center (horizontally) and top "
                  "(vertically) of cell.  This is not the vertical ALE transport, but is "
                  "Eulerian (fixed-frame) in the vertical, and computed from the continuity "
                  "equation from the horizontal GM Bolus velocity.")
    PUT_ATTR_DECOMP (six, 3, dimids_D[5])

    /* double vertAleTransportTop(Time, nCells, nVertLevelsP1) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("vertAleTransportTop", MPI_DOUBLE, ndims, rec_D[5])
    PUT_ATTR_TXT ("units", "m s^{-1}")
    PUT_ATTR_TXT ("long_name",
                  "vertical transport through the layer interface at the top of the cell")
    PUT_ATTR_DECOMP (six, 3, dimids_D[5])

    /* double tendSSH(Time, nCells) */
    DEF_VAR ("tendSSH", MPI_DOUBLE, 2, rec_D[0])
    PUT_ATTR_TXT ("units", "m s^{-1}")
    PUT_ATTR_TXT ("long_name", "time tendency of sea-surface height")
    PUT_ATTR_DECOMP (one, 2, dimids_D[0])

    /* double layerThickness(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("layerThickness", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("units", "m")
    PUT_ATTR_TXT ("long_name", "layer thickness")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double normalVelocity(Time, nEdges, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("normalVelocity", MPI_DOUBLE, ndims, rec_D[3])
    PUT_ATTR_TXT ("units", "m s^{-1}")
    PUT_ATTR_TXT ("long_name", "horizonal velocity, normal component to an edge")
    PUT_ATTR_DECOMP (four, 3, dimids_D[3])

    /* double ssh(Time, nCells) */
    DEF_VAR ("ssh", MPI_DOUBLE, 2, rec_D[0])
    PUT_ATTR_TXT ("units", "m")
    PUT_ATTR_TXT ("long_name", "sea surface height")
    PUT_ATTR_DECOMP (one, 2, dimids_D[0])

    /* int maxLevelEdgeTop(nEdges) */
    DEF_VAR ("maxLevelEdgeTop", MPI_INT, 1, fix_D[1])
    PUT_ATTR_TXT ("units", "unitless")
    PUT_ATTR_TXT ("long_name",
                  "Index to the last edge in a column with active ocean cells on both sides of it.")
    PUT_ATTR_DECOMP (two, 1, dimids_D[1] + 1)

    /* double vertCoordMovementWeights(nVertLevels) */
    DEF_VAR ("vertCoordMovementWeights", MPI_DOUBLE, 1, &dim_nVertLevels)
    PUT_ATTR_TXT ("units", "unitless")
    PUT_ATTR_TXT ("long_name",
                  "Weights used for distribution of sea surface height perturbations through "
                  "multiple vertical levels.")

    /* int edgeMask(nEdges, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 1 : 2;
    DEF_VAR ("edgeMask", MPI_INT, ndims, fix_D[3])
    PUT_ATTR_TXT ("units", "unitless")
    PUT_ATTR_TXT ("long_name",
                  "Mask on edges that determines if computations should be done on edges.")
    PUT_ATTR_DECOMP (four, 2, dimids_D[3] + 1)

    /* int cellMask(nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 1 : 2;
    DEF_VAR ("cellMask", MPI_INT, ndims, fix_D[2])
    PUT_ATTR_TXT ("units", "unitless")
    PUT_ATTR_TXT ("long_name",
                  "Mask on cells that determines if computations should be done on cells.")
    PUT_ATTR_DECOMP (three, 2, dimids_D[2] + 1)

    /* int vertexMask(nVertices, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 1 : 2;
    DEF_VAR ("vertexMask", MPI_INT, ndims, fix_D[4])
    PUT_ATTR_TXT ("units", "unitless")
    PUT_ATTR_TXT ("long_name",
                  "Mask on vertices that determines if computations should be done on vertices.")
    PUT_ATTR_DECOMP (five, 2, dimids_D[4] + 1)

    /* double refZMid(nVertLevels) */
    DEF_VAR ("refZMid", MPI_DOUBLE, 1, &dim_nVertLevels)
    PUT_ATTR_TXT ("units", "m")
    PUT_ATTR_TXT ("long_name",
                  "Reference mid z-coordinate of ocean for each vertical level. This has a "
                  "negative value.")

    /* double refLayerThickness(nVertLevels) */
    DEF_VAR ("refLayerThickness", MPI_DOUBLE, 1, &dim_nVertLevels)
    PUT_ATTR_TXT ("units", "m")
    PUT_ATTR_TXT ("long_name", "Reference layerThickness of ocean for each vertical level.")

    /* char xtime(Time, StrLen) */
    dimids[0] = dim_Time;
    dimids[1] = dim_StrLen;
    DEF_VAR ("xtime", MPI_CHAR, 2, dimids)
    PUT_ATTR_TXT ("units", "unitless")
    PUT_ATTR_TXT ("long_name", "model time, with format \'YYYY-MM-DD_HH:MM:SS\'")

    /* double kineticEnergyCell(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("kineticEnergyCell", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("units", "m^2 s^{-2}")
    PUT_ATTR_TXT ("long_name", "kinetic energy of horizontal velocity on cells")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double relativeVorticityCell(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("relativeVorticityCell", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("units", "s^{-1}")
    PUT_ATTR_TXT ("long_name",
                  "curl of horizontal velocity, averaged from vertices to cell centers")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double relativeVorticity(Time, nVertices, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("relativeVorticity", MPI_DOUBLE, ndims, rec_D[4])
    PUT_ATTR_TXT ("units", "s^{-1}")
    PUT_ATTR_TXT ("long_name", "curl of horizontal velocity, defined at vertices")
    PUT_ATTR_DECOMP (five, 3, dimids_D[4])

    /* double divergence(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("divergence", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("units", "s^{-1}")
    PUT_ATTR_TXT ("long_name", "divergence of horizontal velocity")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double areaCellGlobal(Time) */
    DEF_VAR ("areaCellGlobal", MPI_DOUBLE, 1, &dim_Time)
    PUT_ATTR_TXT ("units", "m^2")
    PUT_ATTR_TXT ("long_name",
                  "sum of the areaCell variable over the full domain, used to normalize global "
                  "statistics")

    /* double areaEdgeGlobal(Time) */
    DEF_VAR ("areaEdgeGlobal", MPI_DOUBLE, 1, &dim_Time)
    PUT_ATTR_TXT ("units", "m^2")
    PUT_ATTR_TXT ("long_name",
                  "sum of the areaEdge variable over the full domain, used to normalize global "
                  "statistics")

    /* double areaTriangleGlobal(Time) */
    DEF_VAR ("areaTriangleGlobal", MPI_DOUBLE, 1, &dim_Time)
    PUT_ATTR_TXT ("units", "m^2")
    PUT_ATTR_TXT ("long_name",
                  "sum of the areaTriangle variable over the full domain, used to normalize "
                  "global statistics")

    /* double volumeCellGlobal(Time) */
    DEF_VAR ("volumeCellGlobal", MPI_DOUBLE, 1, &dim_Time)
    PUT_ATTR_TXT ("units", "m^3")
    PUT_ATTR_TXT ("long_name",
                  "sum of the volumeCell variable over the full domain, used to normalize "
                  "global statistics")

    /* double volumeEdgeGlobal(Time) */
    DEF_VAR ("volumeEdgeGlobal", MPI_DOUBLE, 1, &dim_Time)
    PUT_ATTR_TXT ("units", "m^3")
    PUT_ATTR_TXT ("long_name",
                  "sum of the volumeEdge variable over the full domain, used to normalize "
                  "global statistics")

    /* double CFLNumberGlobal(Time) */
    DEF_VAR ("CFLNumberGlobal", MPI_DOUBLE, 1, &dim_Time)
    PUT_ATTR_TXT ("units", "unitless")
    PUT_ATTR_TXT ("long_name", "maximum CFL number over the full domain")

    /* double BruntVaisalaFreqTop(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("BruntVaisalaFreqTop", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("units", "s^{-2}")
    PUT_ATTR_TXT ("long_name",
                  "Brunt Vaisala frequency defined at the center (horizontally) and top "
                  "(vertically) of cell")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double vertVelocityTop(Time, nCells, nVertLevelsP1) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("vertVelocityTop", MPI_DOUBLE, ndims, rec_D[5])
    PUT_ATTR_TXT ("units", "m s^{-1}")
    PUT_ATTR_TXT ("long_name",
                  "vertical velocity defined at center (horizontally) and top (vertically) of cell")
    PUT_ATTR_DECOMP (six, 3, dimids_D[5])

    /* double velocityZonal(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("velocityZonal", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("units", "m s^{-1}")
    PUT_ATTR_TXT ("long_name", "component of horizontal velocity in the eastward direction")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double velocityMeridional(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("velocityMeridional", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("units", "m s^{-1}")
    PUT_ATTR_TXT ("long_name", "component of horizontal velocity in the northward direction")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double displacedDensity(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("displacedDensity", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("units", "kg m^{-3}")
    PUT_ATTR_TXT ("long_name",
                  "Density displaced adiabatically to the mid-depth one layer deeper.  That "
                  "is, layer k has been displaced to the depth of layer k+1.")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double potentialDensity(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("potentialDensity", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("units", "kg m^{-3}")
    PUT_ATTR_TXT (
        "long_name",
        "potential density: density displaced adiabatically to the mid-depth of top layer")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double pressure(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("pressure", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("units", "N m^{-2}")
    PUT_ATTR_TXT ("long_name", "pressure used in the momentum equation")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double refBottomDepth(nVertLevels) */
    DEF_VAR ("refBottomDepth", MPI_DOUBLE, 1, &dim_nVertLevels)
    PUT_ATTR_TXT ("units", "m")
    PUT_ATTR_TXT (
        "long_name",
        "Reference depth of ocean for each vertical level. Used in \'z-level\' type runs.")

    /* double zMid(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("zMid", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("units", "m")
    PUT_ATTR_TXT ("long_name", "z-coordinate of the mid-depth of the layer")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double bottomDepth(nCells) */
    DEF_VAR ("bottomDepth", MPI_DOUBLE, 1, fix_D[0])
    PUT_ATTR_TXT ("units", "m")
    PUT_ATTR_TXT ("long_name",
                  "Depth of the bottom of the ocean. Given as a positive distance from sea level.")
    PUT_ATTR_DECOMP (one, 1, dimids_D[0] + 1)

    /* int maxLevelCell(nCells) */
    DEF_VAR ("maxLevelCell", MPI_INT, 1, fix_D[0])
    PUT_ATTR_TXT ("units", "unitless")
    PUT_ATTR_TXT ("long_name", "Index to the last active ocean cell in each column.")
    PUT_ATTR_DECOMP (one, 1, dimids_D[0] + 1)

    /* int maxLevelEdgeBot(nEdges) */
    DEF_VAR ("maxLevelEdgeBot", MPI_INT, 1, fix_D[1])
    PUT_ATTR_TXT ("units", "unitless")
    PUT_ATTR_TXT ("long_name",
                  "Index to the last edge in a column with at least one active ocean cell on "
                  "either side of it.")
    PUT_ATTR_DECOMP (two, 1, dimids_D[1] + 1)

    /* double columnIntegratedSpeed(Time, nCells) */
    DEF_VAR ("columnIntegratedSpeed", MPI_DOUBLE, 2, rec_D[0])
    PUT_ATTR_TXT ("units", "m^2 s^{-1}")
    PUT_ATTR_TXT ("long_name",
                  "speed = sum(h*sqrt(2*ke)), where ke is kineticEnergyCell and the sum is "
                  "over the full column at cell centers.")
    PUT_ATTR_DECOMP (one, 2, dimids_D[0])

    /* double temperatureHorizontalAdvectionTendency(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("temperatureHorizontalAdvectionTendency", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("long_name", "potential temperature tendency due to horizontal advection")
    PUT_ATTR_TXT ("units", "degrees Celsius per second")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double salinityHorizontalAdvectionTendency(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("salinityHorizontalAdvectionTendency", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("long_name", "salinity tendency due to horizontal advection")
    PUT_ATTR_TXT ("units", "PSU per second")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double temperatureVerticalAdvectionTendency(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("temperatureVerticalAdvectionTendency", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("long_name", "potential temperature tendency due to vertical advection")
    PUT_ATTR_TXT ("units", "degrees Celsius per second")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double salinityVerticalAdvectionTendency(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("salinityVerticalAdvectionTendency", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("long_name", "salinity tendency due to vertical advection")
    PUT_ATTR_TXT ("units", "PSU per second")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double temperatureVertMixTendency(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("temperatureVertMixTendency", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("long_name", "potential temperature tendency due to vertical mixing")
    PUT_ATTR_TXT ("units", "degrees Celsius per second")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double salinityVertMixTendency(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("salinityVertMixTendency", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("long_name", "salinity tendency due to vertical mixing")
    PUT_ATTR_TXT ("units", "PSU per second")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double temperatureSurfaceFluxTendency(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("temperatureSurfaceFluxTendency", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("long_name", "potential temperature tendency due to surface fluxes")
    PUT_ATTR_TXT ("units", "degrees Celsius per second")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double salinitySurfaceFluxTendency(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("salinitySurfaceFluxTendency", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("long_name", "salinity tendency due to surface fluxes")
    PUT_ATTR_TXT ("units", "PSU per second")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double temperatureShortWaveTendency(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("temperatureShortWaveTendency", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("units", "degrees Celsius per second")
    PUT_ATTR_TXT ("long_name", "potential temperature tendency due to penetrating shortwave")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double temperatureNonLocalTendency(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("temperatureNonLocalTendency", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("long_name", "potential temperature tendency due to kpp non-local flux")
    PUT_ATTR_TXT ("units", "degrees Celsius per second")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double salinityNonLocalTendency(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("salinityNonLocalTendency", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("long_name", "salinity tendency due to kpp non-local flux")
    PUT_ATTR_TXT ("units", "PSU per second")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double temperature(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("temperature", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("long_name", "potential temperature")
    PUT_ATTR_TXT ("units", "degrees Celsius")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    /* double salinity(Time, nCells, nVertLevels) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR ("salinity", MPI_DOUBLE, ndims, rec_D[2])
    PUT_ATTR_TXT ("long_name", "salinity")
    PUT_ATTR_TXT ("units", "grams salt per kilogram seawater")
    PUT_ATTR_DECOMP (three, 3, dimids_D[2])

    assert (varid - varids + 1 == cfg.nvars);

err_out:
    return err;
}
