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
//
#include <mpi.h>
//
#include <e3sm_io.h>
#include <e3sm_io_driver.hpp>
#include <e3sm_io_err.h>

#define GET_ATT_TEXT(F, D, N, S, B)     driver.get_att (F, D, N, (void *)attbuf);
#define GET_ATT_FLOAT(F, D, N, T, S, B) driver.get_att (F, D, N, (float *)attbuf);
#define GET_ATT_INT(F, D, N, T, S, B)   driver.get_att (F, D, N, (int *)attbuf);
#define GET_ATT(F, D, N, T, S, B)       driver.get_att (F, D, N, (void *)attbuf);

#define INQ_VID(F, N, T, S, B, V) driver.inq_var (F, N, V);

static char attbuf[4096];

static int inq_global_attributes (e3sm_io_driver &driver, int ncid) {
    int err, iattr;
    double dattr;
    char txtbuf[1024];

    /* 963 global attributes: */
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "on_a_sphere", 3, "YES");
    CHECK_ERR
    dattr = 6371229.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "sphere_radius", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "model_name", 4, "mpas");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "core_name", 5, "ocean");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "history", 28, "mpirun -n 9600 ./ocean_model");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "source", 4, "MPAS");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "Conventions", 4, "MPAS");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "parent_id", 65,
                        "mt6vdoeok9\nrz1tbn0bed\n555vk5hkh9\n0s1lcmezuy\n7z1uysqc5i\nwj661vqvze");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "mesh_spec", 3, "0.0");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "git_version", 16, "MPAS_GIT_VERSION");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ocean_run_mode", 7, "forward");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_do_restart", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_restart_timestamp_name", 12,
                        "rpointer.ocn");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_start_time", 19, "0001-01-01_00:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_stop_time", 4, "none");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_run_duration", 19, "0001-00-00_00:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_calendar_type", 16, "gregorian_noleap");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_write_output_on_startup", 2, "NO");
    CHECK_ERR
    iattr = 0;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_pio_num_iotasks", MPI_INT, 1, &iattr);
    CHECK_ERR
    iattr = 1;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_pio_stride", MPI_INT, 1, &iattr);
    CHECK_ERR
    iattr = 3;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_num_halos", MPI_INT, 1, &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_block_decomp_file_prefix", 89,
                        "/project/projectdirs/acme/inputdata/ocn/mpas-o/oRRS18to6v3/"
                        "mpas-o.graph.info.170111.part.");
    CHECK_ERR
    iattr = 0;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_number_of_blocks", MPI_INT, 1, &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_explicit_proc_decomp", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_proc_decomp_file_prefix", 16,
                        "graph.info.part.");
    CHECK_ERR
    iattr = -1;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_vert_levels", MPI_INT, 1, &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_init_configuration", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_expand_sphere", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_realistic_coriolis_parameter", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_write_cull_cell_mask", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_vertical_grid", 7, "uniform");
    CHECK_ERR
    dattr = 1.077;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_1dCVTgenerator_stretch1", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 1.0275;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_1dCVTgenerator_stretch2", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 1.2;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_1dCVTgenerator_dzSeed", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iterative_init_variable", 15,
                        "landIcePressure");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_dt", 8, "00:06:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_time_integrator", 14, "split_explicit");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_vert_coord_movement", 18,
                        "uniform_stretching");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_min_max_thickness", 2, "NO");
    CHECK_ERR
    dattr = 1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_min_thickness", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 6.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_max_thickness_factor", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_dzdk_positive", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_freq_filtered_thickness", 2, "NO");
    CHECK_ERR
    dattr = 5.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_thickness_filter_timescale", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_highFreqThick_restore", 2, "NO");
    CHECK_ERR
    dattr = 30.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_highFreqThick_restore_time", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_highFreqThick_del2", 2, "NO");
    CHECK_ERR
    dattr = 100.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_highFreqThick_del2", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_alter_ICs_for_pbcs", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_pbc_alteration_type", 9, "full_cell");
    CHECK_ERR
    dattr = 0.1;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_min_pbc_fraction", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_hmix_scaleWithMesh", 3, "YES");
    CHECK_ERR
    dattr = -1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_maxMeshDensity", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_apvm_scale_factor", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_mom_del2", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_tracer_del2", 2, "NO");
    CHECK_ERR
    dattr = 10.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_mom_del2", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 10.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_tracer_del2", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_mom_del4", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_tracer_del4", 2, "NO");
    CHECK_ERR
    dattr = 3200000000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_mom_del4", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_mom_del4_div_factor", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_tracer_del4", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_Leith_del2", 2, "NO");
    CHECK_ERR
    dattr = 1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_Leith_parameter", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 15000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_Leith_dx", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 2500.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_Leith_visc2_max", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_standardGM", 2, "NO");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_Redi_surface_layer_tapering", 2, "NO");
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_Redi_surface_layer_tapering_extent",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_Redi_bottom_layer_tapering", 2, "NO");
    CHECK_ERR
    dattr = 0.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_Redi_bottom_layer_tapering_depth", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_GM_Bolus_kappa_function", 8, "constant");
    CHECK_ERR
    dattr = 1800.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_standardGM_tracer_kappa", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_GM_Bolus_kappa_min", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 600.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_GM_Bolus_kappa_max", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 20000.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_GM_Bolus_cell_size_min", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 30000.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_GM_Bolus_cell_size_max", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_Redi_kappa", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.3;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_gravWaveSpeed_trunc", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.01;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_max_relative_slope", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_mom_del2_tensor", 2, "NO");
    CHECK_ERR
    dattr = 10.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_mom_del2_tensor", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_mom_del4_tensor", 2, "NO");
    CHECK_ERR
    dattr = 50000000000000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_mom_del4_tensor", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_Rayleigh_friction", 2, "NO");
    CHECK_ERR
    dattr = 0.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_Rayleigh_damping_coeff", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_Rayleigh_bottom_friction", 2, "NO");
    CHECK_ERR
    dattr = 0.0001;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_Rayleigh_bottom_damping_coeff", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_cvmix", 3, "YES");
    CHECK_ERR
    dattr = 1.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_prandtl_number", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_cvmix_background", 3, "YES");
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_background_diffusion", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 0.0001;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_background_viscosity", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_cvmix_convection", 3, "YES");
    CHECK_ERR
    dattr = 1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_convective_diffusion", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_convective_viscosity", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_convective_basedOnBVF", 3, "YES");
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_convective_triggerBVF", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_cvmix_shear", 3, "YES");
    CHECK_ERR
    iattr = 2;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_num_ri_smooth_loops", MPI_INT, 1, &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_use_BLD_smoothing", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_shear_mixing_scheme", 3, "KPP");
    CHECK_ERR
    dattr = 0.005;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_shear_PP_nu_zero", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 5.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_shear_PP_alpha", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 2.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_shear_PP_exp", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.005;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_shear_KPP_nu_zero", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 0.7;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_shear_KPP_Ri_zero", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 3.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_shear_KPP_exp", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_cvmix_tidal_mixing", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_cvmix_double_diffusion", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_cvmix_kpp", 3, "YES");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_cvmix_fixed_boundary_layer", 2, "NO");
    CHECK_ERR
    dattr = 30.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_kpp_boundary_layer_depth", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 0.25;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_kpp_criticalBulkRichardsonNumber",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_kpp_matching", 12, "SimpleShapes");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_kpp_EkmanOBL", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_kpp_MonObOBL", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_kpp_interpolationOMLType", 9,
                        "quadratic");
    CHECK_ERR
    dattr = 0.1;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_kpp_surface_layer_extent", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 5.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_kpp_surface_layer_averaging",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 10.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "configure_cvmix_kpp_minimum_OBL_under_sea_ice",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 100.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_kpp_stop_OBL_search", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_kpp_use_enhanced_diff", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_kpp_nonlocal_with_implicit_mix", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_const_visc", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_const_diff", 2, "NO");
    CHECK_ERR
    dattr = 0.0001;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_vert_visc", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1.e-05;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_vert_diff", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_rich_visc", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_rich_diff", 2, "NO");
    CHECK_ERR
    dattr = 0.0001;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_bkrd_vert_visc", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1.e-05;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_bkrd_vert_diff", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.005;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_rich_mix", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_convective_visc", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_convective_diff", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_tanh_visc", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_tanh_diff", 2, "NO");
    CHECK_ERR
    dattr = 0.25;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_max_visc_tanh", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.0001;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_min_visc_tanh", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.025;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_max_diff_tanh", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1.e-05;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_min_diff_tanh", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -100.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_zMid_tanh", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 100.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_zWidth_tanh", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_bulk_wind_stress", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_bulk_thickness_flux", 3, "YES");
    CHECK_ERR
    dattr = 0.001;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_flux_attenuation_coefficient", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 10.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_flux_attenuation_coefficient_runoff",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 86400.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ssh_grad_relax_timescale", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_remove_AIS_coupler_runoff", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sw_absorption_type", 6, "jerlov");
    CHECK_ERR
    iattr = 3;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_jerlov_water_type", MPI_INT, 1, &iattr);
    CHECK_ERR
    dattr = 1.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_surface_buoyancy_depth", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_frazil_ice_formation", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_frazil_in_open_ocean", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_frazil_under_land_ice", 3, "YES");
    CHECK_ERR
    dattr = 333700.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_frazil_heat_of_fusion", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_frazil_ice_density", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.1;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_frazil_fractional_thickness_limit",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 3996.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_specific_heat_sea_water", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 100.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_frazil_maximum_depth", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 4.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_frazil_sea_ice_reference_salinity",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_frazil_land_ice_reference_salinity",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_frazil_maximum_freezing_temperature",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_frazil_use_surface_pressure", 2, "NO");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_land_ice_flux_mode", 13, "pressure_only");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_land_ice_flux_formulation", 7, "Jenkins");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_land_ice_flux_useHollandJenkinsAdvDiff",
                        2, "NO");
    CHECK_ERR
    dattr = 10.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_land_ice_flux_attenuation_coefficient",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 10.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_land_ice_flux_boundaryLayerThickness",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_land_ice_flux_boundaryLayerNeighborWeight",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 2009.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_land_ice_flux_cp_ice", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 918.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_land_ice_flux_rho_ice", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.0025;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_land_ice_flux_topDragCoeff", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 0.0001;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_land_ice_flux_ISOMIP_gammaT", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 0.05;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_land_ice_flux_rms_tidal_velocity", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 0.011;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR,
                   "config_land_ice_flux_jenkins_heat_transfer_coefficient", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.00031;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR,
                   "config_land_ice_flux_jenkins_salt_transfer_coefficient", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_vert_tracer_adv", 7, "stencil");
    CHECK_ERR
    iattr = 3;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_vert_tracer_adv_order", MPI_INT, 1, &iattr);
    CHECK_ERR
    iattr = 3;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_horiz_tracer_adv_order", MPI_INT, 1, &iattr);
    CHECK_ERR
    dattr = 0.25;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_coef_3rd_order", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_monotonic", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_implicit_bottom_drag", 3, "YES");
    CHECK_ERR
    dattr = 0.001;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_implicit_bottom_drag_coeff", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_explicit_bottom_drag", 2, "NO");
    CHECK_ERR
    dattr = 0.001;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_explicit_bottom_drag_coeff", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 1026.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_density0", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_pressure_gradient_type", 16,
                        "Jacobian_from_TS");
    CHECK_ERR
    dattr = 0.5;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_common_level_weight", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_eos_type", 2, "jm");
    CHECK_ERR
    dattr = -1.8;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_open_ocean_freezing_temperature_coeff_0",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_open_ocean_freezing_temperature_coeff_S",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_open_ocean_freezing_temperature_coeff_p",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_open_ocean_freezing_temperature_coeff_pS",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR,
                   "config_open_ocean_freezing_temperature_reference_pressure", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 0.0622;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_land_ice_cavity_freezing_temperature_coeff_0",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -0.0563;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_land_ice_cavity_freezing_temperature_coeff_S",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -7.43e-08;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_land_ice_cavity_freezing_temperature_coeff_p",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -1.74e-10;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR,
                   "config_land_ice_cavity_freezing_temperature_coeff_pS", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR,
                   "config_land_ice_cavity_freezing_temperature_reference_pressure", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 0.2;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_eos_linear_alpha", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.8;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_eos_linear_beta", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 5.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_eos_linear_Tref", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 35.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_eos_linear_Sref", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1000.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_eos_linear_densityref", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    iattr = 2;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_n_ts_iter", MPI_INT, 1, &iattr);
    CHECK_ERR
    iattr = 1;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_n_bcl_iter_beg", MPI_INT, 1, &iattr);
    CHECK_ERR
    iattr = 2;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_n_bcl_iter_mid", MPI_INT, 1, &iattr);
    CHECK_ERR
    iattr = 2;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_n_bcl_iter_end", MPI_INT, 1, &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_btr_dt", 13, "0000_00:00:12");
    CHECK_ERR
    iattr = 2;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_n_btr_cor_iter", MPI_INT, 1, &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_vel_correction", 3, "YES");
    CHECK_ERR
    iattr = 2;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_btr_subcycle_loop_factor", MPI_INT, 1, &iattr);
    CHECK_ERR
    dattr = 0.5;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_btr_gam1_velWt1", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_btr_gam2_SSHWt1", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_btr_gam3_velWt2", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_btr_solve_SSH2", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_conduct_tests", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_test_tensors", 2, "NO");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_tensor_test_function", 11, "sph_uCosCos");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_redi_k33", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_redi_horizontal_term1", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_redi_horizontal_term2", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_redi_horizontal_term3", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_check_zlevel_consistency", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_check_ssh_consistency", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_filter_btr_mode", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_prescribe_velocity", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_prescribe_thickness", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_include_KE_vertex", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_check_tracer_monotonicity", 2, "NO");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_compute_active_tracer_budgets", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_thick_all_tend", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_thick_hadv", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_thick_vadv", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_thick_sflux", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_vel_all_tend", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_vel_coriolis", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_vel_pgrad", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_vel_hmix", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_vel_surface_stress", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_vel_explicit_bottom_drag", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_vel_vmix", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_vel_vadv", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_tr_all_tend", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_tr_adv", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_tr_hmix", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_tr_vmix", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_tr_sflux", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_disable_tr_nonlocalflux", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_read_nearest_restart", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_rx1_constraint", 2, "NO");
    CHECK_ERR
    iattr = 20;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_rx1_outer_iter_count", MPI_INT, 1, &iattr);
    CHECK_ERR
    iattr = 10;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_rx1_inner_iter_count", MPI_INT, 1, &iattr);
    CHECK_ERR
    dattr = 0.1;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_rx1_init_inner_weight", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 5.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_rx1_max", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_rx1_horiz_smooth_weight", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 1.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_rx1_vert_smooth_weight", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.1;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_rx1_slope_weight", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_rx1_zstar_weight", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    iattr = 20;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_rx1_horiz_smooth_open_ocean_cells", MPI_INT,
                   1, &iattr);
    CHECK_ERR
    iattr = 3;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_rx1_min_levels", MPI_INT, 1, &iattr);
    CHECK_ERR
    dattr = 1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_rx1_min_layer_thickness", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    iattr = 20;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_baroclinic_channel_vert_levels", MPI_INT, 1,
                   &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_baroclinic_channel_use_distances", 2,
                        "NO");
    CHECK_ERR
    dattr = 13.1000003814697;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_baroclinic_channel_surface_temperature",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 10.1000003814697;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_baroclinic_channel_bottom_temperature",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1.20000004768372;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_baroclinic_channel_temperature_difference",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.0799999982118607;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_baroclinic_channel_gradient_width_frac",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 40000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_baroclinic_channel_gradient_width_dist",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_baroclinic_channel_bottom_depth", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 35.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_baroclinic_channel_salinity", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = -0.000119999996968545;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_baroclinic_channel_coriolis_parameter",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    iattr = 20;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_lock_exchange_vert_levels", MPI_INT, 1, &iattr);
    CHECK_ERR
    dattr = 20.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_lock_exchange_bottom_depth", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 5.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_lock_exchange_cold_temperature", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 30.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_lock_exchange_warm_temperature", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_lock_exchange_direction", 1, "y");
    CHECK_ERR
    dattr = 35.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_lock_exchange_salinity", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_lock_exchange_layer_type", 7, "z-level");
    CHECK_ERR
    dattr = 0.00999999977648258;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_lock_exchange_isopycnal_min_thickness",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    iattr = 20;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_internal_waves_vert_levels", MPI_INT, 1,
                   &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_internal_waves_use_distances", 2, "NO");
    CHECK_ERR
    dattr = 20.1000003814697;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_internal_waves_surface_temperature",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 10.1000003814697;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_internal_waves_bottom_temperature",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 2.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_internal_waves_temperature_difference",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.330000013113022;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_internal_waves_amplitude_width_frac",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 50000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_internal_waves_amplitude_width_dist",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 500.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_internal_waves_bottom_depth", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 35.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_internal_waves_salinity", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_internal_waves_layer_type", 7, "z-level");
    CHECK_ERR
    dattr = 125.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_internal_waves_isopycnal_displacement",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    iattr = 100;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_overflow_vert_levels", MPI_INT, 1, &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_overflow_use_distances", 2, "NO");
    CHECK_ERR
    dattr = 2000.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_overflow_bottom_depth", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 500.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_overflow_ridge_depth", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 10.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_overflow_plug_temperature", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 20.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_overflow_domain_temperature", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 35.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_overflow_salinity", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.100000001490116;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_overflow_plug_width_frac", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 0.200000002980232;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_overflow_slope_center_frac", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 0.0500000007450581;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_overflow_slope_width_frac", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 20000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_overflow_plug_width_dist", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 40000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_overflow_slope_center_dist", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 7000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_overflow_slope_width_dist", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_overflow_layer_type", 7, "z-level");
    CHECK_ERR
    dattr = 0.00999999977648258;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_overflow_isopycnal_min_thickness", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 15.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_minimum_depth", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_depth_file", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_depth_dimname", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_depth_varname", 4, "none");
    CHECK_ERR
    dattr = 1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_depth_conversion_factor",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_temperature_file", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_salinity_file", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_tracer_nlat_dimname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_tracer_nlon_dimname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_tracer_ndepth_dimname", 4,
                        "none");
    CHECK_ERR
    dattr = 1.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_tracer_depth_conversion_factor",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    iattr = -1;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_tracer_vert_levels", MPI_INT, 1,
                   &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_temperature_varname", 4,
                        "none");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_salinity_varname", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_tracer_latlon_degrees", 3,
                        "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_tracer_lat_varname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_tracer_lon_varname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_tracer_depth_varname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_tracer_method", 22,
                        "bilinear_interpolation");
    CHECK_ERR
    iattr = 0;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_smooth_TS_iterations", MPI_INT,
                   1, &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_swData_file", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_swData_nlat_dimname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_swData_nlon_dimname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_swData_lat_varname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_swData_lon_varname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_swData_latlon_degrees", 3,
                        "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_swData_method", 22,
                        "bilinear_interpolation");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_chlorophyll_varname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_zenithAngle_varname", 4,
                        "none");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_clearSky_varname", 4, "none");
    CHECK_ERR
    dattr = 4.99999987368938e-05;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_piston_velocity", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_interior_restore_rate",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_topography_file", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_topography_nlat_dimname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_topography_nlon_dimname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_topography_latlon_degrees",
                        3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_topography_lat_varname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_topography_lon_varname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_topography_varname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_topography_has_ocean_frac",
                        2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_global_ocean_topography_ocean_frac_varname", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_topography_method", 22,
                        "bilinear_interpolation");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_smooth_topography", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_deepen_critical_passages",
                        3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_depress_by_land_ice", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_land_ice_topo_file", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_land_ice_topo_nlat_dimname",
                        4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_land_ice_topo_nlon_dimname",
                        4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_global_ocean_land_ice_topo_latlon_degrees", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_land_ice_topo_lat_varname",
                        4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_land_ice_topo_lon_varname",
                        4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_global_ocean_land_ice_topo_thickness_varname", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_global_ocean_land_ice_topo_draft_varname", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_global_ocean_land_ice_topo_ice_frac_varname", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_global_ocean_land_ice_topo_grounded_frac_varname", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_global_ocean_use_constant_land_ice_cavity_temperature", 2, "NO");
    CHECK_ERR
    dattr = -1.79999995231628;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR,
                 "config_global_ocean_constant_land_ice_cavity_temperature", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_cull_inland_seas", 3, "YES");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_windstress_file", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_windstress_nlat_dimname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_windstress_nlon_dimname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_windstress_latlon_degrees",
                        3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_windstress_lat_varname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_windstress_lon_varname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_windstress_zonal_varname",
                        4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_global_ocean_windstress_meridional_varname", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_windstress_method", 22,
                        "bilinear_interpolation");
    CHECK_ERR
    dattr = 1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_windstress_conversion_factor",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_ecosys_file", 7, "unknown");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_ecosys_forcing_file", 7,
                        "unknown");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_ecosys_nlat_dimname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_ecosys_nlon_dimname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_ecosys_ndepth_dimname", 4,
                        "none");
    CHECK_ERR
    dattr = 1.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_ecosys_depth_conversion_factor",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    iattr = -1;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_ecosys_vert_levels", MPI_INT, 1,
                   &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_ecosys_lat_varname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_ecosys_lon_varname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_ecosys_depth_varname", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_ecosys_latlon_degrees", 3,
                        "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_ecosys_method", 22,
                        "bilinear_interpolation");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_global_ocean_ecosys_forcing_time_dimname", 4, "none");
    CHECK_ERR
    iattr = 0;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_global_ocean_smooth_ecosys_iterations",
                   MPI_INT, 1, &iattr);
    CHECK_ERR
    iattr = 100;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_vert_levels", MPI_INT, 1, &iattr);
    CHECK_ERR
    dattr = 15.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_surface_temperature", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 35.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_surface_salinity", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 15.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_surface_restoring_temperature",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 35.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_surface_restoring_salinity",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 3.99999998990097e-06;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_temperature_piston_velocity",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 3.99999998990097e-06;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_salinity_piston_velocity",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_sensible_heat_flux", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_latent_heat_flux", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 0.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_shortwave_heat_flux", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_rain_flux", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_evaporation_flux", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 9.99999997475243e-07;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR,
                 "config_cvmix_WSwSBF_interior_temperature_restoring_rate", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 9.99999997475243e-07;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR,
                   "config_cvmix_WSwSBF_interior_salinity_restoring_rate", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.00999999977648258;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_temperature_gradient",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_salinity_gradient", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR,
                   "config_cvmix_WSwSBF_temperature_gradient_mixed_layer", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_salinity_gradient_mixed_layer",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_mixed_layer_depth_temperature",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_mixed_layer_depth_salinity",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_mixed_layer_temperature_change",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_mixed_layer_salinity_change",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_vertical_grid", 7, "uniform");
    CHECK_ERR
    dattr = 400.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_bottom_depth", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 0.100000001490116;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_max_windstress", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 9.99999974737875e-05;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_cvmix_WSwSBF_coriolis_parameter", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    iattr = 100;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_vert_levels", MPI_INT, 1, &iattr);
    CHECK_ERR
    dattr = 4000.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_main_channel_depth", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -50.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_north_wall_lat", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -70.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_south_wall_lat", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_ridge_flag", 3, "YES");
    CHECK_ERR
    dattr = 180.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_ridge_center_lon", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 2000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_ridge_height", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 2000000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_ridge_width", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_plateau_flag", 3, "YES");
    CHECK_ERR
    dattr = 300.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_plateau_center_lon", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -58.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_plateau_center_lat", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 2000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_plateau_height", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 200000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_plateau_radius", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1000000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_plateau_slope_width", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_shelf_flag", 3, "YES");
    CHECK_ERR
    dattr = 500.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_shelf_depth", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 120000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_shelf_width", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_cont_slope_flag", 3, "YES");
    CHECK_ERR
    dattr = 0.00999999977648258;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_max_cont_slope", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_embayment_flag", 3, "YES");
    CHECK_ERR
    dattr = 60.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_embayment_center_lon", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = -71.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_embayment_center_lat", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 500000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_embayment_radius", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 2000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_embayment_depth", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_depression_flag", 3, "YES");
    CHECK_ERR
    dattr = 60.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_depression_center_lon", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = -72.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_depression_south_lat", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = -65.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_depression_north_lat", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 480000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_depression_width", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 800.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_depression_depth", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 35.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_salinity", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.00999999977648258;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_wind_stress_max", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.200000002980232;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_acc_wind", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -0.0500000007450581;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_asf_wind", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -65.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_wind_trans", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -5.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_heat_flux_south", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 10.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_heat_flux_middle", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -5.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_heat_flux_north", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -70.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_heat_flux_lat_ss", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -65.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_heat_flux_lat_sm", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -53.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_heat_flux_lat_mn", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 60.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_region1_center_lon", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -75.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_region1_center_lat", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 150.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_region2_center_lon", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -71.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_region2_center_lat", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 240.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_region3_center_lon", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -71.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_region3_center_lat", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 330.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_region4_center_lon", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -71.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_region4_center_lat", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_heat_flux_region1_flag", 2, "NO");
    CHECK_ERR
    dattr = -5.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_heat_flux_region1", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 300000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_heat_flux_region1_radius", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_heat_flux_region2_flag", 2, "NO");
    CHECK_ERR
    dattr = -5.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_heat_flux_region2", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 240000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_heat_flux_region2_radius", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 5.80000014451798e-05;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_surface_temperature_piston_velocity",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 3.5;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_initial_temp_t1", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 4.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_initial_temp_t2", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1200.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_initial_temp_h0", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 500.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_initial_temp_h1", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 7.50000035623088e-05;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_initial_temp_mt", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -75.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_initial_temp_latS", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -50.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_initial_temp_latN", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 10.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_sponge_t1", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 1000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_sponge_h1", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 120000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_sponge_l1", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 10.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_sponge_tau1", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_restore_region1_flag", 3,
                        "YES");
    CHECK_ERR
    dattr = -1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_restore_t1", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 600000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_restore_lcx1", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 600000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_restore_lcy1", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_restore_region2_flag", 3,
                        "YES");
    CHECK_ERR
    dattr = -1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_restore_t2", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 600000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_restore_lcx2", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 250000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_restore_lcy2", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_restore_region3_flag", 3,
                        "YES");
    CHECK_ERR
    dattr = -1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_restore_t3", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 600000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_restore_lcx3", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 250000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_restore_lcy3", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_restore_region4_flag", 3,
                        "YES");
    CHECK_ERR
    dattr = -1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_restore_t4", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 600000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_restore_lcx4", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 250000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_iso_temperature_restore_lcy4", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    iattr = 100;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_soma_vert_levels", MPI_INT, 1, &iattr);
    CHECK_ERR
    dattr = 1250000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_soma_domain_width", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 35.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_soma_center_latitude", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_soma_center_longitude", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.100000001490116;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_soma_phi", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 2500.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_soma_bottom_depth", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -0.400000005960464;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_soma_shelf_width", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 100.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_soma_shelf_depth", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_soma_ref_density", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 4.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_soma_density_difference", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 300.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_soma_thermocline_depth", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.0500000007450581;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_soma_density_difference_linear", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 20.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_soma_surface_temperature", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 33.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_soma_surface_salinity", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_soma_use_surface_temp_restoring", 2, "NO");
    CHECK_ERR
    dattr = 7.5;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR,
                   "config_soma_surface_temp_restoring_at_center_latitude", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.5;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR,
                   "config_soma_surface_temp_restoring_latitude_gradient", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 9.99999974737875e-06;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_soma_restoring_temp_piston_vel", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    iattr = 100;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_vert_levels", MPI_INT, 1, &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_add_easterly_wind_stress_ASF", 2,
                        "NO");
    CHECK_ERR
    dattr = 800000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_wind_transition_position", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 600000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_antarctic_shelf_front_width", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = -0.0500000007450581;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_wind_stress_shelf_front_max", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_use_slopping_bathymetry", 2, "NO");
    CHECK_ERR
    dattr = 2000000.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_meridional_extent", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1000000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_zonal_extent", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 2500.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_bottom_depth", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 500.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_shelf_depth", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 100000.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_slope_half_width", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 500000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_slope_center_position", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = -9.99999974737875e-05;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_reference_coriolis", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 0.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_coriolis_gradient", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.200000002980232;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_wind_stress_max", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 3.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_mean_restoring_temp", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 2.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_restoring_temp_dev_ta", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 2.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_restoring_temp_dev_tb", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 30.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_restoring_temp_tau", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 1.89999991562217e-05;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_restoring_temp_piston_vel", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 1250.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_restoring_temp_ze", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 80000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_restoring_sponge_l", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 6.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_initial_temp_t1", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 3.59999990463257;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_initial_temp_t2", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 300.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_initial_temp_h1", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 7.50000035623088e-05;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_initial_temp_mt", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_frazil_enable", 2, "NO");
    CHECK_ERR
    dattr = -3.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ziso_frazil_temperature_anomaly", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    iattr = 20;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sub_ice_shelf_2D_vert_levels", MPI_INT, 1,
                   &iattr);
    CHECK_ERR
    dattr = 2000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sub_ice_shelf_2D_bottom_depth", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 25.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sub_ice_shelf_2D_cavity_thickness",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 500.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sub_ice_shelf_2D_slope_height", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 15000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sub_ice_shelf_2D_edge_width", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 30000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sub_ice_shelf_2D_y1", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 60000.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sub_ice_shelf_2D_y2", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sub_ice_shelf_2D_temperature", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 34.5;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sub_ice_shelf_2D_surface_salinity",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 34.7000007629395;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sub_ice_shelf_2D_bottom_salinity", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    iattr = 100;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_periodic_planar_vert_levels", MPI_INT, 1,
                   &iattr);
    CHECK_ERR
    dattr = 2500.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_periodic_planar_bottom_depth", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 1.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_periodic_planar_velocity_strength",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    iattr = 100;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ecosys_column_vert_levels", MPI_INT, 1, &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ecosys_column_vertical_grid", 14,
                        "100layerACMEv1");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ecosys_column_TS_filename", 7, "unknown");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ecosys_column_ecosys_filename", 7,
                        "unknown");
    CHECK_ERR
    dattr = 6000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ecosys_column_bottom_depth", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    iattr = 10;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_vert_levels", MPI_INT, 1, &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_layer_type", 5, "sigma");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_stratification_type", 11,
                        "exponential");
    CHECK_ERR
    dattr = 1024.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_density_coef_linear", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 1028.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_density_coef_exp", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 0.100000001490116;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_density_gradient_linear",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 3.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_density_gradient_exp", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 4500.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_density_depth_linear", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 500.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_density_depth_exp", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 1028.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_density_ref", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 5.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_density_Tref", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.200000002980232;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_density_alpha", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 5000.;
    err =
        GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_bottom_depth", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 4500.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_height", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 10000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_radius", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 40000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_width", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 35.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_salinity", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -9.99999974737875e-05;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_sea_mount_coriolis_parameter", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    iattr = 30;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_vert_levels", MPI_INT, 1, &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_vertical_level_distribution", 8,
                        "constant");
    CHECK_ERR
    dattr = -900.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_bottom_depth", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -1.89999997615814;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_temperature", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 34.4000015258789;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_salinity", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -1.89999997615814;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_restoring_temperature", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 1.20000004244503e-05;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_temperature_piston_velocity",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 34.4000015258789;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_restoring_salinity", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 1.20000004244503e-05;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_salinity_piston_velocity", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = -0.00014000000373926;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_coriolis_parameter", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_southern_boundary", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 1000000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_northern_boundary", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_western_boundary", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 500000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_eastern_boundary", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 0.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_y1", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -700.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_z1", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_ice_fraction1", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 400000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_y2", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -200.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_z2", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_ice_fraction2", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1000000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_y3", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -200.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_z3", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 1.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_ice_fraction3", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    iattr = 36;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_vert_levels", MPI_INT, 1, &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_vertical_level_distribution",
                        8, "constant");
    CHECK_ERR
    dattr = -720.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_max_bottom_depth", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    iattr = 3;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_minimum_levels", MPI_INT, 1,
                   &iattr);
    CHECK_ERR
    dattr = 10.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_min_column_thickness", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 0.5;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_min_ocean_fraction", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_topography_file", 27,
                        "input_geometry_processed.nc");
    CHECK_ERR
    dattr = -1.89999997615814;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_init_top_temp", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = -1.89999997615814;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_init_bot_temp", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 33.7999992370605;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_init_top_sal", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 34.5;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_init_bot_sal", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = -1.89999997615814;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_restore_top_temp", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 1.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_restore_bot_temp", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 33.7999992370605;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_restore_top_sal", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 34.7000007629395;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_restore_bot_sal", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 10.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_restore_rate", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 200.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_restore_evap_rate", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 790000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_restore_xMin", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 800000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_restore_xMax", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = -0.000140999996801838;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_coriolis_parameter", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 1026.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_isomip_plus_effective_density", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_activeTracers", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_activeTracers_surface_bulk_forcing",
                        3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_activeTracers_surface_restoring", 3,
                        "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_activeTracers_interior_restoring", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_activeTracers_exponential_decay", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_activeTracers_idealAge_forcing", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_activeTracers_ttd_forcing", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_surface_salinity_monthly_restoring",
                        3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_surface_salinity_monthly_restoring_compute_interval", 19,
                        "0000-00-01_00:00:00");
    CHECK_ERR
    dattr = 1.585e-06;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_salinity_restoring_constant_piston_velocity",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 100;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_salinity_restoring_max_difference",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_salinity_restoring_under_sea_ice", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_debugTracers", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_debugTracers_surface_bulk_forcing",
                        2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_debugTracers_surface_restoring", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_debugTracers_interior_restoring", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_debugTracers_exponential_decay", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_debugTracers_idealAge_forcing", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_debugTracers_ttd_forcing", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_ecosysTracers", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ecosys_atm_co2_option", 8, "constant");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ecosys_atm_alt_co2_option", 8, "constant");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ecosys_atm_alt_co2_use_eco", 2, "NO");
    CHECK_ERR
    dattr = 379.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ecosys_atm_co2_constant_value", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_ecosysTracers_surface_bulk_forcing",
                        2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_ecosysTracers_surface_restoring", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_ecosysTracers_interior_restoring", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_ecosysTracers_exponential_decay", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_ecosysTracers_idealAge_forcing", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_ecosysTracers_ttd_forcing", 2, "NO");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_ecosysTracers_surface_value", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_ecosysTracers_sea_ice_coupling", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ecosysTracers_diagnostic_fields_level1",
                        2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ecosysTracers_diagnostic_fields_level2",
                        2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ecosysTracers_diagnostic_fields_level3",
                        2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ecosysTracers_diagnostic_fields_level4",
                        2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_ecosysTracers_diagnostic_fields_level5",
                        2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_DMSTracers", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_DMSTracers_surface_bulk_forcing", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_DMSTracers_surface_restoring", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_DMSTracers_interior_restoring", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_DMSTracers_exponential_decay", 2,
                        "NO");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_DMSTracers_idealAge_forcing", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_DMSTracers_ttd_forcing", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_DMSTracers_surface_value", 2, "NO");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_DMSTracers_sea_ice_coupling", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_MacroMoleculesTracers", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_use_MacroMoleculesTracers_surface_bulk_forcing", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_use_MacroMoleculesTracers_surface_restoring", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_use_MacroMoleculesTracers_interior_restoring", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_use_MacroMoleculesTracers_exponential_decay", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_use_MacroMoleculesTracers_idealAge_forcing", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_MacroMoleculesTracers_ttd_forcing",
                        2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_use_MacroMoleculesTracers_surface_value",
                        2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_use_MacroMoleculesTracers_sea_ice_coupling", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_globalStats_enable", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_globalStats_compute_interval", 15,
                        "output_interval");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_globalStats_compute_on_startup", 3,
                        "YES");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_globalStats_write_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_globalStats_text_file", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_globalStats_directory", 16,
                        "analysis_members");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_globalStats_output_stream", 17,
                        "globalStatsOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_surfaceAreaWeightedAverages_enable",
                        3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_surfaceAreaWeightedAverages_compute_on_startup", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_surfaceAreaWeightedAverages_write_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_surfaceAreaWeightedAverages_compute_interval", 19,
                        "0000-00-00_01:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_surfaceAreaWeightedAverages_output_stream", 33,
                        "surfaceAreaWeightedAveragesOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_waterMassCensus_enable", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_waterMassCensus_compute_interval", 19,
                        "0000-00-00_01:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_waterMassCensus_output_stream", 21,
                        "waterMassCensusOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_waterMassCensus_compute_on_startup",
                        3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_waterMassCensus_write_on_startup", 2,
                        "NO");
    CHECK_ERR
    dattr = -2.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_waterMassCensus_minTemperature",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 30.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_waterMassCensus_maxTemperature",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 32.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_waterMassCensus_minSalinity", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 37.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_waterMassCensus_maxSalinity", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_waterMassCensus_compute_predefined_regions", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_waterMassCensus_region_group", 0, "");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_layerVolumeWeightedAverage_enable", 3,
                        "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_layerVolumeWeightedAverage_compute_interval", 19,
                        "0000-00-00_01:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_layerVolumeWeightedAverage_compute_on_startup", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_layerVolumeWeightedAverage_write_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_layerVolumeWeightedAverage_output_stream", 32,
                        "layerVolumeWeightedAverageOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_zonalMean_enable", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_zonalMean_compute_on_startup", 3,
                        "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_zonalMean_write_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_zonalMean_compute_interval", 19,
                        "0000-00-00_01:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_zonalMean_output_stream", 15,
                        "zonalMeanOutput");
    CHECK_ERR
    iattr = 180;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_zonalMean_num_bins", MPI_INT, 1, &iattr);
    CHECK_ERR
    dattr = -1.e+34;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_zonalMean_min_bin", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -1.e+34;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_zonalMean_max_bin", MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_okuboWeiss_enable", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_okuboWeiss_compute_on_startup", 3,
                        "YES");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_okuboWeiss_write_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_okuboWeiss_compute_interval", 19,
                        "0000-00-00_01:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_okuboWeiss_output_stream", 16,
                        "okuboWeissOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_okuboWeiss_directory", 16,
                        "analysis_members");
    CHECK_ERR
    dattr = -0.2;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_okuboWeiss_threshold_value", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 1.e-10;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_okuboWeiss_normalization", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = 1.e-10;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_okuboWeiss_lambda2_normalization",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_okuboWeiss_use_lat_lon_coords", 3,
                        "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_okuboWeiss_compute_eddy_census", 3,
                        "YES");
    CHECK_ERR
    iattr = 20;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_okuboWeiss_eddy_min_cells", MPI_INT, 1,
                   &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_meridionalHeatTransport_enable", 3,
                        "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_meridionalHeatTransport_compute_interval", 19,
                        "0000-00-00_01:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_meridionalHeatTransport_compute_on_startup", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_meridionalHeatTransport_write_on_startup", 2, "NO");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_meridionalHeatTransport_output_stream",
                      29, "meridionalHeatTransportOutput");
    CHECK_ERR
    iattr = 180;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_meridionalHeatTransport_num_bins", MPI_INT,
                   1, &iattr);
    CHECK_ERR
    dattr = -1.e+34;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_meridionalHeatTransport_min_bin",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = -1.e+34;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_meridionalHeatTransport_max_bin",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_meridionalHeatTransport_region_group",
                        0, "");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_testComputeInterval_enable", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_testComputeInterval_compute_interval",
                        17, "00-00-01_00:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_testComputeInterval_compute_on_startup", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_testComputeInterval_write_on_startup",
                        2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_testComputeInterval_output_stream",
                        25, "testComputeIntervalOutput");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_highFrequencyOutput_enable", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_highFrequencyOutput_compute_interval",
                        15, "output_interval");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_highFrequencyOutput_output_stream",
                        19, "highFrequencyOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_highFrequencyOutput_compute_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_highFrequencyOutput_write_on_startup",
                        2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeFilters_enable", 2, "NO");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeFilters_compute_interval", 2, "dt");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeFilters_output_stream", 17,
                        "timeFiltersOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeFilters_restart_stream", 18,
                        "timeFiltersRestart");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeFilters_compute_on_startup", 3,
                        "YES");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeFilters_write_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeFilters_initialize_filters", 3,
                        "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeFilters_tau", 11, "90_00:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeFilters_compute_cell_centered_values", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_lagrPartTrack_enable", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_lagrPartTrack_compute_interval", 2,
                        "dt");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_lagrPartTrack_compute_on_startup", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_lagrPartTrack_output_stream", 19,
                        "lagrPartTrackOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_lagrPartTrack_restart_stream", 20,
                        "lagrPartTrackRestart");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_lagrPartTrack_input_stream", 18,
                        "lagrPartTrackInput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_lagrPartTrack_write_on_startup", 2,
                        "NO");
    CHECK_ERR
    iattr = 0;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_lagrPartTrack_filter_number", MPI_INT, 1,
                   &iattr);
    CHECK_ERR
    iattr = 2;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_lagrPartTrack_timeIntegration", MPI_INT, 1,
                   &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_lagrPartTrack_reset_criteria", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_lagrPartTrack_reset_global_timestamp",
                        13, "0000_00:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_lagrPartTrack_region_stream", 20,
                        "lagrPartTrackRegions");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_lagrPartTrack_reset_if_outside_region", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_lagrPartTrack_reset_if_inside_region",
                        2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_eliassenPalm_enable", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_eliassenPalm_compute_interval", 15,
                        "output_interval");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_eliassenPalm_output_stream", 18,
                        "eliassenPalmOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_eliassenPalm_restart_stream", 19,
                        "eliassenPalmRestart");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_eliassenPalm_compute_on_startup", 3,
                        "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_eliassenPalm_write_on_startup", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_eliassenPalm_debug", 2, "NO");
    CHECK_ERR
    iattr = 45;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_eliassenPalm_nBuoyancyLayers", MPI_INT, 1,
                   &iattr);
    CHECK_ERR
    dattr = 900.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_eliassenPalm_rhomin_buoycoor", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    dattr = 1080.;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_eliassenPalm_rhomax_buoycoor", MPI_DOUBLE,
                   1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mixedLayerDepths_enable", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mixedLayerDepths_compute_interval",
                        19, "0000-00-00_01:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mixedLayerDepths_output_stream", 22,
                        "mixedLayerDepthsOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mixedLayerDepths_write_on_startup", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mixedLayerDepths_compute_on_startup",
                        3, "YES");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mixedLayerDepths_Tthreshold", 3, "YES");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mixedLayerDepths_Dthreshold", 3, "YES");
    CHECK_ERR
    dattr = 0.2;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mixedLayerDepths_crit_temp_threshold",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 0.03;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mixedLayerDepths_crit_dens_threshold",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 100000.;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mixedLayerDepths_reference_pressure",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mixedLayerDepths_Tgradient", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mixedLayerDepths_Dgradient", 2, "NO");
    CHECK_ERR
    dattr = 5.e-07;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mixedLayerDepths_temp_gradient_threshold",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    dattr = 5.e-08;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mixedLayerDepths_den_gradient_threshold",
                   MPI_DOUBLE, 1, &dattr);
    CHECK_ERR
    iattr = 1;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mixedLayerDepths_interp_method", MPI_INT,
                   1, &iattr);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsDaily_enable", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsDaily_compute_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsDaily_write_on_startup",
                        2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsDaily_compute_interval",
                        15, "output_interval");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsDaily_output_stream", 24,
                        "regionalStatsDailyOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsDaily_restart_stream",
                        18, "regionalMasksInput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsDaily_input_stream", 18,
                        "regionalMasksInput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsDaily_operation", 3,
                        "avg");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsDaily_region_type", 4,
                        "cell");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsDaily_region_group", 3,
                        "all");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsDaily_1d_weighting_function", 3, "mul");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsDaily_2d_weighting_function", 3, "mul");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsDaily_1d_weighting_field", 8, "areaCell");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsDaily_2d_weighting_field", 10, "volumeCell");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsDaily_vertical_mask", 8,
                        "cellMask");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsDaily_vertical_dimension", 11, "nVertLevels");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsWeekly_enable", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsWeekly_compute_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsWeekly_write_on_startup",
                        2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsWeekly_compute_interval",
                        15, "output_interval");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsWeekly_output_stream",
                        25, "regionalStatsWeeklyOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsWeekly_restart_stream",
                        18, "regionalMasksInput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsWeekly_input_stream", 18,
                        "regionalMasksInput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsWeekly_operation", 3,
                        "avg");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsWeekly_region_type", 4,
                        "cell");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsWeekly_region_group", 3,
                        "all");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsWeekly_1d_weighting_function", 3, "mul");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsWeekly_2d_weighting_function", 3, "mul");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsWeekly_1d_weighting_field", 8, "areaCell");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsWeekly_2d_weighting_field", 10, "volumeCell");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsWeekly_vertical_mask", 8,
                        "cellMask");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsWeekly_vertical_dimension", 11, "nVertLevels");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsMonthly_enable", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsMonthly_compute_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsMonthly_write_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsMonthly_compute_interval", 15, "output_interval");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsMonthly_output_stream",
                        26, "regionalStatsMonthlyOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsMonthly_restart_stream",
                        18, "regionalMasksInput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsMonthly_input_stream",
                        18, "regionalMasksInput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsMonthly_operation", 3,
                        "avg");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsMonthly_region_type", 4,
                        "cell");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsMonthly_region_group", 3,
                        "all");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsMonthly_1d_weighting_function", 3, "mul");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsMonthly_2d_weighting_function", 3, "mul");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsMonthly_1d_weighting_field", 8, "areaCell");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsMonthly_2d_weighting_field", 10, "volumeCell");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsMonthly_vertical_mask",
                        8, "cellMask");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsMonthly_vertical_dimension", 11, "nVertLevels");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsCustom_enable", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsCustom_compute_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsCustom_write_on_startup",
                        2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsCustom_compute_interval",
                        15, "output_interval");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsCustom_output_stream",
                        25, "regionalStatsCustomOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsCustom_restart_stream",
                        18, "regionalMasksInput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsCustom_input_stream", 18,
                        "regionalMasksInput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsCustom_operation", 3,
                        "avg");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsCustom_region_type", 4,
                        "cell");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsCustom_region_group", 3,
                        "all");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsCustom_1d_weighting_function", 3, "mul");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsCustom_2d_weighting_function", 3, "mul");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsCustom_1d_weighting_field", 8, "areaCell");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsCustom_2d_weighting_field", 10, "volumeCell");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_regionalStatsCustom_vertical_mask", 8,
                        "cellMask");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_regionalStatsCustom_vertical_dimension", 11, "nVertLevels");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsDaily_enable", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsDaily_compute_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsDaily_write_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsDaily_compute_interval", 17, "00-00-00_01:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsDaily_output_stream",
                        26, "timeSeriesStatsDailyOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsDaily_restart_stream",
                        4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsDaily_operation", 3,
                        "avg");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsDaily_reference_times",
                        12, "initial_time");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsDaily_duration_intervals", 15, "repeat_interval");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsDaily_repeat_intervals", 14, "reset_interval");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsDaily_reset_intervals",
                        17, "00-00-01_00:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsDaily_backward_output_offset", 17,
                        "00-00-01_00:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsMonthly_enable", 3,
                        "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsMonthly_compute_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsMonthly_write_on_startup", 2, "NO");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                      "config_AM_timeSeriesStatsMonthly_compute_interval", 17, "00-00-00_01:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsMonthly_output_stream",
                        28, "timeSeriesStatsMonthlyOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsMonthly_restart_stream", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsMonthly_operation", 3,
                        "avg");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsMonthly_reference_times", 12, "initial_time");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                      "config_AM_timeSeriesStatsMonthly_duration_intervals", 15, "repeat_interval");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsMonthly_repeat_intervals", 14, "reset_interval");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsMonthly_reset_intervals",
                      17, "00-01-00_00:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsMonthly_backward_output_offset", 17,
                        "00-01-00_00:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsClimatology_enable", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsClimatology_compute_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsClimatology_write_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsClimatology_compute_interval", 17,
                        "00-00-00_01:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsClimatology_output_stream", 32,
                        "timeSeriesStatsClimatologyOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsClimatology_restart_stream", 33,
                        "timeSeriesStatsClimatologyRestart");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsClimatology_operation",
                        3, "avg");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsClimatology_reference_times", 71,
                        "00-03-01_00:00:00;00-06-01_00:00:00;00-09-01_00:00:00;00-12-01_00:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsClimatology_duration_intervals", 71,
                        "00-03-00_00:00:00;00-03-00_00:00:00;00-03-00_00:00:00;00-03-00_00:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsClimatology_repeat_intervals", 71,
                        "01-00-00_00:00:00;01-00-00_00:00:00;01-00-00_00:00:00;01-00-00_00:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (
        ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsClimatology_reset_intervals", 79,
        "1000-00-00_00:00:00;1000-00-00_00:00:00;1000-00-00_00:00:00;1000-00-00_00:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsClimatology_backward_output_offset", 17,
                        "00-03-00_00:00:00");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsCustom_enable", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsCustom_compute_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsCustom_write_on_startup", 2, "NO");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsCustom_compute_interval",
                      17, "00-00-00_01:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsCustom_output_stream",
                        27, "timeSeriesStatsCustomOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsCustom_restart_stream",
                        28, "timeSeriesStatsCustomRestart");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_timeSeriesStatsCustom_operation", 3,
                        "avg");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsCustom_reference_times", 12, "initial_time");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                      "config_AM_timeSeriesStatsCustom_duration_intervals", 15, "repeat_interval");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsCustom_repeat_intervals", 14, "reset_interval");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsCustom_reset_intervals", 17, "00-00-07_00:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_timeSeriesStatsCustom_backward_output_offset", 17,
                        "00-00-01_00:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_pointwiseStats_enable", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_pointwiseStats_compute_interval", 15,
                        "output_interval");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_pointwiseStats_output_stream", 20,
                        "pointwiseStatsOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_pointwiseStats_compute_on_startup", 3,
                        "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_pointwiseStats_write_on_startup", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_debugDiagnostics_enable", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_debugDiagnostics_compute_interval", 2,
                        "dt");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_debugDiagnostics_output_stream", 22,
                        "debugDiagnosticsOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_debugDiagnostics_compute_on_startup",
                        3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_debugDiagnostics_write_on_startup", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_debugDiagnostics_check_state", 3,
                        "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_enable", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_compute_on_startup", 3,
                        "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_write_on_startup", 2,
                        "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_compute_interval", 19,
                        "0010-00-00_00:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_output_stream", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_variable_a", 14,
                        "layerThickness");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_variable_b", 8,
                        "areaCell");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_variable_c", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_variable_d", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_variable_e", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_variable_f", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_variable_g", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_variable_h", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_expression_1", 5,
                        "a b *");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_expression_2", 4, "none");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_expression_3", 4, "none");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_expression_4", 4, "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_output_name_1", 10,
                        "volumeCell");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_output_name_2", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_output_name_3", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_rpnCalculator_output_name_4", 4,
                        "none");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_transectTransport_enable", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_transectTransport_compute_interval",
                        15, "output_interval");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_transectTransport_output_stream", 23,
                        "transectTransportOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_transectTransport_compute_on_startup",
                        3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_transectTransport_write_on_startup",
                        2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_transectTransport_transect_group", 3,
                        "all");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_eddyProductVariables_enable", 3, "YES");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_eddyProductVariables_compute_interval",
                      19, "0000-00-00_01:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_eddyProductVariables_output_stream",
                        26, "eddyProductVariablesOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_eddyProductVariables_compute_on_startup", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_eddyProductVariables_write_on_startup", 2, "NO");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mocStreamfunction_enable", 3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mocStreamfunction_compute_interval",
                        19, "0000-00-00_01:00:00");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mocStreamfunction_output_stream", 23,
                        "mocStreamfunctionOutput");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mocStreamfunction_compute_on_startup",
                        3, "YES");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mocStreamfunction_write_on_startup",
                        2, "NO");
    CHECK_ERR
    dattr = -1.e+34;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mocStreamfunction_min_bin", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    dattr = -1.e+34;
    err = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mocStreamfunction_max_bin", MPI_DOUBLE, 1,
                   &dattr);
    CHECK_ERR
    iattr = 180;
    err   = GET_ATT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mocStreamfunction_num_bins", MPI_INT, 1,
                   &iattr);
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                      "config_AM_mocStreamfunction_vertical_velocity_value", 15, "vertVelocityTop");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR,
                        "config_AM_mocStreamfunction_normal_velocity_value", 14, "normalVelocity");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mocStreamfunction_region_group", 3,
                        "all");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "config_AM_mocStreamfunction_transect_group", 3,
                        "all");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, E3SM_IO_GLOBAL_ATTR, "file_id", 10, "3q7jdjett8");
    CHECK_ERR

err_out:
    return err;
}

/*----< def_G_case_h0() >----------------------------------------------------*/
int inq_G_case_h0 (e3sm_io_driver &driver,
                   int ncid,              /* file ID */
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

    int i, err, dimids[3];
    int dim_nVertLevelsP1, dim_nCells, dim_Time, dim_nVertLevels, dim_nEdges, dim_nVertices,
        dim_StrLen;

    err = inq_global_attributes (driver, ncid);
    CHECK_ERR

    /* define dimensions */
    /*
    err = driver.def_dim(ncid, "nCells", dims_D1[0], &dim_nCells); CHECK_ERR
    err = driver.def_dim(ncid, "Time", NC_UNLIMITED, &dim_Time); CHECK_ERR
    err = driver.def_dim(ncid, "nVertLevelsP1", dims_D6[1], &dim_nVertLevelsP1); CHECK_ERR
    err = driver.def_dim(ncid, "nVertLevels", dims_D3[1], &dim_nVertLevels); CHECK_ERR
    err = driver.def_dim(ncid, "nEdges", dims_D2[0], &dim_nEdges); CHECK_ERR
    err = driver.def_dim(ncid, "nVertices", dims_D5[0], &dim_nVertices); CHECK_ERR
    err = driver.def_dim(ncid, "StrLen", 64, &dim_StrLen); CHECK_ERR
    */

    err = driver.inq_dim (ncid, "nCells", &dim_nCells);
    CHECK_ERR
    err = driver.inq_dim (ncid, "nVertLevelsP1", &dim_nVertLevelsP1);
    CHECK_ERR
    err = driver.inq_dim (ncid, "nVertLevels", &dim_nVertLevels);
    CHECK_ERR
    err = driver.inq_dim (ncid, "nEdges", &dim_nEdges);
    CHECK_ERR
    err = driver.inq_dim (ncid, "nVertices", &dim_nVertices);
    CHECK_ERR

    err = driver.inq_dimlen (ncid, dim_nCells, &(((MPI_Offset *)dims_D1)[0]));
    CHECK_ERR
    err = driver.inq_dimlen (ncid, dim_nVertLevelsP1, &(((MPI_Offset *)dims_D6)[1]));
    CHECK_ERR
    err = driver.inq_dimlen (ncid, dim_nVertLevels, &(((MPI_Offset *)dims_D3)[1]));
    CHECK_ERR
    err = driver.inq_dimlen (ncid, dim_nEdges, &(((MPI_Offset *)dims_D2)[0]));
    CHECK_ERR
    err = driver.inq_dimlen (ncid, dim_nVertices, &(((MPI_Offset *)dims_D5)[0]));
    CHECK_ERR

    i = 0;

    /* define variables */
    /* 1 double (Time, nCells) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;

    err = INQ_VID (ncid, "salinitySurfaceRestoringTendency", MPI_DOUBLE, 2, dimids,
                   &salinitySurfaceRestoringTendency);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, salinitySurfaceRestoringTendency, "units", 7, "m PSU/s");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, salinitySurfaceRestoringTendency, "long_name", 42,
                        "salinity tendency due to surface restoring");
    CHECK_ERR
    varids[i++] = salinitySurfaceRestoringTendency;

    /* 3 double (Time, nCells, nVertLevelsP1) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevelsP1;

    err = INQ_VID (ncid, "vertTransportVelocityTop", MPI_DOUBLE, 3, dimids,
                   &vertTransportVelocityTop);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, vertTransportVelocityTop, "units", 8, "m s^{-1}");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, vertTransportVelocityTop, "long_name", 280,
                        "vertical tracer-transport velocity defined at center (horizontally) "
                        "and top (vertically) of cell.  This is not the vertical ALE transport, "
                        "but is Eulerian (fixed-frame) in the vertical, and computed from the "
                        "continuity equation from the horizontal total tracer-transport velocity.");
    CHECK_ERR
    varids[i++] = vertTransportVelocityTop;

    err = INQ_VID (ncid, "vertGMBolusVelocityTop", MPI_DOUBLE, 3, dimids, &vertGMBolusVelocityTop);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, vertGMBolusVelocityTop, "units", 8, "m s^{-1}");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, vertGMBolusVelocityTop, "long_name", 266,
                        "vertical tracer-transport velocity defined at center (horizontally) "
                        "and top (vertically) of cell.  This is not the vertical ALE transport, "
                        "but is Eulerian (fixed-frame) in the vertical, and computed from the "
                        "continuity equation from the horizontal GM Bolus velocity.");
    CHECK_ERR
    varids[i++] = vertGMBolusVelocityTop;

    err = INQ_VID (ncid, "vertAleTransportTop", MPI_DOUBLE, 3, dimids, &vertAleTransportTop);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, vertAleTransportTop, "units", 8, "m s^{-1}");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, vertAleTransportTop, "long_name", 69,
                        "vertical transport through "
                        "the layer interface at the top of the cell");
    CHECK_ERR
    varids[i++] = vertAleTransportTop;

    /* 1 double (Time, nCells) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;

    err = INQ_VID (ncid, "tendSSH", MPI_DOUBLE, 2, dimids, &tendSSH);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, tendSSH, "units", 8, "m s^{-1}");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, tendSSH, "long_name", 35, "time tendency of sea-surface height");
    CHECK_ERR
    varids[i++] = tendSSH;

    /* 1 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "layerThickness", MPI_DOUBLE, 3, dimids, &layerThickness);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, layerThickness, "units", 1, "m");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, layerThickness, "long_name", 15, "layer thickness");
    CHECK_ERR
    varids[i++] = layerThickness;

    /* 1 double (Time, nEdges, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nEdges;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "normalVelocity", MPI_DOUBLE, 3, dimids, &normalVelocity);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, normalVelocity, "units", 8, "m s^{-1}");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, normalVelocity, "long_name", 47,
                        "horizonal velocity, "
                        "normal component to an edge");
    CHECK_ERR
    varids[i++] = normalVelocity;

    /* 1 double (Time, nCells) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;

    err = INQ_VID (ncid, "ssh", MPI_DOUBLE, 2, dimids, &ssh);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, ssh, "units", 1, "m");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, ssh, "long_name", 18, "sea surface height");
    CHECK_ERR
    varids[i++] = ssh;

    /* 1 int (nEdges) */
    dimids[0] = dim_nEdges;

    err = INQ_VID (ncid, "maxLevelEdgeTop", MPI_INT, 1, dimids, &maxLevelEdgeTop);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, maxLevelEdgeTop, "units", 8, "unitless");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, maxLevelEdgeTop, "long_name", 79,
                        "Index to the last edge "
                        "in a column with active ocean cells on both sides of it.");
    CHECK_ERR
    varids[i++] = maxLevelEdgeTop;

    /* 1 double (nVertLevels) */
    dimids[0] = dim_nVertLevels;

    err = INQ_VID (ncid, "vertCoordMovementWeights", MPI_DOUBLE, 1, dimids,
                   &vertCoordMovementWeights);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, vertCoordMovementWeights, "units", 8, "unitless");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, vertCoordMovementWeights, "long_name", 98,
                        "Weights used "
                        "for distribution of sea surface heigh purturbations through "
                        "multiple vertical levels.");
    CHECK_ERR
    varids[i++] = vertCoordMovementWeights;

    /* 1 int (nEdges, nVertLevels) */
    dimids[0] = dim_nEdges;
    dimids[1] = dim_nVertLevels;

    err = INQ_VID (ncid, "edgeMask", MPI_INT, 2, dimids, &edgeMask);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, edgeMask, "units", 8, "unitless");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, edgeMask, "long_name", 69,
                        "Mask on edges that determines "
                        "if computations should be done on edge.");
    CHECK_ERR
    varids[i++] = edgeMask;

    /* 1 int (nCells, nVertLevels) */
    dimids[0] = dim_nCells;
    dimids[1] = dim_nVertLevels;

    err = INQ_VID (ncid, "cellMask", MPI_INT, 2, dimids, &cellMask);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, cellMask, "units", 8, "unitless");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, cellMask, "long_name", 69,
                        "Mask on cells that determines "
                        "if computations should be done on cell.");
    CHECK_ERR
    varids[i++] = cellMask;

    /* 1 int (nVertices, nVertLevels) */
    dimids[0] = dim_nVertices;
    dimids[1] = dim_nVertLevels;

    err = INQ_VID (ncid, "vertexMask", MPI_INT, 2, dimids, &vertexMask);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, vertexMask, "units", 8, "unitless");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, vertexMask, "long_name", 75,
                        "Mask on vertices that determines "
                        "if computations should be done on vertice.");
    CHECK_ERR
    varids[i++] = vertexMask;

    /* 2 double (nVertLevels) */
    dimids[0] = dim_nVertLevels;

    err = INQ_VID (ncid, "refZMid", MPI_DOUBLE, 1, dimids, &refZMid);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, refZMid, "units", 1, "m");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, refZMid, "long_name", 87,
                        "Reference mid z-coordinate of ocean "
                        "for each vertical level. This has a negative value.");
    CHECK_ERR
    varids[i++] = refZMid;

    err = INQ_VID (ncid, "refLayerThickness", MPI_DOUBLE, 1, dimids, &refLayerThickness);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, refLayerThickness, "units", 1, "m");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, refLayerThickness, "long_name", 58,
                        "Reference layerThickness "
                        "of ocean for each vertical level.");
    CHECK_ERR
    varids[i++] = refLayerThickness;

    /* 1 char (Time, StrLen) */
    dimids[0] = dim_Time;
    dimids[1] = dim_StrLen;

    err = INQ_VID (ncid, "xtime", MPI_CHAR, 2, dimids, &xtime);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, xtime, "units", 8, "unitless");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, xtime, "long_name", 45,
                        "model time, with format \'YYYY-MM-DD_HH:MM:SS\'");
    CHECK_ERR
    varids[i++] = xtime;

    /* 2 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "kineticEnergyCell", MPI_DOUBLE, 3, dimids, &kineticEnergyCell);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, kineticEnergyCell, "units", 10, "m^2 s^{-2}");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, kineticEnergyCell, "long_name", 45,
                        "kinetic energy of horizonal "
                        "velocity on cells");
    CHECK_ERR
    varids[i++] = kineticEnergyCell;

    err = INQ_VID (ncid, "relativeVorticityCell", MPI_DOUBLE, 3, dimids, &relativeVorticityCell);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, relativeVorticityCell, "units", 6, "s^{-1}");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, relativeVorticityCell, "long_name", 67,
                        "curl of horizontal velocity, "
                        "averaged from vertices to cell centers");
    CHECK_ERR
    varids[i++] = relativeVorticityCell;

    /* 1 double (Time, nVertices, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nVertices;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "relativeVorticity", MPI_DOUBLE, 3, dimids, &relativeVorticity);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, relativeVorticity, "units", 6, "s^{-1}");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, relativeVorticity, "long_name", 48,
                        "curl of horizontal velocity, "
                        "defined at vertices");
    CHECK_ERR
    varids[i++] = relativeVorticity;

    /* 1 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "divergence", MPI_DOUBLE, 3, dimids, &divergence);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, divergence, "units", 6, "s^{-1}");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, divergence, "long_name", 32, "divergence of horizonal velocity");
    CHECK_ERR
    varids[i++] = divergence;

    /* 6 double (Time) */
    dimids[0] = dim_Time;

    err = INQ_VID (ncid, "areaCellGlobal", MPI_DOUBLE, 1, dimids, &areaCellGlobal);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, areaCellGlobal, "units", 3, "m^2");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, areaCellGlobal, "long_name", 86,
                        "sum of the areaCell variable over "
                        "the full domain, used to normalize global statistics");
    CHECK_ERR
    varids[i++] = areaCellGlobal;

    err = INQ_VID (ncid, "areaEdgeGlobal", MPI_DOUBLE, 1, dimids, &areaEdgeGlobal);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, areaEdgeGlobal, "units", 3, "m^2");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, areaEdgeGlobal, "long_name", 86,
                        "sum of the areaEdge variable over "
                        "the full domain, used to normalize global statistics");
    CHECK_ERR
    varids[i++] = areaEdgeGlobal;

    err = INQ_VID (ncid, "areaTriangleGlobal", MPI_DOUBLE, 1, dimids, &areaTriangleGlobal);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, areaTriangleGlobal, "units", 3, "m^2");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, areaTriangleGlobal, "long_name", 90,
                        "sum of the areaTriangle variable "
                        "over the full domain, used to normalize global statistics");
    CHECK_ERR
    varids[i++] = areaTriangleGlobal;

    err = INQ_VID (ncid, "volumeCellGlobal", MPI_DOUBLE, 1, dimids, &volumeCellGlobal);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, volumeCellGlobal, "units", 3, "m^3");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, volumeCellGlobal, "long_name", 88,
                        "sum of the volumeCell variable over "
                        "the full domain, used to normalize global statistics");
    CHECK_ERR
    varids[i++] = volumeCellGlobal;

    err = INQ_VID (ncid, "volumeEdgeGlobal", MPI_DOUBLE, 1, dimids, &volumeEdgeGlobal);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, volumeEdgeGlobal, "units", 3, "m^3");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, volumeEdgeGlobal, "long_name", 88,
                        "sum of the volumeEdge variable over "
                        "the full domain, used to normalize global statistics");
    CHECK_ERR
    varids[i++] = volumeEdgeGlobal;

    err = INQ_VID (ncid, "CFLNumberGlobal", MPI_DOUBLE, 1, dimids, &CFLNumberGlobal);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, CFLNumberGlobal, "units", 8, "unitless");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, CFLNumberGlobal, "long_name", 39,
                        "maximum CFL number over the full domain");
    CHECK_ERR
    varids[i++] = CFLNumberGlobal;

    /* 1 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "BruntVaisalaFreqTop", MPI_DOUBLE, 3, dimids, &BruntVaisalaFreqTop);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, BruntVaisalaFreqTop, "units", 6, "s^{-2}");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, BruntVaisalaFreqTop, "long_name", 89,
                        "Brunt Vaisala frequency defined at "
                        "the center (horizontally) and top (vertically) of cell");
    CHECK_ERR
    varids[i++] = BruntVaisalaFreqTop;

    /* 1 double (Time, nCells, nVertLevelsP1) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevelsP1;

    err = INQ_VID (ncid, "vertVelocityTop", MPI_DOUBLE, 3, dimids, &vertVelocityTop);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, vertVelocityTop, "units", 8, "m s^{-1}");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, vertVelocityTop, "long_name", 79,
                        "vertical velocity defined at center "
                        "(horizontally) and top (vertically) of cell");
    CHECK_ERR
    varids[i++] = vertVelocityTop;

    /* 5 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "velocityZonal", MPI_DOUBLE, 3, dimids, &velocityZonal);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, velocityZonal, "units", 8, "m s^{-1}");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, velocityZonal, "long_name", 58,
                        "component of horizontal velocity in "
                        "the eastward direction");
    CHECK_ERR
    varids[i++] = velocityZonal;

    err = INQ_VID (ncid, "velocityMeridional", MPI_DOUBLE, 3, dimids, &velocityMeridional);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, velocityMeridional, "units", 8, "m s^{-1}");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, velocityMeridional, "long_name", 59,
                        "component of horizontal velocity in "
                        "the northward direction");
    CHECK_ERR
    varids[i++] = velocityMeridional;

    err = INQ_VID (ncid, "displacedDensity", MPI_DOUBLE, 3, dimids, &displacedDensity);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, displacedDensity, "units", 9, "kg m^{-3}");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, displacedDensity, "long_name", 130,
                      "Density displaced adiabatically to "
                      "the mid-depth one layer deeper.  That is, layer k has been displaced to the "
                      "depth of layer k+1.");
    CHECK_ERR
    varids[i++] = displacedDensity;

    err = INQ_VID (ncid, "potentialDensity", MPI_DOUBLE, 3, dimids, &potentialDensity);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, potentialDensity, "units", 9, "kg m^{-3}");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, potentialDensity, "long_name", 80,
                        "potential density: density displaced "
                        "adiabatically to the mid-depth of top layer");
    CHECK_ERR
    varids[i++] = potentialDensity;

    err = INQ_VID (ncid, "pressure", MPI_DOUBLE, 3, dimids, &pressure);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, pressure, "units", 8, "N m^{-2}");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, pressure, "long_name", 38, "pressure used in the momentum equation");
    CHECK_ERR
    varids[i++] = pressure;

    /* 1 double (nVertLevels) */
    dimids[0] = dim_nVertLevels;

    err = INQ_VID (ncid, "refBottomDepth", MPI_DOUBLE, 1, dimids, &refBottomDepth);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, refBottomDepth, "units", 1, "m");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, refBottomDepth, "long_name", 78,
                        "Reference depth of ocean for each "
                        "vertical level. Used in \'z-level\' type runs.");
    CHECK_ERR
    varids[i++] = refBottomDepth;

    /* 1 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "zMid", MPI_DOUBLE, 3, dimids, &zMid);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, zMid, "units", 1, "m");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, zMid, "long_name", 42, "z-coordinate of the mid-depth of the layer");
    CHECK_ERR
    varids[i++] = zMid;

    /* 1 double (nCells) */
    dimids[0] = dim_nCells;

    err = INQ_VID (ncid, "bottomDepth", MPI_DOUBLE, 1, dimids, &bottomDepth);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, bottomDepth, "units", 1, "m");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, bottomDepth, "long_name", 78,
                        "Depth of the bottom of the ocean. Given "
                        "as a positive distance from sea level.");
    CHECK_ERR
    varids[i++] = bottomDepth;

    /* 1 int (nCells) */
    dimids[0] = dim_nCells;

    err = INQ_VID (ncid, "maxLevelCell", MPI_INT, 1, dimids, &maxLevelCell);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, maxLevelCell, "units", 8, "unitless");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, maxLevelCell, "long_name", 51,
                        "Index to the last active ocean cell in each column.");
    CHECK_ERR
    varids[i++] = maxLevelCell;

    /* 1 int (nEdges) */
    dimids[0] = dim_nEdges;

    err = INQ_VID (ncid, "maxLevelEdgeBot", MPI_INT, 1, dimids, &maxLevelEdgeBot);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, maxLevelEdgeBot, "units", 8, "unitless");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, maxLevelEdgeBot, "long_name", 92,
                        "Index to the last edge in a column with at "
                        "least one active ocean cell on either side of it.");
    CHECK_ERR
    varids[i++] = maxLevelEdgeBot;

    /* 1 double (Time, nCells) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;

    err = INQ_VID (ncid, "columnIntegratedSpeed", MPI_DOUBLE, 2, dimids, &columnIntegratedSpeed);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, columnIntegratedSpeed, "units", 10, "m^2 s^{-1}");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, columnIntegratedSpeed, "long_name", 109,
                      "speed = sum(h*sqrt(2*ke)), where ke "
                      "is kineticEnergyCell and the sum is over the full column at cell centers.");
    CHECK_ERR
    varids[i++] = columnIntegratedSpeed;

    /* 13 double (Time, nCells, nVertLevels) */
    dimids[0] = dim_Time;
    dimids[1] = dim_nCells;
    dimids[2] = dim_nVertLevels;

    err = INQ_VID (ncid, "temperatureHorizontalAdvectionTendency", MPI_DOUBLE, 3, dimids,
                   &temperatureHorizontalAdvectionTendency);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, temperatureHorizontalAdvectionTendency, "long_name", 58,
                        "potential temperature "
                        "tendency due to horizontal advection");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, temperatureHorizontalAdvectionTendency, "units", 26,
                        "degrees Celsius per second");
    CHECK_ERR
    varids[i++] = temperatureHorizontalAdvectionTendency;

    err = INQ_VID (ncid, "salinityHorizontalAdvectionTendency", MPI_DOUBLE, 3, dimids,
                   &salinityHorizontalAdvectionTendency);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, salinityHorizontalAdvectionTendency, "long_name", 45,
                        "salinity tendency due "
                        "to horizontal advection");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, salinityHorizontalAdvectionTendency, "units", 14, "PSU per second");
    CHECK_ERR
    varids[i++] = salinityHorizontalAdvectionTendency;

    err = INQ_VID (ncid, "temperatureVerticalAdvectionTendency", MPI_DOUBLE, 3, dimids,
                   &temperatureVerticalAdvectionTendency);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, temperatureVerticalAdvectionTendency, "long_name", 56,
                        "potential temperature "
                        "tendency due to vertical advection");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, temperatureVerticalAdvectionTendency, "units", 26,
                        "degrees Celsius per second");
    CHECK_ERR
    varids[i++] = temperatureVerticalAdvectionTendency;

    err = INQ_VID (ncid, "salinityVerticalAdvectionTendency", MPI_DOUBLE, 3, dimids,
                   &salinityVerticalAdvectionTendency);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, salinityVerticalAdvectionTendency, "long_name", 43,
                        "salinity tendency due "
                        "to vertical advection");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, salinityVerticalAdvectionTendency, "units", 14, "PSU per second");
    CHECK_ERR
    varids[i++] = salinityVerticalAdvectionTendency;

    err = INQ_VID (ncid, "temperatureVertMixTendency", MPI_DOUBLE, 3, dimids,
                   &temperatureVertMixTendency);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, temperatureVertMixTendency, "long_name", 53,
                        "potential temperature tendency "
                        "due to vertical mixing");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, temperatureVertMixTendency, "units", 26, "degrees Celsius per second");
    CHECK_ERR
    varids[i++] = temperatureVertMixTendency;

    err =
        INQ_VID (ncid, "salinityVertMixTendency", MPI_DOUBLE, 3, dimids, &salinityVertMixTendency);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, salinityVertMixTendency, "long_name", 40,
                        "salinity tendency due to vertical mixing");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, salinityVertMixTendency, "units", 14, "PSU per second");
    CHECK_ERR
    varids[i++] = salinityVertMixTendency;

    err = INQ_VID (ncid, "temperatureSurfaceFluxTendency", MPI_DOUBLE, 3, dimids,
                   &temperatureSurfaceFluxTendency);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, temperatureSurfaceFluxTendency, "long_name", 52,
                        "potential temperature tendency "
                        "due to surface fluxes");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, temperatureSurfaceFluxTendency, "units", 26,
                        "degrees Celsius per second");
    CHECK_ERR
    varids[i++] = temperatureSurfaceFluxTendency;

    err = INQ_VID (ncid, "salinitySurfaceFluxTendency", MPI_DOUBLE, 3, dimids,
                   &salinitySurfaceFluxTendency);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, salinitySurfaceFluxTendency, "long_name", 39,
                        "salinity tendency due to surface fluxes");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, salinitySurfaceFluxTendency, "units", 14, "PSU per second");
    CHECK_ERR
    varids[i++] = salinitySurfaceFluxTendency;

    err = INQ_VID (ncid, "temperatureShortWaveTendency", MPI_DOUBLE, 3, dimids,
                   &temperatureShortWaveTendency);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, temperatureShortWaveTendency, "units", 26,
                        "degrees Celsius per second");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, temperatureShortWaveTendency, "long_name", 59,
                        "potential temperature tendency due "
                        "to penetrating shortwave");
    CHECK_ERR
    varids[i++] = temperatureShortWaveTendency;

    err = INQ_VID (ncid, "temperatureNonLocalTendency", MPI_DOUBLE, 3, dimids,
                   &temperatureNonLocalTendency);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, temperatureNonLocalTendency, "long_name", 56,
                        "potential temperature tendency due "
                        "to kpp non-local flux");
    CHECK_ERR
    err =
        GET_ATT_TEXT (ncid, temperatureNonLocalTendency, "units", 26, "degrees Celsius per second");
    CHECK_ERR
    varids[i++] = temperatureNonLocalTendency;

    err = INQ_VID (ncid, "salinityNonLocalTendency", MPI_DOUBLE, 3, dimids,
                   &salinityNonLocalTendency);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, salinityNonLocalTendency, "long_name", 43,
                        "salinity tendency due to kpp non-local flux");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, salinityNonLocalTendency, "units", 14, "PSU per second");
    CHECK_ERR
    varids[i++] = salinityNonLocalTendency;

    err = INQ_VID (ncid, "temperature", MPI_DOUBLE, 3, dimids, &temperature);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, temperature, "long_name", 21, "potential temperature");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, temperature, "units", 15, "degrees Celsius");
    CHECK_ERR
    varids[i++] = temperature;

    err = INQ_VID (ncid, "salinity", MPI_DOUBLE, 3, dimids, &salinity);
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, salinity, "long_name", 8, "salinity");
    CHECK_ERR
    err = GET_ATT_TEXT (ncid, salinity, "units", 32, "grams salt per kilogram seawater");
    CHECK_ERR
    varids[i++] = salinity;

    assert (i == nvars);

err_out:
    return err;
}
