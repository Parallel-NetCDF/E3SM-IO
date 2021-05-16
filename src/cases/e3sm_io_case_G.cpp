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
#include <e3sm_io_case_G.hpp>
#include <e3sm_io.h>
#include <e3sm_io_case.hpp>
#include <e3sm_io_err.h>

e3sm_io_case_G::e3sm_io_case_G () {}

e3sm_io_case_G::~e3sm_io_case_G () {
    if (this->D1_rec_dbl_buf != NULL) { free (this->D1_rec_dbl_buf); }
    if (this->D3_rec_dbl_buf != NULL) { free (this->D3_rec_dbl_buf); }
    if (this->D4_rec_dbl_buf != NULL) { free (this->D4_rec_dbl_buf); }
    if (this->D5_rec_dbl_buf != NULL) { free (this->D5_rec_dbl_buf); }
    if (this->D6_rec_dbl_buf != NULL) { free (this->D6_rec_dbl_buf); }
    if (this->D1_fix_dbl_buf != NULL) { free (this->D1_fix_dbl_buf); }
    if (this->D1_fix_int_buf != NULL) { free (this->D1_fix_int_buf); }
    if (this->D2_fix_int_buf != NULL) { free (this->D2_fix_int_buf); }
    if (this->D3_fix_int_buf != NULL) { free (this->D3_fix_int_buf); }
    if (this->D4_fix_int_buf != NULL) { free (this->D4_fix_int_buf); }
    if (this->D5_fix_int_buf != NULL) { free (this->D5_fix_int_buf); }
}

int e3sm_io_case_G::wr_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver) {
    int err, nerrs = 0;

    PRINT_MSG (0, "number of requests for D1=%d D2=%d D3=%d D4=%d D5=%d D6=%d\n",
               decom.contig_nreqs[0], decom.contig_nreqs[1], decom.contig_nreqs[2],
               decom.contig_nreqs[3], decom.contig_nreqs[4], decom.contig_nreqs[5]);

    if ((cfg.verbose >= 0) && (cfg.rank == 0)) {
        printf ("Total number of MPI processes      = %d\n", cfg.np);
        printf ("Number of IO processes             = %d\n", cfg.num_iotasks);
        printf ("Input decomposition file           = %s\n", cfg.cfgpath);
        printf ("Number of decompositions           = %d\n", decom.num_decomp);
        printf ("Output file directory              = %s\n", cfg.targetdir);
        printf ("Variable dimensions (C order)      = %lld x %lld\n", decom.dims[2][0],
                decom.dims[2][1]);
        printf ("Write number of records (time dim) = %d\n", cfg.nrec);
        printf ("Using noncontiguous write buffer   = %s\n", cfg.non_contig_buf ? "yes" : "no");
    }

    /* vard APIs require internal data type matches external one */
    if (cfg.vard) {
#if REC_XTYPE != NC_FLOAT
        RET_ERR ("PnetCDF vard API requires internal and external data types match, skip\n");
#endif
        RET_ERR ("Low level API not supported in g case\n");
    } else {
        PRINT_MSG (0,
                   "\n==== benchmarking G case writing using varn API ========================\n");

        if (cfg.two_buf) {
            PRINT_MSG (0, "Variable written order: 2D variables then 3D variables\n\n");
        } else {
            PRINT_MSG (0, "Variable written order: same as variables are defined\n\n");
        }
        fflush (stdout);

        MPI_Barrier (cfg.io_comm);
        nerrs += run_varn_G_case (cfg, decom, driver, "g_case_hist_varn.nc", this->D1_fix_int_buf,
                                  this->D2_fix_int_buf, this->D3_fix_int_buf, this->D4_fix_int_buf,
                                  this->D5_fix_int_buf, this->D1_rec_dbl_buf, this->D3_rec_dbl_buf,
                                  this->D4_rec_dbl_buf, this->D5_rec_dbl_buf, this->D6_rec_dbl_buf,
                                  this->D1_fix_dbl_buf);
    }

err_out:;
    return nerrs;
}

int e3sm_io_case_G::rd_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver) {
    int err, nerrs = 0;

    PRINT_MSG (0, "number of requests for D1=%d D2=%d D3=%d D4=%d D5=%d D6=%d\n",
               decom.contig_nreqs[0], decom.contig_nreqs[1], decom.contig_nreqs[2],
               decom.contig_nreqs[3], decom.contig_nreqs[4], decom.contig_nreqs[5]);

    if ((cfg.verbose >= 0) && (cfg.rank == 0)) {
        printf ("Total number of MPI processes      = %d\n", cfg.np);
        printf ("Number of IO processes             = %d\n", cfg.num_iotasks);
        printf ("Input decomposition file           = %s\n", cfg.cfgpath);
        printf ("Number of decompositions           = %d\n", decom.num_decomp);
        printf ("Output file directory              = %s\n", cfg.targetdir);
        printf ("Variable dimensions (C order)      = %lld x %lld\n", decom.dims[2][0],
                decom.dims[2][1]);
        printf ("Write number of records (time dim) = %d\n", cfg.nrec);
        printf ("Using noncontiguous write buffer   = %s\n", cfg.non_contig_buf ? "yes" : "no");
    }

    /* vard APIs require internal data type matches external one */
    if (cfg.vard) {
#if REC_XTYPE != NC_FLOAT
        RET_ERR ("PnetCDF vard API requires internal and external data types match, skip\n");
#endif
        RET_ERR ("Low level API not supported in g case\n");
    } else {
        PRINT_MSG (0,
                   "\n==== benchmarking G case writing using varn API ========================\n");

        if (cfg.two_buf) {
            PRINT_MSG (0, "Variable written order: 2D variables then 3D variables\n\n");
        } else {
            PRINT_MSG (0, "Variable written order: same as variables are defined\n\n");
        }
        fflush (stdout);

        MPI_Barrier (cfg.io_comm);
        nerrs += run_varn_G_case_rd (
            cfg, decom, driver, "g_case_hist_varn.nc", &(this->D1_fix_int_buf),
            &(this->D2_fix_int_buf), &(this->D3_fix_int_buf), &(this->D4_fix_int_buf),
            &(this->D5_fix_int_buf), &(this->D1_rec_dbl_buf), &(this->D3_rec_dbl_buf),
            &(this->D4_rec_dbl_buf), &(this->D5_rec_dbl_buf), &(this->D6_rec_dbl_buf),
            &(this->D1_fix_dbl_buf));
    }
err_out:;
    return nerrs;
}

int e3sm_io_case_G::load_data (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver) {
    int err, nerrs = 0;
    int verbose;
    char *tmp;

    verbose     = cfg.verbose;
    cfg.verbose = -1;  // Disable output
    // Swap datadir with targetdir temporarily
    tmp           = cfg.targetdir;
    cfg.targetdir = cfg.datadir;
    cfg.datadir   = tmp;

    // Run dummy G case read for data
    MPI_Barrier (cfg.io_comm);
    nerrs += run_varn_G_case_rd (
        cfg, decom, driver, "g_case_hist_varn.nc", &(this->D1_fix_int_buf), &(this->D2_fix_int_buf),
        &(this->D3_fix_int_buf), &(this->D4_fix_int_buf), &(this->D5_fix_int_buf),
        &(this->D1_rec_dbl_buf), &(this->D3_rec_dbl_buf), &(this->D4_rec_dbl_buf),
        &(this->D5_rec_dbl_buf), &(this->D6_rec_dbl_buf), &(this->D1_fix_dbl_buf));

    cfg.verbose = verbose;
    // Swap datadir and targetdir back
    tmp           = cfg.targetdir;
    cfg.targetdir = cfg.datadir;
    cfg.datadir   = tmp;

err_out:;
    return nerrs;
}
