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

int e3sm_io_case_G::wr_test(e3sm_io_config &cfg,
                            e3sm_io_decom  &decom,
                            e3sm_io_driver &driver)
{
    int err;
    
    PRINT_MSG (0, "number of requests for D1=%d D2=%d D3=%d D4=%d D5=%d D6=%d\n",
               decom.contig_nreqs[0], decom.contig_nreqs[1], decom.contig_nreqs[2],
               decom.contig_nreqs[3], decom.contig_nreqs[4], decom.contig_nreqs[5]);

    if ((cfg.verbose >= 0) && (cfg.rank == 0)) {
        printf ("Total number of MPI processes      = %d\n", cfg.np);
        printf ("Number of IO processes             = %d\n", cfg.num_iotasks);
        printf ("Input decomposition file           = %s\n", cfg.cfg_path);
        printf ("Number of decompositions           = %d\n", decom.num_decomp);
        printf ("Output file/directory              = %s\n", cfg.out_path);
        printf ("Variable dimensions (C order)      = %lld x %lld\n",
                decom.dims[2][0], decom.dims[2][1]);
        printf ("Write number of records (time dim) = %d\n", cfg.nrec);
        printf ("Using noncontiguous write buffer   = %s\n", cfg.non_contig_buf ? "yes" : "no");
        printf("==== Benchmarking G case ====\n");
        if (cfg.strategy == blob)
            printf("Using one-file-per-node blob I/O strategy\n");
        else if (!cfg.vard)
            printf("Using varn API\n");
        printf("Variable written order: same as variables are defined\n");
    }

    if (cfg.api == pnetcdf && cfg.strategy == blob) {
        /* construct metadata for blob I/O strategy */
        err = blob_metadata(&cfg, &decom);
        CHECK_ERR

        /* Use one-file-per-compute-node blob I/O strategy */
        err = pnetcdf_blob_G_case(cfg, decom, driver);
        CHECK_ERR

        if (cfg.sub_comm != MPI_COMM_NULL)
            MPI_Comm_free(&cfg.sub_comm);
    }
    else if (cfg.vard) {
        /* vard APIs require internal data type matches external one */
#if REC_XTYPE != NC_FLOAT
        ERR_OUT ("PnetCDF vard API requires internal and external data types match, skip\n");
#endif
        ERR_OUT ("Low level API not supported in g case\n");
    } else {
        err = run_varn_G_case (cfg, decom, driver, this->D1_fix_int_buf,
                                  this->D2_fix_int_buf, this->D3_fix_int_buf, this->D4_fix_int_buf,
                                  this->D5_fix_int_buf, this->D1_rec_dbl_buf, this->D3_rec_dbl_buf,
                                  this->D4_rec_dbl_buf, this->D5_rec_dbl_buf, this->D6_rec_dbl_buf,
                                  this->D1_fix_dbl_buf);
        CHECK_ERR
    }

err_out:
    return err;
}

int e3sm_io_case_G::rd_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver) {
    int err;

    PRINT_MSG (0, "number of requests for D1=%d D2=%d D3=%d D4=%d D5=%d D6=%d\n",
               decom.contig_nreqs[0], decom.contig_nreqs[1], decom.contig_nreqs[2],
               decom.contig_nreqs[3], decom.contig_nreqs[4], decom.contig_nreqs[5]);

    if ((cfg.verbose >= 0) && (cfg.rank == 0)) {
        printf ("Total number of MPI processes      = %d\n", cfg.np);
        printf ("Number of IO processes             = %d\n", cfg.num_iotasks);
        printf ("Input decomposition file           = %s\n", cfg.cfg_path);
        printf ("Number of decompositions           = %d\n", decom.num_decomp);
        printf ("Input file/directory               = %s\n", cfg.in_path);
        printf ("Variable dimensions (C order)      = %lld x %lld\n", decom.dims[2][0],
                decom.dims[2][1]);
        printf ("Read number of records (time dim)  = %d\n", cfg.nrec);
        printf ("Using noncontiguous read buffer    = %s\n", cfg.non_contig_buf ? "yes" : "no");
    }

    /* vard APIs require internal data type matches external one */
    if (cfg.vard) {
#if REC_XTYPE != NC_FLOAT
        ERR_OUT ("PnetCDF vard API requires internal and external data types match, skip\n");
#endif
        ERR_OUT ("Low level API not supported in g case\n");
    } else {
        PRINT_MSG (0,
                   "\n==== benchmarking G case read using varn API ========================\n");

        if (cfg.two_buf) {
            PRINT_MSG (0, "Variable read order: 2D variables then 3D variables\n\n");
        } else {
            PRINT_MSG (0, "Variable read order: same as variables are defined\n\n");
        }
        fflush (stdout);

        MPI_Barrier (cfg.io_comm);
        err = run_varn_G_case_rd (
            cfg, decom, driver, &(this->D1_fix_int_buf),
            &(this->D2_fix_int_buf), &(this->D3_fix_int_buf), &(this->D4_fix_int_buf),
            &(this->D5_fix_int_buf), &(this->D1_rec_dbl_buf), &(this->D3_rec_dbl_buf),
            &(this->D4_rec_dbl_buf), &(this->D5_rec_dbl_buf), &(this->D6_rec_dbl_buf),
            &(this->D1_fix_dbl_buf));
        CHECK_ERR
    }

err_out:
    return err;
}

int e3sm_io_case_G::load_data (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver) {
    int err, verbose;

    verbose     = cfg.verbose;
    cfg.verbose = -1;  // Disable output

    // Run dummy G case read for data
    MPI_Barrier (cfg.io_comm);
    err = run_varn_G_case_rd (
        cfg, decom, driver, &(this->D1_fix_int_buf), &(this->D2_fix_int_buf),
        &(this->D3_fix_int_buf), &(this->D4_fix_int_buf), &(this->D5_fix_int_buf),
        &(this->D1_rec_dbl_buf), &(this->D3_rec_dbl_buf), &(this->D4_rec_dbl_buf),
        &(this->D5_rec_dbl_buf), &(this->D6_rec_dbl_buf), &(this->D1_fix_dbl_buf));
    CHECK_ERR

    cfg.verbose = verbose;

err_out:
    return err;
}
