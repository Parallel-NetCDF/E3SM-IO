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
#include <e3sm_io_case_F.hpp>
#include <e3sm_io.h>
#include <e3sm_io_case.hpp>
#include <e3sm_io_err.h>

e3sm_io_case_F::e3sm_io_case_F () {}

e3sm_io_case_F::~e3sm_io_case_F () {
    if (this->dbl_buf_h0 != NULL) { free (this->dbl_buf_h0); }
    if (this->dbl_buf_h1 != NULL) { free (this->dbl_buf_h1); }
    if (this->rec_buf_h0 != NULL) { free (this->rec_buf_h0); }
    if (this->rec_buf_h1 != NULL) { free (this->rec_buf_h1); }
}

int e3sm_io_case_F::wr_test(e3sm_io_config &cfg,
                            e3sm_io_decom &decom,
                            e3sm_io_driver &driver)
{
    int err=0, nvar;

    PRINT_MSG (0, "number of requests for D1=%d D2=%d D3=%d\n", decom.contig_nreqs[0],
               decom.contig_nreqs[1], decom.contig_nreqs[2]);

    if ((cfg.verbose >= 0) && (cfg.rank == 0)) {
        printf ("Total number of MPI processes      = %d\n", cfg.np);
        printf ("Number of IO processes             = %d\n", cfg.num_iotasks);
        printf ("Input decomposition file           = %s\n", cfg.cfg_path);
        printf ("Number of decompositions           = %d\n", decom.num_decomp);
        printf ("Output file/directory              = %s\n", cfg.out_path);
        printf ("Variable dimensions (C order)      = %lld x %lld\n", decom.dims[2][0],
                decom.dims[2][1]);
        printf ("Write number of records (time dim) = %d\n", cfg.nrec);
        printf ("Using noncontiguous write buffer   = %s\n", cfg.non_contig_buf ? "yes" : "no");
        printf("==== Benchmarking F case ====\n");
        if (cfg.strategy == blob)
            printf("Using one-file-per-node blob I/O strategy\n");
        else if (cfg.vard)
            printf("Using PnetCDF vard API\n");
        else
            printf("Using varn API\n");
        if (cfg.two_buf)
            printf("Variable written order: 2D variables then 3D variables\n");
        else
            printf("Variable written order: same as variables are defined\n");
    }

    nvar = cfg.nvars;

    if (cfg.api == pnetcdf && cfg.strategy == blob) {
        /* construct metadata for blob I/O strategy */
        err = blob_metadata(&cfg, &decom);
        CHECK_ERR

        /* Use one-file-per-compute-node blob I/O strategy */
        if (cfg.hx == 0 || cfg.hx == -1) {  /* h0 file */
            cfg.nvars = 414;
            err = pnetcdf_blob_F_case(cfg, decom, driver);
            CHECK_ERR

        }

        if (cfg.hx == 1 || cfg.hx == -1) {  /* h1 file */
            cfg.nvars = 51;
            err = pnetcdf_blob_F_case(cfg, decom, driver);
            CHECK_ERR
        }

        if (cfg.sub_comm != MPI_COMM_NULL)
            MPI_Comm_free(&cfg.sub_comm);
    }
    else if (cfg.vard) { /* using PnetCDF vard APIs to write/read */
        /* vard APIs require internal data type matches external one */
#if REC_XTYPE != NC_FLOAT
        ERR_OUT ("PnetCDF vard API requires internal and external data types match, skip\n");
#endif
        if (cfg.hx == 0 || cfg.hx == -1) {
            cfg.nvars = 414;
            err = run_vard_F_case (cfg, decom, driver, this->dbl_buf_h0,
                                      this->rec_buf_h0, this->txt_buf[0], this->int_buf[0]);
        }

        if (cfg.hx == 0 || cfg.hx == -1) {
            cfg.nvars = 51;
            err = run_vard_F_case (cfg, decom, driver, this->dbl_buf_h0,
                                      this->rec_buf_h0, this->txt_buf[0], this->int_buf[0]);
        }
    } else { /* using PnetCDF varn APIs to write/read */
        if (cfg.hx == 0 || cfg.hx == -1) {
            cfg.nvars = 414;
            err = run_varn_F_case (cfg, decom, driver, this->dbl_buf_h0,
                                      this->rec_buf_h0, this->txt_buf[0], this->int_buf[0]);
        }

        if (cfg.hx == 1 || cfg.hx == -1) {
            cfg.nvars = 51;
            err = run_varn_F_case (cfg, decom, driver, this->dbl_buf_h0,
                                      this->rec_buf_h0, this->txt_buf[0], this->int_buf[0]);
        }
    }

    cfg.nvars = nvar;

err_out:
    return err;
}

int e3sm_io_case_F::rd_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver) {
    int err, nvar;

    PRINT_MSG (0, "number of requests for D1=%d D2=%d D3=%d\n", decom.contig_nreqs[0],
               decom.contig_nreqs[1], decom.contig_nreqs[2]);

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

    nvar = cfg.nvars;

    /* vard APIs require internal data type matches external one */
    if (cfg.vard) {
        //#if REC_XTYPE != NC_FLOAT
        ERR_OUT ("Low level API requires internal and external data types match, skip\n");
        //#endif
        ERR_OUT ("Reading not supported for low-level API\n");
    } else {
        PRINT_MSG (
            0, "\n==== benchmarking F case read using low level API ========================\n");

        if (cfg.two_buf) {
            PRINT_MSG (0, "Variable read order: 2D variables then 3D variables\n\n");
        } else {
            PRINT_MSG (0, "Variable read order: same as variables are defined\n\n");
        }

        fflush (stdout);

        if (cfg.hx == 0 || cfg.hx == -1) {
            MPI_Barrier (cfg.io_comm);
            cfg.nvars = 414;
            err = run_varn_F_case_rd(cfg, decom, driver, &(this->dbl_buf_h0),
                                     &(this->rec_buf_h0), this->txt_buf[0], this->int_buf[0]);
        }
        if (cfg.hx == 0 || cfg.hx == -1) {
            MPI_Barrier (cfg.io_comm);
            cfg.nvars = 51;
            err = run_varn_F_case_rd(cfg, decom, driver, &(this->dbl_buf_h0),
                                     &(this->rec_buf_h0), this->txt_buf[0], this->int_buf[0]);
        }
    }

    cfg.nvars = nvar;

err_out:
    return err;
}
