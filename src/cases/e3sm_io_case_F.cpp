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

#include <assert.h>

#include <e3sm_io.h>
#include <e3sm_io_err.h>
#include <e3sm_io_case.hpp>

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

    nvar = cfg.nvars;

    if (cfg.strategy == blob) {
        assert (cfg.api == pnetcdf || cfg.api == hdf5);

        /* construct metadata for blob I/O strategy */
        err = blob_metadata(&cfg, &decom);
        CHECK_ERR

        /* Use one-file-per-compute-node blob I/O strategy */
        if (cfg.hx == 0 || cfg.hx == -1) {  /* h0 file */
            cfg.nvars = 414;
            err = blob_F_case(cfg, decom, driver);
            CHECK_ERR

        }

        if (cfg.hx == 1 || cfg.hx == -1) {  /* h1 file */
            cfg.nvars = 51;
            err = blob_F_case(cfg, decom, driver);
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
            err = run_vard_F_case(cfg, decom, driver,
                                  this->dbl_buf_h0,
                                  this->rec_buf_h0,
                                  this->txt_buf[0],
                                  this->int_buf[0]);
            CHECK_ERR
        }

        if (cfg.hx == 0 || cfg.hx == -1) {
            cfg.nvars = 51;
            err = run_vard_F_case (cfg, decom, driver,
                                  this->dbl_buf_h0,
                                  this->rec_buf_h0,
                                  this->txt_buf[0],
                                  this->int_buf[0]);
            CHECK_ERR
        }
    } else { /* using PnetCDF varn APIs to write/read */
        if (cfg.hx == 0 || cfg.hx == -1) {
            cfg.nvars = 414;
            err = run_varn_F_case(cfg, decom, driver);
            CHECK_ERR
        }

        if (cfg.hx == 1 || cfg.hx == -1) {
            cfg.nvars = 51;
            err = run_varn_F_case(cfg, decom, driver);
            CHECK_ERR
        }
    }

    cfg.nvars = nvar;

err_out:
    return err;
}

int e3sm_io_case_F::rd_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver) {
    int err, nvar;

    nvar = cfg.nvars;

    /* vard APIs require internal data type matches external one */
    if (cfg.vard) {
        //#if REC_XTYPE != NC_FLOAT
        ERR_OUT ("Low level API requires internal and external data types match, skip\n");
        //#endif
        ERR_OUT ("Reading not supported for low-level API\n");
    } else {
        if (cfg.hx == 0 || cfg.hx == -1) {
            MPI_Barrier (cfg.io_comm);
            cfg.nvars = 414;
            err = run_varn_F_case_rd(cfg, decom, driver,
                                     &(this->dbl_buf_h0),
                                     &(this->rec_buf_h0),
                                     this->txt_buf[0],
                                     this->int_buf[0]);
            CHECK_ERR
        }
        if (cfg.hx == 0 || cfg.hx == -1) {
            MPI_Barrier (cfg.io_comm);
            cfg.nvars = 51;
            err = run_varn_F_case_rd(cfg, decom, driver,
                                     &(this->dbl_buf_h0),
                                     &(this->rec_buf_h0),
                                     this->txt_buf[0],
                                     this->int_buf[0]);
            CHECK_ERR
        }
    }

    cfg.nvars = nvar;

err_out:
    return err;
}
