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
#include <stdlib.h>
//
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
    int err=0, nvars, nrecs;

    nvars = cfg.nvars;
    nrecs = cfg.nrecs;

    /* construct I/O metadata */
    err = calc_metadata(&cfg, &decom);
    CHECK_ERR

    if (cfg.hx == 0 || cfg.hx == -1) {  /* h0 file */
        cfg.nvars = 414;
        cfg.nrecs = 1;
        err = var_wr_F_case(cfg, decom, driver);
        CHECK_ERR
    }

    if (cfg.hx == 1 || cfg.hx == -1) {  /* h1 file */
        cfg.nvars = 51;
        cfg.nrecs = nrecs;
        err = var_wr_F_case(cfg, decom, driver);
        CHECK_ERR
    }

    cfg.nvars = nvars;
    cfg.nrecs = nrecs;

err_out:
    if (cfg.sub_comm != MPI_COMM_NULL) {
        MPI_Comm_free(&cfg.sub_comm);
        cfg.sub_comm = MPI_COMM_NULL;
    }

    return err;
}

int e3sm_io_case_F::rd_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver) {
    int err, nvars, nrecs;

    nvars = cfg.nvars;
    nrecs = cfg.nrecs;

    if (cfg.hx == 0 || cfg.hx == -1) {
        MPI_Barrier (cfg.io_comm);
        cfg.nvars = 414;
        cfg.nrecs = 1;
        err = run_varn_F_case_rd(cfg, decom, driver,
                                 &(this->dbl_buf_h0),
                                 &(this->rec_buf_h0),
                                 this->txt_buf,
                                 this->int_buf);
        CHECK_ERR
    }
    if (cfg.hx == 0 || cfg.hx == -1) {
        MPI_Barrier (cfg.io_comm);
        cfg.nvars = 51;
        cfg.nrecs = nrecs;
        err = run_varn_F_case_rd(cfg, decom, driver,
                                 &(this->dbl_buf_h0),
                                 &(this->rec_buf_h0),
                                 this->txt_buf,
                                 this->int_buf);
        CHECK_ERR
    }

    cfg.nvars = nvars;
    cfg.nrecs = nrecs;

err_out:
    return err;
}
