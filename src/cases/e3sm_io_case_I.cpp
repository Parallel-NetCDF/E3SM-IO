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
#include <stdlib.h>

#include <e3sm_io.h>
#include <e3sm_io_err.h>
#include <e3sm_io_case.hpp>

e3sm_io_case_I::e3sm_io_case_I () {}

e3sm_io_case_I::~e3sm_io_case_I () {}

int e3sm_io_case_I::wr_test(e3sm_io_config &cfg,
                            e3sm_io_decom  &decom,
                            e3sm_io_driver &driver)
{
    int err=0, nvars, nrecs;
    
    nvars = cfg.nvars;
    nrecs = cfg.nrecs;
    
    /* construct I/O metadata */
    err = calc_metadata(&cfg, &decom);
    CHECK_ERR
    
    if (cfg.hx == 0 || cfg.hx == -1) {  /* h0 file */
        cfg.nvars = 560;
        cfg.nrecs = nrecs;
        err = var_wr_I_case(cfg, decom, driver);
        CHECK_ERR
    }
    
    if (cfg.hx == 1 || cfg.hx == -1) {  /* h1 file */
        cfg.nvars = 552;
        cfg.nrecs = 1;
        err = var_wr_I_case(cfg, decom, driver);
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


int e3sm_io_case_I::rd_test(e3sm_io_config &cfg,
                            e3sm_io_decom &decom,
                            e3sm_io_driver &driver)
{
    return 0;
}
