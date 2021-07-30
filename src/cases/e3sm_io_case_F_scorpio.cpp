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
#include <e3sm_io.h>
#include <e3sm_io_err.h>

#include <e3sm_io_case_scorpio.hpp>

e3sm_io_case_F_scorpio::e3sm_io_case_F_scorpio () {}
e3sm_io_case_F_scorpio::~e3sm_io_case_F_scorpio () {}

int e3sm_io_case_F_scorpio::wr_test(e3sm_io_config &cfg,
                                    e3sm_io_decom &decom,
                                    e3sm_io_driver &driver)
{
    int err=0, nvars, nrecs;

    nvars = cfg.nvars;
    nrecs = cfg.nrecs;

    if (cfg.hx == 0 || cfg.hx == -1) {
        cfg.nvars = 414;
        cfg.nrecs = 1;
        err = run_varn_F_case_scorpio(cfg, decom, driver,
                                      this->dbl_buf_h0,
                                      this->rec_buf_h0,
                                      this->txt_buf,
                                      this->int_buf);
        CHECK_ERR
    }

    if (cfg.hx == 1 || cfg.hx == -1) {
        cfg.nvars = 51;
        cfg.nrecs = nrecs;
        err = run_varn_F_case_scorpio(cfg, decom, driver,
                                      this->dbl_buf_h0,
                                      this->rec_buf_h0,
                                      this->txt_buf,
                                      this->int_buf);
        CHECK_ERR
    }

    cfg.nvars = nvars;
    cfg.nrecs = nrecs;

err_out:
    return err;
}

int e3sm_io_case_F_scorpio::rd_test(e3sm_io_config &cfg E3SM_IO_UNUSED,
                                    e3sm_io_decom &decom E3SM_IO_UNUSED,
                                    e3sm_io_driver &driver E3SM_IO_UNUSED)
{
    int err=0;
    ERR_OUT ("PIO case does not support reading")
err_out:
    return err;
}
