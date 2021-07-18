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

e3sm_io_case_G_scorpio::e3sm_io_case_G_scorpio () {}
e3sm_io_case_G_scorpio::~e3sm_io_case_G_scorpio () {}

int e3sm_io_case_G_scorpio::wr_test(e3sm_io_config &cfg,
                                    e3sm_io_decom &decom,
                                    e3sm_io_driver &driver)
{
    int err=0;

    if (cfg.vard) {
        /* vard APIs require internal data type matches external one */
#if REC_XTYPE != NC_FLOAT
        ERR_OUT ("PnetCDF vard API requires internal and external data types match, skip\n");
#endif
        ERR_OUT ("Low level API not supported in g case\n");
    } else {
        err = run_varn_G_case_scorpio(cfg, decom, driver,
                                      this->D1_fix_int_buf,
                                      this->D2_fix_int_buf,
                                      this->D3_fix_int_buf,
                                      this->D4_fix_int_buf,
                                      this->D5_fix_int_buf,
                                      this->D1_rec_dbl_buf,
                                      this->D3_rec_dbl_buf,
                                      this->D4_rec_dbl_buf,
                                      this->D5_rec_dbl_buf,
                                      this->D6_rec_dbl_buf,
                                      this->D1_fix_dbl_buf);
        CHECK_ERR
    }

err_out:
    return err;
}

int e3sm_io_case_G_scorpio::rd_test(e3sm_io_config &cfg,
                                    e3sm_io_decom &decom,
                                    e3sm_io_driver &driver)
{
    int err=0;
    ERR_OUT ("PIO case does not support reading")
err_out:
    return err;
}
