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
    int err=0;
    char *outfile, *ext;

    ext = strrchr(cfg.out_path, '.');

    if (cfg.hx == 0 || cfg.hx == -1) {
        outfile = cfg.F_case_h0.outfile;

        /* construct file name */
        if (ext == NULL || (strcmp(ext, ".nc") && strcmp(ext, ".h5"))) {
            sprintf(outfile, "%s_h0", cfg.out_path);
        }
        else {
            sprintf(outfile, "%s", cfg.out_path);
            sprintf(outfile + (ext - cfg.out_path), "_h0%s", ext);
        }

        cfg.hist  = h0;;
        cfg.nvars = NVARS_F_CASE_H0;
        err = run_varn_F_case_scorpio(cfg, decom, driver, outfile,
                                      this->dbl_buf_h0,
                                      this->rec_buf_h0,
                                      this->txt_buf,
                                      this->int_buf);
        CHECK_ERR
    }

    if (cfg.hx == 1 || cfg.hx == -1) {
        outfile = cfg.F_case_h1.outfile;

        /* construct file name */
        if (ext == NULL || (strcmp(ext, ".nc") && strcmp(ext, ".h5"))) {
            sprintf(outfile, "%s_h1", cfg.out_path);
        }
        else {
            sprintf(outfile, "%s", cfg.out_path);
            sprintf(outfile + (ext - cfg.out_path), "_h1%s", ext);
        }

        cfg.hist  = h1;;
        cfg.nvars = NVARS_F_CASE_H1;
        err = run_varn_F_case_scorpio(cfg, decom, driver, outfile,
                                      this->dbl_buf_h0,
                                      this->rec_buf_h0,
                                      this->txt_buf,
                                      this->int_buf);
        CHECK_ERR
    }

err_out:
    return err;
}

int e3sm_io_case_F_scorpio::rd_test(e3sm_io_config &cfg,
                                    e3sm_io_decom &decom,
                                    e3sm_io_driver &driver)
{
    int err=0;
    ERR_OUT ("PIO case does not support reading")
err_out:
    return err;
}
