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

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcopy() */
#include <assert.h>

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
    char base_name[1040], outfile[1056], *ext;

    nvars = cfg.nvars;
    nrecs = cfg.nrecs;

    /* construct I/O metadata */
    err = calc_metadata(&cfg, &decom);
    CHECK_ERR

    ext = strrchr(cfg.out_path, '.');

    if (cfg.hx == 0 || cfg.hx == -1) {  /* h0 file */
        /* construct file name */
        if (ext == NULL || (strcmp(ext, ".nc") && strcmp(ext, ".h5"))) {
            sprintf(base_name, "%s_h0", cfg.out_path);
        }
        else {
            sprintf(base_name, "%s", cfg.out_path);
            sprintf(base_name + (ext - cfg.out_path), "_h0%s", ext);
        }
        if (cfg.strategy == blob) /* set output subfile name */
            sprintf(outfile, "%s.%04d", base_name, cfg.subfile_ID);
        else
            strcpy(outfile, base_name);

        cfg.hist  = h0;
        cfg.nvars = 560;
        cfg.nrecs = nrecs;
        err = var_wr_I_case(cfg, decom, driver, outfile);
        CHECK_ERR

        /* report timing breakdowns */
        report_timing_WR(&cfg, &driver, &decom, base_name);
    }

    if (cfg.hx == 1 || cfg.hx == -1) {  /* h1 file */
        /* construct file name */
        if (ext == NULL || (strcmp(ext, ".nc") && strcmp(ext, ".h5"))) {
            sprintf(base_name, "%s_h1", cfg.out_path);
        }
        else {
            sprintf(base_name, "%s", cfg.out_path);
            sprintf(base_name + (ext - cfg.out_path), "_h1%s", ext);
        }
        if (cfg.strategy == blob) /* set output subfile name */
            sprintf(outfile, "%s.%04d", base_name, cfg.subfile_ID);
        else
            strcpy(outfile, base_name);

        cfg.hist  = h1;
        cfg.nvars = 552;
        cfg.nrecs = 1;
        err = var_wr_I_case(cfg, decom, driver, outfile);
        CHECK_ERR

        /* report timing breakdowns */
        report_timing_WR(&cfg, &driver, &decom, base_name);
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
