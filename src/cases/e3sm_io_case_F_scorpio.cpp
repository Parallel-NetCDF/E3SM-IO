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

#include <e3sm_io_case.hpp>
#include <e3sm_io_case_F_scorpio.hpp>

e3sm_io_case_F_scorpio::e3sm_io_case_F_scorpio () {}
e3sm_io_case_F_scorpio::~e3sm_io_case_F_scorpio () {}
int e3sm_io_case_F_scorpio::load_data (e3sm_io_config &cfg,
                                   e3sm_io_decom &decom,
                                   e3sm_io_driver &driver) {
    int err=0;
    ERR_OUT ("PIO case does not support reading")
err_out:;
    return err;
}
int e3sm_io_case_F_scorpio::wr_test (e3sm_io_config &cfg,
                                 e3sm_io_decom &decom,
                                 e3sm_io_driver &driver) {
    int err=0;
    int nvar;

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
    }

    nvar = cfg.nvars;

    /* vard APIs require internal data type matches external one */
    if (cfg.vard) {
        ERR_OUT ("PIO case does not support vard API")
    } else {
        PRINT_MSG (0,
                   "\n==== benchmarking PIO F case writing using varn API ========================\n");

        if (cfg.two_buf) {
            PRINT_MSG (0, "Variable written order: 2D variables then 3D variables\n\n");
        } else {
            PRINT_MSG (0, "Variable written order: same as variables are defined\n\n");
        }

        fflush (stdout);

        if (cfg.hx == 0 || cfg.hx == -1) {
            MPI_Barrier (cfg.io_comm);
            cfg.nvars = 414;
            err = run_varn_F_case_scorpio (cfg, decom, driver, this->dbl_buf_h0,
                                          this->rec_buf_h0, this->txt_buf[0], this->int_buf[0]);
            CHECK_ERR
        }

        if (cfg.hx == 1 || cfg.hx == -1) {
            MPI_Barrier (cfg.io_comm);
            cfg.nvars = 51;
            err = run_varn_F_case_scorpio (cfg, decom, driver, this->dbl_buf_h0,
                                          this->rec_buf_h0, this->txt_buf[0], this->int_buf[0]);
            CHECK_ERR
        }
    }

    cfg.nvars = nvar;
err_out:
    return err;
}

int e3sm_io_case_F_scorpio::rd_test (e3sm_io_config &cfg,
                                 e3sm_io_decom &decom,
                                 e3sm_io_driver &driver) {
    int err=0;
    ERR_OUT ("PIO case does not support reading")
err_out:
    return err;
}
