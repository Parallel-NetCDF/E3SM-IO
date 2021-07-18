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

    if (cfg.strategy == blob) {
        assert (cfg.api == pnetcdf || cfg.api == hdf5);

        /* construct metadata for blob I/O strategy */
        err = blob_metadata(&cfg, &decom);
        CHECK_ERR

        /* Use one-file-per-compute-node blob I/O strategy */
        err = blob_G_case(cfg, decom, driver);
        CHECK_ERR

        if (cfg.sub_comm != MPI_COMM_NULL)
            MPI_Comm_free(&cfg.sub_comm);
    }
    else if (cfg.vard) { /* using PnetCDF vard APIs to write/read */
        /* vard APIs require internal data type matches external one */
#if REC_XTYPE != NC_FLOAT
        ERR_OUT ("PnetCDF vard API requires internal and external data types match, skip\n");
#endif
        ERR_OUT ("Low level API not supported in g case\n");
    } else { /* using PnetCDF varn APIs to write/read */
        err = run_varn_G_case(cfg, decom, driver,
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

int e3sm_io_case_G::rd_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver) {
    int err;

    /* vard APIs require internal data type matches external one */
    if (cfg.vard) {
#if REC_XTYPE != NC_FLOAT
        ERR_OUT ("PnetCDF vard API requires internal and external data types match, skip\n");
#endif
        ERR_OUT ("Low level API not supported in g case\n");
    } else {
        err = run_varn_G_case_rd(cfg, decom, driver,
                                 &(this->D1_fix_int_buf),
                                 &(this->D2_fix_int_buf),
                                 &(this->D3_fix_int_buf),
                                 &(this->D4_fix_int_buf),
                                 &(this->D5_fix_int_buf),
                                 &(this->D1_rec_dbl_buf),
                                 &(this->D3_rec_dbl_buf),
                                 &(this->D4_rec_dbl_buf),
                                 &(this->D5_rec_dbl_buf),
                                 &(this->D6_rec_dbl_buf),
                                 &(this->D1_fix_dbl_buf));
        CHECK_ERR
    }

err_out:
    return err;
}
