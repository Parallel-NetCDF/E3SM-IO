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

    /* construct I/O metadata */
    err = calc_metadata(&cfg, &decom);
    CHECK_ERR

    err = var_wr_G_case(cfg, decom, driver);
    CHECK_ERR

err_out:
    if (cfg.sub_comm != MPI_COMM_NULL) {
        MPI_Comm_free(&cfg.sub_comm);
        cfg.sub_comm = MPI_COMM_NULL;
    }

    return err;
}

int e3sm_io_case_G::rd_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver) {
    int err;

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

err_out:
    return err;
}
