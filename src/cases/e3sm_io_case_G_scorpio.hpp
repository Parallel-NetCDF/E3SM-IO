/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 *********************************************************************/
//
#pragma once
//
#include <string>
#include <vector>
//
#include <mpi.h>
//
#include <e3sm_io.h>

#include <e3sm_io_case_scorpio.hpp>
#include <e3sm_io_driver.hpp>

int def_G_case_scorpio (e3sm_io_config &cfg,
                        e3sm_io_decom &decom,
                        e3sm_io_driver &driver,
                        int ncid, /* file ID */
                        std::vector<int> &decomids,
                        e3sm_io_scorpio_var *varids, /* variable IDs */
                        int *scorpiovars);

int run_varn_G_case_scorpio (e3sm_io_config &cfg,
                            e3sm_io_decom &decom,
                            e3sm_io_driver &driver,
                            int *D1_fix_int_bufp,     /* D1 fix int buffer */
                            int *D2_fix_int_bufp,     /* D2 fix int buffer */
                            int *D3_fix_int_bufp,     /* D3 fix int buffer */
                            int *D4_fix_int_bufp,     /* D4 fix int buffer */
                            int *D5_fix_int_bufp,     /* D5 fix int buffer */
                            double *D1_rec_dbl_bufp,  /* D1 rec double buffer */
                            double *D3_rec_dbl_bufp,  /* D3 rec double buffer */
                            double *D4_rec_dbl_bufp,  /* D4 rec double buffer */
                            double *D5_rec_dbl_bufp,  /* D5 rec double buffer */
                            double *D6_rec_dbl_bufp,  /* D6 rec double buffer */
                            double *D1_fix_dbl_bufp); /* D1 fix double buffer */
