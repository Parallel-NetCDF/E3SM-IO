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

#include <e3sm_io_case_pio.hpp>
#include <e3sm_io_driver.hpp>

int def_F_case_h0_pio (e3sm_io_driver &driver,
                       e3sm_io_config &cfg,
                       e3sm_io_decom &decom,
                       int ncid,                 /* file ID */
                       const MPI_Offset dims[2], /* dimension sizes */
                       int nvars,                /* number of variables */
                       std::vector<int> &decomids,
                       e3sm_io_pio_var *varids,
                       int *piovars);

int def_F_case_h1_pio (e3sm_io_driver &driver,
                       e3sm_io_config &cfg,            
                       e3sm_io_decom &decom,
                       int ncid,                 /* file ID */
                       const MPI_Offset dims[2], /* dimension sizes */
                       int nvars,                /* number of variables */
                       std::vector<int> &decomids,
                       e3sm_io_pio_var *varids,
                       int *piovars);

int run_varn_F_case_pio (e3sm_io_config &cfg,
                         e3sm_io_decom &decom,
                         e3sm_io_driver &driver,
                         std::string outfile, /* output file name */
                         double *dbl_bufp,    /* buffer for fixed size double var */
                         itype *rec_bufp,     /* buffer for rec floating point var */
                         char *txt_buf,       /* buffer for char var */
                         int *int_buf);       /* buffer for int var */