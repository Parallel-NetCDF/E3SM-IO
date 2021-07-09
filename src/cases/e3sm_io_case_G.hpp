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
//
#include <mpi.h>
//
#include <e3sm_io.h>
#include <e3sm_io_driver.hpp>

extern int def_G_case_h0 (e3sm_io_driver &driver,
                          int ncid,                    /* file ID */
                          const MPI_Offset dims_D1[1], /* dimension sizes of decomposition 1 */
                          const MPI_Offset dims_D2[1], /* dimension sizes of decomposition 2 */
                          const MPI_Offset dims_D3[2], /* dimension sizes of decomposition 3 */
                          const MPI_Offset dims_D4[2], /* dimension sizes of decomposition 4 */
                          const MPI_Offset dims_D5[2], /* dimension sizes of decomposition 5 */
                          const MPI_Offset dims_D6[2], /* dimension sizes of decomposition 6 */
                          int nvars,                   /* number of variables */
                          int *varids);                /* variable IDs */

extern int inq_G_case_h0 (e3sm_io_driver &driver,
                          int ncid,              /* file ID */
                          MPI_Offset dims_D1[1], /* dimension sizes of decomposition 1 */
                          MPI_Offset dims_D2[1], /* dimension sizes of decomposition 2 */
                          MPI_Offset dims_D3[2], /* dimension sizes of decomposition 3 */
                          MPI_Offset dims_D4[2], /* dimension sizes of decomposition 4 */
                          MPI_Offset dims_D5[2], /* dimension sizes of decomposition 5 */
                          MPI_Offset dims_D6[2], /* dimension sizes of decomposition 6 */
                          int nvars,             /* number of variables */
                          int *varids);          /* variable IDs */

extern int run_varn_G_case (e3sm_io_config &cfg,
                            e3sm_io_decom &decom,
                            e3sm_io_driver &driver,
                            std::string outfile,
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

extern int run_varn_G_case_rd (e3sm_io_config &cfg,
                               e3sm_io_decom &decom,
                               e3sm_io_driver &driver,
                               std::string outfile,
                               int **D1_fix_int_bufp,     /* D1 fix int buffer */
                               int **D2_fix_int_bufp,     /* D2 fix int buffer */
                               int **D3_fix_int_bufp,     /* D3 fix int buffer */
                               int **D4_fix_int_bufp,     /* D4 fix int buffer */
                               int **D5_fix_int_bufp,     /* D5 fix int buffer */
                               double **D1_rec_dbl_bufp,  /* D1 rec double buffer */
                               double **D3_rec_dbl_bufp,  /* D3 rec double buffer */
                               double **D4_rec_dbl_bufp,  /* D4 rec double buffer */
                               double **D5_rec_dbl_bufp,  /* D5 rec double buffer */
                               double **D6_rec_dbl_bufp,  /* D6 rec double buffer */
                               double **D1_fix_dbl_bufp); /* D1 fix double buffer */

extern int
pnetcdf_blob_G_case(e3sm_io_config &cfg,
                    e3sm_io_decom &decom,
                    e3sm_io_driver &driver,
                    char *outfile);

extern int
blob_def_G_case(e3sm_io_config &cfg,
                e3sm_io_decom  &decom,
                e3sm_io_driver &driver,
                int ncid,
                int *varids);

