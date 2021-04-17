#pragma once
#include <mpi.h>

#include <string>

#include "e3sm_io.hpp"
#include "e3sm_io_driver.hpp"

extern int def_F_case_h0 (
    e3sm_io_driver &driver, int ncid, const MPI_Offset dims[2], int nvars, int *varids);
extern int inq_F_case_h0 (e3sm_io_driver &driver,
                          int ncid,           /* file ID */
                          MPI_Offset dims[2], /* dimension sizes */
                          int nvars,          /* number of variables */
                          int *varids);       /* variable IDs */

extern int def_F_case_h1 (
    e3sm_io_driver &driver, int ncid, const MPI_Offset dims[2], int nvars, int *varids);
extern int inq_F_case_h1 (e3sm_io_driver &driver,
                          int ncid,           /* file ID */
                          MPI_Offset dims[2], /* dimension sizes */
                          int nvars,          /* number of variables */
                          int *varids);       /* variable IDs */

extern int run_vard_F_case (e3sm_io_config &cfg,
                            e3sm_io_decom &decom,
                            e3sm_io_driver &driver,
                            std::string outfile, /* output file name */
                            double *dbl_bufp,    /* buffer for fixed size double var */
                            itype *rec_bufp,     /* buffer for rec floating point var */
                            char *txt_buf,       /* buffer for char var */
                            int *int_buf);       /* request's block lengths */

extern int run_varn_F_case (e3sm_io_config &cfg,
                            e3sm_io_decom &decom,
                            e3sm_io_driver &driver,
                            std::string outfile, /* output file name */
                            double *dbl_bufp,    /* buffer for fixed size double var */
                            itype *rec_bufp,     /* buffer for rec floating point var */
                            char *txt_buf,       /* buffer for char var */
                            int *int_buf);       /* buffer for int var */
extern int run_varn_F_case_rd (e3sm_io_config &cfg,
                               e3sm_io_decom &decom,
                               e3sm_io_driver &driver,
                               std::string outfile, /* output file name */
                               double **dbl_bufp,    /* buffer for fixed size double var */
                               itype **rec_bufp,     /* buffer for rec floating point var */
                               char *txt_buf,       /* buffer for char var */
                               int *int_buf);       /* buffer for int var */