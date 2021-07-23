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

#include <vector>
#include <e3sm_io.h>
#include <e3sm_io_driver.hpp>

class e3sm_io_case {
   public:
    virtual ~e3sm_io_case () {};
    virtual int wr_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver)   = 0;
    virtual int rd_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver)   = 0;
};

class e3sm_io_case_F : public e3sm_io_case {
   protected:
    double *dbl_buf_h0 = NULL, *dbl_buf_h1 = NULL;
    itype *rec_buf_h0 = NULL, *rec_buf_h1 = NULL;
    char *txt_buf = NULL;
    int *int_buf = NULL;

   public:
    e3sm_io_case_F ();
    ~e3sm_io_case_F ();
    int wr_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver);
    int rd_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver);
};

class e3sm_io_case_G : public e3sm_io_case {
   protected:
    double *D1_rec_dbl_buf = NULL, *D3_rec_dbl_buf = NULL, *D4_rec_dbl_buf = NULL,
           *D5_rec_dbl_buf = NULL, *D6_rec_dbl_buf = NULL, *D1_fix_dbl_buf = NULL;
    int *D1_fix_int_buf = NULL, *D2_fix_int_buf = NULL, *D3_fix_int_buf = NULL,
        *D4_fix_int_buf = NULL, *D5_fix_int_buf = NULL;

   public:
    e3sm_io_case_G ();
    ~e3sm_io_case_G ();
    int wr_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver);
    int rd_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver);
};

#ifdef ENABLE_ADIOS2
class e3sm_io_case_F_scorpio : public e3sm_io_case_F {
   public:
    e3sm_io_case_F_scorpio ();
    ~e3sm_io_case_F_scorpio ();
    int wr_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver);
    int rd_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver);
};

class e3sm_io_case_G_scorpio : public e3sm_io_case_G {
   public:
    e3sm_io_case_G_scorpio ();
    ~e3sm_io_case_G_scorpio ();
    int wr_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver);
    int rd_test (e3sm_io_config &cfg, e3sm_io_decom &decom, e3sm_io_driver &driver);
};
#endif

/*---- functions for F case -------------------------------------------------*/
extern int
def_F_case_h0(e3sm_io_driver &driver,
              int ncid,
              const MPI_Offset dims[2],
              int nvars,
              int *varids);

extern int
inq_F_case_h0(e3sm_io_driver &driver,
              int ncid,           /* file ID */
              MPI_Offset dims[2], /* dimension sizes */
              int nvars,          /* number of variables */
              int *varids);       /* variable IDs */

extern int
def_F_case_h1(e3sm_io_driver &driver,
              int ncid,
              const MPI_Offset dims[2],
              int nvars,
              int *varids);

extern int
inq_F_case_h1(e3sm_io_driver &driver,
              int ncid,           /* file ID */
              MPI_Offset dims[2], /* dimension sizes */
              int nvars,          /* number of variables */
              int *varids);       /* variable IDs */

extern int
run_vard_F_case(e3sm_io_config &cfg,
                e3sm_io_decom &decom,
                e3sm_io_driver &driver,
                double *dbl_bufp,    /* buffer for fixed size double var */
                itype *rec_bufp,     /* buffer for rec floating point var */
                char *txt_buf,       /* buffer for char var */
                int *int_buf);       /* request's block lengths */

extern int
run_varn_F_case(e3sm_io_config &cfg,
                e3sm_io_decom &decom,
                e3sm_io_driver &driver);

extern int
run_varn_F_case_rd(e3sm_io_config &cfg,
                   e3sm_io_decom &decom,
                   e3sm_io_driver &driver,
                   double **dbl_bufp,   /* buffer for fixed size double var */
                   itype **rec_bufp,    /* buffer for rec floating point var */
                   char *txt_buf,       /* buffer for char var */
                   int *int_buf);       /* buffer for int var */

extern int
def_F_case_h0(e3sm_io_config &cfg,
              e3sm_io_decom  &decom,
              e3sm_io_driver &driver,
              int ncid,
              int *varids);

extern int
def_F_case_h1(e3sm_io_config &cfg,
              e3sm_io_decom  &decom,
              e3sm_io_driver &driver,
              int ncid,
              int *varids);

extern int
blob_F_case(e3sm_io_config &cfg,
            e3sm_io_decom &decom,
            e3sm_io_driver &driver);


/*---- functions for G case -------------------------------------------------*/
extern int
def_G_case_h0(e3sm_io_driver &driver,
              int ncid,                    /* file ID */
              const MPI_Offset dims_D1[1], /* dimension sizes of decomposition 1 */
              const MPI_Offset dims_D2[1], /* dimension sizes of decomposition 2 */
              const MPI_Offset dims_D3[2], /* dimension sizes of decomposition 3 */
              const MPI_Offset dims_D4[2], /* dimension sizes of decomposition 4 */
              const MPI_Offset dims_D5[2], /* dimension sizes of decomposition 5 */
              const MPI_Offset dims_D6[2], /* dimension sizes of decomposition 6 */
              int nvars,                   /* number of variables */
              int *varids);                /* variable IDs */

extern int
inq_G_case_h0(e3sm_io_driver &driver,
              int ncid,              /* file ID */
              MPI_Offset dims_D1[1], /* dimension sizes of decomposition 1 */
              MPI_Offset dims_D2[1], /* dimension sizes of decomposition 2 */
              MPI_Offset dims_D3[2], /* dimension sizes of decomposition 3 */
              MPI_Offset dims_D4[2], /* dimension sizes of decomposition 4 */
              MPI_Offset dims_D5[2], /* dimension sizes of decomposition 5 */
              MPI_Offset dims_D6[2], /* dimension sizes of decomposition 6 */
              int nvars,             /* number of variables */
              int *varids);          /* variable IDs */

extern int
run_varn_G_case(e3sm_io_config &cfg,
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

extern int
run_varn_G_case_rd(e3sm_io_config &cfg,
                   e3sm_io_decom &decom,
                   e3sm_io_driver &driver,
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
blob_G_case(e3sm_io_config &cfg,
            e3sm_io_decom &decom,
            e3sm_io_driver &driver);

extern int
def_G_case(e3sm_io_config &cfg,
           e3sm_io_decom  &decom,
           e3sm_io_driver &driver,
           int ncid,
           int *varids);

extern int
check_malloc(e3sm_io_config *cfg,
             e3sm_io_driver *driver);

extern int
report_timing_WR(e3sm_io_config *cfg,
                 e3sm_io_driver *driver,
                 const char *outfile);

