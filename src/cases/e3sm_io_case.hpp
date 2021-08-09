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
   protected:
       double *dbl_buf_h0 = NULL, *dbl_buf_h1 = NULL;
       vtype  *rec_buf_h0 = NULL, *rec_buf_h1 = NULL;
       char   *txt_buf = NULL;
       int    *int_buf = NULL;
       double *D1_rec_dbl_buf = NULL, *D3_rec_dbl_buf = NULL, *D4_rec_dbl_buf = NULL,
              *D5_rec_dbl_buf = NULL, *D6_rec_dbl_buf = NULL, *D1_fix_dbl_buf = NULL;
       int    *D1_fix_int_buf = NULL, *D2_fix_int_buf = NULL, *D3_fix_int_buf = NULL,
              *D4_fix_int_buf = NULL, *D5_fix_int_buf = NULL;

    public:
        virtual ~e3sm_io_case () {};
        virtual int wr_test(e3sm_io_config &cfg,
                            e3sm_io_decom &decom,
                            e3sm_io_driver &driver) = 0;
        virtual int rd_test(e3sm_io_config &cfg,
                            e3sm_io_decom &decom,
                            e3sm_io_driver &driver) = 0;
};

class e3sm_io_all_cases : public e3sm_io_case {
   protected:
   public:
        e3sm_io_all_cases ();
        ~e3sm_io_all_cases ();
        int wr_test(e3sm_io_config &cfg,
                    e3sm_io_decom &decom,
                    e3sm_io_driver &driver);
        int rd_test(e3sm_io_config &cfg,
                    e3sm_io_decom &decom,
                    e3sm_io_driver &driver);
};

#ifdef ENABLE_ADIOS2
class e3sm_io_case_F_scorpio : public e3sm_io_case {
   protected:
   public:
        e3sm_io_case_F_scorpio ();
        ~e3sm_io_case_F_scorpio ();
        int wr_test(e3sm_io_config &cfg,
                    e3sm_io_decom &decom,
                    e3sm_io_driver &driver);
        int rd_test(e3sm_io_config &cfg,
                    e3sm_io_decom &decom,
                    e3sm_io_driver &driver);
};

class e3sm_io_case_G_scorpio : public e3sm_io_case {
   protected:
   public:
        e3sm_io_case_G_scorpio ();
        ~e3sm_io_case_G_scorpio ();
        int wr_test(e3sm_io_config &cfg,
                    e3sm_io_decom &decom,
                    e3sm_io_driver &driver);
        int rd_test(e3sm_io_config &cfg,
                    e3sm_io_decom &decom,
                    e3sm_io_driver &driver);
};
#endif

#define BUF_GAP 10

typedef struct {
    int vid;        /* variable ID */
    int decomp_id;  /* decomposition map ID */
    int isRecVar;   /* whether is a record variable */
    size_t vlen;    /* length to be written by this rank */
    MPI_Datatype itype;
} var_meta;

typedef struct {
    size_t  gap;

    /* buffers for fixed-size variables */
    size_t  fix_txt_buflen;
    char   *fix_txt_buf;
    size_t  fix_int_buflen;
    int    *fix_int_buf;
    size_t  fix_dbl_buflen;
    double *fix_dbl_buf;
    size_t  fix_buflen;
    vtype  *fix_buf;

    /* buffers for record variables */
    size_t  rec_txt_buflen;
    char   *rec_txt_buf;
    size_t  rec_int_buflen;
    int    *rec_int_buf;
    size_t  rec_dbl_buflen;
    double *rec_dbl_buf;
    size_t  rec_buflen;
    vtype  *rec_buf;
} io_buffers;

extern
int var_wr_all_cases(e3sm_io_config &cfg,
                     e3sm_io_decom  &decom,
                     e3sm_io_driver &driver,
                     case_meta      *cmeta);

/*---- functions for F case -------------------------------------------------*/
extern
int def_F_case(e3sm_io_config &cfg,
               e3sm_io_decom  &decom,
               e3sm_io_driver &driver,
               int             ncid,
               var_meta       *vars,
               io_buffers     *wr_buf);

extern
int inq_F_case_h0(e3sm_io_driver &driver,
                  int ncid,           /* file ID */
                  MPI_Offset dims[2], /* dimension sizes */
                  int nvars,          /* number of variables */
                  int *varids);       /* variable IDs */

extern
int inq_F_case_h1(e3sm_io_driver &driver,
                  int ncid,           /* file ID */
                  MPI_Offset dims[2], /* dimension sizes */
                  int nvars,          /* number of variables */
                  int *varids);       /* variable IDs */

extern int
run_varn_F_case_rd(e3sm_io_config &cfg,
                   e3sm_io_decom &decom,
                   e3sm_io_driver &driver,
                   double **dbl_bufp,   /* buffer for fixed size double var */
                   vtype **rec_bufp,    /* buffer for rec floating point var */
                   char *txt_buf,       /* buffer for char var */
                   int *int_buf);       /* buffer for int var */


/*---- functions for G case -------------------------------------------------*/
extern int
inq_G_case(e3sm_io_driver &driver,
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

extern
int def_G_case(e3sm_io_config &cfg,
               e3sm_io_decom  &decom,
               e3sm_io_driver &driver,
               int             ncid,
               var_meta       *vars,
               io_buffers     *wr_buf);

extern
int def_I_case(e3sm_io_config &cfg,
               e3sm_io_decom  &decom,
               e3sm_io_driver &driver,
               int             ncid,
               var_meta       *vars,
               io_buffers     *wr_buf);

extern
int check_malloc(e3sm_io_config *cfg,
                 e3sm_io_driver *driver);

extern
void wr_buf_init(io_buffers &buf, int gap);

extern
int wr_buf_malloc(e3sm_io_config &cfg, int one_flush, io_buffers &buf);

extern
void wr_buf_free(io_buffers &buf);

