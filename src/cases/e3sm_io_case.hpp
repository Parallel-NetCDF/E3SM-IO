/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 *********************************************************************/

#pragma once

#include <vector>
#include <map>

#include <e3sm_io.h>
#include <e3sm_io_driver.hpp>

#define BUF_GAP 10

typedef struct {
    size_t  gap;

    /* buffers for fixed-size variables */
    size_t  fix_txt_buflen; char   *fix_txt_buf;
    size_t  fix_int_buflen; int    *fix_int_buf;
    size_t  fix_dbl_buflen; double *fix_dbl_buf;
    size_t  fix_buflen;     vtype  *fix_buf;

    /* buffers for record variables */
    size_t  rec_txt_buflen; char   *rec_txt_buf;
    size_t  rec_int_buflen; int    *rec_int_buf;
    size_t  rec_dbl_buflen; double *rec_dbl_buf;
    size_t  rec_buflen;     vtype  *rec_buf;
} io_buffers;

typedef struct {
    int vid;         /* variable ID, returned from the driver */

    int frame_id;    /* frame variable ID returned from adios driver */
    int fillval_id;  /* fillval variable ID returned from adios driver */
    int decom_id;    /* decomposition map variable ID returned from adios driver */
    int piodecomid;  /* map IDs used on Scorpio starting at 512 */
    int64_t dims[3]; /* dimension sizes */
    int ndims;       /* number of dimensions */

    int decomp_id;      /* decomposition map ID, e.g. 0, 1, 2, ... */
    int isRecVar;       /* whether is a record variable */
    size_t vlen;        /* length to be written by this rank */
    MPI_Datatype itype; /* memory buffer of internal data type */
} var_meta;

class e3sm_io_case {
    protected:
        /* TODO: merge below buffer pointers for read test with wr_buf */
        double *dbl_buf_h0 = NULL, *dbl_buf_h1 = NULL;
        vtype  *rec_buf_h0 = NULL, *rec_buf_h1 = NULL;
        char   *txt_buf = NULL;
        int    *int_buf = NULL;
        double *D1_rec_dbl_buf = NULL,
               *D3_rec_dbl_buf = NULL,
               *D4_rec_dbl_buf = NULL,
               *D5_rec_dbl_buf = NULL,
               *D6_rec_dbl_buf = NULL,
               *D1_fix_dbl_buf = NULL;
        int    *D1_fix_int_buf = NULL,
               *D2_fix_int_buf = NULL,
               *D3_fix_int_buf = NULL,
               *D4_fix_int_buf = NULL,
               *D5_fix_int_buf = NULL;

        /* dimension IDs used to define blob variables */
        int fix_dimids[MAX_NUM_DECOMP];    /* fixed-size variables are 1D */
        int rec_dimids[MAX_NUM_DECOMP][2]; /* record     variables are 2D */

        io_buffers  wr_buf;  /* write buffers and their length metadata */
        var_meta   *vars;    /* variable metadata */

        int var_wr_case(e3sm_io_config &cfg,
                        e3sm_io_decom  &decom,
                        e3sm_io_driver &driver,
                        case_meta      *cmeta);

        int  def_F_case(e3sm_io_config &cfg,
                        e3sm_io_decom  &decom,
                        e3sm_io_driver &driver,
                        int             ncid);

        int  def_G_case(e3sm_io_config &cfg,
                        e3sm_io_decom  &decom,
                        e3sm_io_driver &driver,
                        int             ncid);

        int  def_I_case(e3sm_io_config &cfg,
                        e3sm_io_decom  &decom,
                        e3sm_io_driver &driver,
                        int             ncid);

        int check_malloc(e3sm_io_config *cfg,
                         e3sm_io_driver *driver);

        void wr_buf_init(int gap);

        int wr_buf_malloc(e3sm_io_config &cfg,
                          int ffreq);

        void wr_buf_free(void);

        int def_var_decomp(e3sm_io_config &cfg,
                           e3sm_io_decom  &decom,
                           e3sm_io_driver &driver,
                           int             ncid,
                           int             dim_time,
                           int             dim_nblobs,
                           int             dim_max_nreqs[MAX_NUM_DECOMP],
                           int             g_dimids[MAX_NUM_DECOMP][4]);

    public:
         e3sm_io_case();
        ~e3sm_io_case();

        int wr_test(e3sm_io_config &cfg,
                    e3sm_io_decom  &decom,
                    e3sm_io_driver &driver);

        int rd_test(e3sm_io_config &cfg,
                    e3sm_io_decom  &decom,
                    e3sm_io_driver &driver);
};

/*---- functions for F case -------------------------------------------------*/
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

/*---- wrapper functions for adios driver -----------------------------------*/
extern
int e3sm_io_scorpio_define_var(e3sm_io_config &cfg,
                               e3sm_io_decom &decom,
                               e3sm_io_driver &driver,
                               std::map<int, std::string> &dnames,
                               int decomid,
                               int fid,
                               std::string name,
                               nc_type xtype,
                               int ndims,
                               int *dimids,
                               var_meta *var);

extern
int e3sm_io_scorpio_write_var(e3sm_io_driver &driver,
                              int frameid,
                              int fid,
                              var_meta &var,
                              MPI_Datatype itype,
                              void *buf,
                              e3sm_io_op_mode mode);

/*---- MACROS used by header define functions -------------------------------*/

#define CHECK_VAR_ERR(varid) {                                                \
    if (err != 0) {                                                           \
        char var_name[64];                                                    \
        driver.inq_var_name(ncid, varid, var_name);                           \
        printf("Error in %s:%d: %s() var %s\n",                               \
               __FILE__, __LINE__, __func__, var_name);                       \
        goto err_out;                                                         \
    }                                                                         \
}
#define DEF_DIM(name, num, dimid) {                                           \
    err = driver.def_dim(ncid, name, num, dimid);                             \
    CHECK_ERR                                                                 \
    if (cfg.api == adios) dnames[*dimid] = name;                              \
}
#define PUT_GATTR_TXT(name, buf) {                                            \
    err = driver.put_att(ncid, NC_GLOBAL, prefix name, NC_CHAR, strlen(buf),  \
                         buf);                                                \
    CHECK_ERR                                                                 \
}
#define PUT_GATTR_INT(name, val) {                                            \
    int buf = val;                                                            \
    err = driver.put_att(ncid, NC_GLOBAL, prefix name, NC_INT, 1, &buf);      \
    CHECK_ERR                                                                 \
}
#define PUT_GATTR_DBL(name, val) {                                            \
    double buf = val;                                                         \
    err = driver.put_att(ncid, NC_GLOBAL, prefix name, NC_DOUBLE, 1, &buf);   \
    CHECK_ERR                                                                 \
}
#define PUT_ATTR_TXT(name, buf) {                                             \
    err = driver.put_att(ncid, varp->vid, name, NC_CHAR, strlen(buf), buf);   \
    CHECK_VAR_ERR(varp->vid)                                                  \
}
#define PUT_ATTR_INT1(name, val) {                                            \
    int buf = val;                                                            \
    err = driver.put_att(ncid, varp->vid, name, NC_INT, 1, &buf);             \
    CHECK_VAR_ERR(varp->vid)                                                  \
}
#define PUT_ATTR_INT(name, num, buf) {                                        \
    err = driver.put_att(ncid, varp->vid, name, NC_INT, num, buf);            \
    CHECK_VAR_ERR(varp->vid)                                                  \
}
#define PUT_ATTR_FLT(name, val) {                                             \
    float buf = val;                                                          \
    err = driver.put_att(ncid, varp->vid, name, NC_FLOAT, 1, &buf);           \
    CHECK_VAR_ERR(varp->vid)                                                  \
}
#define PUT_ATTR_INT64(name, num, buf) {                                      \
    err = driver.put_att(ncid, varp->vid, name, NC_INT64, num, buf);          \
    CHECK_VAR_ERR(varp->vid)                                                  \
}
#define SET_BUF_META(dtype, buflen, varlen) {                                 \
    int id = varp->decomp_id;                                                 \
    size_t vlen = (cfg.api == adios && id >= 0) ? decom.raw_nreqs[id]         \
                                                : varlen;                     \
    varp->itype    = dtype; /* internal data type of write buffer */          \
    varp->vlen     = vlen;  /* length of this write request */                \
    wr_buf.buflen += vlen + wr_buf.gap;                                       \
}
#define DEF_VAR(name, xtype, ndims, dimids, decomid) {                        \
    /* ndims and dimids are canonical dimensions */                           \
    int *_dimids = dimids;                                                    \
    varp++;                                                                   \
    varp->decomp_id = decomid;  /* decomposition map ID */                    \
    varp->isRecVar  = (ndims != 0 && *_dimids == dim_time);                   \
    if (cfg.api == adios) {                                                   \
        err = e3sm_io_scorpio_define_var(cfg, decom, driver, dnames, decomid, \
                                         ncid, name, xtype, ndims, dimids,    \
                                         varp);                               \
    } else if (cfg.strategy == blob && decomid >= 0) {                        \
        /* use blob dimensions to define blob variables */                    \
        int ival, _ndims;                                                     \
        if (varp->isRecVar) {                                                 \
            _ndims = 2;  /* all blob record variables are 2D */               \
            _dimids = rec_dimids[decomid];                                    \
        } else {                                                              \
            _ndims = 1;  /* all blob fixed-size variables are 1D */           \
            _dimids = &fix_dimids[decomid];                                   \
        }                                                                     \
        err = driver.def_var(ncid, name, xtype, _ndims, _dimids, &varp->vid); \
        /* save the canonical domensions as attributes */                     \
        ival = decomid + 1;                                                   \
        PUT_ATTR_INT("decomposition_ID", 1, &ival)                            \
        PUT_ATTR_INT("global_dimids", ndims, dimids)                          \
    } else { /* cfg.strategy == canonical or log or decomid == -1 */          \
        err = driver.def_var(ncid, name, xtype, ndims, dimids, &varp->vid);   \
    }                                                                         \
    if (err != 0) {                                                           \
        printf("Error in %s line %d: def_var %s\n", __FILE__, __LINE__,       \
               name);                                                         \
        goto err_out;                                                         \
    }                                                                         \
}
