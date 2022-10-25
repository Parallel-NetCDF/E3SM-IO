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
    size_t  fix_txt_buflen; char      *fix_txt_buf;
    size_t  fix_int_buflen; int       *fix_int_buf;
    size_t  fix_flt_buflen; float     *fix_flt_buf;
    size_t  fix_dbl_buflen; double    *fix_dbl_buf;
    size_t  fix_lld_buflen; long long *fix_lld_buf;

    /* buffers for record variables */
    size_t  rec_txt_buflen; char      *rec_txt_buf;
    size_t  rec_int_buflen; int       *rec_int_buf;
    size_t  rec_flt_buflen; float     *rec_flt_buf;
    size_t  rec_dbl_buflen; double    *rec_dbl_buf;
    size_t  rec_lld_buflen; long long *rec_lld_buf;
} io_buffers;

typedef struct {
    int vid;         /* variable ID, returned from the driver */

    char *_name;     /* name of variable */
    int fill_id;     /* fill value variable ID returned from adios driver */
    int frame_id;    /* frame variable ID returned from adios driver */
    int decom_id;    /* decomposition map variable ID returned from adios driver */
    int num_data_block_writers; /* number of data block_writers returned from adios driver */
    int piodecomid;  /* map IDs used on Scorpio starting at 512 */
    int ndims;                  /* number of dimensions */
    MPI_Offset dims[MAX_NDIMS]; /* dimension sizes */
    nc_type xType;

    int decomp_id;      /* decomposition map ID, e.g. 0, 1, 2, ... */
    int isRecVar;       /* whether is a record variable */
    size_t vlen;        /* length to be written by this rank */
    MPI_Datatype iType; /* memory buffer of internal data type */
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
        int fix_dimids[MAX_NUM_DECOMP];    /* fixed-size blob variables are 1D */
        int rec_dimids[MAX_NUM_DECOMP][2]; /* record     blob variables are 2D */

        io_buffers  wr_buf;  /* write buffers and their length metadata */
        var_meta   *vars;    /* variable metadata */

        int var_wr_case(e3sm_io_config &cfg,
                        e3sm_io_decom  &decom,
                        e3sm_io_driver &driver,
                        case_meta      *cmeta);
        
        int var_rd_case(e3sm_io_config &cfg,
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

        int  inq_F_case(e3sm_io_config &cfg,
                        e3sm_io_decom  &decom,
                        e3sm_io_driver &driver,
                        int             ncid);

        int  inq_G_case(e3sm_io_config &cfg,
                        e3sm_io_decom  &decom,
                        e3sm_io_driver &driver,
                        int             ncid);

        int  inq_I_case(e3sm_io_config &cfg,
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
                           case_meta      *cmeta,
                           int ncid,
                           int dim_time,
                           int dim_nblobs,
                           int dim_max_nreqs[MAX_NUM_DECOMP],
                           int g_dimids[MAX_NUM_DECOMP][MAX_NDIMS]);

        int inq_var_decomp(e3sm_io_config &cfg,
                           e3sm_io_decom  &decom,
                           e3sm_io_driver &driver,
                           case_meta      *cmeta,
                           int ncid,
                           int dim_time,
                           int dim_nblobs,
                           int dim_max_nreqs[MAX_NUM_DECOMP],
                           int g_dimids[MAX_NUM_DECOMP][MAX_NDIMS]);

        int def_var(e3sm_io_config             &cfg,
                    e3sm_io_decom              &decom,
                    e3sm_io_driver             &driver,
                    case_meta                  *cmeta,
                    int                         ncid,
                    std::string                 name,
                    std::map<int, std::string> &dnames,
                    int                         xtype,
                    int                         nDims,
                    int                         dim_time,
                    int                        *dimids,
                    MPI_Datatype                itype,
                    int                         decomid,
                    var_meta                   *varp);

        int inq_var(e3sm_io_config             &cfg,
                    e3sm_io_decom              &decom,
                    e3sm_io_driver             &driver,
                    case_meta                  *cmeta,
                    int                         ncid,
                    const char                 *name,
                    int                         dim_time,
                    int                        *dimids,
                    MPI_Datatype                itype,
                    int                         decomid,
                    var_meta                   *varp);

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
int scorpio_define_var(e3sm_io_config &cfg,
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
int scorpio_write_var(e3sm_io_driver &driver,
                      int frameid,
                      int fid,
                      var_meta &var,
                      MPI_Datatype itype,
                      void *buf,
                      e3sm_io_op_mode mode);

/*---- MACROS used by header define functions -------------------------------*/

#define CHECK_VAR_ERR(var_name) {                                             \
    if (err != 0) {                                                           \
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
    err = driver.put_att(ncid, NC_GLOBAL, prefix+name, NC_CHAR, strlen(buf),  \
                         buf);                                                \
    CHECK_ERR                                                                 \
    cmeta->num_attrs++;                                                       \
}
#define PUT_GATTR_INT(name, val) {                                            \
    int buf = val;                                                            \
    err = driver.put_att(ncid, NC_GLOBAL, prefix+name, NC_INT, 1, &buf);      \
    CHECK_ERR                                                                 \
    cmeta->num_attrs++;                                                       \
}
#define PUT_GATTR_DBL(name, val) {                                            \
    double buf = val;                                                         \
    err = driver.put_att(ncid, NC_GLOBAL, prefix+name, NC_DOUBLE, 1, &buf);   \
    CHECK_ERR                                                                 \
    cmeta->num_attrs++;                                                       \
}
#define PUT_ATTR_TXT(name, buf) {                                             \
    err = driver.put_att(ncid, varp->vid, name, NC_CHAR, strlen(buf), buf);   \
    CHECK_VAR_ERR(varp->_name)                                                \
    cmeta->num_attrs++;                                                       \
}
#define PUT_ATTR_INT1(name, val) {                                            \
    int buf = val;                                                            \
    err = driver.put_att(ncid, varp->vid, name, NC_INT, 1, &buf);             \
    CHECK_VAR_ERR(varp->_name)                                                \
    cmeta->num_attrs++;                                                       \
}
#define PUT_ATTR_INT(name, num, buf) {                                        \
    err = driver.put_att(ncid, varp->vid, name, NC_INT, num, buf);            \
    CHECK_VAR_ERR(varp->_name)                                                \
    cmeta->num_attrs++;                                                       \
}
#define PUT_ATTR_FLT1(name, val) {                                            \
    float buf = val;                                                          \
    err = driver.put_att(ncid, varp->vid, name, NC_FLOAT, 1, &buf);           \
    CHECK_VAR_ERR(varp->_name)                                                \
    cmeta->num_attrs++;                                                       \
}
#define PUT_ATTR_DBL1(name, val) {                                            \
    double buf = val;                                                         \
    err = driver.put_att(ncid, varp->vid, name, NC_DOUBLE, 1, &buf);          \
    CHECK_VAR_ERR(varp->_name)                                                \
}
#define PUT_ATTR_FILL(val) {                                                  \
    if (varp->xType == NC_FLOAT) {                                            \
        float buf = (float)val;                                               \
        err = driver.put_att(ncid,varp->vid,_FillValue,varp->xType, 1, &buf); \
    }                                                                         \
    else if (varp->xType == NC_INT) {                                         \
        int buf = (int)val;                                                   \
        err = driver.put_att(ncid,varp->vid,_FillValue,varp->xType, 1, &buf); \
    }                                                                         \
    else if (varp->xType == NC_DOUBLE) {                                      \
        double buf = (double)val;                                             \
        err = driver.put_att(ncid,varp->vid,_FillValue,varp->xType, 1, &buf); \
    }                                                                         \
    CHECK_VAR_ERR(varp->_name)                                                \
    cmeta->num_attrs++;                                                       \
}
#define PUT_ATTR_INT64(name, num, buf) {                                      \
    err = driver.put_att(ncid, varp->vid, name, NC_INT64, num, buf);          \
    CHECK_VAR_ERR(varp->_name)                                                \
    cmeta->num_attrs++;                                                       \
}
#define DEF_VAR(name, xtype, nDims, dimids, itype, decomid) {                 \
    varp++;                                                                   \
    err = e3sm_io_case::def_var(cfg, decom, driver, cmeta, ncid, name,        \
                                dnames, xtype,  nDims, dim_time, dimids,      \
                                itype, decomid, varp);                        \
    if (err != 0) goto err_out;                                               \
}
#if 0
#define DEF_VAR(name, xtype, nDims, dimids, itype, decomid) {                 \
    /* nDims and dimids are canonical dimensions */                           \
    int _i, *_dimids = dimids;                                                \
    varp++;                                                                   \
    varp->_name     = strdup(name);                                           \
    varp->ndims     = nDims;   /* number of dimensions */                     \
    varp->iType     = itype;   /* internal data type of write buffer */       \
    varp->xType     = xtype;   /* external data type of variable in file */   \
    varp->decomp_id = decomid; /* decomposition map ID */                     \
    varp->isRecVar  = (nDims != 0 && *_dimids == dim_time);                   \
    /* calculate variable size */                                             \
    for (varp->vlen=1, _i=0; _i<nDims; _i++) {                                \
        err = driver.inq_dimlen(ncid, _dimids[_i], &varp->dims[_i]);          \
        CHECK_ERR                                                             \
        if (_i == 0 && varp->isRecVar) varp->dims[_i] = 1;                    \
        varp->vlen *= varp->dims[_i];                                         \
    }                                                                         \
    /* define a new variable */                                               \
    if (cfg.api == adios) {                                                   \
        err = scorpio_define_var(cfg, decom, driver, dnames, decomid, ncid,   \
                                 name, xtype, nDims, dimids, varp);           \
        if (decomid >= 0) varp->vlen = decom.raw_nreqs[decomid];              \
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
        /* save the canonical dimensions as attributes */                     \
        ival = decomid + 1;                                                   \
        PUT_ATTR_INT("decomposition_ID", 1, &ival)                            \
        PUT_ATTR_INT("global_dimids", nDims, dimids)                          \
        varp->vlen = decom.count[decomid];                                    \
    } else { /* cfg.strategy == canonical or log or decomid == -1 */          \
        err = driver.def_var(ncid, name, xtype, nDims, dimids, &varp->vid);   \
        if (decomid >= 0) varp->vlen = decom.count[decomid];                  \
    }                                                                         \
    if (err != 0) {                                                           \
        printf("Error in %s line %d: def_var %s\n", __FILE__, __LINE__,name); \
        goto err_out;                                                         \
    }                                                                         \
    /* increment I/O buffer sizes */                                          \
    if (varp->isRecVar) {                                                     \
        if (varp->iType == MPI_DOUBLE)                                        \
            wr_buf.rec_dbl_buflen   += varp->vlen + wr_buf.gap;               \
        else if (varp->iType == MPI_INT)                                      \
            wr_buf.rec_int_buflen   += varp->vlen + wr_buf.gap;               \
        else if (varp->iType == MPI_CHAR)                                     \
            wr_buf.rec_txt_buflen   += varp->vlen + wr_buf.gap;               \
        else if (varp->iType == MPI_FLOAT)                                    \
            wr_buf.rec_flt_buflen   += varp->vlen + wr_buf.gap;               \
        else if (varp->iType == MPI_LONG_LONG)                                \
            wr_buf.rec_lld_buflen += varp->vlen + wr_buf.gap;                 \
        else assert(0)                                                        \
    } else {                                                                  \
        if (varp->iType == MPI_DOUBLE)                                        \
            wr_buf.fix_dbl_buflen   += varp->vlen + wr_buf.gap;               \
        else if (varp->iType == MPI_INT)                                      \
            wr_buf.fix_int_buflen   += varp->vlen + wr_buf.gap;               \
        else if (varp->iType == MPI_CHAR)                                     \
            wr_buf.fix_txt_buflen   += varp->vlen + wr_buf.gap;               \
        else if (varp->iType == MPI_FLOAT)                                    \
            wr_buf.fix_flt_buflen   += varp->vlen + wr_buf.gap;               \
        else if (varp->iType == MPI_LONG_LONG)                                \
            wr_buf.fix_lld_buflen += varp->vlen + wr_buf.gap;                 \
        else assert(0)                                                        \
    }                                                                         \
}
#endif
#define INQ_DIM(name, num, dimid) {                                           \
    err = driver.inq_dim(ncid, name, dimid);                                  \
    CHECK_ERR                                                                 \
    if (cfg.api == adios) dnames[*dimid] = name;                              \
}
#define GET_GATTR_TXT(name, buf) {                                            \
    err = driver.get_att(ncid, NC_GLOBAL, prefix+name, buf);                  \
    CHECK_ERR                                                                 \
    cmeta->num_attrs++;                                                       \
}
#define GET_GATTR_INT(name, val) {                                            \
    err = driver.get_att(ncid, NC_GLOBAL, prefix+name,  val);                 \
    CHECK_ERR                                                                 \
    cmeta->num_attrs++;                                                       \
}
#define GET_GATTR_DBL(name, val) {                                            \
    double buf = val;                                                         \
    err = driver.get_att(ncid, NC_GLOBAL, prefix+name, &buf);                 \
    CHECK_ERR                                                                 \
    cmeta->num_attrs++;                                                       \
}
#define GET_ATTR_TXT(name, dbuf) {                                            \
    err = driver.get_att(ncid, varp->vid, name, dbuf);                        \
    CHECK_VAR_ERR(varp->_name)                                                \
    cmeta->num_attrs++;                                                       \
}
#define GET_ATTR_INT1(name, val) {                                            \
    err = driver.get_att(ncid, varp->vid, name, val);                         \
    CHECK_VAR_ERR(varp->_name)                                                \
    cmeta->num_attrs++;                                                       \
}
#define GET_ATTR_INT(name, num, dbuf) {                                       \
    err = driver.get_att(ncid, varp->vid, name, dbuf);                        \
    CHECK_VAR_ERR(varp->_name)                                                \
    cmeta->num_attrs++;                                                       \
}
#define GET_ATTR_FLT1(name, val) {                                            \
    err = driver.get_att(ncid, varp->vid, name, val);                         \
    CHECK_VAR_ERR(varp->_name)                                                \
    cmeta->num_attrs++;                                                       \
}
#define GET_ATTR_DBL1(name, val) {                                            \
    err = driver.get_att(ncid, varp->vid, name, val);                         \
    CHECK_VAR_ERR(varp->_name)                                                \
}
#define GET_ATTR_FILL(val) {                                                  \
    if (varp->xType == NC_FLOAT) {                                            \
        err = driver.get_att(ncid,varp->vid,_FillValue,&val);                 \
    }                                                                         \
    else if (varp->xType == NC_INT) {                                         \
        int buf;                                                              \
        err = driver.get_att(ncid,varp->vid,_FillValue,&buf);                 \
        val = (float)buf;                                                     \
    }                                                                         \
    else if (varp->xType == NC_DOUBLE) {                                      \
        double buf;                                                           \
        err = driver.get_att(ncid,varp->vid,_FillValue,&buf);                 \
        val = (float)buf;                                                     \
    }                                                                         \
    CHECK_VAR_ERR(varp->_name)                                                \
    cmeta->num_attrs++;                                                       \
}
#define GET_ATTR_INT64(name, num, dbuf) {                                     \
    err = driver.get_att(ncid, varp->vid, name, dbuf);                        \
    CHECK_VAR_ERR(varp->_name)                                                \
    cmeta->num_attrs++;                                                       \
}
#define INQ_VAR(name, xtype, nDims, dimids, itype, decomid) {                 \
    varp++;                                                                   \
    err = e3sm_io_case::inq_var(cfg, decom, driver, cmeta, ncid, (char*)name, \
                                dim_time, dimids, itype, decomid, varp);      \
    if (err != 0) goto err_out;                                               \
}
#if 0
#define INQ_VAR(name, xtype, nDims, dimids, itype, decomid) {                 \
    /* nDims and dimids are canonical dimensions */                           \
    int _i, *_dimids = dimids;                                                \
    varp++;                                                                   \
    varp->_name     = strdup(name);                                           \
    varp->ndims     = nDims;   /* number of dimensions */                     \
    varp->iType     = itype;   /* internal data type of write buffer */       \
    varp->xType     = xtype;   /* external data type of variable in file */   \
    varp->decomp_id = decomid; /* decomposition map ID */                     \
    varp->isRecVar  = (nDims != 0 && *_dimids == dim_time);                   \
    /* calculate variable size */                                             \
    for (varp->vlen=1, _i=0; _i<nDims; _i++) {                                \
        err = driver.inq_dimlen(ncid, _dimids[_i], &varp->dims[_i]);          \
        CHECK_ERR                                                             \
        if (_i == 0 && varp->isRecVar) varp->dims[_i] = 1;                    \
        varp->vlen *= varp->dims[_i];                                         \
    }                                                                         \
    /* define a new variable */                                               \
    err = driver.inq_var(ncid, name, &varp->vid);                             \
    if (decomid >= 0) varp->vlen = decom.count[decomid];                      \
                                                                              \
    if (err != 0) {                                                           \
        printf("Error in %s line %d: def_var %s\n", __FILE__, __LINE__,name); \
        goto err_out;                                                         \
    }                                                                         \
    /* increment I/O buffer sizes */                                          \
    if (varp->isRecVar) {                                                     \
        if (varp->iType == MPI_DOUBLE)                                        \
            wr_buf.rec_dbl_buflen   += varp->vlen + wr_buf.gap;               \
        else if (varp->iType == MPI_INT)                                      \
            wr_buf.rec_int_buflen   += varp->vlen + wr_buf.gap;               \
        else if (varp->iType == MPI_CHAR)                                     \
            wr_buf.rec_txt_buflen   += varp->vlen + wr_buf.gap;               \
        else if (varp->iType == MPI_FLOAT)                                    \
            wr_buf.rec_flt_buflen   += varp->vlen + wr_buf.gap;               \
        else if (varp->iType == MPI_LONG_LONG)                                \
            wr_buf.rec_lld_buflen += varp->vlen + wr_buf.gap;                 \
        else assert(0)                                                        \
    } else {                                                                  \
        if (varp->iType == MPI_DOUBLE)                                        \
            wr_buf.fix_dbl_buflen   += varp->vlen + wr_buf.gap;               \
        else if (varp->iType == MPI_INT)                                      \
            wr_buf.fix_int_buflen   += varp->vlen + wr_buf.gap;               \
        else if (varp->iType == MPI_CHAR)                                     \
            wr_buf.fix_txt_buflen   += varp->vlen + wr_buf.gap;               \
        else if (varp->iType == MPI_FLOAT)                                    \
            wr_buf.fix_flt_buflen   += varp->vlen + wr_buf.gap;               \
        else if (varp->iType == MPI_LONG_LONG)                                \
            wr_buf.fix_lld_buflen += varp->vlen + wr_buf.gap;                 \
        else assert(0)                                                        \
    }                                                                         \
}
#endif

