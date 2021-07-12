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
#include <mpi.h>
//
#include <string>
//
#include <e3sm_io.h>
#include <e3sm_io_driver.hpp>

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
                            double *dbl_bufp,    /* buffer for fixed size double var */
                            itype *rec_bufp,     /* buffer for rec floating point var */
                            char *txt_buf,       /* buffer for char var */
                            int *int_buf);       /* request's block lengths */

extern int run_varn_F_case (e3sm_io_config &cfg,
                            e3sm_io_decom &decom,
                            e3sm_io_driver &driver,
                            double *dbl_bufp,    /* buffer for fixed size double var */
                            itype *rec_bufp,     /* buffer for rec floating point var */
                            char *txt_buf,       /* buffer for char var */
                            int *int_buf);       /* buffer for int var */
extern int run_varn_F_case_rd (e3sm_io_config &cfg,
                               e3sm_io_decom &decom,
                               e3sm_io_driver &driver,
                               double **dbl_bufp,   /* buffer for fixed size double var */
                               itype **rec_bufp,    /* buffer for rec floating point var */
                               char *txt_buf,       /* buffer for char var */
                               int *int_buf);       /* buffer for int var */

extern int
pnetcdf_blob_F_case(e3sm_io_config &cfg,
                    e3sm_io_decom &decom,
                    e3sm_io_driver &driver);

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

#include <e3sm_io_err.h>
#include <vector>
#include <map>

#include <e3sm_io_driver_adios2.hpp>
#include <e3sm_io_driver_pnc.hpp>

typedef struct e3sm_io_scorpio_var {
    int data;
    int frame_id;
    int decomp_id;
    int fillval_id;
    MPI_Datatype type;
    int decomid;
} e3sm_io_scorpio_var;

int def_F_case_h0_scorpio (e3sm_io_driver &driver,
                       e3sm_io_config &cfg,
                       e3sm_io_decom &decom,
                       int ncid,                 /* file ID */
                       const MPI_Offset dims[2], /* dimension sizes */
                       int nvars,                /* number of variables */
                       std::vector<int> &decomids,
                       e3sm_io_scorpio_var *varids,
                       int *scorpiovars);

int def_F_case_h1_scorpio (e3sm_io_driver &driver,
                       e3sm_io_config &cfg,
                       e3sm_io_decom &decom,
                       int ncid,                 /* file ID */
                       const MPI_Offset dims[2], /* dimension sizes */
                       int nvars,                /* number of variables */
                       std::vector<int> &decomids,
                       e3sm_io_scorpio_var *varids,
                       int *scorpiovars);

int run_varn_F_case_scorpio (e3sm_io_config &cfg,
                         e3sm_io_decom &decom,
                         e3sm_io_driver &driver,
                         double *dbl_bufp,    /* buffer for fixed size double var */
                         itype *rec_bufp,     /* buffer for rec floating point var */
                         char *txt_buf,       /* buffer for char var */
                         int *int_buf);       /* buffer for int var */

inline int e3sm_io_scorpio_define_dim (e3sm_io_driver &driver,
                                   int fid,
                                   std::string name,
                                   MPI_Offset size,
                                   std::map<int, std::string> &dnames,
                                   int *did) {
    int err, nerrs = 0;

    err = driver.def_dim (fid, name, size, did);
    CHECK_ERR

    dnames[*did] = name;

err_out:
    return nerrs;
}

inline int e3sm_io_scorpio_define_var (e3sm_io_driver &driver,
                                   e3sm_io_config &cfg,
                                   std::map<int, std::string> &dnames,
                                   e3sm_io_decom &decom,
                                   int decomid,
                                   int fid,
                                   std::string name,
                                   MPI_Datatype type,
                                   int ndim,
                                   int *dimids,
                                   e3sm_io_scorpio_var *var) {
    int err, nerrs = 0;
    int i, j;
    char cbuf[64], *cbufp;
    int ibuf;
    int ret;

    var->type    = type;
    var->decomid = decomid;

    if (decomid >= 0) {
        MPI_Offset one = 1;

        err = driver.def_local_var (fid, name, type, 1, &(decom.raw_nreqs[decomid]), &(var->data));
        CHECK_ERR

        if (ndim > 0) {
            err =
                driver.def_local_var (fid, "frame_id/" + name, MPI_INT, 1, &one, &(var->frame_id));
            CHECK_ERR

            err = driver.def_local_var (fid, "decomp_id/" + name, MPI_INT, 1, &one,
                                        &(var->decomp_id));
            CHECK_ERR

            // Only double vars have fillval_id
            if (type == MPI_DOUBLE) {
                err = driver.def_local_var (fid, "fillval_id/" + name, type, 1, &one,
                                            &(var->fillval_id));
                CHECK_ERR
            } else {
                var->fillval_id = -1;
            }

            // Attributes
            // err = driver.put_att (fid, var->data, "_FillValue", var->type, 1, cbuf);
            // CHECK_ERR

            if (cfg.rank == 0) {
                ibuf = var->decomp_id + 512;
                err  = driver.put_att (fid, var->data, "__pio__/decomp", MPI_INT, 1, &ibuf);
                CHECK_ERR

                cbufp = cbuf;
                ret   = sprintf (cbufp, "{%s", dnames[dimids[0]].c_str ());
                cbufp += ret;
                for (i = 1; i < ndim; i++) {
                    ret = sprintf (cbufp, ", %s", dnames[dimids[i]].c_str ());
                    cbufp += ret;
                }
                err =
                    driver.put_att (fid, var->data, "__pio__/dims", MPI_CHAR, strlen (cbuf), cbuf);
                CHECK_ERR

                err =
                    driver.put_att (fid, var->data, "__pio__/ncop", MPI_CHAR, 7, (void *)"darray");
                CHECK_ERR

                ibuf = 5;
                err  = driver.put_att (fid, var->data, "__pio__/nctype", MPI_INT, 1, &ibuf);
                CHECK_ERR

                err = driver.put_att (fid, var->data, "__pio__/ndims", MPI_INT, 1, &ndim);
                CHECK_ERR
            }
        }
    } else {
        std::vector<MPI_Offset> dsize (ndim);
        MPI_Offset vsize = 1;
        int esize;

        for (i = j = 0; i < ndim; i++) {
            err = driver.inq_dimlen (fid, dimids[i], &(dsize[j]));
            CHECK_ERR

            // Skip time dim
            if (dsize[j] > 0) { j++; }
        }

        // flatten into 1 dim
        for (i = 0; i < j; i++) { vsize *= dsize[i]; }
        // Convert into byte array
        err = MPI_Type_size(type, &esize);
        CHECK_MPIERR
        vsize *= esize;
        vsize += 8 * 2 * ndim; // Include start and count array
        err = driver.def_local_var (fid, name, MPI_BYTE, 1, &vsize, &(var->data));
        CHECK_ERR

        // Attributes for non-constant small vars
        if (cfg.rank == 0) {
            if (ndim > 0) {
                ibuf = (int)mpi_type_to_adios2_type (type);
                err  = driver.put_att (fid, var->data, "__pio__/adiostype", MPI_INT, 1, &ibuf);
                CHECK_ERR

                cbufp = cbuf;
                ret   = sprintf (cbufp, "{%s", dnames[dimids[0]].c_str ());
                cbufp += ret;
                for (i = 1; i < ndim; i++) {
                    ret = sprintf (cbufp, ", %s", dnames[dimids[i]].c_str ());
                    cbufp += ret;
                }
                err =
                    driver.put_att (fid, var->data, "__pio__/dims", MPI_CHAR, strlen (cbuf), cbuf);
                CHECK_ERR

                err =
                    driver.put_att (fid, var->data, "__pio__/ncop", MPI_CHAR, 7, (void *)"put_var");
                CHECK_ERR

                ibuf = (int)mpitype2nctype (type);
                err  = driver.put_att (fid, var->data, "__pio__/nctype", MPI_INT, 1, &ibuf);
                CHECK_ERR

                err = driver.put_att (fid, var->data, "__pio__/ndims", MPI_INT, 1, &ndim);
                CHECK_ERR
            }
        }

        var->decomp_id  = -1;
        var->frame_id   = -1;
        var->fillval_id = -1;
    }
err_out:
    return nerrs;
}

inline int e3sm_io_scorpio_write_var (e3sm_io_driver &driver,
                                  int frameid,
                                  int fid,
                                  e3sm_io_scorpio_var &var,
                                  MPI_Datatype type,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err, nerrs = 0;

    err = driver.put_varl (fid, var.data, type, buf, mode);
    CHECK_ERR

    if (var.frame_id >= 0) {
        err = driver.put_varl (fid, var.frame_id, MPI_INT, &frameid, mode);
        CHECK_ERR

        err = driver.put_varl (fid, var.decomp_id, MPI_INT, &(var.decomid), mode);
        CHECK_ERR

        if (var.fillval_id >= 0) {
            double fbuf = 1e+20;

            err = driver.put_varl (fid, var.fillval_id, var.type, &fbuf, mode);
            CHECK_ERR
        }
    }

err_out:
    return nerrs;
}

inline int e3sm_io_scorpio_put_att (e3sm_io_driver &driver,
                                int fid,
                                int vid,
                                std::string name,
                                MPI_Datatype type,
                                MPI_Offset size,
                                void *buf) {
    return driver.put_att (fid, vid, name, type, size, buf);
}

inline int e3sm_io_scorpio_put_att (e3sm_io_driver &driver,
                                int fid,
                                e3sm_io_scorpio_var &var,
                                std::string name,
                                MPI_Datatype type,
                                MPI_Offset size,
                                void *buf) {
    return driver.put_att (fid, var.data, name, type, size, buf);
}
