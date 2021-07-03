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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//
#include <cstring>
#include <map>
//
#include <e3sm_io.h>
#include <e3sm_io_err.h>

#include <e3sm_io_case.hpp>
#include <e3sm_io_case_F.hpp>
#include <e3sm_io_driver.hpp>
#include <e3sm_io_driver_adios2.hpp>
#include <e3sm_io_driver_pnc.hpp>

typedef struct e3sm_io_pio_var {
    int data;
    int frame_id;
    int decomp_id;
    int fillval_id;
    MPI_Datatype type;
    int decomid;
} e3sm_io_pio_var;

inline int e3sm_io_pio_define_dim (e3sm_io_driver &driver,
                                   int fid,
                                   std::string name,
                                   MPI_Offset size,
                                   std::map<int, std::string> &dnames,
                                   int *did) {
    int err, nerrs = 0;

    err = driver.def_dim (fid, name, size, did);
    CHECK_ERR

    dnames[*did] = name;

err_out:;
    return nerrs;
}

inline int e3sm_io_pio_define_var (e3sm_io_driver &driver,
                                   std::map<int, std::string> &dnames,
                                   e3sm_io_decom &decom,
                                   int decomid,
                                   int fid,
                                   std::string name,
                                   MPI_Datatype type,
                                   int ndim,
                                   int *dimids,
                                   e3sm_io_pio_var *var) {
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

            ibuf = var->decomp_id + 512;
            err  = driver.put_att (fid, var->data, "__e3sm_io__/decomp", MPI_INT, 1, &ibuf, master);
            CHECK_ERR

            cbufp = cbuf;
            ret   = sprintf (cbufp, "{%s", dnames[dimids[0]].c_str ());
            cbufp += ret;
            for (i = 1; i < ndim; i++) {
                ret = sprintf (cbufp, ", %s", dnames[dimids[i]].c_str ());
                cbufp += ret;
            }
            err =
                driver.put_att (fid, var->data, "__e3sm_io__/dims", MPI_CHAR, strlen (cbuf), cbuf, master);
            CHECK_ERR

            err =
                driver.put_att (fid, var->data, "__e3sm_io__/ncop", MPI_CHAR, 7, (void *)"darray", master);
            CHECK_ERR

            ibuf = 5;
            err  = driver.put_att (fid, var->data, "__e3sm_io__/nctype", MPI_INT, 1, &ibuf, master);
            CHECK_ERR

            err = driver.put_att (fid, var->data, "__e3sm_io__/ndims", MPI_INT, 1, &ndim, master);
            CHECK_ERR
        }
    } else {
        std::vector<MPI_Offset> dsize (ndim);
        MPI_Offset vsize = 1;

        for (i = j = 0; i < ndim; i++) {
            err = driver.inq_dimlen (fid, dimids[i], &(dsize[j]));
            CHECK_ERR

            // Skip time dim
            if (dsize[j] > 0) { j++; }
        }

        // flatten into 1 dim
        for (i = 0; i < j; i++) { vsize *= dsize[i]; }
        err = driver.def_local_var (fid, name, type, 1, &vsize, &(var->data));
        CHECK_ERR

        // Attributes for non-constant small vars
        if (ndim > 0) {
            ibuf = (int)mpi_type_to_adios2_type (type);
            err  = driver.put_att (fid, var->data, "__e3sm_io__/adiostype", MPI_INT, 1, &ibuf, master);
            CHECK_ERR

            cbufp = cbuf;
            ret   = sprintf (cbufp, "{%s", dnames[dimids[0]].c_str ());
            cbufp += ret;
            for (i = 1; i < ndim; i++) {
                ret = sprintf (cbufp, ", %s", dnames[dimids[i]].c_str ());
                cbufp += ret;
            }
            err =
                driver.put_att (fid, var->data, "__e3sm_io__/dims", MPI_CHAR, strlen (cbuf), cbuf, master);
            CHECK_ERR

            err =
                driver.put_att (fid, var->data, "__e3sm_io__/ncop", MPI_CHAR, 7, (void *)"put_var", master);
            CHECK_ERR

            ibuf = (int)mpitype2nctype (type);
            err  = driver.put_att (fid, var->data, "__e3sm_io__/nctype", MPI_INT, 1, &ibuf, master);
            CHECK_ERR

            err = driver.put_att (fid, var->data, "__e3sm_io__/ndims", MPI_INT, 1, &ndim, master);
            CHECK_ERR
        }

        var->decomp_id  = -1;
        var->frame_id   = -1;
        var->fillval_id = -1;
    }
err_out:;
    return nerrs;
}

inline int e3sm_io_pio_write_var (e3sm_io_driver &driver,
                                  int frameid,
                                  int fid,
                                  e3sm_io_pio_var &var,
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

err_out:;
    return nerrs;
}

inline int e3sm_io_pio_put_att (e3sm_io_driver &driver,
                                int fid,
                                int vid,
                                std::string name,
                                MPI_Datatype type,
                                MPI_Offset size,
                                void *buf) {
    return driver.put_att (fid, vid, name, type, size, buf, coll);
}

inline int e3sm_io_pio_put_att (e3sm_io_driver &driver,
                                int fid,
                                e3sm_io_pio_var &var,
                                std::string name,
                                MPI_Datatype type,
                                MPI_Offset size,
                                void *buf) {
    return driver.put_att (fid, var.data, name, type, size, buf, coll);
}