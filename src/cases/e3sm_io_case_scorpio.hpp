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
#include <e3sm_io_driver.hpp>
#include <e3sm_io_driver_adios2.hpp>
#include <e3sm_io_driver_pnc.hpp>

typedef struct e3sm_io_scorpio_var {
    int data;
    int frame_id;
    int decomp_id;
    int fillval_id;
    int piodecomid;
    int64_t bsize[3];
    int ndim;
} e3sm_io_scorpio_var;

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

err_out:;
    return nerrs;
}

inline int e3sm_io_scorpio_define_var (e3sm_io_driver &driver,
                                   e3sm_io_config &cfg,
                                   std::map<int, std::string> &dnames,
                                   e3sm_io_decom &decom,
                                   int decomid,
                                   int pio_decomid,
                                   int fid,
                                   std::string name,
                                   nc_type xtype,
                                   int ndim,
                                   int *dimids,
                                   e3sm_io_scorpio_var *var) {
    int err, nerrs = 0;
    int i, j;
    char cbuf[64];
    int ibuf;
    std::vector<const char*> dnames_array (ndim);

    var->piodecomid = pio_decomid + 512;
    var->ndim = 0;

    for(i = 0; i < ndim; i++){
        dnames_array[i] = dnames[dimids[i]].c_str();
    }

    // If there is a decomposition map associated with the variable,
    // created 2 associated scalar variables frame_id (timesteps) and decom_id (decomposition map)
    if (decomid >= 0) {
        MPI_Offset one = 1;
        
        err = driver.def_local_var (fid, name, xtype, 1, &(decom.raw_nreqs[decomid]), &(var->data));
        CHECK_ERR

        if (ndim > 0) {
            err =
                driver.def_local_var (fid, "frame_id/" + name, NC_INT, 1, &one, &(var->frame_id));
            CHECK_ERR

            err = driver.def_local_var (fid, "decomp_id/" + name, NC_INT, 1, &one,
                                        &(var->decomp_id));
            CHECK_ERR

            // Double vars have an additional fillval_id
            if (xtype == NC_DOUBLE) {
                err = driver.def_local_var (fid, "fillval_id/" + name, xtype, 1, &one,
                                            &(var->fillval_id));
                CHECK_ERR
            } else {
                var->fillval_id = -1;
            }

            // Attributes
            // err = driver.put_att (fid, var->data, "_FillValue", var->type, 1, cbuf);
            // CHECK_ERR

            // Scorpio attributes are only written by rank 0
            if (cfg.rank == 0) {
                // Decomposition map
                sprintf(cbuf, "%d", var->piodecomid);
                err  = driver.put_att (fid, var->data, "__pio__/decomp", NC_CHAR, strlen(cbuf), &cbuf);
                CHECK_ERR

                err =
                    driver.put_att (fid, var->data, "__pio__/dims", NC_STRING, ndim, dnames_array.data());
                CHECK_ERR

                // Type of NetCDF API called
                err =
                    driver.put_att (fid, var->data, "__pio__/ncop", NC_CHAR, 7, (void *)"darray");
                CHECK_ERR

                // NetCDF type enum
                ibuf = xtype;
                err  = driver.put_att (fid, var->data, "__pio__/nctype", NC_INT, 1, &ibuf);
                CHECK_ERR

                // Number of dimensions
                err = driver.put_att (fid, var->data, "__pio__/ndims", NC_INT, 1, &ndim);
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

        // flatten into 1 dim only apply to non-scalar variables
        if (ndim){
            // Record block size to be written into the data later
            var->ndim = ndim;
            for (i = 0; i < ndim; i++) { 
                var->bsize[i] = dsize[ndim - i - 1]; 
                // Time dim is always 0, but block size should be 1
                if (var->bsize[i] == 0){
                    var->bsize[i] = 1;
                }
            }
            // Convert into byte array
            for (i = 0; i < j; i++) { vsize *= dsize[i]; }
            err = e3sm_io_xlen_nc_type(xtype, &esize);
            CHECK_MPIERR
            vsize *= esize;
            vsize += 8 * 2 * ndim; // Include start and count array
            err = driver.def_local_var (fid, name, NC_BYTE, 1, &vsize, &(var->data));
            CHECK_ERR
        }
        else{
            err = driver.def_local_var (fid, name, xtype, ndim, NULL, &(var->data));
            CHECK_ERR
        }

        // Attributes for non-constant small vars
        // Scorpio attributes are only written by rank 0
        if (cfg.rank == 0) {
            // ADIOS type enum, variables without decomposition map are stored as byte array
            ibuf = (int)e3sm_io_type_nc2adios (xtype);
            err  = driver.put_att (fid, var->data, "__pio__/adiostype", NC_INT, 1, &ibuf);
            CHECK_ERR

            // Scalar var does not have dims
            if (ndim > 0) {
                err =
                    driver.put_att (fid, var->data, "__pio__/dims", NC_STRING, ndim, dnames_array.data());
                CHECK_ERR
            }
            
            // Type of NetCDF API called
            err =
                driver.put_att (fid, var->data, "__pio__/ncop", NC_CHAR, 7, (void *)"put_var");
            CHECK_ERR

            // NetCDF type enum
            ibuf = xtype;
            err  = driver.put_att (fid, var->data, "__pio__/nctype", NC_INT, 1, &ibuf);
            CHECK_ERR

            // Number of dimensions
            err = driver.put_att (fid, var->data, "__pio__/ndims", NC_INT, 1, &ndim);
            CHECK_ERR
        }

        var->decomp_id  = -1;
        var->frame_id   = -1;
        var->fillval_id = -1;
    }
err_out:;
    return nerrs;
}

inline int e3sm_io_scorpio_write_var (e3sm_io_driver &driver,
                                  int frameid,
                                  int fid,
                                  e3sm_io_scorpio_var &var,
                                  MPI_Datatype itype,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err, nerrs = 0;
    int decomid;

    // Attach start and count before the data for small, non-scalar variables
    if (var.ndim){
        memset(buf , 0, var.ndim * sizeof(int64_t));
        memcpy((char*)buf + var.ndim * sizeof(int64_t), var.bsize, var.ndim * sizeof(int64_t));
    }

    // Small variables are stored as byte array
    if ((var.frame_id < 0) && var.ndim){
        itype = MPI_BYTE;
    }

    err = driver.put_varl (fid, var.data, itype, buf, nbe);
    CHECK_ERR

    if (var.frame_id >= 0) {

        err = driver.put_varl (fid, var.frame_id, MPI_INT, &frameid, nbe);
        CHECK_ERR

        decomid = var.piodecomid ;
        if (var.fillval_id < 0) {
            decomid *= -1;
        }
        err = driver.put_varl (fid, var.decomp_id, MPI_INT, &(decomid), nbe);
        CHECK_ERR

        if (var.fillval_id >= 0) {
            double fbuf = 1e+20;

            err = driver.put_varl (fid, var.fillval_id, MPI_DOUBLE, &fbuf, nbe);
            CHECK_ERR
        }
    }

err_out:;
    return nerrs;
}

inline int e3sm_io_scorpio_put_att (e3sm_io_driver &driver,
                                int fid,
                                int vid,
                                std::string name,
                                nc_type xtype,
                                MPI_Offset size,
                                void *buf) {
    return driver.put_att (fid, vid, name, xtype, size, buf);
}

inline int e3sm_io_scorpio_put_att (e3sm_io_driver &driver,
                                int fid,
                                e3sm_io_scorpio_var &var,
                                std::string name,
                                nc_type xtype,
                                MPI_Offset size,
                                void *buf) {
    return driver.put_att (fid, var.data, name, xtype, size, buf);
}

extern int
def_F_case_h0_scorpio(e3sm_io_driver &driver,
                      e3sm_io_config &cfg,
                      e3sm_io_decom &decom,
                      int ncid,                 /* file ID */
                      const MPI_Offset dims[2], /* dimension sizes */
                      int nvars,                /* number of variables */
                      std::vector<int> &decomids,
                      e3sm_io_scorpio_var *varids,
                      int *scorpiovars);

extern int
def_F_case_h1_scorpio(e3sm_io_driver &driver,
                      e3sm_io_config &cfg,
                      e3sm_io_decom &decom,
                      int ncid,                 /* file ID */
                      const MPI_Offset dims[2], /* dimension sizes */
                      int nvars,                /* number of variables */
                      std::vector<int> &decomids,
                      e3sm_io_scorpio_var *varids,
                      int *scorpiovars);

extern int
run_varn_F_case_scorpio(e3sm_io_config &cfg,
                        e3sm_io_decom &decom,
                        e3sm_io_driver &driver,
                        double *dbl_bufp,    /* buffer for fixed size double var */
                        vtype *rec_bufp,     /* buffer for rec floating point var */
                        char *txt_buf,       /* buffer for char var */
                        int *int_buf);       /* buffer for int var */

extern int
def_G_case_scorpio(e3sm_io_config &cfg,
                   e3sm_io_decom &decom,
                   e3sm_io_driver &driver,
                   int ncid, /* file ID */
                   std::vector<int> &decomids,
                   e3sm_io_scorpio_var *varids, /* variable IDs */
                   int *scorpiovars);

extern int
run_varn_G_case_scorpio(e3sm_io_config &cfg,
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

