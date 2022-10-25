/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 *********************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>

#include <e3sm_io.h>
#include <e3sm_io_err.h>
#include <e3sm_io_case.hpp>
#include <e3sm_io_driver.hpp>

#ifdef ENABLE_ADIOS2
#include <e3sm_io_driver_adios2.hpp>
#endif

int scorpio_define_var(e3sm_io_config &cfg,
                       e3sm_io_decom &decom,
                       e3sm_io_driver &driver,
                       std::map<int, std::string> &dnames,
                       int decomp_id,
                       int fid,
                       std::string name,
                       nc_type xtype,
                       int ndims,
                       int *dimids,
                       var_meta *var)
{
#ifdef ENABLE_ADIOS2
    int i, err, ibuf;
    char cbuf[64];
    std::vector<const char*> dnames_array (ndims);

    var->piodecomid = decomp_id + 512;
    var->ndims = ndims;

    for (i = 0; i < ndims; i++)
        dnames_array[i] = dnames[dimids[i]].c_str();

    // If there is a decomposition map associated with the variable,
    if (decomp_id >= 0) {
        MPI_Offset one = 1;

        err = driver.def_local_var(fid, "/__pio__/var/" + name, xtype, 1,
                                   &(decom.raw_nreqs[decomp_id]), &var->vid);
        CHECK_ERR

        assert(ndims > 0);

        /* add an associated scalar variable frame_id (i.e. timesteps) */
        if (cfg.rank == 0) {
            err = driver.def_local_var(fid, "/__pio__/track/frame_id/" + name, NC_INT, 1,
                                       &one, &var->frame_id);
            CHECK_ERR
        }

        /* add an associated scalar variable decom_id (decomposition map ID) */
        if (cfg.rank == 0) {
            err = driver.def_local_var(fid, "/__pio__/track/decomp_id/" + name, NC_INT, 1,
                                       &one, &var->decom_id);
            CHECK_ERR
        }

        /* add an associated scalar variable fillval_id (fill value ID) for double variables */
        if (xtype == NC_DOUBLE){
            if (cfg.rank == 0) {
                err = driver.def_local_var(fid, "/__pio__/track/fillval_id/" + name, NC_DOUBLE, 1,
                                           &one, &var->fill_id);
                CHECK_ERR
            }
        }
        else {
            var->fill_id = -1;
        }

        /* add an associated variable num_data_block_writers (number of data block writers) */
        err = driver.def_local_var(fid, "/__pio__/track/num_data_block_writers/" + name, NC_INT, 1,
                                   &one, &var->num_data_block_writers);
        CHECK_ERR

        /* Add attributes */

        /* In Scorpio, variable's attributes are only written by rank 0 */
        if (cfg.rank == 0) {
            /* Decomposition map ID number */
            sprintf(cbuf, "%d", var->piodecomid);
            err = driver.put_att(fid, NC_GLOBAL, "/__pio__/var/" + name + "/def/decomp", NC_CHAR,
                                 strlen(cbuf), &cbuf);
            CHECK_ERR

            /* dimension names */
            err = driver.put_att(fid, NC_GLOBAL, "/__pio__/var/" + name + "/def/dims", NC_STRING,
                                 ndims, dnames_array.data());
            CHECK_ERR

            /* Type of NetCDF API called */
            err = driver.put_att(fid, NC_GLOBAL, "/__pio__/var/" + name + "/def/ncop", NC_CHAR, 7,
                                 (void *)"darray");
            CHECK_ERR

            /* NetCDF external data type */
            ibuf = xtype;
            err = driver.put_att(fid, NC_GLOBAL, "/__pio__/var/" + name + "/def/nctype", NC_INT, 1,
                                 &ibuf);
            CHECK_ERR

            /* Number of dimensions */
            err = driver.put_att(fid, NC_GLOBAL, "/__pio__/var/" + name + "/def/ndims", NC_INT, 1,
                                 &ndims);
            CHECK_ERR
        }
    } else { /* this variable is not partitioned */
        // flatten into 1 dim only apply to non-scalar variables
        if (ndims){
            int esize;
            MPI_Offset vsize;

            // Convert into byte array
            err = e3sm_io_xlen_nc_type(xtype, &esize);
            CHECK_MPIERR

            vsize = 2 * sizeof(int64_t) * ndims  // Include start and count arrays
                  + var->vlen * esize;
            err = driver.def_local_var(fid, "/__pio__/var/" + name, NC_UBYTE, 1, &vsize,
                                       &var->vid);
            CHECK_ERR
        }
        else {
            err = driver.def_local_var(fid, "/__pio__/var/" + name, xtype, ndims, NULL,
                                       &var->vid);
            CHECK_ERR
        }

        // Attributes for non-constant small vars
        // Scorpio attributes are only written by rank 0
        if (cfg.rank == 0) {
            // ADIOS type enum, variables without decomposition map are stored as byte array
            ibuf = (int)e3sm_io_type_nc2adios (xtype);
            err  = driver.put_att(fid, NC_GLOBAL, "/__pio__/var/" + name + "/def/adiostype", NC_INT,
                                  1, &ibuf);
            CHECK_ERR

            // Scalar var does not have dims
            if (ndims > 0) {
                err = driver.put_att(fid, NC_GLOBAL, "/__pio__/var/" + name + "/def/dims", NC_STRING,
                                     ndims, dnames_array.data());
                CHECK_ERR
            }

            // Type of NetCDF API called
            err = driver.put_att(fid, NC_GLOBAL, "/__pio__/var/" + name + "/def/ncop", NC_CHAR, 7,
                                 (void *)"put_var");
            CHECK_ERR

            // NetCDF type enum
            ibuf = xtype;
            err = driver.put_att(fid, NC_GLOBAL, "/__pio__/var/" + name + "/def/nctype", NC_INT, 1,
                                 &ibuf);
            CHECK_ERR

            // Number of dimensions
            err = driver.put_att(fid, NC_GLOBAL, "/__pio__/var/" + name + "/def/ndims", NC_INT, 1,
                                 &ndims);
            CHECK_ERR
        }

        var->decom_id   = -1;
        var->frame_id   = -1;
        var->fill_id    = -1;
    }
err_out:
    return err;
#else
    return -1;
#endif
}

int scorpio_write_var(e3sm_io_driver &driver,
                      int frameid,
                      int fid,
                      var_meta &var,
                      MPI_Datatype itype,
                      void *buf,
                      e3sm_io_op_mode mode)
{
#ifdef ENABLE_ADIOS2
    int err;
    void *wbuf=buf;

    /* prepend start and count to write buffer for small, not-partitioned,
     * non-scalar variables This must be done by allocating another buffer,
     * add start/count, and copy over user write buffer.
     * These variables are stored as arrays of type byte.
     */
    if (var.decomp_id < 0 && var.ndims) {
        int i, esize;
        int64_t *ptr;

        MPI_Type_size(itype, &esize);
        wbuf = (void*) malloc(2 * var.ndims * sizeof(int64_t) + var.vlen * esize);
        ptr = (int64_t*)wbuf;
        *(ptr++) = (var.isRecVar) ? frameid : 0; /* prepend start of 1st dim */
        for (i=1; i<var.ndims; i++) *(ptr++) = 0; /* prepend start of remaining dims */
        for (i=0; i<var.ndims; i++) *(ptr++) = var.dims[i]; /* prepend count */
        memcpy(ptr, buf, var.vlen * esize);

        /* not partitioned variables are stored as byte type */
        itype = MPI_BYTE;
    }

    err = driver.put_varl(fid, var.vid, itype, wbuf, nbe);
    CHECK_ERR

    if (var.decomp_id >= 0) {
        int num_data_block_writers = 1;
        err = driver.put_varl(fid, var.num_data_block_writers, MPI_INT, &num_data_block_writers, nbe);
        CHECK_ERR
    }

    if (var.frame_id >= 0) {
        /* Scorpio allows variables to be written with different numbers of
         * time steps.
         */
        if (driver.cfg->rank == 0) {
            err = driver.put_varl(fid, var.frame_id, MPI_INT, &frameid, nbe);
            CHECK_ERR
        }

        /* Scorpio allows variables to switch decomposition at different time
         * steps.
         */
        if (driver.cfg->rank == 0) {
            if (var.fill_id >= 0) {
                err = driver.put_varl(fid, var.decom_id, MPI_INT, &var.piodecomid, nbe);
            }
            else {
                int tmp = -(var.piodecomid);
                err = driver.put_varl(fid, var.decom_id, MPI_INT, &tmp, nbe);
            }
            CHECK_ERR
        }
    }

    if (var.fill_id >= 0) {
        if (driver.cfg->rank == 0) {
            double tmp = 0;
            err = driver.put_varl(fid, var.fill_id, MPI_DOUBLE, &tmp, nbe);
            CHECK_ERR
        }
    }

err_out:
    if (wbuf != buf) free(wbuf);

    return err;
#else
    return -1;
#endif
}

