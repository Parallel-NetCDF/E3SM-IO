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
    var->name = strdup(name.c_str());

    for (i = 0; i < ndims; i++)
        dnames_array[i] = dnames[dimids[i]].c_str();

    // If there is a decomposition map associated with the variable,
    if (decomp_id >= 0) {
        MPI_Offset one = 1;

        err = driver.def_local_var(fid, name, xtype, 1,
                                   &(decom.raw_nreqs[decomp_id]), &var->vid);
        CHECK_ERR

        assert(ndims > 0);

        /* add an associated scalar variable frame_id (i.e. timesteps) */
        err = driver.def_local_var(fid, "frame_id/" + name, NC_INT, 1,
                                   &one, &var->frame_id);
        CHECK_ERR

        /* add an associated scalar variable decom_id (decomposition map ID) */
        err = driver.def_local_var(fid, "decomp_id/" + name, NC_INT, 1,
                                   &one, &var->decom_id);
        CHECK_ERR

        /* Add attributes */

        /* In Scorpio, variable's attributes are only written by rank 0 */
        if (cfg.rank == 0) {
            /* Decomposition map ID number */
            sprintf(cbuf, "%d", var->piodecomid);
            err = driver.put_att(fid, var->vid, "__pio__/decomp", NC_CHAR,
                                 strlen(cbuf), &cbuf);
            CHECK_ERR

            /* dimension names */
            err = driver.put_att(fid, var->vid, "__pio__/dims", NC_STRING,
                                 ndims, dnames_array.data());
            CHECK_ERR

            /* Type of NetCDF API called */
            err = driver.put_att(fid, var->vid, "__pio__/ncop", NC_CHAR, 7,
                                 (void *)"darray");
            CHECK_ERR

            /* NetCDF external data type */
            ibuf = xtype;
            err = driver.put_att(fid, var->vid, "__pio__/nctype", NC_INT, 1,
                                 &ibuf);
            CHECK_ERR

            /* Number of dimensions */
            err = driver.put_att(fid, var->vid, "__pio__/ndims", NC_INT, 1,
                                 &ndims);
            CHECK_ERR
        }
    } else { /* this variable is not partitioned */
        std::vector<MPI_Offset> dsize (ndims);
        MPI_Offset vsize = 1;
        int esize;

        for (i = 0; i < ndims; i++) {
            err = driver.inq_dimlen (fid, dimids[i], &dsize[i]);
            CHECK_ERR

            // Time dim is always 0, but block size should be 1
            if (dsize[i] == 0) dsize[i] = 1;
            var->dims[i] = dsize[i];
            vsize *= dsize[i];
        }

        // flatten into 1 dim only apply to non-scalar variables
        if (ndims){
            // Convert into byte array
            err = e3sm_io_xlen_nc_type(xtype, &esize);
            CHECK_MPIERR
            vsize *= esize;
            vsize += 8 * 2 * ndims; // Include start and count array
            err = driver.def_local_var(fid, name, NC_UBYTE, 1, &vsize,
                                       &var->vid);
            CHECK_ERR
        }
        else {
            err = driver.def_local_var(fid, name, xtype, ndims, NULL,
                                       &var->vid);
            CHECK_ERR
        }

        // Attributes for non-constant small vars
        // Scorpio attributes are only written by rank 0
        if (cfg.rank == 0) {
            // ADIOS type enum, variables without decomposition map are stored as byte array
            ibuf = (int)e3sm_io_type_nc2adios (xtype);
            err  = driver.put_att(fid, var->vid, "__pio__/adiostype", NC_INT,
                                  1, &ibuf);
            CHECK_ERR

            // Scalar var does not have dims
            if (ndims > 0) {
                err = driver.put_att(fid, var->vid, "__pio__/dims", NC_STRING,
                                     ndims, dnames_array.data());
                CHECK_ERR
            }

            // Type of NetCDF API called
            err = driver.put_att(fid, var->vid, "__pio__/ncop", NC_CHAR, 7,
                                 (void *)"put_var");
            CHECK_ERR

            // NetCDF type enum
            ibuf = xtype;
            err = driver.put_att(fid, var->vid, "__pio__/nctype", NC_INT, 1,
                                 &ibuf);
            CHECK_ERR

            // Number of dimensions
            err = driver.put_att(fid, var->vid, "__pio__/ndims", NC_INT, 1,
                                 &ndims);
            CHECK_ERR
        }

        var->decom_id   = -1;
        var->frame_id   = -1;
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
    int err, decomid;
    void *wbuf;

    if (var.isRecVar) frameid = -1;

    /* prepend start and count to write buffer for small, not-partitioned,
     * non-scalar variables This must be done by allocating another buffer,
     * add start/count, and copy over user write buffer.
     * These variables are stored as arrays of type byte.
     */
    if (var.decomp_id < 0 && var.ndims) {
        int esize;
        size_t cp_len;

#ifdef DEBUG
        size_t i, vlen=1;
        for (i=0; i<var.ndims; i++) vlen *= var.dims[i];
        assert(vlen == var.vlen);
#endif

        MPI_Type_size(itype, &esize);
        cp_len = var.ndims * sizeof(int64_t);
        wbuf = (void*) malloc(2 * cp_len + var.vlen * esize);
        memset(wbuf, 0, cp_len);
        memcpy((char*)wbuf+cp_len,   var.dims, cp_len);
        memcpy((char*)wbuf+cp_len*2, buf,      var.vlen * esize);
    }
    else
        wbuf = buf;

    if (var.frame_id < 0 && var.ndims)
        itype = MPI_BYTE;

    err = driver.put_varl (fid, var.vid, itype, wbuf, nbe);
    CHECK_ERR

    if (var.frame_id >= 0) {

        err = driver.put_varl (fid, var.frame_id, MPI_INT, &frameid, nbe);
        CHECK_ERR

        decomid = var.piodecomid;
        if (var.decomp_id < 0) decomid *= -1;

        err = driver.put_varl (fid, var.decom_id, MPI_INT, &decomid, nbe);
        CHECK_ERR
    }

err_out:
    if (wbuf != buf) free(wbuf);

    return err;
#else
    return -1;
#endif
}

int scorpio_put_fill_att(e3sm_io_driver &driver,
                         int             fid,
                         void           *buf,
                         var_meta       *varp)
{
#ifdef ENABLE_ADIOS2
    int err, varid;
    std::string name("fillval_id/");
    MPI_Offset one = 1;
    MPI_Datatype xType;

    /* define a local variable to store fill value */
    name.append(varp->name);
    err = driver.def_local_var(fid, name, varp->xType, 1, &one, &varid);
    CHECK_ERR

         if (varp->xType == NC_FLOAT)  xType = MPI_FLOAT;
    else if (varp->xType == NC_INT)    xType = MPI_INT;
    else if (varp->xType == NC_DOUBLE) xType = MPI_DOUBLE;
    else assert(0);

    err = driver.put_varl(fid, varid, xType, buf, nbe);
    CHECK_ERR

    /* also add attribute _FillValue */
    err = driver.put_att(fid, varp->vid, _FillValue, varp->xType, 1, buf);
    CHECK_ERR

err_out:
    return err;
#else
    return -1;
#endif
}

