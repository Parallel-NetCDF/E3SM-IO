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
//
#include <sys/stat.h>
#include <string.h> /* strdup() */
//
#include <e3sm_io.h>
#include <e3sm_io_err.h>

#include <e3sm_io_driver_pnc.hpp>

#define CHECK_NCERR {                                                       \
    if (err != NC_NOERR) {                                                  \
        printf ("Error in %s line %d function %s:\n", __FILE__, __LINE__,   \
                __func__);                                                  \
        printf ("\t(%s) %s\n", ncmpi_strerrno (err), ncmpi_strerror (err)); \
        DEBUG_ABORT;                                                        \
        goto err_out;                                                       \
    }                                                                       \
}

inline nc_type mpitype2nctype (MPI_Datatype type) {
    if (type == MPI_DOUBLE){
        return NC_DOUBLE;
    }
    else if (type == MPI_FLOAT){
        return NC_FLOAT;
    }
    else if (type == MPI_INT){
        return NC_INT;
    }
    else if (type == MPI_LONG_LONG){
        return NC_INT64;
    }
    else if (type == MPI_CHAR){
        return NC_CHAR;
    }
    else if (type == MPI_BYTE){
        return NC_BYTE;
    }

    throw "Unsupported datatype";
}

e3sm_io_driver_pnc::e3sm_io_driver_pnc (e3sm_io_config *cfg) : e3sm_io_driver (cfg) {
    if ((cfg->chunksize != 0) && (cfg->filter != none)) {
        throw "Fitler requries chunking in PnetCDF";
    }
}

e3sm_io_driver_pnc::~e3sm_io_driver_pnc () {}

int e3sm_io_driver_pnc::create (std::string path, MPI_Comm comm, MPI_Info info, int *fid) {
    int err;
    MPI_Offset put_buffer_size_limit, size;

    // Use the zip driver for chunked I/O
    if (cfg->chunksize > 0) {
        MPI_Info_set (info, "nc_compression", "enable");
        MPI_Info_set (info, "nc_zip_comm_unit", "chunk");
    }

    err = ncmpi_create (comm, path.c_str (), NC_CLOBBER | NC_64BIT_DATA, info, fid);
    CHECK_NCERR

    put_buffer_size_limit = 10485760;
    err = ncmpi_buffer_attach (*fid, put_buffer_size_limit);
    CHECK_NCERR

    err = ncmpi_inq_get_size(*fid, &size); CHECK_NCERR
    this->amount_RD += size;
    err = ncmpi_inq_put_size(*fid, &size); CHECK_NCERR
    this->amount_WR += size;

err_out:
    return err;
}

int e3sm_io_driver_pnc::open (std::string path, MPI_Comm comm, MPI_Info info, int *fid) {
    int err;
    MPI_Offset size;
    MPI_Offset put_buffer_size_limit;

    err = ncmpi_open (comm, path.c_str (), NC_64BIT_DATA, info, fid);
    CHECK_NCERR

    put_buffer_size_limit = 10485760;
    err = ncmpi_buffer_attach (*fid, put_buffer_size_limit);
    CHECK_NCERR

    err = ncmpi_inq_get_size(*fid, &size); CHECK_NCERR
    this->amount_RD += size;
    err = ncmpi_inq_put_size(*fid, &size); CHECK_NCERR
    this->amount_WR += size;

err_out:
    return err;
}

int e3sm_io_driver_pnc::close (int fid) {
    int err;

    err = ncmpi_buffer_detach (fid);
    CHECK_NCERR

    err = ncmpi_close (fid);
    CHECK_NCERR

err_out:
    return err;
}

int e3sm_io_driver_pnc::inq_file_info (int fid, MPI_Info *info) {
    int err;

    err = ncmpi_inq_file_info (fid, info);
    CHECK_NCERR

err_out:
    return err;
}

int e3sm_io_driver_pnc::inq_file_size (std::string path, MPI_Offset *size) {
    int err;
    struct stat file_stat;

    err = stat (path.c_str (), &file_stat);
    CHECK_ERR

    *size = (MPI_Offset) (file_stat.st_size);

err_out:
    return err;
}

int e3sm_io_driver_pnc::inq_put_size (MPI_Offset *size) {
    *size = this->amount_WR;
    return 0;
}

int e3sm_io_driver_pnc::inq_get_size (MPI_Offset *size) {
    *size = this->amount_RD;
    return 0;
}

int e3sm_io_driver_pnc::inq_malloc_size (MPI_Offset *size) {
    int err;

    err = ncmpi_inq_malloc_size (size);
    CHECK_NCERR

err_out:
    return err;
}
int e3sm_io_driver_pnc::inq_malloc_max_size (MPI_Offset *size) {
    int err;

    err = ncmpi_inq_malloc_max_size (size);
    CHECK_NCERR

err_out:
    return err;
}

int e3sm_io_driver_pnc::inq_rec_size (int fid, MPI_Offset *size) {
    int err;

    err = ncmpi_inq_recsize (fid, size);
    CHECK_NCERR

err_out:
    return err;
}

int e3sm_io_driver_pnc::expand_rec_size (int fid, MPI_Offset size) {
    return 0;
}

int e3sm_io_driver_pnc::def_var (
    int fid, std::string name, nc_type xtype, int ndim, int *dimids, int *varid) {
    int err;

    err = ncmpi_def_var (fid, name.c_str (), xtype, ndim, dimids, varid);
    CHECK_NCERR

    /* skip the followings if compression is not enabled */
    if (cfg->chunksize == 0) return err;

    if (ndim) {
        int i, tsize, cdim[E3SM_IO_DRIVER_MAX_RANK];
        int csize = 0;

        if ((cfg->chunksize > 0) || (this->dim_lens[dimids[0]] == NC_UNLIMITED)) {
            e3sm_io_xlen_nc_type(xtype, &csize);
            for (i = 0; i < ndim; i++) {
                if ((size_t)csize < cfg->chunksize) {
                    cdim[i] = this->dim_lens[dimids[i]];
                    csize *= cdim[i];
                } else {
                    cdim[i] = 1;
                }
            }
            // Chunk size along rec dim is always 1
            if (this->dim_lens[dimids[0]] == NC_UNLIMITED) { cdim[0] = 1; }
        }

        if (csize > 0) {
            err = ncmpi_put_att_int (fid, *varid, "_chunkdim", NC_INT, ndim, cdim);
            CHECK_NCERR

            switch (cfg->filter) {
                case none:
                    break;
                case deflate:
                    tsize = 2;  // TODO: Use formal PnetCDF filter ID
                    err   = ncmpi_put_att_int (fid, *varid, "_zipdriver", NC_INT, 1, &tsize);
                    CHECK_NCERR
                    break;
                default:
                    ERR_OUT ("Unknown filter")
            }
        }
    }

    /* variable IDs returned by ncmpi_def_var(), if successful, are always
     * non-negative.
     */
    if (this->var_nelems.size () <= (size_t)*varid) {
        this->var_nelems.resize (*varid + 1);
        this->var_ndims.resize (*varid + 1);
    }
    this->var_ndims[*varid]  = ndim;

err_out:
    return err;
}

int e3sm_io_driver_pnc::def_local_var (
    int fid, std::string name, nc_type xtype, int ndim, MPI_Offset *dsize, int *varid) {
    int err = 0;

    ERR_OUT ("PNC does not support local variables")

err_out:
    return err;
}

int e3sm_io_driver_pnc::inq_varid (int fid, const char *name, int *varid) {

    // inq_var is used to check whether a variable exist so error is expected
    return ncmpi_inq_varid (fid, name, varid);
}

int e3sm_io_driver_pnc::inq_var (int fid, int varid, char *name, nc_type *xtypep,
                                 int *ndimsp, int *dimids, int *nattsp)
{
    int err;
    err = ncmpi_inq_var(fid, varid, name, xtypep, ndimsp, dimids, nattsp);
    return err;
}

int e3sm_io_driver_pnc::inq_var_name(int   ncid,
                                     int   varid,
                                     char *name)
{
    int err = ncmpi_inq_varname(ncid, varid, name);
    CHECK_NCERR

err_out:
    return err;
}

int e3sm_io_driver_pnc::inq_var_off (int fid, int vid, MPI_Offset *off) {
    int err;

    err = ncmpi_inq_varoffset (fid, vid, off);
    CHECK_NCERR

err_out:
    return err;
}

int e3sm_io_driver_pnc::def_dim (int fid, std::string name, MPI_Offset size, int *dimid) {
    int err;

    err = ncmpi_def_dim (fid, name.c_str (), size, dimid);
    CHECK_NCERR

    if (cfg->chunksize == 0) return err;

    /* dimension IDs returned by ncmpi_def_dim(), if successful, are always
     * non-negative.
     */
    if (this->dim_lens.size () <= (size_t)*dimid)
        this->dim_lens.resize (*dimid + 1);
    this->dim_lens[*dimid] = size;

err_out:
    return err;
}

int e3sm_io_driver_pnc::inq_dim (int fid, std::string name, int *dimid) {
    int err;
    MPI_Offset size;

    err = ncmpi_inq_dimid (fid, name.c_str (), dimid);
    if (err == NC_EBADDIM){
        err = -1;
        goto err_out;
    }
    CHECK_NCERR

    if (cfg->chunksize == 0) return err;

    err = ncmpi_inq_dim (fid, *dimid, NULL, &size);
    CHECK_NCERR

    /* dimension IDs in PnetCDF are always non-negative */
    if (this->dim_lens.size () <= (size_t)*dimid)
        this->dim_lens.resize (*dimid + 1);
    this->dim_lens[*dimid] = size;

err_out:
    return err;
}

int e3sm_io_driver_pnc::inq_dimlen (int fid, int dimid, MPI_Offset *size) {
    int err=0;

    if (cfg->chunksize == 0) {
        err = ncmpi_inq_dimlen(fid, dimid, size);
        CHECK_NCERR
    }
    else
        *size = this->dim_lens[dimid];

err_out:
    return err;
}

int e3sm_io_driver_pnc::enddef (int fid) {
    int err;
    MPI_Offset size, prev_WR, prev_RD;

    err = ncmpi_inq_get_size(fid, &prev_RD); CHECK_NCERR
    err = ncmpi_inq_put_size(fid, &prev_WR); CHECK_NCERR

    err = ncmpi_enddef (fid); CHECK_NCERR

    err = ncmpi_inq_get_size(fid, &size); CHECK_NCERR
    this->amount_RD += size - prev_RD;
    err = ncmpi_inq_put_size(fid, &size); CHECK_NCERR
    this->amount_WR += size - prev_WR;

err_out:
    return err;
}

int e3sm_io_driver_pnc::redef (int fid) {
    int err;
    MPI_Offset size, prev_WR, prev_RD;

    err = ncmpi_inq_get_size(fid, &prev_RD); CHECK_NCERR
    err = ncmpi_inq_put_size(fid, &prev_WR); CHECK_NCERR

    err = ncmpi_redef (fid); CHECK_NCERR

    err = ncmpi_inq_get_size(fid, &size); CHECK_NCERR
    this->amount_RD += size - prev_RD;
    err = ncmpi_inq_put_size(fid, &size); CHECK_NCERR
    this->amount_WR += size - prev_WR;

err_out:
    return err;
}

int e3sm_io_driver_pnc::wait (int fid) {
    int err;
    MPI_Offset size, prev_WR, prev_RD;

    err = ncmpi_inq_get_size(fid, &prev_RD); CHECK_NCERR
    err = ncmpi_inq_put_size(fid, &prev_WR); CHECK_NCERR

    err = ncmpi_wait_all(fid, NC_REQ_ALL, NULL, NULL); CHECK_NCERR

    err = ncmpi_inq_get_size(fid, &size); CHECK_NCERR
    this->amount_RD += size - prev_RD;
    err = ncmpi_inq_put_size(fid, &size); CHECK_NCERR
    this->amount_WR += size - prev_WR;

err_out:
    return err;
}

int e3sm_io_driver_pnc::put_att (
    int fid, int vid, std::string name, nc_type xtype, MPI_Offset size, const void *buf) {
    int err;

    err = ncmpi_put_att (fid, vid, name.c_str (), xtype, size, buf);
    CHECK_NCERR

err_out:
    return err;
}

int e3sm_io_driver_pnc::get_att (int fid, int vid, std::string name, void *buf) {
    int err;

    err = ncmpi_get_att (fid, vid, name.c_str (), buf);
    CHECK_NCERR

err_out:
    return err;
}

int e3sm_io_driver_pnc::inq_att (int fid, int vid, std::string name, MPI_Offset *size){
    int err;

    err = ncmpi_inq_attlen (fid, vid, name.c_str (), size);
    CHECK_NCERR

err_out:
    return err;
}

int e3sm_io_driver_pnc::put_varl (
    int fid, int vid, MPI_Datatype itype, void *buf, e3sm_io_op_mode mode) {
    int err = 0;

    ERR_OUT ("PNC does not support local variables")

err_out:
    return err;
}

int e3sm_io_driver_pnc::put_vara (int fid,
                                  int vid,
                                  MPI_Datatype itype,
                                  MPI_Offset *start,
                                  MPI_Offset *count,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err=NC_NOERR;
    MPI_Offset size, prev_WR;

    if (mode == coll || mode == indep) {
        err = ncmpi_inq_put_size(fid, &prev_WR); CHECK_NCERR
    }

    if (start) {
        if (count) {
            switch (mode) {
                case nb:
                    err = ncmpi_iput_vara (fid, vid, start, count, buf, -1, itype, NULL);
                    break;
                case nbe:
                    err = ncmpi_bput_vara (fid, vid, start, count, buf, -1, itype, NULL);
                    break;
                case coll:
                    err = ncmpi_put_vara_all (fid, vid, start, count, buf, -1, itype);
                    break;
                case indep:
                    err = ncmpi_put_vara (fid, vid, start, count, buf, -1, itype);
                    break;
                default:
                    throw "Unrecognized mode";
            }
        } else {
            switch (mode) {
                case nb:
                    err = ncmpi_iput_var1 (fid, vid, start, buf, -1, itype, NULL);
                    break;
                case nbe:
                    err = ncmpi_bput_var1 (fid, vid, start, buf, -1, itype, NULL);
                    break;
                case coll:
                    err = ncmpi_put_var1_all (fid, vid, start, buf, -1, itype);
                    break;
                case indep:
                    err = ncmpi_put_var1 (fid, vid, start, buf, -1, itype);
                    break;
                default:
                    throw "Unrecognized mode";
            }
        }
    } else {
        switch (mode) {
            case nb:
                err = ncmpi_iput_var (fid, vid, buf, -1, itype, NULL);
                break;
            case nbe:
                err = ncmpi_bput_var (fid, vid, buf, -1, itype, NULL);
                break;
            case coll:
                err = ncmpi_put_var_all (fid, vid, buf, -1, itype);
                break;
            case indep:
                err = ncmpi_put_var (fid, vid, buf, -1, itype);
                break;
            default:
                throw "Unrecognized mode";
        }
    }
    CHECK_NCERR

    if (mode == coll || mode == indep) {
        err = ncmpi_inq_put_size (fid, &size); CHECK_NCERR
        this->amount_WR += size - prev_WR;
    }

err_out:
    return err;
}

int e3sm_io_driver_pnc::put_varn (int fid,
                                  int vid,
                                  MPI_Datatype itype,
                                  int nreq,
                                  MPI_Offset **starts,
                                  MPI_Offset **counts,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err = NC_NOERR;
    MPI_Offset size, prev_WR;

    if (mode == coll || mode == indep) {
        err = ncmpi_inq_put_size(fid, &prev_WR); CHECK_NCERR
    }

    switch (mode) {
        case nb:
            err = ncmpi_iput_varn (fid, vid, nreq, starts, counts, buf, -1, itype, NULL);
            break;
        case nbe:
            err = ncmpi_bput_varn (fid, vid, nreq, starts, counts, buf, -1, itype, NULL);
            break;
        case coll:
            err = ncmpi_put_varn_all (fid, vid, nreq, starts, counts, buf, -1, itype);
            break;
        case indep:
            err = ncmpi_put_varn (fid, vid, nreq, starts, counts, buf, -1, itype);
            break;
        default:
            throw "Unrecognized mode";
    }
    CHECK_NCERR

    if (mode == coll || mode == indep) {
        err = ncmpi_inq_put_size (fid, &size); CHECK_NCERR
        this->amount_WR += size - prev_WR;
    }

err_out:
    return err;
}

int e3sm_io_driver_pnc::get_vara (int fid,
                                  int vid,
                                  MPI_Datatype itype,
                                  MPI_Offset *start,
                                  MPI_Offset *count,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err=NC_NOERR;
    MPI_Offset size, prev_RD;

    if (mode == coll || mode == indep) {
        err = ncmpi_inq_get_size(fid, &prev_RD); CHECK_NCERR
    }

    if (start) {
        if (count) {
            switch (mode) {
                case nb:
                    err = ncmpi_iget_vara (fid, vid, start, count, buf, -1, itype, NULL);
                    break;
                case coll:
                    err = ncmpi_get_vara_all (fid, vid, start, count, buf, -1, itype);
                    break;
                case indep:
                    err = ncmpi_get_vara (fid, vid, start, count, buf, -1, itype);
                    break;
                default:
                    throw "Unrecognized mode";
            }
        } else {
            switch (mode) {
                case nb:
                    err = ncmpi_iget_var1 (fid, vid, start, buf, -1, itype, NULL);
                    break;
                case coll:
                    err = ncmpi_get_var1_all (fid, vid, start, buf, -1, itype);
                    break;
                case indep:
                    err = ncmpi_get_var1 (fid, vid, start, buf, -1, itype);
                    break;
                default:
                    throw "Unrecognized mode";
            }
        }
    } else {
        switch (mode) {
            case nb:
                err = ncmpi_iget_var (fid, vid, buf, -1, itype, NULL);
                break;
            case coll:
                err = ncmpi_get_var_all (fid, vid, buf, -1, itype);
                break;
            case indep:
                err = ncmpi_get_var (fid, vid, buf, -1, itype);
                break;
            default:
                throw "Unrecognized mode";
        }
    }
    CHECK_NCERR

    if (mode == coll || mode == indep) {
        err = ncmpi_inq_get_size (fid, &size); CHECK_NCERR
        this->amount_RD += size - prev_RD;
    }

err_out:
    return err;
}

int e3sm_io_driver_pnc::get_varn (int fid,
                                  int vid,
                                  MPI_Datatype itype,
                                  int nreq,
                                  MPI_Offset **starts,
                                  MPI_Offset **counts,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err=NC_NOERR;
    MPI_Offset size, prev_RD;

    if (mode == coll || mode == indep) {
        err = ncmpi_inq_get_size(fid, &prev_RD); CHECK_NCERR
    }

    switch (mode) {
        case nb:
            err = ncmpi_iget_varn (fid, vid, nreq, starts, counts, buf, -1, itype, NULL);
            break;
        case coll:
            err = ncmpi_get_varn_all (fid, vid, nreq, starts, counts, buf, -1, itype);
            break;
        case indep:
            err = ncmpi_get_varn (fid, vid, nreq, starts, counts, buf, -1, itype);
            break;
        default:
            throw "Unrecognized mode";
    }
    CHECK_NCERR

    if (mode == coll || mode == indep) {
        err = ncmpi_inq_get_size (fid, &size); CHECK_NCERR
        this->amount_RD += size - prev_RD;
    }

err_out:
    return err;
}

