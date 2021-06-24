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
#include <pnetcdf.h>
//
#include <e3sm_io.h>
#include <e3sm_io_err.h>

#include <e3sm_io_driver_pnc.hpp>

#define CHECK_NCERR                                                             \
    {                                                                           \
        if (err != NC_NOERR) {                                                  \
            printf ("Error at line %d in %s:\n", __LINE__, __FILE__);           \
            printf ("\t(%s) %s\n", ncmpi_strerrno (err), ncmpi_strerror (err)); \
            err = -1;                                                           \
            DEBUG_ABORT;                                                        \
            goto err_out;                                                       \
        }                                                                       \
    }

e3sm_io_driver_pnc::e3sm_io_driver_pnc (e3sm_io_config *cfg) : e3sm_io_driver (cfg) {
    if ((cfg->chunksize != 0) && (cfg->filter != none)) {
        throw "Fitler requries chunking in PnetCDF";
    }
}

e3sm_io_driver_pnc::~e3sm_io_driver_pnc () {}

int e3sm_io_driver_pnc::create (std::string path, MPI_Comm comm, MPI_Info info, int *fid) {
    int err, nerrs = 0;
    MPI_Offset put_buffer_size_limit;

    // Use the zip driver for chunked I/O
    if (cfg->chunksize > 0) {
        MPI_Info_set (info, "nc_compression", "enable");
        MPI_Info_set (info, "nc_zip_comm_unit", "chunk");
    }

    if (cfg->num_group > 1) {
        char ng[32];
        sprintf (ng, "%d", cfg->num_group);
        MPI_Info_set (info, "nc_num_subfiles", ng);
    }

    err = ncmpi_create (comm, path.c_str (), NC_CLOBBER | NC_64BIT_DATA, info, fid);
    CHECK_NCERR

    put_buffer_size_limit = 10485760;
    err                   = ncmpi_buffer_attach (*fid, put_buffer_size_limit);
    CHECK_ERR

err_out:;
    return nerrs;
}

int e3sm_io_driver_pnc::open (std::string path, MPI_Comm comm, MPI_Info info, int *fid) {
    int err, nerrs = 0;

    err = ncmpi_open (comm, path.c_str (), NC_64BIT_DATA, info, fid);
    CHECK_NCERR

err_out:;
    return nerrs;
}

int e3sm_io_driver_pnc::close (int fid) {
    int err, nerrs = 0;

    err = ncmpi_buffer_detach (fid);
    CHECK_ERR

    err = ncmpi_close (fid);
    CHECK_NCERR

err_out:;
    return nerrs;
}

int e3sm_io_driver_pnc::inq_file_info (int fid, MPI_Info *info) {
    int err, nerrs = 0;

    err = ncmpi_inq_file_info (fid, info);
    CHECK_NCERR

err_out:;
    return nerrs;
}

int e3sm_io_driver_pnc::inq_put_size (int fid, MPI_Offset *size) {
    int err, nerrs = 0;

    err = ncmpi_inq_put_size (fid, size);
    CHECK_NCERR

err_out:;
    return nerrs;
}

int e3sm_io_driver_pnc::inq_get_size (int fid, MPI_Offset *size) {
    int err, nerrs = 0;
    err = ncmpi_inq_get_size (fid, size);
    CHECK_NCERR

err_out:;
    return nerrs;
}

int e3sm_io_driver_pnc::inq_malloc_size (MPI_Offset *size) {
    int err, nerrs = 0;

    err = ncmpi_inq_malloc_size (size);
    CHECK_NCERR

err_out:;
    return nerrs;
}
int e3sm_io_driver_pnc::inq_malloc_max_size (MPI_Offset *size) {
    int err, nerrs = 0;

    err = ncmpi_inq_malloc_max_size (size);
    CHECK_NCERR

err_out:;
    return nerrs;
}

int e3sm_io_driver_pnc::inq_rec_size (int fid, MPI_Offset *size) {
    int err, nerrs = 0;

    err = ncmpi_inq_recsize (fid, size);
    CHECK_NCERR

err_out:;
    return nerrs;
}

int e3sm_io_driver_pnc::def_var (
    int fid, std::string name, MPI_Datatype type, int ndim, int *dimids, int *did) {
    int err, nerrs = 0;
    int i;
    MPI_Offset bufcounts;
    int cdim[E3SM_IO_DRIVER_MAX_RANK];
    int tsize;
    size_t csize = 0;

    err = ncmpi_def_var (fid, name.c_str (), mpitype2nctype (type), ndim, dimids, did);
    CHECK_NCERR

    if (ndim) {
        if ((cfg->chunksize > 0) || (this->dim_lens[dimids[0]] == NC_UNLIMITED)) {
            csize = MPI_Type_size (type, &tsize);
            for (i = 0; i < ndim; i++) {
                if (csize < cfg->chunksize) {
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
            err = ncmpi_put_att_int (fid, *did, "_chunkdim", NC_INT, ndim, cdim);
            CHECK_ERR

            switch (cfg->filter) {
                case none:
                    break;
                case deflate:
                    tsize = 2;  // TODO: Use formal PnetCDF filter ID
                    err   = ncmpi_put_att_int (fid, *did, "_zipdriver", NC_INT, 1, &tsize);
                    CHECK_ERR
                    break;
                default:
                    RET_ERR ("Unknown filter")
            }
        }
    }

    if (this->var_nelems.size () <= *did) {
        this->var_nelems.resize (*did + 1);
        this->var_ndims.resize (*did + 1);
    }
    bufcounts = 1;
    for (i = 0; i < ndim; i++) { bufcounts *= this->dim_lens[dimids[i]]; }
    this->var_nelems[*did] = bufcounts;
    this->var_ndims[*did]  = ndim;

err_out:;
    return nerrs;
}
int e3sm_io_driver_pnc::def_local_var (
    int fid, std::string name, MPI_Datatype type, int ndim, MPI_Offset *dsize, int *did) {
    int nerrs = 0;

    RET_ERR ("PNC does not support local variables")

err_out:;
    return nerrs;
}

int e3sm_io_driver_pnc::inq_var (int fid, std::string name, int *did) {
    int err, nerrs = 0;

    err = ncmpi_inq_varid (fid, name.c_str (), did);
    CHECK_NCERR

err_out:;
    return nerrs;
}
int e3sm_io_driver_pnc::inq_var_off (int fid, int vid, MPI_Offset *off) {
    int err, nerrs = 0;

    err = ncmpi_inq_varoffset (fid, vid, off);
    CHECK_NCERR

err_out:;
    return nerrs;
}
int e3sm_io_driver_pnc::def_dim (int fid, std::string name, MPI_Offset size, int *dimid) {
    int err, nerrs = 0;

    err = ncmpi_def_dim (fid, name.c_str (), size, dimid);
    CHECK_NCERR

    if (this->dim_lens.size () <= *dimid) { this->dim_lens.resize (*dimid + 1); }
    this->dim_lens[*dimid] = size;

err_out:;
    return nerrs;
}
int e3sm_io_driver_pnc::inq_dim (int fid, std::string name, int *dimid) {
    int err, nerrs = 0;
    MPI_Offset size;

    err = ncmpi_inq_dimid (fid, name.c_str (), dimid);
    CHECK_NCERR

    err = ncmpi_inq_dim (fid, *dimid, NULL, &size);
    CHECK_NCERR

    if (this->dim_lens.size () <= *dimid) { this->dim_lens.resize (*dimid + 1); }
    this->dim_lens[*dimid] = size;

err_out:;
    return nerrs;
}

int e3sm_io_driver_pnc::inq_dimlen (int fid, int dimid, MPI_Offset *size) {
    *size = this->dim_lens[dimid];
    return 0;
}
int e3sm_io_driver_pnc::enddef (int fid) {
    int err, nerrs = 0;

    err = ncmpi_enddef (fid);
    CHECK_NCERR

err_out:;
    return nerrs;
}
int e3sm_io_driver_pnc::redef (int fid) {
    int err, nerrs = 0;

    err = ncmpi_redef (fid);
    CHECK_NCERR

err_out:;
    return nerrs;
}
int e3sm_io_driver_pnc::wait (int fid) {
    int err, nerrs = 0;

    err = ncmpi_wait_all (fid, NC_REQ_ALL, NULL, NULL);
    CHECK_NCERR

err_out:;
    return nerrs;
}

int e3sm_io_driver_pnc::put_att (
    int fid, int vid, std::string name, MPI_Datatype type, MPI_Offset size, void *buf) {
    int err, nerrs = 0;

    if (vid == E3SM_IO_GLOBAL_ATTR) { vid = NC_GLOBAL; }

    err = ncmpi_put_att (fid, vid, name.c_str (), mpitype2nctype (type), size, buf);
    CHECK_NCERR

err_out:;
    return nerrs;
}

int e3sm_io_driver_pnc::get_att (int fid, int vid, std::string name, void *buf) {
    int err, nerrs = 0;

    if (vid == E3SM_IO_GLOBAL_ATTR) { vid = NC_GLOBAL; }

    err = ncmpi_get_att (fid, vid, name.c_str (), buf);
    CHECK_NCERR

err_out:;
    return nerrs;
}

int e3sm_io_driver_pnc::put_varl (
    int fid, int vid, MPI_Datatype type, void *buf, e3sm_io_op_mode mode) {
    int nerrs = 0;

    RET_ERR ("PNC does not support local variables")

err_out:;
    return nerrs;
}

int e3sm_io_driver_pnc::put_vara (int fid,
                                  int vid,
                                  MPI_Datatype type,
                                  MPI_Offset *start,
                                  MPI_Offset *count,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err, nerrs = 0;
    int i;
    MPI_Offset bufcount;

    if (start) {
        if (count) {
            bufcount = 1;
            for (i = 0; i < this->var_ndims[vid]; i++) { bufcount *= count[i]; }
            switch (mode) {
                case indep: {
                    err = ncmpi_put_vara (fid, vid, start, count, buf, bufcount, type);
                    break;
                }
                case coll: {
                    err = ncmpi_put_vara_all (fid, vid, start, count, buf, bufcount, type);
                    break;
                }
                case nb: {
                    err = ncmpi_iput_vara (fid, vid, start, count, buf, bufcount, type, NULL);
                    break;
                }
                case nbe: {
                    err = ncmpi_bput_vara (fid, vid, start, count, buf, bufcount, type, NULL);
                    break;
                }
                default:
                    throw "Unrecognized mode";
            }
        } else {
            bufcount = 1;
            switch (mode) {
                case indep: {
                    err = ncmpi_put_var1 (fid, vid, start, buf, bufcount, type);
                    break;
                }
                case coll: {
                    err = ncmpi_put_var1_all (fid, vid, start, buf, bufcount, type);
                    break;
                }
                case nb: {
                    err = ncmpi_iput_var1 (fid, vid, start, buf, bufcount, type, NULL);
                    break;
                }
                case nbe: {
                    err = ncmpi_bput_var1 (fid, vid, start, buf, bufcount, type, NULL);
                    break;
                }
                default:
                    throw "Unrecognized mode";
            }
        }
    } else {
        bufcount = this->var_nelems[vid];
        switch (mode) {
            case indep: {
                err = ncmpi_put_var (fid, vid, buf, bufcount, type);
                break;
            }
            case coll: {
                err = ncmpi_put_var_all (fid, vid, buf, bufcount, type);
                break;
            }
            case nb: {
                err = ncmpi_iput_var (fid, vid, buf, bufcount, type, NULL);
                break;
            }
            case nbe: {
                err = ncmpi_bput_var (fid, vid, buf, bufcount, type, NULL);
                break;
            }
            default:
                throw "Unrecognized mode";
        }
    }
    CHECK_NCERR

err_out:;
    return nerrs;
}
int e3sm_io_driver_pnc::put_vars (int fid,
                                  int vid,
                                  MPI_Datatype type,
                                  MPI_Offset *start,
                                  MPI_Offset *count,
                                  MPI_Offset *stride,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err, nerrs = 0;
    int i;
    MPI_Offset bufcount;

    bufcount = 1;
    for (i = 0; i < this->var_ndims[vid]; i++) { bufcount *= count[i]; }
    switch (mode) {
        case indep: {
            err = ncmpi_put_vars (fid, vid, start, count, stride, buf, bufcount, type);
            break;
        }
        case coll: {
            err = ncmpi_put_vars_all (fid, vid, start, count, stride, buf, bufcount, type);
            break;
        }
        case nb: {
            err = ncmpi_iput_vars (fid, vid, start, count, stride, buf, bufcount, type, NULL);
            break;
        }
        case nbe: {
            err = ncmpi_bput_vars (fid, vid, start, count, stride, buf, bufcount, type, NULL);
            break;
        }
        default:
            throw "Unrecognized mode";
    }
    CHECK_NCERR

err_out:;
    return nerrs;
}
int e3sm_io_driver_pnc::put_varn (int fid,
                                  int vid,
                                  MPI_Datatype type,
                                  int nreq,
                                  MPI_Offset **starts,
                                  MPI_Offset **counts,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err = NC_NOERR, nerrs = 0;
    int i, j;
    MPI_Offset bufcount;
    MPI_Offset blockcount;

    bufcount = 0;
    for (i = 0; i < nreq; i++) {
        blockcount = 1;
        for (j = 0; j < this->var_ndims[vid]; j++) { blockcount *= counts[i][j]; }
        bufcount += blockcount;
    }
    switch (mode) {
        case indep: {
            err = ncmpi_put_varn (fid, vid, nreq, starts, counts, buf, bufcount, type);
            break;
        }
        case coll: {
            err = ncmpi_put_varn_all (fid, vid, nreq, starts, counts, buf, bufcount, type);
            break;
        }
        case nb: {
            err = ncmpi_iput_varn (fid, vid, nreq, starts, counts, buf, bufcount, type, NULL);
            break;
        }
        case nbe: {
            err = ncmpi_bput_varn (fid, vid, nreq, starts, counts, buf, bufcount, type, NULL);
            break;
        }
        default:
            throw "Unrecognized mode";
    }
    CHECK_NCERR

err_out:;
    return nerrs;
}

int e3sm_io_driver_pnc::put_vard (int fid,
                                  int vid,
                                  MPI_Datatype type,
                                  MPI_Offset bufcount,
                                  MPI_Datatype ftype,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err, nerrs = 0;

    switch (mode) {
        case indep: {
            err = ncmpi_put_vard (fid, vid, ftype, buf, bufcount, type);
            break;
        }
        case coll: {
            err = ncmpi_put_vard_all (fid, vid, ftype, buf, bufcount, type);
            break;
        }
        default:
            throw "Unrecognized mode";
    }
    CHECK_NCERR

err_out:;
    return nerrs;
}

int e3sm_io_driver_pnc::get_vara (int fid,
                                  int vid,
                                  MPI_Datatype type,
                                  MPI_Offset *start,
                                  MPI_Offset *count,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err, nerrs = 0;
    int i;
    MPI_Offset bufcount;

    if (start) {
        if (count) {
            bufcount = 1;
            for (i = 0; i < this->var_ndims[vid]; i++) { bufcount *= count[i]; }
            switch (mode) {
                case indep: {
                    err = ncmpi_get_vara (fid, vid, start, count, buf, bufcount, type);
                    break;
                }
                case coll: {
                    err = ncmpi_get_vara_all (fid, vid, start, count, buf, bufcount, type);
                    break;
                }
                case nb: {
                    err = ncmpi_iget_vara (fid, vid, start, count, buf, bufcount, type, NULL);
                    break;
                }
                default:
                    throw "Unrecognized mode";
            }
        } else {
            bufcount = 1;
            switch (mode) {
                case indep: {
                    err = ncmpi_get_var1 (fid, vid, start, buf, bufcount, type);
                    break;
                }
                case coll: {
                    err = ncmpi_get_var1_all (fid, vid, start, buf, bufcount, type);
                    break;
                }
                case nb: {
                    err = ncmpi_iget_var1 (fid, vid, start, buf, bufcount, type, NULL);
                    break;
                }
                default:
                    throw "Unrecognized mode";
            }
        }
    } else {
        bufcount = this->var_nelems[vid];
        switch (mode) {
            case indep: {
                err = ncmpi_get_var (fid, vid, buf, bufcount, type);
                break;
            }
            case coll: {
                err = ncmpi_get_var_all (fid, vid, buf, bufcount, type);
                break;
            }
            case nb: {
                err = ncmpi_iget_var (fid, vid, buf, bufcount, type, NULL);
                break;
            }
            default:
                throw "Unrecognized mode";
        }
    }
    CHECK_NCERR

err_out:;
    return nerrs;
}
int e3sm_io_driver_pnc::get_vars (int fid,
                                  int vid,
                                  MPI_Datatype type,
                                  MPI_Offset *start,
                                  MPI_Offset *count,
                                  MPI_Offset *stride,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err, nerrs = 0;
    int i;
    MPI_Offset bufcount;

    bufcount = 1;
    for (i = 0; i < this->var_ndims[vid]; i++) { bufcount *= count[i]; }
    switch (mode) {
        case indep: {
            err = ncmpi_get_vars (fid, vid, start, count, stride, buf, bufcount, type);
            break;
        }
        case coll: {
            err = ncmpi_get_vars_all (fid, vid, start, count, stride, buf, bufcount, type);
            break;
        }
        case nb: {
            err = ncmpi_iget_vars (fid, vid, start, count, stride, buf, bufcount, type, NULL);
            break;
        }
        default:
            throw "Unrecognized mode";
    }
    CHECK_NCERR

err_out:;
    return nerrs;
}
int e3sm_io_driver_pnc::get_varn (int fid,
                                  int vid,
                                  MPI_Datatype type,
                                  int nreq,
                                  MPI_Offset **starts,
                                  MPI_Offset **counts,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err, nerrs = 0;
    int i, j;
    MPI_Offset bufcount;
    MPI_Offset blockcount;

    bufcount = 0;
    for (i = 0; i < nreq; i++) {
        blockcount = 1;
        for (j = 0; j < this->var_ndims[vid]; j++) { blockcount *= counts[i][j]; }
        bufcount += blockcount;
    }
    switch (mode) {
        case indep: {
            err = ncmpi_get_varn (fid, vid, nreq, starts, counts, buf, bufcount, type);
            break;
        }
        case coll: {
            err = ncmpi_get_varn_all (fid, vid, nreq, starts, counts, buf, bufcount, type);
            break;
        }
        case nb: {
            err = ncmpi_iget_varn (fid, vid, nreq, starts, counts, buf, bufcount, type, NULL);
            break;
        }
        default:
            throw "Unrecognized mode";
    }
    CHECK_NCERR

err_out:;
    return nerrs;
}

int e3sm_io_driver_pnc::get_vard (int fid,
                                  int vid,
                                  MPI_Datatype type,
                                  MPI_Offset bufcount,
                                  MPI_Datatype ftype,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err, nerrs = 0;

    switch (mode) {
        case indep: {
            err = ncmpi_get_vard (fid, vid, ftype, buf, bufcount, type);
            break;
        }
        case coll: {
            err = ncmpi_get_vard_all (fid, vid, ftype, buf, bufcount, type);
            break;
        }
        default:
            throw "Unrecognized mode";
    }
    CHECK_NCERR

err_out:;
    return nerrs;
}
