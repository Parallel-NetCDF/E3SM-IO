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
#include <cstring>
//
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
//
#include <adios2_c.h>
#include <mpi.h>
//
#include <e3sm_io.h>
#include <e3sm_io_err.h>

#include <e3sm_io_driver_adios2.hpp>

#define CHECK_AERR                                                    \
    {                                                                 \
        if (aerr != adios2_error_none) {                              \
            printf ("Error at line %d in %s:\n", __LINE__, __FILE__); \
            nerrs++;                                                  \
            DEBUG_ABORT                                               \
            goto err_out;                                             \
        }                                                             \
    }

#define CHECK_APTR(P)                                              \
    {                                                              \
        if (P == NULL) {                                           \
            printf ("Error in %s:%d: NULL\n", __FILE__, __LINE__); \
            nerrs++;                                               \
            DEBUG_ABORT                                            \
            goto err_out;                                          \
        }                                                          \
    }

size_t adios2_type_size (adios2_type type) {
    switch (type) {
        case adios2_type_int32_t:
            return 4;
        case adios2_type_float:
            return 4;
        case adios2_type_double:
            return 9;
        case adios2_type_uint8_t:
            return 1;
        default:
            printf ("Error at line %d in %s: Unknown type %d\n", __LINE__, __FILE__, type);
    }

    return 0;
}

e3sm_io_driver_adios2::e3sm_io_driver_adios2 (e3sm_io_config *cfg) : e3sm_io_driver (cfg) {
    twrite = tsel = text = 0;
}

e3sm_io_driver_adios2::~e3sm_io_driver_adios2 () {
    int nerrs = 0;
    int rank;
    double tsel_all, twrite_all, text_all;

    // printf("adios2 destructor\n");

    MPI_Allreduce (&twrite, &twrite_all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce (&tsel, &tsel_all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce (&text, &text_all, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        printf ("#%%$: H5Dwrite_time_max: %lf\n", twrite_all);
        printf ("#%%$: H5Sselect_hyperslab_time_max: %lf\n", tsel_all);
        printf ("#%%$: H5Dset_extent_time_max: %lf\n", text_all);
    }
}

int e3sm_io_driver_adios2::create (std::string path, MPI_Comm comm, MPI_Info info, int *fid) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp;
    char ng[32];

    fp       = new adios2_file ();
    fp->path = std::string (path);

    fp->adp = adios2_init (comm, "");
    CHECK_APTR (fp->adp)

    switch (cfg->filter) {
        case none:
            fp->op = NULL;
            break;
        case bzip2:
            fp->op = adios2_define_operator (fp->adp, "testOp", "bzip2");
            CHECK_APTR (fp->op)
            break;
        default:
            throw "Unsupported filter";
    }

    fp->iop = adios2_declare_io (fp->adp, "e3sm_wrap");
    CHECK_APTR (fp->iop)
    if (cfg->api == adios2_bp3) {
        aerr = adios2_set_engine (fp->iop, "BP3");
    } else {
        aerr = adios2_set_engine (fp->iop, "BP4");
    }
    CHECK_AERR

    sprintf (ng, "%d", cfg->num_group);
    aerr = adios2_set_parameter (fp->iop, "substreams", ng);
    CHECK_AERR
    aerr = adios2_set_parameter (fp->iop, "CollectiveMetadata", "OFF");
    CHECK_AERR
    fp->ep   = NULL;
    fp->read = 0;

    *fid = this->files.size ();
    this->files.push_back (fp);

err_out:;
    return nerrs;
}

int e3sm_io_driver_adios2::open (std::string path, MPI_Comm comm, MPI_Info info, int *fid) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp;
    adios2_step_status stat;

    fp       = new adios2_file ();
    fp->path = std::string (path);

    fp->adp = adios2_init (comm, "");
    CHECK_APTR (fp->adp)

    fp->iop = adios2_declare_io (fp->adp, "e3sm_wrap");
    CHECK_APTR (fp->iop)
    if (cfg->api == adios2_bp3) {
        aerr = adios2_set_engine (fp->iop, "BP3");
    } else {
        aerr = adios2_set_engine (fp->iop, "BP4");
    }
    aerr = adios2_set_parameter (fp->iop, "substreams", "1");
    CHECK_AERR
    aerr = adios2_set_parameter (fp->iop, "CollectiveMetadata", "OFF");
    CHECK_AERR
    fp->ep   = NULL;
    fp->read = 1;

    // Open into data mode
    fp->ep = adios2_open (fp->iop, fp->path.c_str (), adios2_mode_read);
    CHECK_APTR (fp->ep)
    aerr = adios2_begin_step (fp->ep, adios2_step_mode_read, -1, &stat);
    CHECK_AERR

    *fid = this->files.size ();
    this->files.push_back (fp);

err_out:;
    return nerrs;
}

int e3sm_io_driver_adios2::close (int fid) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    adios2_bool result;

    if (fp->ep) {
        aerr = adios2_end_step (fp->ep);
        CHECK_AERR

        aerr = adios2_close (fp->ep);
        CHECK_AERR
        fp->ep = NULL;
    }
    aerr = adios2_remove_io (&result, fp->adp, "e3sm_wrap");
    CHECK_AERR
    fp->iop = NULL;
    adios2_finalize (fp->adp);

    delete fp;

err_out:;
    return nerrs;
}

int e3sm_io_driver_adios2::inq_file_info (int fid, MPI_Info *info) {
    *info = MPI_INFO_NULL;
    return 0;
}

static MPI_Offset get_dir_size (std::string path) {
    DIR *dir;
    struct dirent *dit;
    struct stat st;
    long size             = 0;
    MPI_Offset total_size = 0;
    std::string subpath;

    dir = opendir (path.c_str ());
    if (dir) {
        while ((dit = readdir (dir)) != NULL) {
            if ((strcmp (dit->d_name, ".") == 0) || (strcmp (dit->d_name, "..") == 0)) continue;

            subpath = path + '/' + std::string (dit->d_name);
            if (lstat (subpath.c_str (), &st) != 0) continue;
            size = st.st_size;

            if (S_ISDIR (st.st_mode)) {
                total_size += get_dir_size (subpath.c_str ()) + size;
            } else {
                total_size += (MPI_Offset)size;
            }
        }
        closedir (dir);
    } else {  // Try again as file
        struct stat file_stat;
        stat (path.c_str (), &file_stat);
        total_size = (MPI_Offset) (file_stat.st_size);
    }

    return total_size;
}

int e3sm_io_driver_adios2::inq_put_size (int fid, MPI_Offset *size) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];

    *size = get_dir_size (fp->path);

err_out:;
    return nerrs;
}
int e3sm_io_driver_adios2::inq_get_size (int fid, MPI_Offset *size) {
    return this->inq_put_size (fid, size);
}
int e3sm_io_driver_adios2::inq_malloc_size (MPI_Offset *size) {
    *size = 0;
    return 0;
}
int e3sm_io_driver_adios2::inq_malloc_max_size (MPI_Offset *size) {
    *size = 0;
    return 0;
}
int e3sm_io_driver_adios2::inq_rec_size (int fid, MPI_Offset *size) {
    adios2_file *fp = this->files[fid];

    *size = (MPI_Offset) (fp->recsize);

    return 0;
}

int e3sm_io_driver_adios2::def_var (
    int fid, std::string name, MPI_Datatype type, int ndim, int *dimids, int *did) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    int i;
    adios2_variable *dp;
    size_t dims[E3SM_IO_DRIVER_MAX_RANK], zero[E3SM_IO_DRIVER_MAX_RANK];
    size_t opidx;

    for (i = 0; i < ndim; i++) {
        dims[i] = fp->dsizes[dimids[i]];
        if (dims[i] == NC_UNLIMITED) { dims[i] = 1; }
        zero[i] = 0;
    }

    dp = adios2_define_variable (fp->iop, name.c_str (), mpi_type_to_adios2_type (type),
                                 (size_t)ndim, dims, zero, zero, adios2_constant_dims_false);
    CHECK_APTR (dp)

    if (fp->op) {
        aerr = adios2_add_operation (&opidx, dp, fp->op, NULL, NULL);
        CHECK_AERR
    }

    *did = fp->dids.size ();
    fp->dids.push_back (dp);
    fp->ndims.push_back (ndim);

err_out:;
    return nerrs;
}
int e3sm_io_driver_adios2::def_local_var (
    int fid, std::string name, MPI_Datatype type, int ndim, MPI_Offset *dsize, int *did) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    int i;
    adios2_variable *dp;
    size_t dims[E3SM_IO_DRIVER_MAX_RANK];
    size_t opidx;

    for (i = 0; i < ndim; i++) {
        dims[i] = (size_t)dsize[i];
        // Unlimited dimension is simulated by timesteps
        if (dims[i] == NC_UNLIMITED) { dims[i] = 1; }
    }

    dp = adios2_define_variable (fp->iop, name.c_str (), mpi_type_to_adios2_type (type),
                                 (size_t)ndim, NULL, NULL, dims, adios2_constant_dims_false);
    CHECK_APTR (dp)

    if (fp->op) {
        aerr = adios2_add_operation (&opidx, dp, fp->op, NULL, NULL);
        CHECK_AERR
    }

    *did = fp->dids.size ();
    fp->dids.push_back (dp);
    fp->ndims.push_back (ndim);

err_out:;
    return nerrs;
}
int e3sm_io_driver_adios2::inq_var (int fid, std::string name, int *did) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    adios2_variable *dp;
    size_t ndim;

    dp = adios2_inquire_variable (fp->iop, name.c_str ());
    CHECK_APTR (dp)

    aerr = adios2_variable_ndims (&ndim, dp);
    CHECK_AERR

    *did = fp->dids.size ();
    fp->dids.push_back (dp);
    fp->ndims.push_back (ndim);

err_out:;
    return nerrs;
}

int e3sm_io_driver_adios2::inq_var_off (int fid, int vid, MPI_Offset *off) {
    throw "Function not supported";
    return 1;
}
int e3sm_io_driver_adios2::def_dim (int fid, std::string name, MPI_Offset size, int *dimid) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    adios2_variable *dp;

    dp =
        adios2_define_variable (fp->iop, ("/__e3sm_io__/dim/" + name).c_str (), adios2_type_int64_t,
                                0, NULL, NULL, NULL, adios2_constant_dims_true);
    CHECK_APTR (dp)

    *dimid = fp->dsizes.size ();
    fp->dsizes.push_back ((size_t)size);
    fp->ddids.push_back (dp);
err_out:;
    return nerrs;
}
int e3sm_io_driver_adios2::inq_dim (int fid, std::string name, int *dimid) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    MPI_Offset size;
    adios2_variable *dp;

    dp = adios2_inquire_variable (fp->iop, ("/__e3sm_io__/dim/" + name).c_str ());

    // Read must happen in data mode
    if (fp->ep == NULL) {
        nerrs += this->enddef (fid);
        CHECK_NERR
        aerr = adios2_get (fp->ep, dp, &size, adios2_mode_sync);
        CHECK_AERR
        nerrs += this->redef (fid);
        CHECK_NERR
    } else {
        aerr = adios2_get (fp->ep, dp, &size, adios2_mode_sync);
        CHECK_AERR
    }

    *dimid = fp->dsizes.size ();
    fp->dsizes.push_back ((size_t)size);
    fp->ddids.push_back (dp);
err_out:;
    return nerrs;
}

int e3sm_io_driver_adios2::inq_dimlen (int fid, int dimid, MPI_Offset *size) {
    adios2_file *fp = this->files[fid];

    *size = (MPI_Offset) (fp->dsizes[dimid]);

    return 0;
}

int e3sm_io_driver_adios2::enddef (int fid) {
    int nerrs = 0;
    int i;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    adios2_step_status stat;

    if (fp->read) {
        fp->ep = adios2_open (fp->iop, fp->path.c_str (), adios2_mode_read);
        CHECK_APTR (fp->ep)
        aerr = adios2_begin_step (fp->ep, adios2_step_mode_read, -1, &stat);
        CHECK_AERR
    } else {
        fp->ep = adios2_open (fp->iop, fp->path.c_str (), adios2_mode_write);
        CHECK_APTR (fp->ep)
        aerr = adios2_begin_step (fp->ep, adios2_step_mode_append, -1, &stat);
        CHECK_AERR

        // Write dim variables
        if (this->cfg->rank == 0) {
            for (i = 0; i < fp->ddids.size (); i++) {
                adios2_put (fp->ep, fp->ddids[i], &(fp->dsizes[i]), adios2_mode_sync);
            }
        }
    }

err_out:;
    return nerrs;
}
int e3sm_io_driver_adios2::redef (int fid) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];

    aerr = adios2_close (fp->ep);
    CHECK_AERR
    fp->ep = NULL;

err_out:;
    return nerrs;
}
int e3sm_io_driver_adios2::wait (int fid) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    adios2_step_status stat;

    aerr = adios2_end_step (fp->ep);
    CHECK_AERR

    if (fp->read) {
        aerr = adios2_perform_gets (fp->ep);
        CHECK_AERR
        aerr = adios2_begin_step (fp->ep, adios2_step_mode_read, -1, &stat);
    } else {
        aerr = adios2_perform_puts (fp->ep);
        CHECK_AERR
        aerr = adios2_begin_step (fp->ep, adios2_step_mode_update, -1, &stat);
    }
    CHECK_AERR

err_out:;
    return nerrs;
}

int e3sm_io_driver_adios2::put_att (
    int fid, int vid, std::string name, MPI_Datatype type, MPI_Offset size, void *buf) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    adios2_variable *did;
    adios2_attribute *aid;

    // adios2 can't create empty attr
    if (size == 0) goto err_out;

    if (vid == E3SM_IO_GLOBAL_ATTR) {
        aid = adios2_define_attribute_array (fp->iop, name.c_str (), mpi_type_to_adios2_type (type),
                                             buf, (size_t)size);
        CHECK_APTR (aid);
    } else {
        char vname[1024];
        size_t namesize = 1024;

        did  = fp->dids[vid];
        aerr = adios2_variable_name (vname, &namesize, did);
        CHECK_AERR
        vname[namesize] = '\0';

        aid = adios2_define_variable_attribute_array (
            fp->iop, name.c_str (), mpi_type_to_adios2_type (type), buf, (size_t)size, vname, "/");
        CHECK_APTR (aid)
    }

err_out:;
    return nerrs;
}

int e3sm_io_driver_adios2::get_att (int fid, int vid, std::string name, void *buf) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    adios2_variable *did;
    adios2_attribute *aid;
    size_t asize;

    if (vid == E3SM_IO_GLOBAL_ATTR) {
        aid = adios2_inquire_attribute (fp->iop, name.c_str ());
        CHECK_APTR (aid)
    } else {
        char vname[1024];
        size_t tmp = 1024;

        did  = fp->dids[vid];
        aerr = adios2_variable_name (vname, &tmp, did);
        CHECK_AERR

        aid = adios2_inquire_variable_attribute (fp->iop, name.c_str (), vname, "/");
        CHECK_APTR (aid)
    }

    aerr = adios2_attribute_size (&asize, aid);
    CHECK_AERR
    aerr = adios2_attribute_data (buf, &asize, aid);
    CHECK_AERR

err_out:;
    return nerrs;
}

int e3sm_io_driver_adios2::put_varl (
    int fid, int vid, MPI_Datatype type, void *buf, e3sm_io_op_mode mode) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    int i;
    size_t ndim;
    adios2_variable *did;
    size_t astart[E3SM_IO_DRIVER_MAX_RANK], ablock[E3SM_IO_DRIVER_MAX_RANK];
    adios2_mode iomode;
    double ts, te;

    did = fp->dids[vid];

    te = MPI_Wtime ();

    switch (mode) {
        case indep:
        case coll: {
            iomode = adios2_mode_sync;
            break;
        }
        case nb: {
            iomode = adios2_mode_deferred;
            break;
        }
        case nbe: {
            iomode = adios2_mode_sync;
            break;
        }
        default:
            throw "Unrecognized mode";
    }

    aerr = adios2_put (fp->ep, did, buf, iomode);
    CHECK_AERR
    this->twrite += MPI_Wtime () - te;

err_out:;
    return nerrs;
}

int e3sm_io_driver_adios2::put_vara (int fid,
                                     int vid,
                                     MPI_Datatype type,
                                     MPI_Offset *start,
                                     MPI_Offset *count,
                                     void *buf,
                                     e3sm_io_op_mode mode) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    int i;
    size_t ndim;
    adios2_variable *did;
    size_t astart[E3SM_IO_DRIVER_MAX_RANK], ablock[E3SM_IO_DRIVER_MAX_RANK];
    adios2_mode iomode;
    double ts, te;

    did  = fp->dids[vid];
    ndim = fp->ndims[vid];

    ts = MPI_Wtime ();
    if (start) {
        if (count) {
            for (i = 0; i < ndim; i++) {
                astart[i] = (size_t)start[i];
                ablock[i] = (size_t)count[i];
            }
        } else {
            for (i = 0; i < ndim; i++) {
                astart[i] = (size_t)start[i];
                ablock[i] = 1;
            }
        }
    } else {
        for (i = 0; i < ndim; i++) { astart[i] = 0; }
        aerr = adios2_variable_shape (ablock, did);
        CHECK_AERR
    }

    ts = MPI_Wtime ();
    if (ndim > 0) {
        aerr = adios2_set_selection (did, ndim, astart, ablock);
        CHECK_AERR
    }
    te = MPI_Wtime ();
    this->tsel += te - ts;

    switch (mode) {
        case indep:
        case coll: {
            iomode = adios2_mode_sync;
            break;
        }
        case nb: {
            iomode = adios2_mode_deferred;
            break;
        }
        case nbe: {
            iomode = adios2_mode_sync;
            break;
        }
        default:
            throw "Unrecognized mode";
    }

    aerr = adios2_put (fp->ep, did, buf, iomode);
    CHECK_AERR
    this->twrite += MPI_Wtime () - te;

err_out:;
    return nerrs;
}

int e3sm_io_driver_adios2::put_vars (int fid,
                                     int vid,
                                     MPI_Datatype type,
                                     MPI_Offset *start,
                                     MPI_Offset *count,
                                     MPI_Offset *stride,
                                     void *buf,
                                     e3sm_io_op_mode mode) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    int i;
    size_t esize;
    size_t ndim;
    adios2_type atype;
    adios2_variable *did;
    size_t astart[E3SM_IO_DRIVER_MAX_RANK], aend[E3SM_IO_DRIVER_MAX_RANK],
        ablock[E3SM_IO_DRIVER_MAX_RANK], dims[E3SM_IO_DRIVER_MAX_RANK];
    char *bufp;
    adios2_mode iomode;
    double ts, te;

    did  = fp->dids[vid];
    ndim = fp->ndims[vid];

    aerr = adios2_variable_type (&atype, did);
    CHECK_AERR
    esize = adios2_type_size (atype);

    // Extend rec dim
    ts = MPI_Wtime ();

    aerr = adios2_variable_shape (dims, did);
    CHECK_AERR

    text = MPI_Wtime () - ts;

    for (i = 0; i < ndim; i++) {
        astart[i] = start[i];
        aend[i]   = astart[i] + count[i] * stride[i];
        ablock[i] = 1;
    }
    if (aend[0] > dims[0]) {
        astart[0] = 0;
        aend[0]   = 1;
    }

    switch (mode) {
        case indep:
        case coll: {
            iomode = adios2_mode_sync;
            break;
        }
        case nb: {
            iomode = adios2_mode_deferred;
            break;
        }
        case nbe: {
            iomode = adios2_mode_sync;
            break;
        }
        default:
            throw "Unrecognized mode";
    }

    bufp = (char *)buf;
    while (astart[0] < aend[0]) {
        // Record dim replaced by time step
        if (astart[0] + ablock[0] > dims[0]) { astart[0] = 0; }

        ts = MPI_Wtime ();

        aerr = adios2_set_selection (did, ndim, astart, ablock);
        CHECK_AERR

        te = MPI_Wtime ();
        tsel += te - ts;

        aerr = adios2_put (fp->ep, did, bufp, iomode);
        CHECK_AERR

        twrite += MPI_Wtime () - te;

        bufp += esize;

        astart[ndim - 1] += stride[ndim - 1];
        for (i = ndim - 1; i > 0; i--) {
            if (astart[i] >= aend[i]) {
                astart[i] = (size_t) (start[i]);
                astart[i - 1] += stride[i - 1];
            }
        }
    }

err_out:;
    return nerrs;
}
int e3sm_io_driver_adios2::put_varn (int fid,
                                     int vid,
                                     MPI_Datatype type,
                                     int nreq,
                                     MPI_Offset **starts,
                                     MPI_Offset **counts,
                                     void *buf,
                                     e3sm_io_op_mode mode) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    int i, j;
    double ts, te;
    size_t esize, rsize, rsize_old = 0;
    int ndim;
    adios2_type atype;
    adios2_variable *did;
    char *bufp;
    size_t start[E3SM_IO_DRIVER_MAX_RANK], block[E3SM_IO_DRIVER_MAX_RANK];
    size_t dims[E3SM_IO_DRIVER_MAX_RANK];
    adios2_mode iomode;

    did  = fp->dids[vid];
    ndim = fp->ndims[vid];

    aerr = adios2_variable_type (&atype, did);
    CHECK_AERR
    esize = adios2_type_size (atype);

    // Extend rec dim if needed
    ts = MPI_Wtime ();

    aerr = adios2_variable_shape (dims, did);
    CHECK_AERR

    text += MPI_Wtime () - ts;

    switch (mode) {
        case indep:
        case coll: {
            iomode = adios2_mode_sync;
            break;
        }
        case nb: {
            iomode = adios2_mode_deferred;
            break;
        }
        case nbe: {
            iomode = adios2_mode_sync;
            break;
        }
        default:
            throw "Unrecognized mode";
    }

    // Call adios2_put_vara
    bufp = (char *)buf;
    for (i = 0; i < nreq; i++) {
        rsize = esize;
        for (j = 0; j < ndim; j++) { rsize *= counts[i][j]; }

        if (rsize) {
            for (j = 0; j < ndim; j++) {
                start[j] = (size_t)starts[i][j];
                block[j] = (size_t)counts[i][j];
            }
            // Record dim replaced by time step
            if (start[0] + block[0] > dims[0]) { start[0] = 0; }

            ts = MPI_Wtime ();

            aerr = adios2_set_selection (did, ndim, start, block);
            CHECK_AERR

            te = MPI_Wtime ();
            tsel += te - ts;

            aerr = adios2_put (fp->ep, did, bufp, iomode);
            CHECK_AERR

            twrite += MPI_Wtime () - te;
            bufp += rsize;
        }
    }

err_out:;
    return nerrs;
}

int e3sm_io_driver_adios2::put_vard (int fid,
                                     int vid,
                                     MPI_Datatype type,
                                     MPI_Offset nelem,
                                     MPI_Datatype ftype,
                                     void *buf,
                                     e3sm_io_op_mode mode) {
    throw "vard API not supported for HDF5 driver";
    return 1;
}

int e3sm_io_driver_adios2::get_vara (int fid,
                                     int vid,
                                     MPI_Datatype type,
                                     MPI_Offset *start,
                                     MPI_Offset *count,
                                     void *buf,
                                     e3sm_io_op_mode mode) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    int i;
    size_t ndim;
    adios2_variable *did;
    size_t astart[E3SM_IO_DRIVER_MAX_RANK], ablock[E3SM_IO_DRIVER_MAX_RANK];
    adios2_mode iomode;
    double ts, te;

    did  = fp->dids[vid];
    ndim = fp->ndims[vid];

    ts = MPI_Wtime ();
    if (start) {
        if (count) {
            for (i = 0; i < ndim; i++) {
                astart[i] = (size_t)start[i];
                ablock[i] = (size_t)count[i];
            }
        } else {
            for (i = 0; i < ndim; i++) {
                astart[i] = (size_t)start[i];
                ablock[i] = 1;
            }
        }
    } else {
        for (i = 0; i < ndim; i++) { astart[i] = 0; }
        aerr = adios2_variable_shape (ablock, did);
        CHECK_AERR
    }

    ts = MPI_Wtime ();
    if (ndim > 0) {
        aerr = adios2_set_selection (did, ndim, astart, ablock);
        CHECK_AERR
    }
    te = MPI_Wtime ();
    this->tsel += te - ts;

    switch (mode) {
        case indep:
        case coll: {
            iomode = adios2_mode_sync;
            break;
        }
        case nb: {
            iomode = adios2_mode_deferred;
            break;
        }
        default:
            throw "Unrecognized mode";
    }

    aerr = adios2_get (fp->ep, did, buf, iomode);
    CHECK_AERR
    this->twrite += MPI_Wtime () - te;

err_out:;
    return nerrs;
}
int e3sm_io_driver_adios2::get_vars (int fid,
                                     int vid,
                                     MPI_Datatype type,
                                     MPI_Offset *start,
                                     MPI_Offset *count,
                                     MPI_Offset *stride,
                                     void *buf,
                                     e3sm_io_op_mode mode) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    int i;
    size_t esize;
    size_t ndim;
    adios2_type atype;
    adios2_variable *did;
    size_t astart[E3SM_IO_DRIVER_MAX_RANK], aend[E3SM_IO_DRIVER_MAX_RANK],
        ablock[E3SM_IO_DRIVER_MAX_RANK], dims[E3SM_IO_DRIVER_MAX_RANK];
    char *bufp;
    adios2_mode iomode;
    double ts, te;

    did  = fp->dids[vid];
    ndim = fp->ndims[vid];

    aerr = adios2_variable_type (&atype, did);
    CHECK_AERR
    esize = adios2_type_size (atype);

    // Extend rec dim
    ts = MPI_Wtime ();

    aerr = adios2_variable_shape (dims, did);
    CHECK_AERR

    text = MPI_Wtime () - ts;

    for (i = 0; i < ndim; i++) {
        astart[i] = start[i];
        aend[i]   = astart[i] + count[i] * stride[i];
        ablock[i] = 1;
    }
    if (aend[0] > dims[0]) {
        astart[0] = 0;
        aend[0]   = 1;
    }

    switch (mode) {
        case indep:
        case coll: {
            iomode = adios2_mode_sync;
            break;
        }
        case nb: {
            iomode = adios2_mode_deferred;
            break;
        }
        default:
            throw "Unrecognized mode";
    }

    bufp = (char *)buf;
    while (astart[0] < aend[0]) {
        // Record dim replaced by time step
        if (astart[0] + ablock[0] > dims[0]) { astart[0] = 0; }

        ts = MPI_Wtime ();

        aerr = adios2_set_selection (did, ndim, astart, ablock);
        CHECK_AERR

        te = MPI_Wtime ();
        tsel += te - ts;

        aerr = adios2_get (fp->ep, did, bufp, iomode);
        CHECK_AERR

        twrite += MPI_Wtime () - te;

        bufp += esize;

        astart[ndim - 1] += stride[ndim - 1];
        for (i = ndim - 1; i > 0; i--) {
            if (astart[i] >= aend[i]) {
                astart[i] = (size_t) (start[i]);
                astart[i - 1] += stride[i - 1];
            }
        }
    }

err_out:;
    return nerrs;
}
int e3sm_io_driver_adios2::get_varn (int fid,
                                     int vid,
                                     MPI_Datatype type,
                                     int nreq,
                                     MPI_Offset **starts,
                                     MPI_Offset **counts,
                                     void *buf,
                                     e3sm_io_op_mode mode) {
    int nerrs = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    int i, j;
    double ts, te;
    size_t esize, rsize, rsize_old = 0;
    int ndim;
    adios2_type atype;
    adios2_variable *did;
    char *bufp;
    size_t start[E3SM_IO_DRIVER_MAX_RANK], block[E3SM_IO_DRIVER_MAX_RANK];
    size_t dims[E3SM_IO_DRIVER_MAX_RANK];
    adios2_mode iomode;

    did  = fp->dids[vid];
    ndim = fp->ndims[vid];

    aerr = adios2_variable_type (&atype, did);
    CHECK_AERR
    esize = adios2_type_size (atype);

    // Extend rec dim if needed
    ts = MPI_Wtime ();

    aerr = adios2_variable_shape (dims, did);
    CHECK_AERR

    text += MPI_Wtime () - ts;

    switch (mode) {
        case indep:
        case coll: {
            iomode = adios2_mode_sync;
            break;
        }
        case nb: {
            iomode = adios2_mode_deferred;
            break;
        }
        default:
            throw "Unrecognized mode";
    }

    // Call adios2_put_vara
    bufp = (char *)buf;
    for (i = 0; i < nreq; i++) {
        rsize = esize;
        for (j = 0; j < ndim; j++) { rsize *= counts[i][j]; }

        if (rsize) {
            for (j = 0; j < ndim; j++) {
                start[j] = (size_t)starts[i][j];
                block[j] = (size_t)counts[i][j];
            }
            // Record dim replaced by time step
            if (start[0] + block[0] > dims[0]) { start[0] = 0; }

            ts = MPI_Wtime ();

            aerr = adios2_set_selection (did, ndim, start, block);
            CHECK_AERR

            te = MPI_Wtime ();
            tsel += te - ts;

            aerr = adios2_get (fp->ep, did, bufp, iomode);
            CHECK_AERR

            twrite += MPI_Wtime () - te;
            bufp += rsize;
        }
    }

err_out:;
    return nerrs;
}

int e3sm_io_driver_adios2::get_vard (int fid,
                                     int vid,
                                     MPI_Datatype type,
                                     MPI_Offset nelem,
                                     MPI_Datatype ftype,
                                     void *buf,
                                     e3sm_io_op_mode mode) {
    throw "vard API not supported for HDF5 driver";
    return 1;
}
