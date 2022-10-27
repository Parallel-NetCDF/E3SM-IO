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
#include <cstdlib>
#include <cstring>
#include <iostream>
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
#include <e3sm_io_profile.hpp>

#define CHECK_AERR                                                    \
    {                                                                 \
        if (aerr != adios2_error_none) {                              \
            printf ("Error at line %d in %s:\n", __LINE__, __FILE__); \
            err = -1;                                                 \
            DEBUG_ABORT                                               \
            goto err_out;                                             \
        }                                                             \
    }

#define CHECK_APTR(P)                                              \
    {                                                              \
        if (P == NULL) {                                           \
            printf ("Error in %s:%d: NULL\n", __FILE__, __LINE__); \
            err = -1;                                              \
            DEBUG_ABORT                                            \
            goto err_out;                                          \
        }                                                          \
    }

size_t adios2_type_size (adios2_type type) {
    switch (type) {
        case adios2_type_int32_t:
            return 4;
        case adios2_type_int64_t:
            return 4;
        case adios2_type_float:
            return 4;
        case adios2_type_double:
            return 8;
        case adios2_type_uint8_t:
            return 1;
        case adios2_type_int8_t:
            return 1;
        case adios2_type_string:
            return 1;
        default:
            printf ("Error at line %d in %s: Unknown type %d\n", __LINE__, __FILE__, type);
    }

    return 0;
}

inline int adios2_type_convert (
    adios2_type in_type, size_t len, void *inbuf, adios2_type out_type, void *outbuf) {
    int err = 0;

    if (in_type == adios2_type_double) {
        if (out_type == adios2_type_float) {
            double *inptr = (double *)inbuf;
            float *outptr = (float *)outbuf;

            for (; inptr < ((double *)inbuf) + len; inptr++, outptr++) {
                *outptr = (float)(*inptr);
            }
        } else {
            ERR_OUT ("Output type not supproted")
        }
    } else if (in_type == adios2_type_float) {
        if (out_type == adios2_type_double) {
            float *inptr   = (float *)inbuf;
            double *outptr = (double *)outbuf;

            for (; inptr < ((float *)inbuf) + len; inptr++, outptr++) {
                *outptr = (double)(*inptr);
            }
        } else {
            ERR_OUT ("Output type not supproted")
        }
    } else {
        ERR_OUT ("Input type not supproted")
    }

err_out:;
    return err;
}

bool e3sm_io_driver_adios2::compatible (std::string path) {
    int err = 0;
    adios2_error aerr;
    adios2_adios *adp = NULL;
    adios2_io *iop    = NULL;
    adios2_engine *ep = NULL;
    adios2_bool result;
    bool ret = false;

    adp = adios2_init (MPI_COMM_SELF, "");
    CHECK_APTR (adp)

    iop = adios2_declare_io (adp, "e3sm_check");
    CHECK_APTR (iop)

    aerr = adios2_set_engine (iop, "BP3");
    CHECK_AERR

    std::cerr.setstate(std::ios_base::failbit);
    ep = adios2_open (iop, path.c_str(), adios2_mode_read);
    std::cerr.clear();
    if (ep) { 
        ret = true; 
        adios2_close (ep);
    }

    adios2_remove_io (&result, adp, "e3sm_check");
    adios2_finalize (adp);

err_out:;
    return (err == -1) ? false : ret;
}


e3sm_io_driver_adios2::e3sm_io_driver_adios2 (e3sm_io_config *cfg) : e3sm_io_driver (cfg) {}

e3sm_io_driver_adios2::~e3sm_io_driver_adios2 () {}

int e3sm_io_driver_adios2::create (std::string path, MPI_Comm comm, MPI_Info info, int *fid) {
    int err = 0;
    adios2_error aerr;
    adios2_file *fp;
    char ng[32];
    adios2_step_status stat;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    fp       = new adios2_file ();
    fp->path = std::string (path); /* BP5 files directory */
    err      = MPI_Comm_rank (comm, &(fp->rank));
    CHECK_MPIERR

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
    aerr = adios2_set_engine (fp->iop, "BP5");
    CHECK_AERR

    sprintf (ng, "%d", cfg->num_subfiles);
    aerr = adios2_set_parameter (fp->iop, "substreams", ng);
    CHECK_AERR
    aerr = adios2_set_parameter (fp->iop, "CollectiveMetadata", "OFF");
    CHECK_AERR
    fp->ep   = NULL;
    fp->read = 0;

    fp->ep = adios2_open (fp->iop, fp->path.c_str (), adios2_mode_write);
    CHECK_APTR (fp->ep)
    aerr = adios2_begin_step (fp->ep, adios2_step_mode_append, -1, &stat);
    CHECK_AERR

    fp->wr = true;

    *fid = this->files.size ();
    this->files.push_back (fp);

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    return err;
}

int e3sm_io_driver_adios2::open (std::string path, MPI_Comm comm, MPI_Info info, int *fid) {
    int err = 0;
    adios2_error aerr;
    adios2_file *fp;
    adios2_step_status stat;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    fp       = new adios2_file ();
    fp->path = std::string (path);
    err      = MPI_Comm_rank (comm, &(fp->rank));
    CHECK_MPIERR

    fp->adp = adios2_init (comm, "");
    CHECK_APTR (fp->adp)

    fp->iop = adios2_declare_io (fp->adp, "e3sm_wrap");
    CHECK_APTR (fp->iop)

    aerr = adios2_set_engine (fp->iop, "BP3");
    CHECK_AERR
    aerr = adios2_set_parameter (fp->iop, "substreams", "1");
    CHECK_AERR
    aerr = adios2_set_parameter (fp->iop, "CollectiveMetadata", "OFF");
    CHECK_AERR
    fp->ep   = NULL;
    fp->read = 1;

    // Open into data mode

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_OPEN)
    fp->ep = adios2_open (fp->iop, fp->path.c_str (), adios2_mode_read);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_OPEN)
    CHECK_APTR (fp->ep)
    aerr = adios2_begin_step (fp->ep, adios2_step_mode_read, -1, &stat);
    CHECK_AERR

    fp->wr = false;

    *fid = this->files.size ();
    this->files.push_back (fp);

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    return err;
}

int e3sm_io_driver_adios2::close (int fid) {
    int err = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    adios2_bool result;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    if (fp->ep) {
        aerr = adios2_end_step (fp->ep);
        CHECK_AERR

        E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_CLOSE)
        aerr = adios2_close (fp->ep);
        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_CLOSE)
        CHECK_AERR
        fp->ep = NULL;
    }

    aerr = adios2_remove_io (&result, fp->adp, "e3sm_wrap");
    CHECK_AERR
    fp->iop = NULL;
    adios2_finalize (fp->adp);

    delete fp;

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    return err;
}

int e3sm_io_driver_adios2::inq_file_info (int fid, MPI_Info *info) {
    *info = MPI_INFO_NULL;
    return 0;
}

int e3sm_io_driver_adios2::inq_put_size (MPI_Offset *size) {
    int err = 0;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    *size = this->amount_WR;

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    return err;
}
int e3sm_io_driver_adios2::inq_get_size (MPI_Offset *size) {
    int err = 0;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    *size = this->amount_RD;

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    return err;
}

inline MPI_Offset get_dir_size (std::string path) {
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

int e3sm_io_driver_adios2::inq_file_size (std::string path, MPI_Offset *size) {
    int err = 0;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_GET_FSIZE)
    *size = get_dir_size (path);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_GET_FSIZE)

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    return err;
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

int e3sm_io_driver_adios2::expand_rec_size (int fid, MPI_Offset size) {
    return 0;
}

int e3sm_io_driver_adios2::def_var (
    int fid, std::string name, nc_type xtype, int ndim, int *dimids, int *did) {
    int err = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    int i;
    adios2_variable *dp;
    size_t dims[E3SM_IO_DRIVER_MAX_RANK], zero[E3SM_IO_DRIVER_MAX_RANK];
    size_t opidx;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    for (i = 0; i < ndim; i++) {
        dims[i] = fp->dsizes[dimids[i]];
        if (dims[i] == NC_UNLIMITED) { dims[i] = 1; }
        zero[i] = 0;
    }
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_DEF_VAR)
    dp = adios2_define_variable (fp->iop, name.c_str (), e3sm_io_type_nc2adios (xtype),
                                 (size_t)ndim, dims, zero, zero, adios2_constant_dims_false);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_DEF_VAR)
    CHECK_APTR (dp)

    if (fp->op) {
        aerr = adios2_add_operation (&opidx, dp, fp->op, NULL, NULL);
        CHECK_AERR
    }

    *did = fp->dids.size ();
    fp->dids.push_back (dp);
    fp->ndims.push_back (ndim);

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    return err;
}
int e3sm_io_driver_adios2::def_local_var (
    int fid, std::string name, nc_type xtype, int ndim, MPI_Offset *dsize, int *did) {
    int err = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    int i;
    adios2_variable *dp;
    size_t dims[E3SM_IO_DRIVER_MAX_RANK];
    size_t *ldim;
    size_t opidx;
    adios2_type dtype = e3sm_io_type_nc2adios (xtype);

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    for (i = 0; i < ndim; i++) {
        // Unlimited dimension is simulated by timesteps
        if (dsize[i] == -1)
            dims[i] = 1; 
        else
            dims[i] = (size_t)dsize[i];
    }

    if (dtype == adios2_type_string) {
        ldim = NULL;
    } else {
        ldim = dims;
    }
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_DEF_VAR)
    dp = adios2_define_variable (fp->iop, name.c_str (), dtype, (size_t)ndim, NULL, NULL, ldim,
                                 adios2_constant_dims_false);
    CHECK_APTR (dp)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_DEF_VAR)

    if (fp->op) {
        aerr = adios2_add_operation (&opidx, dp, fp->op, NULL, NULL);
        CHECK_AERR
    }

    *did = fp->dids.size ();
    fp->dids.push_back (dp);
    fp->ndims.push_back (ndim);

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    return err;
}

int e3sm_io_driver_adios2::inq_varid(int fid, const char *name, int *did)
{
    int err = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    adios2_variable *dp;
    size_t ndim;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_INQ_VAR)
    dp = adios2_inquire_variable (fp->iop, name);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_INQ_VAR)
    // inq_var is used to check whether a variable exist so error is expected
    if (!dp) {
        err = -1;
        goto err_out;
    }

    aerr = adios2_variable_ndims (&ndim, dp);
    CHECK_AERR

    *did = fp->dids.size ();
    fp->dids.push_back (dp);
    fp->ndims.push_back (ndim);

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    return err;
}

int e3sm_io_driver_adios2::inq_var (int fid, int varid, char *name, nc_type *xtypep,
                                    int *ndimsp, int *dimids, int *nattsp)
{
    int err=0;
    return err;
}

int e3sm_io_driver_adios2::inq_var_name (int fid, int vid, char *name) {
    name[0] = '\0';
    printf ("inq_var_name is not yet implementaed\n");
    return -1;
}

int e3sm_io_driver_adios2::inq_var_off (int fid, int vid, MPI_Offset *off) {
    throw "Function not supported";
    return -1;
}
int e3sm_io_driver_adios2::def_dim (int fid, std::string name, MPI_Offset size, int *dimid) {
    int err = 0;
    adios2_file *fp = this->files[fid];
    adios2_variable *dp;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_DEF_VAR)
    dp = adios2_define_variable (fp->iop, ("/__pio__/dim/" + name).c_str (), adios2_type_uint64_t,
                                 0, NULL, NULL, NULL, adios2_constant_dims_true);

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_DEF_VAR)
    CHECK_APTR (dp)

    *dimid = fp->dsizes.size ();
    fp->dsizes.push_back ((size_t)size);
    fp->ddids.push_back (dp);
    fp->dimmap[name] = *dimid;
err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    return err;
}
int e3sm_io_driver_adios2::inq_dim (int fid, std::string name, int *dimid) {
    int err = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    MPI_Offset size;
    adios2_variable *dp;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    if (fp->dimmap.find(name) != fp->dimmap.end()) {
        *dimid = fp->dimmap[name];
    }
    else {
        E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_INQ_VAR)
        dp = adios2_inquire_variable (fp->iop, ("/__pio__/dim/" + name).c_str ());
        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_INQ_VAR)

        // Read must happen in data mode
        if (fp->ep == NULL) {
            E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
            err = this->enddef (fid);
            E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)
            CHECK_ERR
            E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_GET_VAR)
            aerr = adios2_get (fp->ep, dp, &size, adios2_mode_sync);
            E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_GET_VAR)
            CHECK_AERR
            E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
            err = this->redef (fid);
            E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)
            CHECK_ERR
        } else {
            E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_GET_VAR)
            aerr = adios2_get (fp->ep, dp, &size, adios2_mode_sync);
            E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_GET_VAR)
            CHECK_AERR
        }

        *dimid = fp->dsizes.size ();
        fp->dsizes.push_back ((size_t)size);
        fp->ddids.push_back (dp);
    }
err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    return err;
}

int e3sm_io_driver_adios2::inq_dimlen (int fid, int dimid, MPI_Offset *size) {
    adios2_file *fp = this->files[fid];

    *size = (MPI_Offset) (fp->dsizes[dimid]);

    return 0;
}

int e3sm_io_driver_adios2::enddef (int fid) {
    int err = 0;
    size_t i;
    adios2_file *fp = this->files[fid];

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_ENDDEF)

    if (!fp->read) {
        // Write dim variables
        E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_PUT_VAR)
        if (this->cfg->rank == 0) {
            for (i = 0; i < fp->ddids.size (); i++) {
                adios2_put (fp->ep, fp->ddids[i], &(fp->dsizes[i]), adios2_mode_sync);
            }
            this->amount_WR += fp->ddids.size() * 8;
        }
        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_PUT_VAR)
    }

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_ENDDEF)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    return err;
}
int e3sm_io_driver_adios2::redef (int fid) {
#if 0
    int err = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_CLOSE)
    // aerr = adios2_close (fp->ep);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_CLOSE)
    // CHECK_AERR
    // fp->ep = NULL;

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    return err;
#else
    return 0;
#endif
}
int e3sm_io_driver_adios2::wait (int fid) {
    int err = 0;
//    adios2_error aerr;
//   adios2_file *fp = this->files[fid];
//   adios2_step_status stat;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)
    /*
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
    */

// err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    return err;
}

int e3sm_io_driver_adios2::put_att (
    int fid, int vid, std::string name, nc_type xtype, MPI_Offset size, const void *buf) {
    int err = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    adios2_variable *did;
    adios2_attribute *aid;
    adios2_type atype = e3sm_io_type_nc2adios (xtype);
    size_t esize;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    if (vid == NC_GLOBAL) {
        // ADIOS2 have no char type, we translate char array into unit sized string
        // MPI has not string type, we use WCHAR to represent string
        E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_PUT_ATT)
        if (xtype == NC_CHAR) {
            aid = adios2_define_attribute (fp->iop, name.c_str (), adios2_type_string, buf);
        } else { /* NC_STRING */
            aid = adios2_define_attribute_array (fp->iop, name.c_str (), atype, buf, (size_t)size);
        }
        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_PUT_ATT)
        CHECK_APTR (aid);
    } else {
        char vname[1024];
        size_t namesize = 1024;

        did  = fp->dids[vid];
        aerr = adios2_variable_name (vname, &namesize, did);
        CHECK_AERR
        vname[namesize] = '\0';
        // ADIOS2 have no char type, we translate char array into unit sized string
        // MPI has not string type, we use WCHAR to represent string
        E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_PUT_ATT)
        if (xtype == NC_CHAR) {
            aid = adios2_define_variable_attribute (fp->iop, name.c_str (), adios2_type_string, buf,
                                                    vname, "/");
        } else { /* NC_STRING */
            aid = adios2_define_variable_attribute_array (fp->iop, name.c_str (), atype, buf,
                                                          (size_t)size, vname, "/");
        }
        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_PUT_ATT)
        CHECK_APTR (aid)
    }

    // ADIOS2 write attributes in each subfile
    // Instead of finding the exact process writing the attribute 
    // we count them in the first num_subfiles processes
    if (fp->rank < cfg->num_subfiles) {
        esize = adios2_type_size (atype);
        this->amount_WR += size * esize;
    }

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    return err;
}

int e3sm_io_driver_adios2::get_att (int fid, int vid, std::string name, void *buf) {
    int err = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    adios2_variable *did;
    adios2_attribute *aid;
    adios2_type atype;
    size_t asize, esize;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    if (vid == NC_GLOBAL) {
        E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_GET_ATT)
        aid = adios2_inquire_attribute (fp->iop, name.c_str ());
        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_GET_ATT)
        CHECK_APTR (aid)
    } else {
        char vname[1024];
        size_t tmp = 1024;

        did  = fp->dids[vid];
        aerr = adios2_variable_name (vname, &tmp, did);
        CHECK_AERR

        E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_GET_ATT)
        aid = adios2_inquire_variable_attribute (fp->iop, name.c_str (), vname, "/");
        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_GET_ATT)
        CHECK_APTR (aid)
    }

    aerr = adios2_attribute_size (&asize, aid);
    CHECK_AERR
    aerr = adios2_attribute_data (buf, &asize, aid);
    CHECK_AERR

    aerr = adios2_attribute_type (&atype, aid);
    CHECK_AERR
    esize = adios2_type_size (atype);

    this->amount_RD += asize * esize;

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    return err;
}

int e3sm_io_driver_adios2::inq_att (int fid, int vid, std::string name, MPI_Offset *size){
    int err = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    adios2_variable *did;
    adios2_attribute *aid;
    size_t asize;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    if (vid == NC_GLOBAL) {
        E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_GET_ATT)
        aid = adios2_inquire_attribute (fp->iop, name.c_str ());
        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_GET_ATT)
        CHECK_APTR (aid)
    } else {
        char vname[1024];
        size_t tmp = 1024;

        did  = fp->dids[vid];
        aerr = adios2_variable_name (vname, &tmp, did);
        CHECK_AERR

        E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_GET_ATT)
        aid = adios2_inquire_variable_attribute (fp->iop, name.c_str (), vname, "/");
        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_GET_ATT)
        CHECK_APTR (aid)
    }

    aerr = adios2_attribute_size (&asize, aid);
    CHECK_AERR
   
    *size = (MPI_Offset) asize;

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    return err;
}

int e3sm_io_driver_adios2::put_varl (
    int fid, int vid, MPI_Datatype itype, void *buf, e3sm_io_op_mode mode) {
    int err = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    adios2_variable *did;
    adios2_mode iomode;
    size_t putsize, esize;
    adios2_type dtype, mtype;
    void *xbuf = NULL;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    did = fp->dids[vid];

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_CONVERT)
    // Size of selection
    aerr = adios2_selection_size (&putsize, did);
    CHECK_AERR

    // Memory tpye
    mtype = mpi_type_to_adios2_type (itype);
    // Variable type
    aerr = adios2_variable_type (&dtype, did);
    CHECK_AERR
    esize = adios2_type_size (dtype);
    // Need convert
    if (mtype != dtype) {
        xbuf = malloc (esize * putsize);
        CHECK_PTR (xbuf)

        err = adios2_type_convert (mtype, putsize, buf, dtype, xbuf);
        CHECK_ERR
    } else {
        xbuf = buf;
    }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_CONVERT)

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

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_PUT_VAR)
    aerr = adios2_put (fp->ep, did, xbuf, iomode);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_PUT_VAR)
    CHECK_AERR

    this->amount_WR += putsize * esize;

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    if (xbuf && (xbuf != buf)) { free (xbuf); }
    return err;
}

int e3sm_io_driver_adios2::put_vara (int fid,
                                     int vid,
                                     MPI_Datatype itype,
                                     MPI_Offset *start,
                                     MPI_Offset *count,
                                     void *buf,
                                     e3sm_io_op_mode mode) {
    int err = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    size_t i, ndim;
    adios2_variable *did;
    size_t astart[E3SM_IO_DRIVER_MAX_RANK], ablock[E3SM_IO_DRIVER_MAX_RANK];
    adios2_mode iomode;
    size_t putsize, esize;
    adios2_type dtype, mtype;
    void *xbuf = NULL;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    did  = fp->dids[vid];
    ndim = fp->ndims[vid];

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_SEL)
    if (start) {
        if (count) {
            putsize = 1;
            for (i = 0; i < ndim; i++) {
                astart[i] = (size_t)start[i];
                ablock[i] = (size_t)count[i];
                putsize *= count[i];
            }
            this->amount_WR += putsize;
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
    if (ndim > 0) {
        aerr = adios2_set_selection (did, ndim, astart, ablock);
        CHECK_AERR
    }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_SEL)

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_CONVERT)
    // Size of selection
    aerr = adios2_selection_size (&putsize, did);
    CHECK_AERR
    // Memory tpye
    mtype = mpi_type_to_adios2_type (itype);
    // Variable type
    aerr = adios2_variable_type (&dtype, did);
    CHECK_AERR
    esize = adios2_type_size (dtype);
    // Need convert
    if (mtype != dtype) {
        xbuf = malloc (esize * putsize);
        CHECK_PTR (xbuf)

        err = adios2_type_convert (mtype, putsize, buf, dtype, xbuf);
        CHECK_ERR
    } else {
        xbuf = buf;
    }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_CONVERT)

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

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_PUT_VAR)
    aerr = adios2_put (fp->ep, did, xbuf, iomode);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_PUT_VAR)
    CHECK_AERR

    this->amount_WR += putsize * esize;

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    if (xbuf && (xbuf != buf)) { free (xbuf); }
    return err;
}

int e3sm_io_driver_adios2::put_varn (int fid,
                                     int vid,
                                     MPI_Datatype itype,
                                     int nreq,
                                     MPI_Offset **starts,
                                     MPI_Offset **counts,
                                     void *buf,
                                     e3sm_io_op_mode mode) {
    int err = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    int i, j;
    size_t esize, rsize;
    size_t putsize, reqsize;
    int ndim;
    adios2_type dtype, mtype;
    adios2_variable *did;
    char *bufp;
    size_t start[E3SM_IO_DRIVER_MAX_RANK], block[E3SM_IO_DRIVER_MAX_RANK];
    size_t dims[E3SM_IO_DRIVER_MAX_RANK];
    adios2_mode iomode;
    void *xbuf = NULL;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    did  = fp->dids[vid];
    ndim = fp->ndims[vid];

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_SEL)
    aerr = adios2_variable_shape (dims, did);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_SEL)
    CHECK_AERR

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_CONVERT)
    // Size of selection
    putsize = 0;
    for (j = 0; j < nreq; j++) {
        reqsize = 1;
        for (i = 0; i < ndim; i++) { reqsize *= counts[j][i]; }
        putsize += reqsize;
    }
    // Memory tpye
    mtype = mpi_type_to_adios2_type (itype);
    // Variable type
    aerr = adios2_variable_type (&dtype, did);
    CHECK_AERR
    esize = adios2_type_size (dtype);
    // Need convert
    if (mtype != dtype) {
        xbuf = malloc (esize * putsize);
        CHECK_PTR (xbuf)

        err = adios2_type_convert (mtype, putsize, buf, dtype, xbuf);
        CHECK_ERR
    } else {
        xbuf = buf;
    }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_CONVERT)

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
    bufp = (char *)xbuf;
    for (i = 0; i < nreq; i++) {
        rsize = esize;
        for (j = 0; j < ndim; j++) { rsize *= counts[i][j]; }

        if (rsize) {
            E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_SEL)
            for (j = 0; j < ndim; j++) {
                start[j] = (size_t)starts[i][j];
                block[j] = (size_t)counts[i][j];
            }
            // Record dim replaced by time step
            if (start[0] + block[0] > dims[0]) { start[0] = 0; }

            aerr = adios2_set_selection (did, ndim, start, block);
            CHECK_AERR
            E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_SEL)

            E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_PUT_VAR)
            aerr = adios2_put (fp->ep, did, bufp, iomode);
            CHECK_AERR
            E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_PUT_VAR)

            bufp += rsize;
        }
    }

    this->amount_WR += putsize * esize;

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    if (xbuf && (xbuf != buf)) { free (xbuf); }
    return err;
}

int e3sm_io_driver_adios2::get_vara (int fid,
                                     int vid,
                                     MPI_Datatype itype,
                                     MPI_Offset *start,
                                     MPI_Offset *count,
                                     void *buf,
                                     e3sm_io_op_mode mode) {
    int err = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    size_t i, ndim;
    adios2_variable *did;
    size_t astart[E3SM_IO_DRIVER_MAX_RANK], ablock[E3SM_IO_DRIVER_MAX_RANK];
    adios2_mode iomode;
    size_t getsize, esize;
    adios2_type dtype, mtype;
    void *xbuf = NULL;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    did  = fp->dids[vid];
    ndim = fp->ndims[vid];

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_SEL)
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

    if (ndim > 0) {
        aerr = adios2_set_selection (did, ndim, astart, ablock);
        CHECK_AERR
    }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_SEL)

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_CONVERT)
    // Size of selection
    aerr = adios2_selection_size (&getsize, did);
    CHECK_AERR
    // Memory tpye
    mtype = mpi_type_to_adios2_type (itype);
    // Variable type
    aerr = adios2_variable_type (&dtype, did);
    CHECK_AERR
    esize = adios2_type_size (dtype);
    // Need convert
    if (mtype != dtype) {
        xbuf = malloc (esize * getsize);
        CHECK_PTR (xbuf)
    } else {
        xbuf = buf;
    }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_CONVERT)

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

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_GET_VAR)
    aerr = adios2_get (fp->ep, did, xbuf, iomode);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_GET_VAR)
    CHECK_AERR

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_CONVERT)
    if (mtype != dtype) {
        err = adios2_type_convert (dtype, getsize, xbuf, mtype, buf);
        CHECK_ERR
    }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_CONVERT)

    this->amount_RD += getsize * esize;

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    if (xbuf && (xbuf != buf)) { free (xbuf); }
    return err;
}
int e3sm_io_driver_adios2::get_varn (int fid,
                                     int vid,
                                     MPI_Datatype itype,
                                     int nreq,
                                     MPI_Offset **starts,
                                     MPI_Offset **counts,
                                     void *buf,
                                     e3sm_io_op_mode mode) {
    int err = 0;
    adios2_error aerr;
    adios2_file *fp = this->files[fid];
    int i, j;
    size_t esize, rsize;
    size_t getsize, reqsize;
    int ndim;
    adios2_type dtype, mtype;
    void *xbuf = NULL;
    adios2_variable *did;
    char *bufp;
    size_t start[E3SM_IO_DRIVER_MAX_RANK], block[E3SM_IO_DRIVER_MAX_RANK];
    size_t dims[E3SM_IO_DRIVER_MAX_RANK];
    adios2_mode iomode;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2)

    did  = fp->dids[vid];
    ndim = fp->ndims[vid];

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_SEL)
    aerr = adios2_variable_shape (dims, did);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_SEL)
    CHECK_AERR

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_CONVERT)
    // Size of selection
    getsize = 0;
    for (j = 0; j < nreq; j++) {
        reqsize = 1;
        for (i = 0; i < ndim; i++) { reqsize *= counts[j][i]; }
        getsize += reqsize;
    }
    // Memory tpye
    mtype = mpi_type_to_adios2_type (itype);
    // Variable type
    aerr = adios2_variable_type (&dtype, did);
    CHECK_AERR
    esize = adios2_type_size (dtype);
    // Need convert
    if (mtype != dtype) {
        xbuf = malloc (esize * getsize);
        CHECK_PTR (xbuf)
    } else {
        xbuf = buf;
    }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_CONVERT)

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
    bufp = (char *)xbuf;
    for (i = 0; i < nreq; i++) {
        rsize = esize;
        for (j = 0; j < ndim; j++) { rsize *= counts[i][j]; }

        if (rsize) {
            E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_SEL)
            for (j = 0; j < ndim; j++) {
                start[j] = (size_t)starts[i][j];
                block[j] = (size_t)counts[i][j];
            }
            // Record dim replaced by time step
            if (start[0] + block[0] > dims[0]) { start[0] = 0; }

            aerr = adios2_set_selection (did, ndim, start, block);
            CHECK_AERR
            E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_SEL)

            E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_GET_VAR)
            aerr = adios2_get (fp->ep, did, bufp, iomode);
            E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_GET_VAR)
            CHECK_AERR

            bufp += rsize;
        }
    }

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_ADIOS2_CONVERT)
    if (mtype != dtype) {
        err = adios2_type_convert (dtype, getsize, xbuf, mtype, buf);
        CHECK_ERR
    }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2_CONVERT)

    this->amount_RD += getsize * esize;

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_ADIOS2)
    if (xbuf && (xbuf != buf)) { free (xbuf); }
    return err;
}
