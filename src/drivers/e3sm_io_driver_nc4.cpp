/*********************************************************************
 *
 * Copyright (C) 2022, Northwestern University
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

#include <e3sm_io_driver_nc4.hpp>
#include <e3sm_io_profile.hpp>

#define CHECK_NCERR                                                                             \
    {                                                                                           \
        if (err != NC_NOERR) {                                                                  \
            printf ("Error in %s line %d function %s: %s (%d)\n", __FILE__, __LINE__, __func__, \
                    nc_strerror (err), err);                                                    \
            DEBUG_ABORT;                                                                        \
            goto err_out;                                                                       \
        }                                                                                       \
    }

inline nc_type mpitype2nctype (MPI_Datatype type) {
    if (type == MPI_DOUBLE) {
        return NC_DOUBLE;
    } else if (type == MPI_FLOAT) {
        return NC_FLOAT;
    } else if (type == MPI_INT) {
        return NC_INT;
    } else if (type == MPI_LONG_LONG) {
        return NC_INT64;
    } else if (type == MPI_CHAR) {
        return NC_CHAR;
    } else if (type == MPI_BYTE) {
        return NC_BYTE;
    }

    throw "Unsupported datatype";
}

bool e3sm_io_driver_nc4::compatible (std::string path) {
    int err;
    int ncid;
    err = nc_open (path.data (), NC_NOWRITE, &ncid);
    if (err == NC_NOERR) {
        nc_close (ncid);
        return true;
    }
    return false;
}

e3sm_io_driver_nc4::e3sm_io_driver_nc4 (e3sm_io_config *cfg) : e3sm_io_driver (cfg) {
    if ((cfg->chunksize != 0) || (cfg->filter != none)) {
        throw "Fitler and chunking is not supported by the NetCDF 4 driver";
    }
    this->amount_RD = 0;
    this->amount_WR = 0;
}

e3sm_io_driver_nc4::~e3sm_io_driver_nc4 () {}

int e3sm_io_driver_nc4::create (std::string path, MPI_Comm comm, MPI_Info info, int *fid) {
    int err, cmode;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4_OPEN)

    cmode = NC_CLOBBER | NC_NETCDF4;
#if defined HAVE_DECL_NC_NODIMSCALE_ATTACH && HAVE_DECL_NC_NODIMSCALE_ATTACH
    cmode |= NC_NODIMSCALE_ATTACH;
#endif

    err = nc_create_par (path.c_str (), cmode, comm, info, fid);
    CHECK_NCERR

    /* turn off fill mode for the entire file */
    err = nc_set_fill(*fid, NC_NOFILL, NULL);
    CHECK_NCERR

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4_OPEN)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)

err_out:
    return err;
}

int e3sm_io_driver_nc4::open (std::string path, MPI_Comm comm, MPI_Info info, int *fid) {
    int err;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4_OPEN)

    err = nc_open_par (path.c_str (), NC_NOWRITE, comm, info, fid);
    CHECK_NCERR

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4_OPEN)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)

err_out:
    return err;
}

int e3sm_io_driver_nc4::close (int fid) {
    int err;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4_CLOSE)

    err = nc_close (fid);
    CHECK_NCERR

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4_CLOSE)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)

err_out:
    return err;
}

int e3sm_io_driver_nc4::inq_file_info (int fid, MPI_Info *info) {
    *info = MPI_INFO_NULL;
    return NC_NOERR;
}

int e3sm_io_driver_nc4::inq_file_size (std::string path, MPI_Offset *size) {
    int err;
    struct stat file_stat;

    err = stat (path.c_str (), &file_stat);
    CHECK_ERR

    *size = (MPI_Offset) (file_stat.st_size);

err_out:
    return err;
}

int e3sm_io_driver_nc4::inq_put_size (MPI_Offset *size) {
    *size = this->amount_WR;
    return 0;
}

int e3sm_io_driver_nc4::inq_get_size (MPI_Offset *size) {
    *size = this->amount_RD;
    return 0;
}

int e3sm_io_driver_nc4::inq_malloc_size (MPI_Offset *size) {
    *size = 0;
    return NC_NOERR;
}
int e3sm_io_driver_nc4::inq_malloc_max_size (MPI_Offset *size) {
    *size = 0;
    return NC_NOERR;
}

int e3sm_io_driver_nc4::inq_rec_size (int fid, MPI_Offset *size) {
    int err;
    int dimid;
    size_t len;

    err = nc_inq_unlimdim (fid, &dimid);
    CHECK_NCERR

    err = nc_inq_dimlen (fid, dimid, &len);
    CHECK_NCERR

    *size = (MPI_Offset)len;

err_out:
    return err;
}

int e3sm_io_driver_nc4::expand_rec_size (int fid, MPI_Offset size) {
    int err, i, j, ndim, udimid, nvar;
    std::vector<int> dimids;    // Var dimids
    std::vector<size_t> start;  // Starting coordinate to write in var
    char buf[16];               // Dummy buffer

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4_RESIZE)

    err = nc_inq_unlimdim (fid, &udimid);
    CHECK_NCERR

    err = nc_inq_nvars (fid, &nvar);
    CHECK_NCERR

    for (i = 0; i < nvar; i++) {
        err = nc_inq_varndims (fid, i, &ndim);
        CHECK_NCERR

        if (ndim > 0) {
            if (dimids.size () < (size_t)ndim) {
                dimids.resize (ndim);
                start.resize (ndim);
            }
            err = nc_inq_vardimid (fid, i, dimids.data ());
            CHECK_NCERR

            if (dimids[0] == udimid) {
                start[0] = size - 1;
                for (j = 1; j < ndim; j++) { start[j] = 0; }

                // Write a cell to extend the variable
                err = nc_var_par_access (fid, i, NC_COLLECTIVE);
                CHECK_NCERR
                err = nc_put_var1 (fid, i, start.data (), buf);
                CHECK_NCERR
                err = nc_var_par_access (fid, i, NC_INDEPENDENT);
                CHECK_NCERR
            }
        }
    }

err_out:
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4_RESIZE)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)
    return err;
}

int e3sm_io_driver_nc4::def_var (
    int fid, std::string name, nc_type xtype, int ndim, int *dimids, int *varid) {
    int err;
    int i;
    int esize;           // Size of var type element
    int udimid;          // Record dim ID
    size_t *dlen = NULL; // Dim len
    MPI_Offset vsize;    // Var size

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4_DEF_VAR)

    err = nc_def_var (fid, name.c_str (), xtype, ndim, dimids, varid);
    CHECK_NCERR

    err = e3sm_io_xlen_nc_type (xtype, &esize);
    CHECK_ERR
    
    dlen = (size_t*)malloc(sizeof(size_t) * ndim);
    CHECK_PTR(dlen)
    vsize = esize;
    for (i = 0; i < ndim; i++) {
        nc_inq_dimlen (fid, dimids[i], dlen + i);
        vsize *= dlen[i];
    }
    var_size[{fid, *varid}] = vsize;


    err = nc_inq_unlimdim (fid, &udimid);
    CHECK_NCERR

    if (ndim && dimids[0] == udimid) {
        dlen[0] = 1;
        err = nc_def_var_chunking(fid, *varid, NC_CHUNKED, dlen);
        CHECK_ERR
    }

err_out:
    if(dlen) { free(dlen); }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4_DEF_VAR)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)
    return err;
}

int e3sm_io_driver_nc4::def_local_var (
    int fid, std::string name, nc_type xtype, int ndim, MPI_Offset *dsize, int *varid) {
    int err = 0;

    ERR_OUT ("NC4 does not support local variables")

err_out:
    return err;
}

int e3sm_io_driver_nc4::inq_varid (int fid, const char *name, int *varid) {
    int err;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4_INQ_VAR)

    err = nc_inq_varid (fid, name, varid);

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4_INQ_VAR)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)
    return err;
}

int e3sm_io_driver_nc4::inq_var (int fid, int varid, char *name, nc_type *xtypep,
                                 int *ndimsp, int *dimids, int *nattsp)
{
    int err;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4_INQ_VAR)

    err = nc_inq_var (fid, varid, name, xtypep, ndimsp, dimids, nattsp);

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4_INQ_VAR)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)
    return err;
}

int e3sm_io_driver_nc4::inq_var_name (int ncid, int varid, char *name) {
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4_INQ_VAR)

    int err = nc_inq_varname (ncid, varid, name);
    CHECK_NCERR

err_out:
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4_INQ_VAR)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)
    return err;
}

int e3sm_io_driver_nc4::inq_var_off (int fid, int vid, MPI_Offset *off) {
    int err;

    ERR_OUT ("NC4 driver does not support local variables")

err_out:
    return err;
}

int e3sm_io_driver_nc4::def_dim (int fid, std::string name, MPI_Offset size, int *dimid) {
    int err;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4_DEF_DIM)

    err = nc_def_dim (fid, name.c_str (), size, dimid);
    CHECK_NCERR

err_out:
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4_DEF_DIM)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)
    return err;
}

int e3sm_io_driver_nc4::inq_dim (int fid, std::string name, int *dimid) {
    int err;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4_INQ_DIM)

    err = nc_inq_dimid (fid, name.c_str (), dimid);
    CHECK_NCERR

err_out:
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4_INQ_DIM)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)
    return err;
}

int e3sm_io_driver_nc4::inq_dimlen (int fid, int dimid, MPI_Offset *size) {
    int err = 0;
    size_t len;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4_INQ_DIM)

    err = nc_inq_dimlen (fid, dimid, &len);
    CHECK_NCERR

    *size = len;

err_out:
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4_INQ_DIM)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)
    return err;
}

int e3sm_io_driver_nc4::enddef (int fid) {
    // NetCDF 4 enddef is automatic
    return 0;
}

int e3sm_io_driver_nc4::redef (int fid) {
    // NetCDF 4 redef is automatic
    return 0;
}

int e3sm_io_driver_nc4::wait (int fid) { return 0; }

int e3sm_io_driver_nc4::put_att (
    int fid, int vid, std::string name, nc_type xtype, MPI_Offset size, const void *buf) {
    int err;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4_PUT_ATT)

    err = nc_put_att (fid, vid, name.c_str (), xtype, size, buf);
    CHECK_NCERR

err_out:
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4_PUT_ATT)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)
    return err;
}

int e3sm_io_driver_nc4::get_att (int fid, int vid, std::string name, void *buf) {
    int err;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4_GET_ATT)

    err = nc_get_att (fid, vid, name.c_str (), buf);
    CHECK_NCERR

err_out:
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4_GET_ATT)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)
    return err;
}

int e3sm_io_driver_nc4::inq_att (int fid, int vid, std::string name, MPI_Offset *size){
    int err;
    size_t len;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)

    err = nc_inq_attlen (fid, vid, name.c_str (), &len);
    CHECK_NCERR
    *size = (MPI_Offset) len;

err_out:
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)
    return err;
}


int e3sm_io_driver_nc4::put_varl (
    int fid, int vid, MPI_Datatype itype, void *buf, e3sm_io_op_mode mode) {
    int err = 0;

    ERR_OUT ("NC4 does not support local variables")

err_out:
    return err;
}

int e3sm_io_driver_nc4::put_vara (int fid,
                                  int vid,
                                  MPI_Datatype itype,
                                  MPI_Offset *start,
                                  MPI_Offset *count,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err = NC_NOERR;
    int i;
    int xesize;        // Element size of variable type
    int ndim;        // Variable #dim
    size_t bsize;    // Block size
    nc_type xtype;    // Type of the variable

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4_PUT_VAR)

    err = nc_inq_vartype (fid, vid, &xtype);
    CHECK_NCERR
    err = e3sm_io_xlen_nc_type (xtype, &xesize);
    CHECK_ERR
    err = nc_inq_varndims (fid, vid, &ndim);
    CHECK_NCERR

    if (start) {
        if (count) {
            if (itype == MPI_FLOAT)
                err = nc_put_vara_float (fid, vid, (size_t *)start, (size_t *)count, (const float *)buf);
            else if (itype == MPI_DOUBLE)
                err = nc_put_vara_double (fid, vid, (size_t *)start, (size_t *)count, (const double *)buf);
            else if (itype == MPI_INT)
                err = nc_put_vara_int (fid, vid, (size_t *)start, (size_t *)count, (const int *)buf);
            else if (itype == MPI_LONG_LONG)
                err = nc_put_vara_longlong (fid, vid, (size_t *)start, (size_t *)count, (const long long *)buf);
            else if (itype == MPI_CHAR)
                err = nc_put_vara_text (fid, vid, (size_t *)start, (size_t *)count, (const char *)buf);
            else if (itype == MPI_UNSIGNED)
                err = nc_put_vara_uint (fid, vid, (size_t *)start, (size_t *)count, (const unsigned int *)buf);
            else if (itype == MPI_UNSIGNED_LONG_LONG)
                err = nc_put_vara_ulonglong (fid, vid, (size_t *)start, (size_t *)count, (const unsigned long long *)buf);
            else
                ERR_OUT ("Uknown type")

            bsize = xesize;
            for (i = 0; i < ndim; i++) { bsize *= count[i]; }
            this->amount_WR += bsize;
        } else {
            if (itype == MPI_FLOAT)
                err = nc_put_var1_float (fid, vid, (size_t *)start, (const float *)buf);
            else if (itype == MPI_DOUBLE)
                err = nc_put_var1_double (fid, vid, (size_t *)start, (const double *)buf);
            else if (itype == MPI_INT)
                err = nc_put_var1_int (fid, vid, (size_t *)start, (const int *)buf);
            else if (itype == MPI_LONG_LONG)
                err = nc_put_var1_longlong (fid, vid, (size_t *)start, (const long long *)buf);
            else if (itype == MPI_CHAR)
                err = nc_put_var1_text (fid, vid, (size_t *)start, (const char *)buf);
            else if (itype == MPI_UNSIGNED)
                err = nc_put_var1_uint (fid, vid, (size_t *)start, (const unsigned int *)buf);
            else if (itype == MPI_UNSIGNED_LONG_LONG)
                err = nc_put_var1_ulonglong (fid, vid, (size_t *)start, (const unsigned long long *)buf);
            else
                ERR_OUT ("Uknown type")

            this->amount_WR += xesize;
        }
    } else {
        if (itype == MPI_FLOAT)
            err = nc_put_var_float (fid, vid, (const float *)buf);
        else if (itype == MPI_DOUBLE)
            err = nc_put_var_double (fid, vid, (const double *)buf);
        else if (itype == MPI_INT)
            err = nc_put_var_int (fid, vid, (const int *)buf);
        else if (itype == MPI_LONG_LONG)
            err = nc_put_var_longlong (fid, vid, (const long long *)buf);
        else if (itype == MPI_CHAR)
            err = nc_put_var_text (fid, vid, (const char *)buf);
        else if (itype == MPI_UNSIGNED)
            err = nc_put_var_uint (fid, vid, (const unsigned int *)buf);
        else if (itype == MPI_UNSIGNED_LONG_LONG)
            err = nc_put_var_ulonglong (fid, vid, (const unsigned long long *)buf);
        else
            ERR_OUT ("Uknown type")

        this->amount_WR += var_size[{fid, vid}];
    }
    CHECK_NCERR

err_out:
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4_PUT_VAR)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)
    return err;
}

int e3sm_io_driver_nc4::put_varn (int fid,
                                  int vid,
                                  MPI_Datatype itype,
                                  int nreq,
                                  MPI_Offset **starts,
                                  MPI_Offset **counts,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err = NC_NOERR;
    int i, j;
    char *bufp = (char *)buf;  // Current block
    size_t bsize;               // Block size
    int esize;                   // Element size of memory type
    int xesize;                   // Element size of variable type
    int ndim;                   // Variable #dim
    nc_type xtype;               // Type of the variable

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4_PUT_VAR)

    err = nc_inq_varndims (fid, vid, &ndim);
    CHECK_NCERR
    err = nc_inq_vartype (fid, vid, &xtype);
    CHECK_NCERR
    err = e3sm_io_xlen_nc_type (xtype, &xesize);
    CHECK_NCERR
    err = MPI_Type_size (itype, &esize);
    CHECK_MPIERR

    for (i = 0; i < nreq; i++) {
        bsize = 1;
        for (j = 0; j < ndim; j++) { bsize *= counts[i][j]; }

        if (itype == MPI_FLOAT)
            err = nc_put_vara_float (fid, vid, (size_t *)(starts[i]), (size_t *)(counts[i]), (const float *)bufp);
        else if (itype == MPI_DOUBLE)
            err = nc_put_vara_double (fid, vid, (size_t *)(starts[i]), (size_t *)(counts[i]), (const double *)bufp);
        else if (itype == MPI_INT)
            err = nc_put_vara_int (fid, vid, (size_t *)(starts[i]), (size_t *)(counts[i]), (const int *)bufp);
        else if (itype == MPI_LONG_LONG)
            err = nc_put_vara_longlong (fid, vid, (size_t *)(starts[i]), (size_t *)(counts[i]), (const long long *)bufp);
        else if (itype == MPI_CHAR)
            err = nc_put_vara_text (fid, vid, (size_t *)(starts[i]), (size_t *)(counts[i]), (const char *)bufp);
        else if (itype == MPI_UNSIGNED)
            err = nc_put_vara_uint (fid, vid, (size_t *)(starts[i]), (size_t *)(counts[i]), (const unsigned int *)bufp);
        else if (itype == MPI_UNSIGNED_LONG_LONG)
            err = nc_put_vara_ulonglong (fid, vid, (size_t *)(starts[i]), (size_t *)(counts[i]), (const unsigned long long *)bufp);
        else
            ERR_OUT ("Uknown type")

        CHECK_NCERR
        bufp += bsize * esize;

        this->amount_WR += bsize * xesize;
    }

err_out:
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4_PUT_VAR)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)
    return err;
}

int e3sm_io_driver_nc4::get_vara (int fid,
                                  int vid,
                                  MPI_Datatype itype,
                                  MPI_Offset *start,
                                  MPI_Offset *count,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err = NC_NOERR;
    int i;
    int xesize;        // Element size of variable type
    int ndim;        // Variable #dim
    size_t bsize;    // Block size
    nc_type xtype;    // Type of the variable

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4_GET_VAR)

    err = nc_inq_vartype (fid, vid, &xtype);
    CHECK_NCERR
    err = e3sm_io_xlen_nc_type (xtype, &xesize);
    CHECK_ERR
    err = nc_inq_varndims (fid, vid, &ndim);
    CHECK_NCERR

    if (start) {
        if (count) {
            if (itype == MPI_FLOAT)
                err = nc_get_vara_float (fid, vid, (size_t *)start, (size_t *)count, (float *)buf);
            else if (itype == MPI_DOUBLE)
                err = nc_get_vara_double (fid, vid, (size_t *)start, (size_t *)count, (double *)buf);
            else if (itype == MPI_INT)
                err = nc_get_vara_int (fid, vid, (size_t *)start, (size_t *)count, (int *)buf);
            else if (itype == MPI_LONG_LONG)
                err = nc_get_vara_longlong (fid, vid, (size_t *)start, (size_t *)count, (long long *)buf);
            else if (itype == MPI_CHAR)
                err = nc_get_vara_text (fid, vid, (size_t *)start, (size_t *)count, (char *)buf);
            else if (itype == MPI_UNSIGNED)
                err = nc_get_vara_uint (fid, vid, (size_t *)start, (size_t *)count, (unsigned int *)buf);
            else if (itype == MPI_UNSIGNED_LONG_LONG)
                err = nc_get_vara_ulonglong (fid, vid, (size_t *)start, (size_t *)count, (unsigned long long *)buf);
            else
                ERR_OUT ("Uknown type")

            bsize = xesize;
            for (i = 0; i < ndim; i++) { bsize *= count[i]; }
            this->amount_RD += bsize;
        } else {
            if (itype == MPI_FLOAT)
                err = nc_get_var1_float (fid, vid, (size_t *)start, (float *)buf);
            else if (itype == MPI_DOUBLE)
                err = nc_get_var1_double (fid, vid, (size_t *)start, (double *)buf);
            else if (itype == MPI_INT)
                err = nc_get_var1_int (fid, vid, (size_t *)start, (int *)buf);
            else if (itype == MPI_LONG_LONG)
                err = nc_get_var1_longlong (fid, vid, (size_t *)start, (long long *)buf);
            else if (itype == MPI_CHAR)
                err = nc_get_var1_text (fid, vid, (size_t *)start, (char *)buf);
            else if (itype == MPI_UNSIGNED)
                err = nc_get_var1_uint (fid, vid, (size_t *)start, (unsigned int *)buf);
            else if (itype == MPI_UNSIGNED_LONG_LONG)
                err = nc_get_var1_ulonglong (fid, vid, (size_t *)start, (unsigned long long *)buf);
            else
                ERR_OUT ("Uknown type")

            this->amount_RD += xesize;
        }
    } else {
        if (itype == MPI_FLOAT)
            err = nc_get_var_float (fid, vid, (float *)buf);
        else if (itype == MPI_DOUBLE)
            err = nc_get_var_double (fid, vid, (double *)buf);
        else if (itype == MPI_INT)
            err = nc_get_var_int (fid, vid, (int *)buf);
        else if (itype == MPI_LONG_LONG)
            err = nc_get_var_longlong (fid, vid, (long long *)buf);
        else if (itype == MPI_CHAR)
            err = nc_get_var_text (fid, vid, (char *)buf);
        else if (itype == MPI_UNSIGNED)
            err = nc_get_var_uint (fid, vid, (unsigned int *)buf);
        else if (itype == MPI_UNSIGNED_LONG_LONG)
            err = nc_get_var_ulonglong (fid, vid, (unsigned long long *)buf);
        else
            ERR_OUT ("Uknown type")

        this->amount_RD += var_size[{fid, vid}];
    }
    CHECK_NCERR

err_out:
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4_GET_VAR)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)
    return err;
}

int e3sm_io_driver_nc4::get_varn (int fid,
                                  int vid,
                                  MPI_Datatype itype,
                                  int nreq,
                                  MPI_Offset **starts,
                                  MPI_Offset **counts,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err = NC_NOERR;
    int i, j;
    char *bufp = (char *)buf;  // Current block
    size_t bsize;               // Block size
    int esize;                   // Element size of memory type
    int xesize;                   // Element size of variable type
    int ndim;                   // Variable #dim
    nc_type xtype;               // Type of the variable

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_NC4_GET_VAR)

    err = nc_inq_varndims (fid, vid, &ndim);
    CHECK_NCERR
    err = nc_inq_vartype (fid, vid, &xtype);
    CHECK_NCERR
    err = e3sm_io_xlen_nc_type (xtype, &xesize);
    CHECK_NCERR
    err = MPI_Type_size (itype, &esize);
    CHECK_MPIERR

    for (i = 0; i < nreq; i++) {
        if (itype == MPI_FLOAT)
            err = nc_get_vara_float (fid, vid, (size_t *)(starts[i]), (size_t *)(counts[i]), (float *)bufp);
        else if (itype == MPI_DOUBLE)
            err = nc_get_vara_double (fid, vid, (size_t *)(starts[i]), (size_t *)(counts[i]), (double *)bufp);
        else if (itype == MPI_INT)
            err = nc_get_vara_int (fid, vid, (size_t *)(starts[i]), (size_t *)(counts[i]), (int *)bufp);
        else if (itype == MPI_LONG_LONG)
            err = nc_get_vara_longlong (fid, vid, (size_t *)(starts[i]), (size_t *)(counts[i]), (long long *)bufp);
        else if (itype == MPI_CHAR)
            err = nc_get_vara_text (fid, vid, (size_t *)(starts[i]), (size_t *)(counts[i]), (char *)bufp);
        else if (itype == MPI_UNSIGNED)
            err = nc_get_vara_uint (fid, vid, (size_t *)(starts[i]), (size_t *)(counts[i]), (unsigned int *)bufp);
        else if (itype == MPI_UNSIGNED_LONG_LONG)
            err = nc_get_vara_ulonglong (fid, vid, (size_t *)(starts[i]), (size_t *)(counts[i]), (unsigned long long *)bufp);
        else
            ERR_OUT ("Uknown type")

        CHECK_NCERR

        bsize = xesize;
        for (j = 0; j < ndim; j++) { bsize *= counts[i][j]; }
        this->amount_RD += bsize;

        bufp += bsize;
    }

err_out:
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4_GET_VAR)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_NC4)
    return err;
}
