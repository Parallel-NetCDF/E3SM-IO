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
//
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
//
#include <H5VL_log.h>
#include <hdf5.h>
//
#include <e3sm_io.h>
#include <e3sm_io_err.h>

#include <e3sm_io_driver_hdf5.hpp>
#include <e3sm_io_driver_hdf5_int.hpp>
#include <e3sm_io_driver_hdf5_log.hpp>
#include <e3sm_io_profile.hpp>

e3sm_io_driver_hdf5_log::e3sm_io_driver_hdf5_log (e3sm_io_config *cfg) : e3sm_io_driver_hdf5 (cfg) {
    int err     = 0;
    herr_t herr = 0;
    char *env;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

	this->dxplid_coll = H5Pcreate (H5P_DATASET_XFER);
	CHECK_HID (this->dxplid_coll)
	herr = H5Pset_dxpl_mpio (this->dxplid_coll, H5FD_MPIO_COLLECTIVE);
	CHECK_HERR
	this->dxplid_coll_nb = H5Pcreate (H5P_DATASET_XFER);
	CHECK_HID (this->dxplid_coll_nb)
	herr = H5Pset_dxpl_mpio (this->dxplid_coll_nb, H5FD_MPIO_COLLECTIVE);
	CHECK_HERR
	this->dxplid_indep_nb = H5Pcreate (H5P_DATASET_XFER);
	CHECK_HID (this->dxplid_indep_nb)

    herr = H5Pset_nonblocking (this->dxplid_coll, H5VL_LOG_REQ_BLOCKING);
    CHECK_HERR
    herr = H5Pset_nonblocking (this->dxplid_indep, H5VL_LOG_REQ_BLOCKING);
    CHECK_HERR
    herr = H5Pset_nonblocking (this->dxplid_coll_nb, H5VL_LOG_REQ_NONBLOCKING);
    CHECK_HERR
    herr = H5Pset_nonblocking (this->dxplid_indep_nb, H5VL_LOG_REQ_NONBLOCKING);
    CHECK_HERR

    env = getenv ("E3SM_IO_HDF5_USE_LOGVOL_WRITEN");
    if (env) {
        if (std::string (env) == "0") { this->use_logvol_varn = false; }
    }

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    if (err < 0) { throw "HDF5 log driver init fail"; }
}

e3sm_io_driver_hdf5_log::~e3sm_io_driver_hdf5_log () {
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    if (dxplid_coll >= 0) H5Pclose (dxplid_coll);
    //if (dxplid_indep >= 0) H5Pclose (dxplid_indep);
    if (dxplid_coll_nb >= 0) H5Pclose (dxplid_coll_nb);
    if (dxplid_indep_nb >= 0) H5Pclose (dxplid_indep_nb);

    // if (dxplid_nb >= 0) H5Pclose (dxplid_nb);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
}

int e3sm_io_driver_hdf5_log::create (std::string path, MPI_Comm comm, MPI_Info info, int *fid) {
    int err = 0;
    herr_t herr;
    hid_t faplid;
    H5AC_cache_config_t mdcc;
    hdf5_file *fp;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    fp = new hdf5_file (*this);

    err = MPI_Comm_dup(comm, &(fp->comm));
    CHECK_MPIERR

    err = MPI_Comm_rank (comm, &(fp->rank));
    CHECK_MPIERR

    faplid = H5Pcreate (H5P_FILE_ACCESS);
    CHECK_HID (faplid)
    herr = H5Pset_fapl_mpio (faplid, comm, info);
    CHECK_HERR
    herr = H5Pset_coll_metadata_write (faplid, true);
    CHECK_HERR
    herr = H5Pset_vol (faplid, this->log_vlid, NULL);
    CHECK_HERR
    // Enlarge metadata cache
    mdcc.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    herr         = H5Pget_mdc_config (faplid, &mdcc);
    CHECK_HERR
    mdcc.max_size         = 128 * 1024 * 1024;
    mdcc.min_size         = 32 * 1024 * 1024;
    mdcc.initial_size     = 32 * 1024 * 1024;
    mdcc.set_initial_size = true;
    herr                  = H5Pset_mdc_config (faplid, &mdcc);
    CHECK_HERR

    fp->id = H5Fcreate (path.c_str (), H5F_ACC_TRUNC, H5P_DEFAULT, faplid);
    CHECK_HID (fp->id)

    *fid = this->files.size ();
    this->files.push_back (fp);

err_out:;
    if (faplid >= 0) { H5Pclose (faplid); }

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

int e3sm_io_driver_hdf5_log::open (std::string path, MPI_Comm comm, MPI_Info info, int *fid) {
    int err = 0;
    herr_t herr;
    hid_t faplid;
    H5AC_cache_config_t mdcc;
    hdf5_file *fp;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    fp = new hdf5_file (*this);

    err = MPI_Comm_dup(comm, &(fp->comm));
    CHECK_MPIERR
    
    err = MPI_Comm_rank (comm, &(fp->rank));
    CHECK_MPIERR

    faplid = H5Pcreate (H5P_FILE_ACCESS);
    CHECK_HID (faplid)
    herr = H5Pset_fapl_mpio (faplid, comm, info);
    CHECK_HERR
    herr = H5Pset_coll_metadata_write (faplid, true);
    CHECK_HERR
    herr = H5Pset_vol (faplid, this->log_vlid, NULL);
    CHECK_HERR
    // Enlarge metadata cache
    mdcc.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    herr         = H5Pget_mdc_config (faplid, &mdcc);
    CHECK_HERR
    mdcc.max_size         = 128 * 1024 * 1024;
    mdcc.min_size         = 32 * 1024 * 1024;
    mdcc.initial_size     = 32 * 1024 * 1024;
    mdcc.set_initial_size = true;
    herr                  = H5Pset_mdc_config (faplid, &mdcc);
    CHECK_HERR

    fp->id = H5Fopen (path.c_str (), H5F_ACC_RDONLY, faplid);
    CHECK_HID (fp->id)

    *fid = this->files.size ();
    this->files.push_back (fp);

err_out:;
    if (faplid >= 0) { H5Pclose (faplid); }

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
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

            if (S_ISREG (st.st_mode)) {
                total_size += (MPI_Offset)size;
            } 
        }
        closedir (dir);
    }

    return total_size;
}
int e3sm_io_driver_hdf5_log::inq_file_size (std::string path, MPI_Offset *size) {
    int err = 0;
    struct stat file_stat;

    err = stat (path.c_str (), &file_stat);
    CHECK_ERR

    *size = (MPI_Offset) (file_stat.st_size) + get_dir_size (path + ".subfiles");

err_out:;
    return err;
}
int e3sm_io_driver_hdf5_log::put_vara (int fid,
                                       int vid,
                                       MPI_Datatype itype,
                                       MPI_Offset *start,
                                       MPI_Offset *count,
                                       void *buf,
                                       e3sm_io_op_mode mode) {
    int err = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    int i;
    int ndim   = -1;
    hid_t dsid = -1;
    hid_t did, h5_itype;
    hid_t dxplid;
    hsize_t hstart[E3SM_IO_DRIVER_MAX_RANK], hblock[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK];
    hid_t tid = -1;
    MPI_Offset putsize;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    did = fp->dids[vid];

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim = H5Sget_simple_extent_dims (dsid, dims, NULL);
    CHECK_HID (ndim)

    if (start) {
        if (count) {
            for (i = 0; i < ndim; i++) {
                hstart[i] = (hsize_t)start[i];
                hblock[i] = (hsize_t)count[i];
            }
        } else {
            for (i = 0; i < ndim; i++) {
                hstart[i] = (hsize_t)start[i];
                hblock[i] = 1;
            }
        }
    } else {
        for (i = 0; i < ndim; i++) {
            hstart[i] = 0;
            hblock[i] = dims[i];
        }
    }

    // Extend rec dim
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_EXT_DIM)
    if (ndim && (dims[0] < hstart[0] + hblock[0])) {
        dims[0] = hstart[0] + hblock[0];
        if (fp->recsize < (MPI_Offset) (dims[0])) { fp->recsize = dims[0]; }

        H5Sclose (dsid);
        herr = H5Dset_extent (did, dims);
        CHECK_HERR
        dsid = H5Dget_space (did);
        CHECK_HID (dsid)
    }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_EXT_DIM)

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)
    if (ndim > 0) {
        herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, hstart, NULL, this->one, hblock);
        CHECK_HERR
    }

    E3SM_IO_TIMER_SWAP (E3SM_IO_TIMER_HDF5_SEL, E3SM_IO_TIMER_HDF5_WR)
    switch (mode) {
        case nb: {
            dxplid = this->dxplid_indep_nb;
            break;
        }
        case nbe: {
            dxplid = this->dxplid_indep;
            break;
        }
        case indep: {
            dxplid = this->dxplid_indep;
            break;
        }
        case coll: {
            dxplid = this->dxplid_coll;
            break;
        }
        default:
            throw "Unrecognized mode";
    }

    h5_itype = mpi_type_to_hdf5_type (itype);

    // Call H5Dwrite
    herr = H5Dwrite (did, h5_itype, H5S_CONTIG, dsid, dxplid, buf);
    CHECK_HERR

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_WR)

    tid     = H5Dget_type (did);
    putsize = H5Tget_size (tid);
    for (i = 0; i < ndim; i++) { putsize *= hblock[i]; }
    this->amount_WR += putsize;

err_out:;
    if (tid >= 0) H5Tclose (tid);
    if (dsid >= 0) H5Sclose (dsid);

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

int e3sm_io_driver_hdf5_log::put_vars (int fid,
                                       int vid,
                                       MPI_Datatype itype,
                                       MPI_Offset *start,
                                       MPI_Offset *count,
                                       MPI_Offset *stride,
                                       void *buf,
                                       e3sm_io_op_mode mode) {
    int err = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    int i;
    int ndim   = -1;
    hid_t dsid = -1;
    hid_t did, h5_itype;
    hid_t dxplid;
    hsize_t hstart[E3SM_IO_DRIVER_MAX_RANK], hblock[E3SM_IO_DRIVER_MAX_RANK],
        hstride[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK];
    hid_t tid = -1;
    MPI_Offset putsize;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    did = fp->dids[vid];

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim = H5Sget_simple_extent_dims (dsid, dims, NULL);
    CHECK_HID (ndim)

    for (i = 0; i < ndim; i++) {
        hstart[i]  = (hsize_t)start[i];
        hblock[i]  = (hsize_t)count[i];
        hstride[i] = (hsize_t)stride[i];
    }

    // Extend rec dim
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_EXT_DIM)
    if (ndim && (dims[0] < hstart[0] + (hblock[0] - 1) * hstride[0] + 1)) {
        dims[0] = hstart[0] + (hblock[0] - 1) * hstride[0] + 1;
        if (fp->recsize < (MPI_Offset) (dims[0])) { fp->recsize = dims[0]; }

        H5Sclose (dsid);
        herr = H5Dset_extent (did, dims);
        CHECK_HERR
        dsid = H5Dget_space (did);
        CHECK_HID (dsid)
    }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_EXT_DIM)

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)
    herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, hstart, hstride, this->one, hblock);
    CHECK_HERR

    E3SM_IO_TIMER_SWAP (E3SM_IO_TIMER_HDF5_SEL, E3SM_IO_TIMER_HDF5_WR)
    switch (mode) {
        case nb: {
            dxplid = this->dxplid_indep_nb;
            break;
        }
        case nbe: {
            dxplid = this->dxplid_indep;
            break;
        }
        case indep: {
            dxplid = this->dxplid_indep;
            break;
        }
        case coll: {
            dxplid = this->dxplid_coll;
            break;
        }
        default:
            throw "Unrecognized mode";
    }

    h5_itype = mpi_type_to_hdf5_type (itype);

    // Call H5Dwrite
    herr = H5Dwrite (did, h5_itype, H5S_CONTIG, dsid, dxplid, buf);
    CHECK_HERR

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_WR)

    tid     = H5Dget_type (did);
    putsize = H5Tget_size (tid);
    for (i = 0; i < ndim; i++) { putsize *= hblock[i]; }
    this->amount_WR += putsize;

err_out:;
    if (tid >= 0) H5Tclose (tid);
    if (dsid >= 0) H5Sclose (dsid);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

int e3sm_io_driver_hdf5_log::put_varn (int fid,
                                       int vid,
                                       MPI_Datatype itype,
                                       int nreq,
                                       MPI_Offset **starts,
                                       MPI_Offset **counts,
                                       void *buf,
                                       e3sm_io_op_mode mode) {
    int err = 0;
    size_t tsize, putsize;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    int i, j;
    hsize_t esize, rsize = 0;
    int ndim;
    hid_t dsid = -1;
    hid_t mtype;
    char *bufp = (char *)buf;
    hid_t did;
    hid_t dxplid;
    hsize_t start[E3SM_IO_DRIVER_MAX_RANK], block[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK], mdims[E3SM_IO_DRIVER_MAX_RANK];
    hid_t tid         = -1;
    hsize_t **hstarts = NULL, **hcounts;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    did = fp->dids[vid];

    mtype = mpi_type_to_hdf5_type (itype);
    esize = (hsize_t)H5Tget_size (mtype);
    if (esize <= 0) { ERR_OUT ("Unknown memory type") }

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim = H5Sget_simple_extent_dims (dsid, dims, mdims);
    CHECK_HID (ndim)

    // Extend rec dim if needed
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_EXT_DIM)
    if (ndim && mdims[0] == H5S_UNLIMITED) {
        MPI_Offset max_rec = 0;
        for (i = 0; i < nreq; i++) {
            if (max_rec < starts[i][0] + counts[i][0]) { max_rec = starts[i][0] + counts[i][0]; }
        }
        if (dims[0] < (hsize_t)max_rec) {
            dims[0] = (hsize_t)max_rec;
            if (fp->recsize < (MPI_Offset) (dims[0])) { fp->recsize = dims[0]; }

            H5Sclose (dsid);
            herr = H5Dset_extent (did, dims);
            CHECK_HERR
            dsid = H5Dget_space (did);
            CHECK_HID (dsid)
        }
    }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_EXT_DIM)

    switch (mode) {
        case nb: {
            dxplid = this->dxplid_indep_nb;
            break;
        }
        case nbe: {
            dxplid = this->dxplid_indep;
            break;
        }
        case indep: {
            dxplid = this->dxplid_indep;
            break;
        }
        case coll: {
            dxplid = this->dxplid_coll;
            break;
        }
        default:
            throw "Unrecognized mode";
    }

    if (this->use_logvol_varn) {
        E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)
        // Convert starts and counts;
        hstarts    = (hsize_t **)malloc (sizeof (hsize_t *) * nreq * 2);
        hcounts    = hstarts + nreq;
        hstarts[0] = (hsize_t *)malloc (sizeof (hsize_t) * nreq * ndim * 2);
        hcounts[0] = hstarts[0] + nreq * ndim;
        for (i = 1; i < nreq; i++) {
            hstarts[i] = hstarts[i - 1] + ndim;
            hcounts[i] = hcounts[i - 1] + ndim;
        }
        for (i = 0; i < nreq; i++) {
            for (j = 0; j < ndim; j++) {
                hstarts[i][j] = (hsize_t)starts[i][j];
                hcounts[i][j] = (hsize_t)counts[i][j];
            }
        }
        E3SM_IO_TIMER_SWAP (E3SM_IO_TIMER_HDF5_SEL, E3SM_IO_TIMER_HDF5_WR)
        herr = H5Dwrite_n (did, mtype, nreq, hstarts, hcounts, dxplid, buf);
        CHECK_HERR
        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_WR)
    } else {
        // Call H5DWrite
        for (i = 0; i < nreq; i++) {
            rsize = esize;
            for (j = 0; j < ndim; j++) { rsize *= counts[i][j]; }

            if (rsize) {
                for (j = 0; j < ndim; j++) {
                    start[j] = (hsize_t)starts[i][j];
                    block[j] = (hsize_t)counts[i][j];
                }

                E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)
                herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, start, NULL, one, block);
                CHECK_HERR

                E3SM_IO_TIMER_SWAP (E3SM_IO_TIMER_HDF5_SEL, E3SM_IO_TIMER_HDF5_WR)

                herr = H5Dwrite (did, mtype, H5S_CONTIG, dsid, dxplid, bufp);
                CHECK_HERR

                E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_WR)

                bufp += rsize;
            }
        }
    }

    tid   = H5Dget_type (did);
    tsize = H5Tget_size (tid);

    for (j = 0; j < nreq; j++) {
        putsize = tsize;
        for (i = 0; i < ndim; i++) putsize *= counts[j][i];
        this->amount_WR += putsize;
    }

err_out:;
    if (tid >= 0) H5Tclose (tid);
    if (dsid >= 0) H5Sclose (dsid);
    if (hstarts) {
        free (hstarts[0]);
        free (hstarts);
    }

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

int e3sm_io_driver_hdf5_log::get_vara (int fid,
                                       int vid,
                                       MPI_Datatype itype,
                                       MPI_Offset *start,
                                       MPI_Offset *count,
                                       void *buf,
                                       e3sm_io_op_mode mode) {
    int err = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    int i;
    int ndim   = -1;
    hid_t dsid = -1, H5S_CONTIG = -1;
    hid_t did, h5_itype;
    hid_t dxplid;
    hsize_t hstart[E3SM_IO_DRIVER_MAX_RANK], hblock[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK];
    hid_t tid = -1;
    MPI_Offset getsize;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    did = fp->dids[vid];

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim = H5Sget_simple_extent_dims (dsid, dims, NULL);
    CHECK_HID (ndim)

    if (start) {
        if (count) {
            for (i = 0; i < ndim; i++) {
                hstart[i] = (hsize_t)start[i];
                hblock[i] = (hsize_t)count[i];
            }
        } else {
            for (i = 0; i < ndim; i++) {
                hstart[i] = (hsize_t)start[i];
                hblock[i] = 1;
            }
        }
    } else {
        for (i = 0; i < ndim; i++) {
            hstart[i] = 0;
            hblock[i] = dims[i];
        }
    }

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)
    herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, hstart, NULL, this->one, hblock);
    CHECK_HERR

    E3SM_IO_TIMER_SWAP (E3SM_IO_TIMER_HDF5_SEL, E3SM_IO_TIMER_HDF5_RD)
    switch (mode) {
        case nb: {
            dxplid = this->dxplid_indep_nb;
            break;
        }
        case nbe: {
            throw "bput not supported for read";
            break;
        }
        case indep: {
            dxplid = this->dxplid_indep;
            break;
        }
        case coll: {
            dxplid = this->dxplid_coll;
            break;
        }
        default:
            throw "Unrecognized mode";
    }

    h5_itype = mpi_type_to_hdf5_type (itype);

    // Call H5Dread
    herr = H5Dread (did, h5_itype, H5S_CONTIG, dsid, dxplid, buf);
    CHECK_HERR

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_RD)

    tid     = H5Dget_type (did);
    getsize = H5Tget_size (tid);
    for (i = 0; i < ndim; i++) { getsize *= hblock[i]; }
    fp->getsize += getsize;

err_out:;
    if (tid >= 0) H5Sclose (tid);
    if (dsid >= 0) H5Sclose (dsid);

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}
int e3sm_io_driver_hdf5_log::get_vars (int fid,
                                       int vid,
                                       MPI_Datatype itype,
                                       MPI_Offset *start,
                                       MPI_Offset *count,
                                       MPI_Offset *stride,
                                       void *buf,
                                       e3sm_io_op_mode mode) {
    int err = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    int i;
    int ndim   = -1;
    hid_t dsid = -1, H5S_CONTIG = -1;
    hid_t did, h5_itype;
    hid_t dxplid;
    hsize_t hstart[E3SM_IO_DRIVER_MAX_RANK], hblock[E3SM_IO_DRIVER_MAX_RANK],
        hstride[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK];
    hid_t tid = -1;
    MPI_Offset getsize;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    did = fp->dids[vid];

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim = H5Sget_simple_extent_dims (dsid, dims, NULL);
    CHECK_HID (ndim)

    for (i = 0; i < ndim; i++) {
        hstart[i]  = (hsize_t)start[i];
        hblock[i]  = (hsize_t)count[i];
        hstride[i] = (hsize_t)stride[i];
    }

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)
    herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, hstart, hstride, this->one, hblock);
    CHECK_HERR

    E3SM_IO_TIMER_SWAP (E3SM_IO_TIMER_HDF5_SEL, E3SM_IO_TIMER_HDF5_RD)
    switch (mode) {
        case nb: {
            dxplid = this->dxplid_indep_nb;
            break;
        }
        case nbe: {
            throw "bput not supported for read";
            break;
        }
        case indep: {
            dxplid = this->dxplid_indep;
            break;
        }
        case coll: {
            dxplid = this->dxplid_coll;
            break;
        }
        default:
            throw "Unrecognized mode";
    }

    h5_itype = mpi_type_to_hdf5_type (itype);

    // Call H5Dread
    herr = H5Dread (did, h5_itype, H5S_CONTIG, dsid, dxplid, buf);
    CHECK_HERR

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_RD)

    tid     = H5Dget_type (did);
    getsize = H5Tget_size (tid);
    for (i = 0; i < ndim; i++) { getsize *= hblock[i]; }
    fp->getsize += getsize;

err_out:;
    if (tid >= 0) H5Sclose (tid);
    if (dsid >= 0) H5Sclose (dsid);

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

int e3sm_io_driver_hdf5_log::get_varn (int fid,
                                       int vid,
                                       MPI_Datatype itype,
                                       int nreq,
                                       MPI_Offset **starts,
                                       MPI_Offset **counts,
                                       void *buf,
                                       e3sm_io_op_mode mode) {
    int err = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    int i, j;
    hsize_t esize, rsize = 0;
    int ndim;
    hid_t dsid = -1;
    hid_t mtype;
    char *bufp = (char *)buf;
    hid_t did;
    hid_t dxplid;
    hsize_t start[E3SM_IO_DRIVER_MAX_RANK], block[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK], mdims[E3SM_IO_DRIVER_MAX_RANK];
    hid_t tid = -1;
    MPI_Offset getsize;
    MPI_Offset **count;
    hsize_t **hstarts = NULL, **hcounts;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    did = fp->dids[vid];

    mtype = mpi_type_to_hdf5_type (itype);
    esize = (hsize_t)H5Tget_size (mtype);
    if (esize <= 0) { ERR_OUT ("Unknown memory type") }

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim = H5Sget_simple_extent_dims (dsid, dims, mdims);
    CHECK_HID (ndim)

    switch (mode) {
        case nb: {
            dxplid = this->dxplid_indep_nb;
            break;
        }
        case nbe: {
            throw "bput not supported in read";
            break;
        }
        case indep: {
            dxplid = this->dxplid_indep;
            break;
        }
        case coll: {
            dxplid = this->dxplid_coll;
            break;
        }
        default:
            throw "Unrecognized mode";
    }

    if (this->use_logvol_varn) {
        E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)
        // Convert starts and counts;
        hstarts    = (hsize_t **)malloc (sizeof (hsize_t *) * nreq * 2);
        hcounts    = hstarts + nreq;
        hstarts[0] = (hsize_t *)malloc (sizeof (hsize_t) * nreq * ndim * 2);
        hcounts[0] = hstarts[0] + nreq * ndim;
        for (i = 1; i < nreq; i++) {
            hstarts[i] = hstarts[i - 1] + ndim;
            hcounts[i] = hcounts[i - 1] + ndim;
        }
        for (i = 1; i < nreq; i++) {
            for (j = 0; j < ndim; j++) {
                hstarts[i][j] = (hsize_t)starts[i][j];
                hcounts[i][j] = (hsize_t)counts[i][j];
            }
        }
        E3SM_IO_TIMER_SWAP (E3SM_IO_TIMER_HDF5_SEL,E3SM_IO_TIMER_HDF5_RD)
        herr = H5Dread_n (did, mtype, nreq, hstarts, hcounts, dxplid, buf);
        CHECK_HERR
        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_RD)
    } else {
        // Call H5DWrite
        for (i = 0; i < nreq; i++) {
            rsize = esize;
            for (j = 0; j < ndim; j++) { rsize *= counts[i][j]; }

            if (rsize) {
                for (j = 0; j < ndim; j++) {
                    start[j] = (hsize_t)starts[i][j];
                    block[j] = (hsize_t)counts[i][j];
                }

                E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)
                herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, start, NULL, one, block);
                CHECK_HERR

                E3SM_IO_TIMER_SWAP (E3SM_IO_TIMER_HDF5_SEL, E3SM_IO_TIMER_HDF5_RD)
                herr = H5Dread (did, mtype, H5S_CONTIG, dsid, dxplid, bufp);
                CHECK_HERR
                E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_RD)

                bufp += rsize;
            }
        }
    }

    tid     = H5Dget_type (did);
    getsize = H5Tget_size (tid);
    for (count = counts; count < counts + nreq; count++) {
        for (i = 0; i < ndim; i++) { getsize *= *count[i]; }
    }
    fp->getsize += getsize;

err_out:;
    if (tid >= 0) H5Tclose (tid);
    if (dsid >= 0) H5Sclose (dsid);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}
