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
#include <cstdlib>
//
#include <sys/stat.h>
//
#include <hdf5.h>
//
#include <e3sm_io.h>
#include <e3sm_io_err.h>

#include <e3sm_io_driver_hdf5.hpp>
#include <e3sm_io_driver_hdf5_int.hpp>
#include <e3sm_io_profile.hpp>

e3sm_io_driver_hdf5::e3sm_io_driver_hdf5 (e3sm_io_config *cfg) : e3sm_io_driver (cfg) {
    int err     = 0;
    herr_t herr = 0;
    int i;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    this->dxplid_coll = H5Pcreate (H5P_DATASET_XFER);
    CHECK_HID (this->dxplid_coll)
    herr = H5Pset_dxpl_mpio (this->dxplid_coll, H5FD_MPIO_COLLECTIVE);
    CHECK_HERR
    this->dxplid_coll_nb = H5Pcreate (H5P_DATASET_XFER);
    CHECK_HID (this->dxplid_coll_nb)
    herr = H5Pset_dxpl_mpio (this->dxplid_coll_nb, H5FD_MPIO_COLLECTIVE);
    CHECK_HERR

    this->dxplid_indep = H5Pcreate (H5P_DATASET_XFER);
    CHECK_HID (this->dxplid_indep)
    this->dxplid_indep_nb = H5Pcreate (H5P_DATASET_XFER);
    CHECK_HID (this->dxplid_indep_nb)
    // dxplid_nb = H5Pcreate (H5P_DATASET_XFER);
    // CHECK_HID (dxplid_nb)

    for (i = 0; i < E3SM_IO_DRIVER_MAX_RANK; i++) { one[i] = 1; }

    if (cfg->api == hdf5_md) {
#ifdef HDF5_HAVE_DWRITE_MULTI
        this->use_dwrite_multi = true;
#else
        throw "The HDF5 used does not support multi-dataset write";
#endif
    }
#ifdef E3SM_IO_DEBUG
    else if (cfg->api == hdf5_ra) {
    }
#endif

    if ((cfg->chunksize != 0) && (cfg->filter != none)) {
        throw "Fitler requries chunking in HDF5";
    }

    if (cfg->num_group != 1) { throw "Subfiling not supported by HDF5 driver"; }

    /*
    env = getenv ("E3SM_IO_HDF5_ENABLE_LOGVOL");
    if (env) {
        if (std::string (env) == "1") {
            this->use_logvol = true;
            env              = getenv ("E3SM_IO_HDF5_USE_LOGVOL_WRITEN");
            if (env) {
                if (std::string (env) == "1") { this->use_logvol_varn = true; }
            }
        }
    }

    env = getenv ("E3SM_IO_HDF5_MERGE_VARN");
    if (env) {
        if (std::string (env) == "1") { this->merge_varn = true; }
    }
    */

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    if (err < 0) { throw "HDF5 driver init fail"; }
}
e3sm_io_driver_hdf5::~e3sm_io_driver_hdf5 () {
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    if (dxplid_coll >= 0) H5Pclose (dxplid_coll);
    if (dxplid_indep >= 0) H5Pclose (dxplid_indep);
    if (dxplid_coll_nb >= 0) H5Pclose (dxplid_coll_nb);
    if (dxplid_indep_nb >= 0) H5Pclose (dxplid_indep_nb);
    // if (dxplid_nb >= 0) H5Pclose (dxplid_nb);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
}

int e3sm_io_driver_hdf5::create (std::string path, MPI_Comm comm, MPI_Info info, int *fid) {
    int err = 0;
    herr_t herr;
    hid_t faplid;
    H5AC_cache_config_t mdcc;
    hdf5_file *fp;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    fp = new hdf5_file (*this);

    err = MPI_Comm_rank (comm, &(fp->rank));
    CHECK_MPIERR

    faplid = H5Pcreate (H5P_FILE_ACCESS);
    CHECK_HID (faplid)
    herr = H5Pset_fapl_mpio (faplid, comm, info);
    CHECK_HERR
    herr = H5Pset_coll_metadata_write (faplid, true);
    CHECK_HERR
    // Enlarge metadata cache
    mdcc.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    herr         = H5Pget_mdc_config (faplid, &mdcc);
    CHECK_HERR
    mdcc.max_size         = 128 * 1024 * 1024;
    mdcc.min_size         = 32 * 1024 * 1024;
    mdcc.initial_size     = 32 * 1024 * 1024;
    mdcc.evictions_enabled = false;
    mdcc.incr_mode = H5C_incr__off;
    mdcc.decr_mode = H5C_decr__off;
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

int e3sm_io_driver_hdf5::open (std::string path, MPI_Comm comm, MPI_Info info, int *fid) {
    int err = 0;
    herr_t herr;
    hid_t faplid;
    H5AC_cache_config_t mdcc;
    hdf5_file *fp;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    fp = new hdf5_file (*this);

    err = MPI_Comm_rank (comm, &(fp->rank));
    CHECK_MPIERR

    faplid = H5Pcreate (H5P_FILE_ACCESS);
    CHECK_HID (faplid)
    herr = H5Pset_fapl_mpio (faplid, comm, info);
    CHECK_HERR
    herr = H5Pset_coll_metadata_write (faplid, true);
    CHECK_HERR
    // Enlarge metadata cache
    mdcc.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    herr         = H5Pget_mdc_config (faplid, &mdcc);
    CHECK_HERR
    mdcc.max_size         = 128 * 1024 * 1024;
    mdcc.min_size         = 32 * 1024 * 1024;
    mdcc.initial_size     = 32 * 1024 * 1024;
    mdcc.evictions_enabled = false;
    mdcc.incr_mode = H5C_incr__off;
    mdcc.decr_mode = H5C_decr__off;
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

int e3sm_io_driver_hdf5::close (int fid) {
    int err = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    for (auto did : fp->dids) {
        herr = H5Dclose (did);
        CHECK_HERR
    }

    herr = H5Fclose (fp->id);
    CHECK_HERR

    delete fp;

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

int e3sm_io_driver_hdf5::inq_file_info (int fid, MPI_Info *info) {
    int err = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    hid_t pid;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    pid = H5Fget_access_plist (fp->id);
    CHECK_HID (pid);
    herr = H5Pget_fapl_mpio (pid, NULL, info);
    CHECK_HERR

err_out:;
    if (pid != -1) H5Pclose (pid);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}
int e3sm_io_driver_hdf5::inq_file_size (std::string path, MPI_Offset *size) {
    int err = 0;
    struct stat file_stat;

    err = stat (path.c_str (), &file_stat);
    CHECK_ERR

    *size = (MPI_Offset) (file_stat.st_size);

err_out:;
    return err;
}
int e3sm_io_driver_hdf5::inq_put_size (MPI_Offset *size) {
    *size = this->amount_WR;
    return 0;
}
int e3sm_io_driver_hdf5::inq_get_size (MPI_Offset *size) {
    *size = this->amount_RD;
    return 0;
}
int e3sm_io_driver_hdf5::inq_malloc_size (MPI_Offset *size) {
    *size = 0;
    return 0;
}
int e3sm_io_driver_hdf5::inq_malloc_max_size (MPI_Offset *size) {
    *size = 0;
    return 0;
}
int e3sm_io_driver_hdf5::inq_rec_size (int fid, MPI_Offset *size) {
    hdf5_file *fp = this->files[fid];

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    *size = (MPI_Offset) (fp->recsize);

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return 0;
}

int e3sm_io_driver_hdf5::def_var (
    int fid, std::string name, nc_type xtype, int ndim, int *dimids, int *did) {
    int err = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    int i;
    hid_t h5did, h5_xtype;
    hid_t sid    = -1;
    hid_t dcplid = -1;
    hsize_t cdim[E3SM_IO_DRIVER_MAX_RANK], dims[E3SM_IO_DRIVER_MAX_RANK],
        mdims[E3SM_IO_DRIVER_MAX_RANK];
    size_t csize = 0;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    dcplid = H5Pcreate (H5P_DATASET_CREATE);
    CHECK_HID (dcplid)

    H5Pset_fill_value (dcplid, 0, NULL);
    H5Pset_fill_time (dcplid, H5D_FILL_TIME_NEVER);
    H5Pset_alloc_time (dcplid, H5D_ALLOC_TIME_DEFAULT);

    h5_xtype = e3sm_io_type_nc2hdf5(xtype);

    for (i = 0; i < ndim; i++) { dims[i] = mdims[i] = fp->dsizes[dimids[i]]; }
    if (ndim) {
        if ((cfg->chunksize > 0) || (dims[0] == H5S_UNLIMITED)) {
            csize = H5Tget_size (h5_xtype);
            for (i = 0; i < ndim; i++) {
                if (csize < cfg->chunksize) {
                    cdim[i] = mdims[i];
                    csize *= cdim[i];
                } else {
                    cdim[i] = 1;
                }
            }
            // Chunk size along rec dim is always 1
            if (dims[0] == H5S_UNLIMITED) {
                cdim[0] = 1;
                dims[0] = 0;
            }
        }

        if (csize > 0) {
            herr = H5Pset_chunk (dcplid, ndim, cdim);
            CHECK_HERR

            switch (cfg->filter) {
                case none:
                    break;
                case deflate:
                    herr = H5Pset_deflate (dcplid, 6);
                    CHECK_HERR
                    break;
                default:
                    ERR_OUT ("Unknown filter")
            }
        }
    }

    sid = H5Screate_simple (ndim, dims, mdims);
    CHECK_HID (sid);

    h5did = H5Dcreate2 (fp->id, name.c_str (), h5_xtype, sid, H5P_DEFAULT,
                        dcplid, H5P_DEFAULT);
    CHECK_HID (h5did)

    *did = fp->dids.size ();
    fp->dids.push_back (h5did);

err_out:;
    if (sid != -1) H5Sclose (sid);
    if (dcplid != -1) H5Pclose (dcplid);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

int e3sm_io_driver_hdf5::def_local_var (
    int fid, std::string name, nc_type xtype, int ndim, MPI_Offset *dsize, int *did) {
    int err = 0;

    ERR_OUT ("HDF5 does not support local variables")

err_out:;
    return err;
}

int e3sm_io_driver_hdf5::inq_var (int fid, std::string name, int *did) {
    int err = 0;
    hdf5_file *fp = this->files[fid];
    hid_t h5did;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    h5did = H5Dopen2 (fp->id, name.c_str (), H5P_DEFAULT);
    CHECK_HID (h5did)

    *did = fp->dids.size ();
    fp->dids.push_back (h5did);

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

int e3sm_io_driver_hdf5::inq_var_name (int fid, int vid, char *name) {
    name[0] = '\0';
    printf("inq_var_name is not yet implementaed\n");
    return -1;
}

int e3sm_io_driver_hdf5::inq_var_off (int fid, int vid, MPI_Offset *off) {
    throw "Function not supported";
    return -1;
}
int e3sm_io_driver_hdf5::def_dim (int fid, std::string name, MPI_Offset size, int *dimid) {
    int err = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    hid_t sid     = -1;
    hid_t aid     = -1;
    hsize_t hsize;
    char aname[128];

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    hsize = (hsize_t)size;
    if (hsize == NC_UNLIMITED) hsize = H5S_UNLIMITED;

    sid = H5Screate (H5S_SCALAR);
    CHECK_HID (sid)

    sprintf (aname, "_DIM_%s", name.c_str ());
    aid = H5Acreate2 (fp->id, aname, H5T_NATIVE_HSIZE, sid, H5P_DEFAULT, H5P_DEFAULT);
    CHECK_HID (aid)

    herr = H5Awrite (aid, H5T_NATIVE_HSIZE, &hsize);
    CHECK_HERR

    *dimid = fp->dsizes.size ();
    fp->dsizes.push_back (hsize);

err_out:;
    if (aid != -1) H5Aclose (aid);
    if (sid != -1) H5Sclose (sid);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

int e3sm_io_driver_hdf5::inq_dim (int fid, std::string name, int *dimid) {
    int err = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    hid_t aid     = -1;
    hsize_t hsize;
    char aname[128];

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    sprintf (aname, "_DIM_%s", name.c_str ());
    aid = H5Aopen (fp->id, aname, H5P_DEFAULT);
    CHECK_HID (aid)

    herr = H5Aread (aid, H5T_NATIVE_HSIZE, &hsize);
    CHECK_HERR

    *dimid = fp->dsizes.size ();
    fp->dsizes.push_back (hsize);

err_out:;
    if (aid != -1) H5Aclose (aid);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

int e3sm_io_driver_hdf5::inq_dimlen (int fid, int dimid, MPI_Offset *size) {
    hdf5_file *fp = this->files[fid];

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    *size = (MPI_Offset) (fp->dsizes[dimid]);

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return 0;
}

int e3sm_io_driver_hdf5::enddef (int fid) { return 0; }
int e3sm_io_driver_hdf5::redef (int fid) { return 0; }
int e3sm_io_driver_hdf5::wait (int fid) {
    int err = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

// #ifdef HDF5_HAVE_DWRITE_MULTI
    err = fp->flush_multidatasets ();
    CHECK_ERR
// #endif

    herr = H5Fflush (fp->id, H5F_SCOPE_GLOBAL);
    CHECK_HERR

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

int e3sm_io_driver_hdf5::put_att (
    int fid, int vid, std::string name, nc_type xtype, MPI_Offset size, const void *buf) {
    int err = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    hid_t asid = -1, aid = -1;
    hid_t did, h5_xtype = -1;
    hsize_t asize;
    htri_t exists;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    asize = (size_t)size;

    if (vid == NC_GLOBAL)
        did = fp->id;
    else
        did = fp->dids[vid];

    if (xtype == NC_CHAR) { /* string attribute */
        asid = H5Screate(H5S_SCALAR);
        CHECK_HID(asid)
        h5_xtype = H5Tcopy(H5T_C_S1);
        CHECK_HID(h5_xtype)
        herr = H5Tset_size(h5_xtype, asize+1);
        CHECK_HERR
        herr = H5Tset_strpad(h5_xtype, H5T_STR_NULLTERM);
        CHECK_HERR
    }
    else {
        asid = H5Screate_simple (1, &asize, &asize);
        CHECK_HID (asid)
        h5_xtype = e3sm_io_type_nc2hdf5(xtype);
    }

    exists = H5Aexists (did, name.c_str ());
    if (!exists) {
        aid = H5Acreate2 (did, name.c_str (), h5_xtype, asid, H5P_DEFAULT,
                          H5P_DEFAULT);
    } else {
        aid = H5Aopen (did, name.c_str (), H5P_DEFAULT);
    }
    CHECK_HID (aid)

    herr = H5Awrite (aid, h5_xtype, buf);
    CHECK_HERR

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)

    if (fp->rank == 0) {
        size_t esize = H5Tget_size(h5_xtype);
        this->amount_WR += asize * esize;
    }


err_out:;
    if (asid >= 0) H5Sclose (asid);
    if (aid >= 0) H5Aclose (aid);
    if (xtype == NC_CHAR && h5_xtype >= 0) H5Tclose(h5_xtype); /* string attribute */

    return err;
}

int e3sm_io_driver_hdf5::get_att (int fid, int vid, std::string name, void *buf) {
    int err = 0;
    herr_t herr;
    int esize;
    hdf5_file *fp = this->files[fid];
    hid_t asid = -1, aid = -1;
    hid_t did;
    hid_t tid;
    hsize_t asize;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    if (vid == NC_GLOBAL)
        did = fp->id;
    else
        did = fp->dids[vid];

    aid = H5Aopen (did, name.c_str (), H5P_DEFAULT);
    CHECK_HID (aid)
    tid = H5Aget_type (aid);
    CHECK_HID (tid)
    herr = H5Aread (aid, tid, buf);
    CHECK_HERR

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)

    asid = H5Aget_space (aid);
    CHECK_HID (asid)
    H5Sget_simple_extent_dims (asid, &asize, NULL);
    esize = H5Tget_size (tid);
    fp->getsize += asize * esize;

err_out:;
    if (asid >= 0) H5Sclose (asid);
    if (aid >= 0) H5Aclose (aid);
    if (tid >= 0) H5Tclose (tid);

    return err;
}

int e3sm_io_driver_hdf5::put_varl (
    int fid, int vid, MPI_Datatype itype, void *buf, e3sm_io_op_mode mode) {
    int err = 0;

    ERR_OUT ("HDF5 does not support local variables")

err_out:;
    return err;
}

int e3sm_io_driver_hdf5::put_vara (int fid,
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
    hid_t dsid = -1, msid = -1;
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
        if (fp->recsize < (MPI_Offset)(dims[0])) { fp->recsize = dims[0]; }

        H5Sclose (dsid);
        herr = H5Dset_extent (did, dims);
        CHECK_HERR
        dsid = H5Dget_space (did);
        CHECK_HID (dsid)
    }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_EXT_DIM)

    // Init memory space
    msid = H5Screate_simple (ndim, hblock, hblock);
    CHECK_HID (msid)

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
            throw "bput not supported";
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

    h5_itype = mpi_type_to_hdf5_type(itype);

    // Call H5Dwrite immediately if blocking or Log VOL is used
    if ((mode == indep) || (mode == coll)) {
        herr = H5Dwrite (did, h5_itype, msid, dsid, dxplid, buf);
        CHECK_HERR
    }
// #ifdef HDF5_HAVE_DWRITE_MULTI
    else {  // Otherwier, queue request in driver
        herr = fp->register_multidataset (buf, did, dsid, msid, h5_itype, 1);
        CHECK_HERR
        // Prevent freeing of dsid and msid, they will be freed after flush
        dsid = msid = -1;
    }
// #endif

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_WR)

    tid     = H5Dget_type (did);
    putsize = H5Tget_size (tid);
    for (i = 0; i < ndim; i++) { putsize *= hblock[i]; }
    this->amount_WR += putsize;

err_out:;
    if (tid >= 0) H5Tclose (tid);
    if (dsid >= 0) H5Sclose (dsid);
    if (msid >= 0) H5Sclose (msid);

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

int e3sm_io_driver_hdf5::put_vars (int fid,
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
    hid_t dsid = -1, msid = -1;
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
        if (fp->recsize < (MPI_Offset)(dims[0])) { fp->recsize = dims[0]; }

        H5Sclose (dsid);
        herr = H5Dset_extent (did, dims);
        CHECK_HERR
        dsid = H5Dget_space (did);
        CHECK_HID (dsid)
    }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_EXT_DIM)

    // Init memory space
    msid = H5Screate_simple (ndim, hblock, hblock);
    CHECK_HID (msid)

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
            throw "bput not supported";
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

    h5_itype = mpi_type_to_hdf5_type(itype);

    // Call H5Dwrite immediately if blocking or Log VOL is used
    if ((mode == indep) || (mode == coll)) {
        herr = H5Dwrite (did, h5_itype, msid, dsid, dxplid, buf);
        CHECK_HERR
    }
// #ifdef HDF5_HAVE_DWRITE_MULTI
    else {  // Otherwier, queue request in driver
        herr = fp->register_multidataset (buf, did, dsid, msid, h5_itype, 1);
        CHECK_HERR
        // Prevent freeing of dsid and msid, they will be freed after flush
        dsid = msid = -1;
    }
// #endif

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_WR)

    tid     = H5Dget_type (did);
    putsize = H5Tget_size (tid);
    for (i = 0; i < ndim; i++) { putsize *= hblock[i]; }
    this->amount_WR += putsize;

err_out:;
    if (tid >= 0) H5Tclose (tid);
    if (dsid >= 0) H5Sclose (dsid);
    if (msid >= 0) H5Sclose (msid);

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}
int e3sm_io_driver_hdf5::put_varn (int fid,
                                   int vid,
                                   MPI_Datatype itype,
                                   int nreq,
                                   MPI_Offset **starts,
                                   MPI_Offset **counts,
                                   void *buf,
                                   e3sm_io_op_mode mode) {
// #ifdef HDF5_HAVE_DWRITE_MULTI
    if (this->merge_varn)
        return this->put_varn_merge (fid, vid, itype, nreq, starts, counts, buf, mode);
    else
// #endif
        return this->put_varn_expand (fid, vid, itype, nreq, starts, counts, buf, mode);
}
int e3sm_io_driver_hdf5::put_varn_expand (int fid,
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
    hsize_t esize, rsize, rsize_old = 0;
    int ndim;
    hid_t dsid = -1, msid = -1;
    hid_t mtype;
    char *bufp = (char *)buf;
    hid_t did;
    hid_t dxplid;
    hsize_t start[E3SM_IO_DRIVER_MAX_RANK], block[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK], mdims[E3SM_IO_DRIVER_MAX_RANK];
    hid_t tid = -1;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    did = fp->dids[vid];

    mtype = mpi_type_to_hdf5_type (itype);
    esize = (hsize_t)H5Tget_size (mtype);
    if (esize <= 0) {
        ERR_OUT("Unknown memory type")
    }

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
            if (fp->recsize < (MPI_Offset)(dims[0])) { fp->recsize = dims[0]; }

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
            throw "bput not supported";

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

    // Call H5Dwrite immediately if blocking or Log VOL is used
    if ((mode == indep) || (mode == coll)) {
        // Call H5DWrite
        for (i = 0; i < nreq; i++) {
            rsize = esize;
            for (j = 0; j < ndim; j++) { rsize *= counts[i][j]; }

            if (rsize) {
                for (j = 0; j < ndim; j++) {
                    start[j] = (hsize_t)starts[i][j];
                    block[j] = (hsize_t)counts[i][j];
                }

                    // Recreate only when size mismatch
                    if (rsize != rsize_old) {
                        if (msid >= 0) H5Sclose (msid);
                        msid = H5Screate_simple (1, &rsize, &rsize);
                        CHECK_HID (msid)

                        rsize_old = rsize;
                    }

                E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)
                herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, start, NULL, one, block);
                CHECK_HERR

                E3SM_IO_TIMER_SWAP (E3SM_IO_TIMER_HDF5_SEL, E3SM_IO_TIMER_HDF5_WR)

                herr = H5Dwrite (did, mtype, msid, dsid, dxplid, bufp);
                CHECK_HERR

                E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_WR)

                bufp += rsize;
            }
        }
    } else {  // Otherwise, queue request in driver
        hid_t dsid_clone;

        // Call H5DWrite
        for (i = 0; i < nreq; i++) {
            rsize = 1;
            for (j = 0; j < ndim; j++) { rsize *= counts[i][j]; }

            if (rsize) {
                for (j = 0; j < ndim; j++) {
                    start[j] = (hsize_t)starts[i][j];
                    block[j] = (hsize_t)counts[i][j];
                }

                msid = H5Screate_simple (1, &rsize, &rsize);
                CHECK_HID (msid)

                E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)
                dsid_clone = H5Scopy (dsid);
                CHECK_HID (dsid_clone)
                herr =
                    H5Sselect_hyperslab (dsid_clone, H5S_SELECT_SET, start, NULL, one, block);
                CHECK_HERR

                E3SM_IO_TIMER_SWAP (E3SM_IO_TIMER_HDF5_SEL, E3SM_IO_TIMER_HDF5_WR)
                herr = fp->register_multidataset (bufp, did, dsid_clone, msid, mtype, 1);
                CHECK_HERR
                E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_WR)

                bufp += rsize * esize;
            }
        }

        msid = -1;
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
    if (msid >= 0) H5Sclose (msid);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

int e3sm_io_driver_hdf5::get_vara (int fid,
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
    hid_t dsid = -1, msid = -1;
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

    // Init memory space
    msid = H5Screate_simple (ndim, hblock, hblock);
    CHECK_HID (msid)

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

    h5_itype = mpi_type_to_hdf5_type(itype);

    // Call H5Dread immediately if blocking or Log VOL is used
    if ((mode == indep) || (mode == coll)) {
        herr = H5Dread (did, h5_itype, msid, dsid, dxplid, buf);
        CHECK_HERR
    }
// #ifdef HDF5_HAVE_DWRITE_MULTI
    else {  // Otherwier, queue request in driver
        herr = fp->register_multidataset (buf, did, dsid, msid, h5_itype, 0);
        CHECK_HERR
        // Prevent freeing of dsid and msid, they will be freed after flush
        dsid = msid = -1;
    }
// #endif
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_RD)

    tid     = H5Dget_type (did);
    getsize = H5Tget_size (tid);
    for (i = 0; i < ndim; i++) { getsize *= hblock[i]; }
    fp->getsize += getsize;

err_out:;
    if (tid >= 0) H5Sclose (tid);
    if (dsid >= 0) H5Sclose (dsid);
    if (msid >= 0) H5Sclose (msid);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}
int e3sm_io_driver_hdf5::get_vars (int fid,
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
    hid_t dsid = -1, msid = -1;
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

    // Init memory space
    msid = H5Screate_simple (ndim, hblock, hblock);
    CHECK_HID (msid)

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

    h5_itype = mpi_type_to_hdf5_type(itype);

    // Call H5Dread immediately if blocking or Log VOL is used
    if ((mode == indep) || (mode == coll)) {
        herr = H5Dread (did, h5_itype, msid, dsid, dxplid, buf);
        CHECK_HERR
    }
// #ifdef HDF5_HAVE_DWRITE_MULTI
    else {  // Otherwier, queue request in driver
        herr = fp->register_multidataset (buf, did, dsid, msid, h5_itype, 0);
        CHECK_HERR
        // Prevent freeing of dsid and msid, they will be freed after flush
        dsid = msid = -1;
    }
// #endif
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_RD)

    tid     = H5Dget_type (did);
    getsize = H5Tget_size (tid);
    for (i = 0; i < ndim; i++) { getsize *= hblock[i]; }
    fp->getsize += getsize;

err_out:;
    if (tid >= 0) H5Sclose (tid);
    if (dsid >= 0) H5Sclose (dsid);
    if (msid >= 0) H5Sclose (msid);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

int e3sm_io_driver_hdf5::get_varn (int fid,
                                   int vid,
                                   MPI_Datatype itype,
                                   int nreq,
                                   MPI_Offset **starts,
                                   MPI_Offset **counts,
                                   void *buf,
                                   e3sm_io_op_mode mode) {
// #ifdef HDF5_HAVE_DWRITE_MULTI
    if (this->merge_varn)
        return this->get_varn_merge (fid, vid, itype, nreq, starts, counts, buf, mode);
    else
// #endif
        return this->get_varn_expand (fid, vid, itype, nreq, starts, counts, buf, mode);
}
int e3sm_io_driver_hdf5::get_varn_expand (int fid,
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
    hsize_t esize, rsize, rsize_old = 0;
    int ndim;
    hid_t dsid = -1, msid = -1;
    hid_t mtype;
    char *bufp = (char *)buf;
    hid_t did;
    hid_t dxplid;
    hsize_t start[E3SM_IO_DRIVER_MAX_RANK], block[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK], mdims[E3SM_IO_DRIVER_MAX_RANK];
    hid_t tid = -1;
    MPI_Offset getsize;
    MPI_Offset **count;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    did = fp->dids[vid];

    mtype = mpi_type_to_hdf5_type (itype);
    esize = (hsize_t)H5Tget_size (mtype);
    if (esize <= 0) {
        ERR_OUT("Unknown memory type")
    }

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

    // Call H5Dread immediately if blocking or Log VOL is used
    if ((mode == indep) || (mode == coll)) {
        // Call H5DWrite
        for (i = 0; i < nreq; i++) {
            rsize = esize;
            for (j = 0; j < ndim; j++) { rsize *= counts[i][j]; }

            if (rsize) {
                for (j = 0; j < ndim; j++) {
                    start[j] = (hsize_t)starts[i][j];
                    block[j] = (hsize_t)counts[i][j];
                }

                // Recreate only when size mismatch
                if (rsize != rsize_old) {
                    if (msid >= 0) H5Sclose (msid);
                    msid = H5Screate_simple (1, &rsize, &rsize);
                    CHECK_HID (msid)

                    rsize_old = rsize;
                }

                E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)
                herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, start, NULL, one, block);
                CHECK_HERR

                E3SM_IO_TIMER_SWAP (E3SM_IO_TIMER_HDF5_SEL, E3SM_IO_TIMER_HDF5_RD)
                herr = H5Dread (did, mtype, msid, dsid, dxplid, bufp);
                CHECK_HERR
                E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_RD)

                bufp += rsize;
            }
        }
    } else {  // Otherwier, queue request in driver
        hid_t dsid_clone;
        
        for (i = 0; i < nreq; i++) {
            rsize = 1;
            for (j = 0; j < ndim; j++) { rsize *= counts[i][j]; }

            if (rsize) {
                for (j = 0; j < ndim; j++) {
                    start[j] = (hsize_t)starts[i][j];
                    block[j] = (hsize_t)counts[i][j];
                }

                msid = H5Screate_simple (1, &rsize, &rsize);
                CHECK_HID (msid)

                E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)
                dsid_clone = H5Scopy (dsid);
                CHECK_HID (dsid_clone)
                herr =
                    H5Sselect_hyperslab (dsid_clone, H5S_SELECT_SET, start, NULL, one, block);
                CHECK_HERR

                E3SM_IO_TIMER_SWAP (E3SM_IO_TIMER_HDF5_SEL, E3SM_IO_TIMER_HDF5_RD)
                herr = fp->register_multidataset (bufp, did, dsid, msid, mtype, 0);
                CHECK_HERR
                E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_RD)

                bufp += rsize * esize;
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
    if (msid >= 0) H5Sclose (msid);

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

