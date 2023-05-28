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
#include <sys/stat.h>
//
#include <hdf5.h>
#ifdef ENABLE_LOGVOL
#include <H5VL_log.h>
#endif
//
#include <e3sm_io.h>
#include <e3sm_io_err.h>

#include <e3sm_io_driver_hdf5.hpp>
#include <e3sm_io_profile.hpp>

/*----< e3sm_io_type_nc2hdf5() >---------------------------------------------*/
hid_t e3sm_io_type_nc2hdf5(nc_type xtype)
{
    switch(xtype) {
        case NC_BYTE:   return H5T_NATIVE_INT8;
        case NC_UBYTE:  return H5T_NATIVE_UINT8;
        case NC_CHAR:   return H5T_NATIVE_CHAR;
        case NC_SHORT:  return H5T_NATIVE_SHORT;
        case NC_USHORT: return H5T_NATIVE_USHORT;
        case NC_INT:    return H5T_NATIVE_INT;
        case NC_UINT:   return H5T_NATIVE_UINT;
        case NC_FLOAT : return H5T_NATIVE_FLOAT;
        case NC_DOUBLE: return H5T_NATIVE_DOUBLE;
        case NC_INT64:  return H5T_NATIVE_INT64;
        case NC_UINT64: return H5T_NATIVE_UINT64;
        default: return -1;
    }
}

/*----< e3sm_io_type_mpi2hdf5() >--------------------------------------------*/
hid_t e3sm_io_type_mpi2hdf5(MPI_Datatype itype)
{
         if (itype == MPI_DOUBLE)    return H5T_NATIVE_DOUBLE;
    else if (itype == MPI_FLOAT)     return H5T_NATIVE_FLOAT;
    else if (itype == MPI_INT)       return H5T_NATIVE_INT;
    else if (itype == MPI_LONG_LONG) return H5T_NATIVE_INT64;
    else if (itype == MPI_CHAR)      return H5T_NATIVE_CHAR;
    else if (itype == MPI_BYTE)      return H5T_NATIVE_INT8;

    return -1;
}

e3sm_io_driver_hdf5::e3sm_io_driver_hdf5 (e3sm_io_config *cfg) : e3sm_io_driver (cfg) {
    int i, err = 0;
    herr_t herr = 0;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    // Register LOG VOL plugin
    this->log_vlid = -1;
#ifdef ENABLE_LOGVOL
    if (cfg->strategy == log && cfg->env_log == 0) {
        /* Set to use Log VOL connector only if the env HDF5_VOL_CONNECTOR has
         * not been set to use Log VOL yet.
         */
        this->log_vlid = H5VLregister_connector (&H5VL_log_g, H5P_DEFAULT);
        CHECK_HID (this->log_vlid)
    }
#endif

    this->nfixVars = 0;
    this->dxplid_coll = H5Pcreate (H5P_DATASET_XFER);
    CHECK_HID (this->dxplid_coll)
    herr = H5Pset_dxpl_mpio (this->dxplid_coll, H5FD_MPIO_COLLECTIVE);
    CHECK_HERR
    /*
    this->dxplid_coll_nb = H5Pcreate (H5P_DATASET_XFER);
    CHECK_HID (this->dxplid_coll_nb)
    herr = H5Pset_dxpl_mpio (this->dxplid_coll_nb, H5FD_MPIO_COLLECTIVE);
    CHECK_HERR
    */

    this->dxplid_indep = H5Pcreate (H5P_DATASET_XFER);
    CHECK_HID (this->dxplid_indep)
    // this->dxplid_indep_nb = H5Pcreate (H5P_DATASET_XFER);
    // CHECK_HID (this->dxplid_indep_nb)
    // dxplid_nb = H5Pcreate (H5P_DATASET_XFER);
    // CHECK_HID (dxplid_nb)

    for (i = 0; i < E3SM_IO_DRIVER_MAX_RANK; i++) { one[i] = 1; }

    if (cfg->api == hdf5_md) {
#ifdef HDF5_HAVE_MULTI_DATASET_API
        this->use_dwrite_multi = true;
#ifdef HDF5_HAVE_SELECTION_IO
        /* enable collective I/O in H5Dwrite_multi(). H5Pset_selection_io is
         * first introduced in HDF5 1.14.1
         */
        herr = H5Pset_selection_io(this->dxplid_coll,
                                   H5D_SELECTION_IO_MODE_ON);
        CHECK_HERR
        herr = H5Pset_dxpl_mpio_collective_opt(this->dxplid_coll,
                                               H5FD_MPIO_COLLECTIVE_IO);
        CHECK_HERR
        /* give HDF5 32MB space for type conversion */
        herr = H5Pset_buffer(this->dxplid_coll, 33554432, NULL, NULL);
        CHECK_HERR
#endif
#else
        throw "The HDF5 used does not support multi-dataset write";
#endif
    }

    if ((cfg->chunksize != 0) && (cfg->filter != none)) {
        throw "Fitler requries chunking in HDF5";
    }

    /* Note cfg->num_subfiles is used by the log-based VOL method and ignored
     * by others.
     */

err_out:
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    if (err < 0) { throw "HDF5 driver init fail"; }
}

e3sm_io_driver_hdf5::~e3sm_io_driver_hdf5 () {
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    if (this->dxplid_coll >= 0) {
        H5Pclose (this->dxplid_coll);
        this->dxplid_coll = -1;
    }
    if (this->dxplid_indep >= 0) {
        H5Pclose (this->dxplid_indep);
        this->dxplid_indep = -1;
    }
    // if (dxplid_coll_nb >= 0) H5Pclose (dxplid_coll_nb);
    // if (dxplid_indep_nb >= 0) H5Pclose (dxplid_indep_nb);
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

    err = MPI_Comm_dup (comm, &(fp->comm));
    CHECK_MPIERR

    err = MPI_Comm_rank (comm, &(fp->rank));
    CHECK_MPIERR

    faplid = H5Pcreate (H5P_FILE_ACCESS);
    CHECK_HID (faplid)
    herr = H5Pset_libver_bounds(faplid, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
    CHECK_HID (faplid)
    herr = H5Pset_fapl_mpio (faplid, comm, info);
    CHECK_HERR
    herr = H5Pset_coll_metadata_write (faplid, true);
    CHECK_HERR

#ifdef ENABLE_LOGVOL
    if (this->log_vlid >= 0) {
        /* Set to use Log VOL connector only if the env HDF5_VOL_CONNECTOR has
         * not been set to use Log VOL yet.
         */
        if (cfg->env_log_info != NULL) {
            /* use VOL connector info string from env HDF5_VOL_CONNECTOR for
             * case of stacking Log VOL on top of other VOLs
             */
            void *log_info;
            herr = H5VLconnector_str_to_info(cfg->env_log_info, this->log_vlid, &log_info);
            CHECK_HERR
            herr = H5Pset_vol(faplid, this->log_vlid, log_info);
            CHECK_HERR
            herr = H5VLfree_connector_info(this->log_vlid, log_info);
            CHECK_HERR
        }
        else {
            /* HDF5_VOL_CONNECTOR is not set, use native VOL connector
             * See https://github.com/HDFGroup/hdf5/issues/2417
             */
            H5VL_pass_through_info_t passthru_info;
            passthru_info.under_vol_id   = H5VL_NATIVE;
            passthru_info.under_vol_info = NULL;
            herr = H5Pset_vol(faplid, this->log_vlid, &passthru_info);
            CHECK_HERR
        }
    }
#endif

    // Enlarge metadata cache
    mdcc.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    herr         = H5Pget_mdc_config (faplid, &mdcc);
    CHECK_HERR
    mdcc.max_size     = 128 * 1024 * 1024;
    mdcc.min_size     = 32 * 1024 * 1024;
    mdcc.initial_size = 32 * 1024 * 1024;
    // mdcc.evictions_enabled = false;
    mdcc.incr_mode        = H5C_incr__off;
    mdcc.decr_mode        = H5C_decr__off;
    mdcc.set_initial_size = true;
    herr                  = H5Pset_mdc_config (faplid, &mdcc);
    CHECK_HERR

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_OPEN)

    fp->id = H5Fcreate (path.c_str (), H5F_ACC_TRUNC, H5P_DEFAULT, faplid);
    CHECK_HID (fp->id)

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_OPEN)

    /* obtain MPI file info right after file create */
    herr = H5Pget_fapl_mpio(faplid, NULL, &fp->info_used);
    CHECK_HERR

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

    err = MPI_Comm_dup (comm, &(fp->comm));
    CHECK_MPIERR

    err = MPI_Comm_rank (comm, &(fp->rank));
    CHECK_MPIERR

    faplid = H5Pcreate (H5P_FILE_ACCESS);
    CHECK_HID (faplid)
    herr = H5Pset_fapl_mpio (faplid, comm, info);
    CHECK_HERR
    herr = H5Pset_coll_metadata_write (faplid, true);
    CHECK_HERR

#ifdef ENABLE_LOGVOL
    if (this->log_vlid >= 0) {
        /* Set to use Log VOL connector only if the env HDF5_VOL_CONNECTOR has
         * not been set to use Log VOL yet.
         */
        if (cfg->env_log_info != NULL) {
            /* use VOL connector info string from env HDF5_VOL_CONNECTOR for
             * case of stacking Log VOL on top of other VOLs
             */
            void *log_info;
            herr = H5VLconnector_str_to_info(cfg->env_log_info, this->log_vlid, &log_info);
            CHECK_HERR
            herr = H5Pset_vol(faplid, this->log_vlid, log_info);
            CHECK_HERR
            herr = H5VLfree_connector_info(this->log_vlid, log_info);
            CHECK_HERR
        }
        else {
            /* HDF5_VOL_CONNECTOR is not set, use native VOL connector
             * See https://github.com/HDFGroup/hdf5/issues/2417
             */
            H5VL_pass_through_info_t passthru_info;
            passthru_info.under_vol_id   = H5VL_NATIVE;
            passthru_info.under_vol_info = NULL;
            herr = H5Pset_vol(faplid, this->log_vlid, &passthru_info);
            CHECK_HERR
        }
    }
#endif

    // Enlarge metadata cache
    mdcc.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    herr         = H5Pget_mdc_config (faplid, &mdcc);
    CHECK_HERR
    mdcc.max_size     = 128 * 1024 * 1024;
    mdcc.min_size     = 32 * 1024 * 1024;
    mdcc.initial_size = 32 * 1024 * 1024;
    // mdcc.evictions_enabled = false;
    mdcc.incr_mode        = H5C_incr__off;
    mdcc.decr_mode        = H5C_decr__off;
    mdcc.set_initial_size = true;
    herr                  = H5Pset_mdc_config (faplid, &mdcc);
    CHECK_HERR

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_OPEN)

    fp->id = H5Fopen (path.c_str (), H5F_ACC_RDONLY, faplid);
    CHECK_HID (fp->id)

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_OPEN)

    /* obtain MPI file info right after file open */
    herr = H5Pget_fapl_mpio(faplid, NULL, &fp->info_used);
    CHECK_HERR

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

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_CLOSE)

    herr = H5Fclose (fp->id);
    CHECK_HERR

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_CLOSE)

    MPI_Info_free(&fp->info_used);
    MPI_Comm_free(&(fp->comm));

    delete fp;

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

int e3sm_io_driver_hdf5::inq_file_info (int fid, MPI_Info *info) {
    MPI_Info_dup(this->files[fid]->info_used, info);
    return 0;
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

    *size = (MPI_Offset) (fp->recsize);

    return 0;
}

int e3sm_io_driver_hdf5::expand_rec_size (int fid, MPI_Offset size) {
    int err = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    hid_t dsid    = -1;           // Dataaset space ID
    int ndim;                     // Dataset #dims
    hsize_t dims[H5S_MAX_RANK];   // Dataspace dsid dimensions
    hsize_t mdims[H5S_MAX_RANK];  // Dataspace dsid maxdimensions

    E3SM_IO_TIMER_START(E3SM_IO_TIMER_HDF5)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_EXT_DIM)

    for (auto did : fp->dids) {
        // Get current dim size
        dsid = H5Dget_space (did);
        CHECK_HID (dsid)
        ndim = H5Sget_simple_extent_dims (dsid, dims, mdims);
        CHECK_HID (ndim)
        H5Sclose (dsid);
        dsid = -1;

        // Extend rec dim if needed
        if (ndim && mdims[0] == H5S_UNLIMITED && dims[0] < (hsize_t)size) {
            dims[0] = size;
            herr    = H5Dset_extent (did, dims);
            CHECK_HERR
        }
    }

err_out:;
    if (dsid >= 0) { H5Sclose (dsid); }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_EXT_DIM)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
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
    bool isRec = false;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    dcplid = H5Pcreate (H5P_DATASET_CREATE);
    CHECK_HID (dcplid)

    herr = H5Pset_fill_value (dcplid, 0, NULL);
    CHECK_HERR
    herr = H5Pset_fill_time (dcplid, H5D_FILL_TIME_NEVER);
    CHECK_HERR
    herr = H5Pset_alloc_time (dcplid, H5D_ALLOC_TIME_DEFAULT);
    CHECK_HERR

    h5_xtype = e3sm_io_type_nc2hdf5 (xtype);

    for (i = 0; i < ndim; i++) { dims[i] = mdims[i] = fp->dsizes[dimids[i]]; }
    if (ndim) {
        if ((cfg->chunksize > 0) || (dims[0] == H5S_UNLIMITED)) {
            csize = H5Tget_size (h5_xtype);
            if (cfg->chunksize > 0) {
                for (i = ndim - 1; i > -1; i--) {
                    if (csize < cfg->chunksize) {
                        cdim[i] = mdims[i];
                        csize *= cdim[i];
                    } else {
                        cdim[i] = 1;
                    }
                }
            } else {
                cdim[0] = 1;
                for (i = 1; i < ndim; i++) { cdim[i] = mdims[i]; }
            }
            // Chunk size along rec dim is always 1
            if (dims[0] == H5S_UNLIMITED) {
                isRec = true;
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

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_DEF_VAR)
    h5did = H5Dcreate2 (fp->id, name.c_str (), h5_xtype, sid, H5P_DEFAULT, dcplid, H5P_DEFAULT);
    CHECK_HID (h5did)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_DEF_VAR)

    *did = fp->dids.size ();
    fp->dids.push_back (h5did);
    fp->dset_isRec.push_back(isRec);

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

int e3sm_io_driver_hdf5::inq_varid (int fid, const char *name, int *did) {
    int err       = 0;
    hdf5_file *fp = this->files[fid];
    hid_t h5did;
    htri_t exist;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    // inq_var is used to check whether a variable exist
    exist = H5Lexists(fp->id, name, H5P_DEFAULT);
    if (exist == false){
        err = -1;
        goto err_out;
    }

    h5did = H5Dopen2 (fp->id, name, H5P_DEFAULT);
    CHECK_HID(h5did)

    *did = fp->dids.size ();
    fp->dids.push_back (h5did);

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

int e3sm_io_driver_hdf5::inq_var (int fid, int varid, char *name, nc_type *xtypep,
                                  int *ndimsp, int *dimids, int *nattsp)
{
    int err=0;
    return err;
}

int e3sm_io_driver_hdf5::inq_var_name (int fid, int vid, char *name) {
    name[0] = '\0';
    printf ("inq_var_name is not yet implemented\n");
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

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_DEF_DIM)

    sprintf (aname, "-DIM_%s", name.c_str ());
    aid = H5Acreate2 (fp->id, aname, H5T_NATIVE_HSIZE, sid, H5P_DEFAULT, H5P_DEFAULT);
    CHECK_HID (aid)

    herr = H5Awrite (aid, H5T_NATIVE_HSIZE, &hsize);
    CHECK_HERR

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_DEF_DIM)

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

    sprintf (aname, "-DIM_%s", name.c_str ());
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

    *size = (MPI_Offset) (fp->dsizes[dimid]);

    return 0;
}

int e3sm_io_driver_hdf5::enddef (int fid) { return 0; }

int e3sm_io_driver_hdf5::redef (int fid) { return 0; }

int e3sm_io_driver_hdf5::wait (int fid) {
    int err = 0;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    hdf5_file *fp = this->files[fid];

    if (this->use_dwrite_multi) {
        err = fp->flush_multidatasets ();
        CHECK_ERR
    }
#ifdef ENABLE_LOGVOL
    else if (this->log_vlid >= 0) {
        herr_t herr = H5Fflush (fp->id, H5F_SCOPE_GLOBAL);
        CHECK_HERR
    }
#endif

    /* reset the number of fixed-size variables post in varn call */
    this->nfixVars = 0;

err_out:;
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

int e3sm_io_driver_hdf5::put_att (
    int fid, int vid, std::string name, nc_type xtype, MPI_Offset size, const void *buf) {
    int err = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    hid_t asid = -1, aid = -1, did, h5_xtype  = -1;
    hsize_t asize, esize;
    htri_t exists;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    asize = (size_t)size;

    if (vid == NC_GLOBAL)
        did = fp->id;
    else
        did = fp->dids[vid];

    if (xtype == NC_CHAR) { /* string attribute */
        asid = H5Screate (H5S_SCALAR);
        CHECK_HID (asid)
        h5_xtype = H5Tcopy (H5T_C_S1);
        CHECK_HID (h5_xtype)
        herr = H5Tset_size (h5_xtype, asize + 1);
        CHECK_HERR
        herr = H5Tset_strpad (h5_xtype, H5T_STR_NULLTERM);
        CHECK_HERR
        esize = 1;
    } else {
        asid = H5Screate_simple (1, &asize, &asize);
        CHECK_HID (asid)
        h5_xtype = e3sm_io_type_nc2hdf5 (xtype);
        esize = H5Tget_size (h5_xtype);
    }

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_PUT_ATT)

    exists = H5Aexists (did, name.c_str ());
    if (!exists) {
        aid = H5Acreate2 (did, name.c_str (), h5_xtype, asid, H5P_DEFAULT, H5P_DEFAULT);
    } else {
        aid = H5Aopen (did, name.c_str (), H5P_DEFAULT);
    }
    CHECK_HID (aid)

    herr = H5Awrite (aid, h5_xtype, buf);
    CHECK_HERR

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_PUT_ATT)

    if (fp->rank == 0)
        this->amount_WR += asize * esize;

err_out:;
    if (asid >= 0) H5Sclose (asid);
    if (aid >= 0) H5Aclose (aid);
    if (xtype == NC_CHAR && h5_xtype >= 0) H5Tclose (h5_xtype); /* string attribute */

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)

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

    /* For string attributes, asize must be 1 as E3SM writes only single string
     * attributes. For other types, attributes are either scalars or 1D arrays.
     */
    asid = H5Aget_space (aid);
    CHECK_HID (asid)
    herr = H5Sget_simple_extent_dims (asid, &asize, NULL);
    CHECK_HERR
    esize = H5Tget_size (tid);
    fp->getsize += asize * esize;

err_out:;
    if (asid >= 0) H5Sclose (asid);
    if (aid >= 0) H5Aclose (aid);
    if (tid >= 0) H5Tclose (tid);

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)

    return err;
}

int e3sm_io_driver_hdf5::inq_att (int fid, int vid, std::string name, MPI_Offset *size){
    int err = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    hid_t asid = -1, aid = -1;
    hid_t did;
    hsize_t asize;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    if (vid == NC_GLOBAL)
        did = fp->id;
    else
        did = fp->dids[vid];

    aid = H5Aopen (did, name.c_str (), H5P_DEFAULT);
    CHECK_HID (aid)

    asid = H5Aget_space (aid);
    CHECK_HID (asid)

    herr = H5Sget_simple_extent_dims (asid, &asize, NULL);
    CHECK_HERR

    *size = (MPI_Offset) asize;

err_out:;
    if (asid >= 0) H5Sclose (asid);
    if (aid >= 0) H5Aclose (aid);

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)

    return err;
}

int e3sm_io_driver_hdf5::put_varl (
    int fid, int vid, MPI_Datatype itype, void *buf, e3sm_io_op_mode mode) {
    int err = 0;

    ERR_OUT ("HDF5 does not support local variables")

err_out:;
    return err;
}

int e3sm_io_driver_hdf5::put_vara (int              fid,
                                   int              vid,
                                   MPI_Datatype     itype,
                                   MPI_Offset      *start,
                                   MPI_Offset      *count,
                                   void            *buf,
                                   e3sm_io_op_mode  mode)
{
    int err = 0, i, ndim;
    herr_t herr;
    hid_t did, dsid;
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK];

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    did = this->files[fid]->dids[vid];

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    /* retrieve the dataset's dimension sizes */
    ndim = H5Sget_simple_extent_dims (dsid, dims, NULL);
    CHECK_HID (ndim)

    if (this->use_dwrite_multi) {
        bool free_start = false;
        bool free_count = false;

        E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)

        if (start == NULL) { /* var API: put the entire variable */
            start = (MPI_Offset*) malloc(2 * ndim * sizeof(MPI_Offset));
            count = start + ndim;
            for (i=0; i<ndim; i++) {
                start[i] = 0;
                count[i] = dims[i];
            }
            free_start = true;
        }
        else if (count == NULL) { /* var1 API: put 1 element */
            count = (MPI_Offset*) malloc(ndim * sizeof(MPI_Offset));
            for (i=0; i<ndim; i++) count[i] = 1;
            free_count = true;
        }

        err = this->post_varn (fid, vid, itype, 1, &start, &count, buf, true);
        CHECK_ERR

        if (free_start) free(start);
        if (free_count) free(count);

        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_SEL)
    }
    else {
        hid_t msid, tid;
        hid_t h5_itype;
        hsize_t hstart[E3SM_IO_DRIVER_MAX_RANK];
        hsize_t hblock[E3SM_IO_DRIVER_MAX_RANK];
        MPI_Offset putsize;

        /* type cast MPI_Offset to hsize_t */
        if (start) {
            if (count) {
                for (i = 0; i < ndim; i++) {
                    hstart[i] = (hsize_t)start[i];
                    hblock[i] = (hsize_t)count[i];
                }
            } else { /* var1 API: put 1 element */
                for (i = 0; i < ndim; i++) {
                    hstart[i] = (hsize_t)start[i];
                    hblock[i] = 1;
                }
            }
        } else { /* var API: put the entire variable */
            for (i = 0; i < ndim; i++) {
                hstart[i] = 0;
                hblock[i] = dims[i];
            }
        }

        /* write buffer is contiguous in memory space */
        msid = H5Screate_simple (ndim, hblock, NULL);
        CHECK_HID (msid)

        h5_itype = e3sm_io_type_mpi2hdf5 (itype);

        E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)
        if (ndim > 0) {
            herr = H5Sselect_hyperslab(dsid, H5S_SELECT_SET, hstart, NULL,
                                       this->one, hblock);
            CHECK_HERR
        }

        E3SM_IO_TIMER_SWAP (E3SM_IO_TIMER_HDF5_SEL, E3SM_IO_TIMER_HDF5_WR)

        /* Call H5Dwrite - use MPI independent I/O to write */
        herr = H5Dwrite (did, h5_itype, msid, dsid, this->dxplid_indep, buf);
        CHECK_HERR

        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_WR)

        tid     = H5Dget_type (did);
        putsize = H5Tget_size (tid);
        for (i = 0; i < ndim; i++) { putsize *= hblock[i]; }
        this->amount_WR += putsize;

        H5Tclose (tid);
        H5Sclose (msid);

    }

err_out:;
    if (dsid >= 0) H5Sclose (dsid);

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)

    return err;
}

int e3sm_io_driver_hdf5::varn_expand(int            fid,
                                     int            vid,
                                     MPI_Datatype   itype,
                                     int            nreq,
                                     MPI_Offset   **starts,
                                     MPI_Offset   **counts,
                                     void          *buf,
                                     bool           isWrite)
{
    int i, j, err = 0, ndim, union_filespaces = 1;
    size_t tsize;
    herr_t herr;
    hid_t dsid = -1, msid = -1, mtype, did, tid = -1;
    hsize_t start[E3SM_IO_DRIVER_MAX_RANK], block[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t rsize, rsize_all = 0;
    hdf5_file *fp = this->files[fid];

    did = fp->dids[vid];
    mtype = e3sm_io_type_mpi2hdf5 (itype);

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim = H5Sget_simple_extent_dims (dsid, dims, NULL);
    CHECK_HID (ndim)

    if (union_filespaces) {
        /* WARNING: H5Sselect_hyperslab's H5S_SELECT_OR will flatten and sort
         * the element offsets of the union into an increasing order. The
         * filespace created may not match with the user's buffer.
         */

        H5S_seloper_t op = H5S_SELECT_SET;

        E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)

        // set filespace
#ifdef CHECK_INCREASING_ORDER
static int print_once=0;
long long off_start, off_end, prev_end, dim_len; prev_end=-1;
#endif
        for (i = 0; i < nreq; i++) {
            rsize = 1;
            for (j = 0; j < ndim; j++) { rsize *= counts[i][j]; }
            if (rsize == 0) continue;
            rsize_all += rsize;

#ifdef CHECK_INCREASING_ORDER
off_start=0; off_end=0; dim_len=1;
for (j=ndim-1; j>=0; j--) {
off_start += starts[i][j] * dim_len;
off_end += (starts[i][j] + counts[i][j]-1) * dim_len;
dim_len *= dims[j];
}
if (ndim==2 && prev_end >= 0 && prev_end > off_start) printf("Error: %s line %d: req=%d of nreq=%d starts[%d][0]=%lld %lld counts[%d][0]=%lld %lld prev_end=%lld > off_start=%lld off_end=%lld start[%d][0]=%lld %lld counts[%d][0]=%lld %lld (rsize=%ld, dims[]=%ld %ld, ndim=%d)\n",__FILE__,__LINE__,i,nreq,i-1,starts[i-1][0],starts[i-1][1],i-1,counts[i-1][0],counts[i-1][1],prev_end,off_start,off_end,i,starts[i][0],starts[i][1],i,counts[i][0],counts[i][1],rsize,dims[0],dims[1],ndim);
else if (prev_end >= 0 && prev_end > off_start) printf("Error: %s line %d: req=%d of nreq=%d starts[%d][0]=%lld counts[%d][0]=%lld prev_end=%lld > off_start=%lld off_end=%lld start[%d][0]=%lld counts[%d][0]=%lld (rsize=%ld, dims[]=%ld %ld, ndim=%d)\n",__FILE__,__LINE__,i,nreq,i-1,starts[i-1][0],i-1,counts[i-1][0],prev_end,off_start,off_end,i,starts[i][0],i,counts[i][0],rsize,dims[0],dims[1],ndim);
prev_end=off_end+1;
#endif

            /* type cast from MPI_Offset to hsize_t */
            for (j = 0; j < ndim; j++) {
                start[j] = (hsize_t)starts[i][j];
                block[j] = (hsize_t)counts[i][j];
            }

            /* union all nreq hyperslabs */
            herr = H5Sselect_hyperslab (dsid, op, start, NULL, this->one, block);
            CHECK_HERR
            op = H5S_SELECT_OR;
        }

        // create memory space
        if (rsize_all > 0) {
            msid = H5Screate_simple (1, &rsize_all, NULL);
            CHECK_HID (msid)
        }
        else { /* this process has nothing to write, but must participate the
                * collective write call */
          /* create a zero-sized memory space */
           msid = H5Screate(H5S_NULL);
           CHECK_HID (msid)
           /* set the selection of dataset's file space to zero size */
           herr = H5Sselect_none (dsid);
           CHECK_ERR
        }
        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_SEL)

        /* collective write/read, one per variable */
        if (isWrite) {
            E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_WR)
            herr = H5Dwrite (did, mtype, msid, dsid, this->dxplid_coll, buf);
            CHECK_HERR
            E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_WR)
        }
        else {
            E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_RD)
            herr = H5Dread (did, mtype, msid, dsid, this->dxplid_coll, buf);
            CHECK_HERR
            E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_RD)
        }
    }
    else {
        char *bufp = (char *)buf;
        hsize_t rsize_old = 0;
        hsize_t esize = (hsize_t)H5Tget_size (mtype);
        if (esize <= 0) { ERR_OUT ("Unknown memory type") }

        // select filespace, set memspace, and independent write
        for (i = 0; i < nreq; i++) {
            rsize = 1;
            for (j = 0; j < ndim; j++) { rsize *= counts[i][j]; }
            if (rsize == 0) continue;

            E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)
            rsize_all += rsize;

            /* type cast from MPI_Offset to hsize_t */
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

            herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, start, NULL, this->one, block);
            CHECK_HERR
            E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_SEL)

            /* independent write/read, as nreq can be different among processes */
            if (isWrite) {
                E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_WR)
                herr = H5Dwrite (did, mtype, msid, dsid, this->dxplid_indep, bufp);
                CHECK_HERR
                E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_WR)
            }
            else {
                E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_RD)
                herr = H5Dread (did, mtype, msid, dsid, this->dxplid_indep, bufp);
                CHECK_HERR
                E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_RD)
            }

            bufp += rsize * esize;
        }
    }

    tid   = H5Dget_type (did);
    tsize = H5Tget_size (tid);
    if (isWrite)
        this->amount_WR += tsize * rsize_all;
    else
        this->amount_RD += tsize * rsize_all;

err_out:;
    if (tid >= 0) H5Tclose (tid);
    if (dsid >= 0) H5Sclose (dsid);
    if (msid >= 0) H5Sclose (msid);

    return err;
}

int e3sm_io_driver_hdf5::put_varn (int               fid,
                                   int               vid,
                                   MPI_Datatype      itype,
                                   int               nreq,
                                   MPI_Offset      **starts,
                                   MPI_Offset      **counts,
                                   void             *buf,
                                   e3sm_io_op_mode   mode)
{
    int err;

    E3SM_IO_TIMER_START(E3SM_IO_TIMER_HDF5)

    if (this->use_dwrite_multi)
        err = this->post_varn(fid, vid, itype, nreq, starts, counts, buf, true);
    else
        err = varn_expand(fid, vid, itype, nreq, starts, counts, buf, true);

    E3SM_IO_TIMER_STOP(E3SM_IO_TIMER_HDF5)

    return err;
}

int e3sm_io_driver_hdf5::get_vara (int              fid,
                                   int              vid,
                                   MPI_Datatype     itype,
                                   MPI_Offset      *start,
                                   MPI_Offset      *count,
                                   void            *buf,
                                   e3sm_io_op_mode  mode)
{
    int err = 0, i, ndim;
    herr_t herr;
    hid_t did, dsid;
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK];

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    did = this->files[fid]->dids[vid];

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    /* retrieve the dataset's dimension sizes */
    ndim = H5Sget_simple_extent_dims (dsid, dims, NULL);
    CHECK_HID (ndim)

    if (this->use_dwrite_multi) {
        bool free_start = false;
        bool free_count = false;

        E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)

        if (start == NULL) { /* var API: get the entire variable */
            start = (MPI_Offset*) malloc(2 * ndim * sizeof(MPI_Offset));
            count = start + ndim;
            for (i=0; i<ndim; i++) {
                start[i] = 0;
                count[i] = dims[i];
            }
            free_start = true;
        }
        else if (count == NULL) { /* var1 API: get 1 element */
            count = (MPI_Offset*) malloc(ndim * sizeof(MPI_Offset));
            for (i=0; i<ndim; i++) count[i] = 1;
            free_count = true;
        }

        err = this->post_varn (fid, vid, itype, 1, &start, &count, buf, false);
        CHECK_ERR

        if (free_start) free(start);
        if (free_count) free(count);

        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_SEL)
    }
    else {
        hid_t msid, tid;
        hid_t h5_itype;
        hsize_t hstart[E3SM_IO_DRIVER_MAX_RANK];
        hsize_t hblock[E3SM_IO_DRIVER_MAX_RANK];
        MPI_Offset getsize;

        E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)

        /* type cast MPI_Offset to hsize_t */
        if (start) {
            if (count) {
                for (i = 0; i < ndim; i++) {
                    hstart[i] = (hsize_t)start[i];
                    hblock[i] = (hsize_t)count[i];
                }
            } else { /* var1 API: put 1 element */
                for (i = 0; i < ndim; i++) {
                    hstart[i] = (hsize_t)start[i];
                    hblock[i] = 1;
                }
            }
        } else { /* var API: put the entire variable */
            for (i = 0; i < ndim; i++) {
                hstart[i] = 0;
                hblock[i] = dims[i];
            }
        }

        /* write buffer is contiguous in memory space */
        msid = H5Screate_simple (ndim, hblock, NULL);
        CHECK_HID (msid)

        h5_itype = e3sm_io_type_mpi2hdf5 (itype);

        if (ndim > 0) {
            herr = H5Sselect_hyperslab(dsid, H5S_SELECT_SET, hstart, NULL,
                                       this->one, hblock);
            CHECK_HERR
        }

        E3SM_IO_TIMER_SWAP (E3SM_IO_TIMER_HDF5_SEL, E3SM_IO_TIMER_HDF5_RD)

        /* Call H5Dread - use MPI independent I/O to write */
        herr = H5Dread (did, h5_itype, msid, dsid, this->dxplid_indep, buf);
        CHECK_HERR

        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_RD)

        tid     = H5Dget_type (did);
        getsize = H5Tget_size (tid);
        for (i = 0; i < ndim; i++) { getsize *= hblock[i]; }
        this->amount_RD += getsize;

        H5Tclose (tid);
        H5Sclose (msid);
    }

err_out:;
    if (dsid >= 0) H5Sclose (dsid);

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)

    return err;
}

int e3sm_io_driver_hdf5::get_varn (int               fid,
                                   int               vid,
                                   MPI_Datatype      itype,
                                   int               nreq,
                                   MPI_Offset      **starts,
                                   MPI_Offset      **counts,
                                   void             *buf,
                                   e3sm_io_op_mode   mode)
{
    int err;

    E3SM_IO_TIMER_START(E3SM_IO_TIMER_HDF5)

    if (this->use_dwrite_multi)
        err = this->post_varn(fid, vid, itype, nreq, starts, counts, buf, false);
    else
        err = varn_expand(fid, vid, itype, nreq, starts, counts, buf, false);

    E3SM_IO_TIMER_STOP(E3SM_IO_TIMER_HDF5)

    return err;
}

