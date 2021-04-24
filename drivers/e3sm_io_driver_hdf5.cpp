#include "e3sm_io.hpp"
#include "e3sm_io_driver_hdf5.hpp"

#include <hdf5.h>
#include <sys/stat.h>

#include "config.h"
#include "e3sm_io_err.hpp"

#ifdef ENABLE_LOGVOL
#include "H5VL_log.h"
#endif

#define CHECK_HERR                                                    \
    {                                                                 \
        if (herr != 0) {                                              \
            printf ("Error at line %d in %s:\n", __LINE__, __FILE__); \
            H5Eprint1 (stdout);                                       \
            DEBUG_ABORT;                                              \
            nerrs++;                                                  \
            goto err_out;                                             \
        }                                                             \
    }

#define CHECK_HID(A)                                                  \
    {                                                                 \
        if (A < 0) {                                                  \
            printf ("Error at line %d in %s:\n", __LINE__, __FILE__); \
            H5Eprint1 (stdout);                                       \
            DEBUG_ABORT;                                              \
            nerrs++;                                                  \
            goto err_out;                                             \
        }                                                             \
    }

static inline hid_t nc_type_to_hdf5_type (nc_type nctype) {
    switch (nctype) {
        case NC_INT:
            return H5T_NATIVE_INT;
        case NC_FLOAT:
            return H5T_NATIVE_FLOAT;
        case NC_DOUBLE:
            return H5T_NATIVE_DOUBLE;
        case NC_CHAR:
            return H5T_NATIVE_CHAR;
        default:
            printf ("Error at line %d in %s: Unknown type %d\n", __LINE__, __FILE__, nctype);
            DEBUG_ABORT
    }

    return -1;
}

static inline hid_t mpi_type_to_hdf5_type (MPI_Datatype mpitype) {
    switch (mpitype) {
        case MPI_INT:
            return H5T_NATIVE_INT;
        case MPI_FLOAT:
            return H5T_NATIVE_FLOAT;
        case MPI_DOUBLE:
            return H5T_NATIVE_DOUBLE;
        case MPI_CHAR:
            return H5T_NATIVE_CHAR;
        default:
            printf ("Error at line %d in %s: Unknown type %d\n", __LINE__, __FILE__, mpitype);
            DEBUG_ABORT
    }

    return -1;
}

e3sm_io_driver_hdf5::e3sm_io_driver_hdf5 () {
    int nerrs   = 0;
    herr_t herr = 0;
    char *env   = NULL;
    int i;

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

#ifdef ENABLE_LOGVOL
    herr = H5Pset_nonblocking (this->dxplid_coll, H5VL_LOG_REQ_BLOCKING);
    CHECK_HERR
    herr = H5Pset_nonblocking (this->dxplid_indep, H5VL_LOG_REQ_BLOCKING);
    CHECK_HERR
    herr = H5Pset_nonblocking (this->dxplid_coll_nb, H5VL_LOG_REQ_NONBLOCKING);
    CHECK_HERR
    herr = H5Pset_nonblocking (this->dxplid_indep_nb, H5VL_LOG_REQ_NONBLOCKING);
    CHECK_HERR

    // Register LOG VOL plugin
    this->log_vlid = H5VLregister_connector (&H5VL_log_g, H5P_DEFAULT);
    CHECK_HID (this->log_vlid)
#endif

    for (i = 0; i < E3SM_IO_DRIVER_MAX_RANK; i++) {
        one[i]  = 1;
        mone[i] = 1;
    }

    this->tsel = this->twrite = this->tread = this->text = 0;

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

err_out:;
    if (nerrs > 0) { throw "HDF5 driver init fail"; }
}
e3sm_io_driver_hdf5::~e3sm_io_driver_hdf5 () {
    int rank;
    int nerrs = 0;
    double tsel_all, twrite_all, text_all;

    if (dxplid_coll >= 0) H5Pclose (dxplid_coll);
    if (dxplid_indep >= 0) H5Pclose (dxplid_indep);
    if (dxplid_coll_nb >= 0) H5Pclose (dxplid_coll_nb);
    if (dxplid_indep_nb >= 0) H5Pclose (dxplid_indep_nb);
    // if (dxplid_nb >= 0) H5Pclose (dxplid_nb);

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

int e3sm_io_driver_hdf5::create (std::string path, MPI_Comm comm, MPI_Info info, int *fid) {
    int nerrs = 0;
    herr_t herr;
    hid_t faplid;
    hdf5_file *fp;

    fp = new hdf5_file ();

    faplid = H5Pcreate (H5P_FILE_ACCESS);
    CHECK_HID (faplid)
    herr = H5Pset_fapl_mpio (faplid, comm, info);
    CHECK_HERR
    herr = H5Pset_coll_metadata_write (faplid, true);
    CHECK_HERR
#ifdef ENABLE_LOGVOL
    if (this->use_logvol) {
        herr = H5Pset_vol (faplid, this->log_vlid, NULL);
        CHECK_HERR
    }
#endif

    fp->id = H5Fcreate (path.c_str (), H5F_ACC_TRUNC, H5P_DEFAULT, faplid);
    CHECK_HID (fp->id)

    *fid = this->files.size ();
    this->files.push_back (fp);

err_out:;

    if (faplid >= 0) { H5Pclose (faplid); }

    return nerrs;
}

int e3sm_io_driver_hdf5::open (std::string path, MPI_Comm comm, MPI_Info info, int *fid) {
    int nerrs = 0;
    herr_t herr;
    hid_t faplid;
    hdf5_file *fp;

    fp = new hdf5_file ();

    faplid = H5Pcreate (H5P_FILE_ACCESS);
    CHECK_HID (faplid)
    herr = H5Pset_fapl_mpio (faplid, comm, info);
    CHECK_HERR
    herr = H5Pset_coll_metadata_write (faplid, true);
    CHECK_HERR
#ifdef ENABLE_LOGVOL
    if (this->use_logvol) {
        herr = H5Pset_vol (faplid, this->log_vlid, NULL);
        CHECK_HERR
    }
#endif

    fp->id = H5Fopen (path.c_str (), H5F_ACC_RDONLY, faplid);
    CHECK_HID (fp->id)

    *fid = this->files.size ();
    this->files.push_back (fp);

err_out:;

    if (faplid >= 0) { H5Pclose (faplid); }

    return nerrs;
}

int e3sm_io_driver_hdf5::close (int fid) {
    int nerrs = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];

    for (auto did : fp->dids) {
        herr = H5Dclose (did);
        CHECK_HERR
    }

    herr = H5Fclose (fp->id);
    CHECK_HERR

    delete fp;

err_out:;
    return nerrs;
}

int e3sm_io_driver_hdf5::inq_file_info (int fid, MPI_Info *info) {
    int nerrs = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    hid_t pid;

    pid = H5Fget_access_plist (fp->id);
    CHECK_HID (pid);
    herr = H5Pget_fapl_mpio (pid, NULL, info);
    CHECK_HERR

err_out:;
    if (pid != -1) H5Pclose (pid);
    return nerrs;
}
int e3sm_io_driver_hdf5::inq_put_size (int fid, MPI_Offset *size) {
    int nerrs = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    ssize_t namelen;
    struct stat file_stat;
    char fname[1024];

    namelen = H5Fget_name (fp->id, fname, 1024);
    CHECK_HID (namelen)

    stat (fname, &file_stat);
    *size = (MPI_Offset) (file_stat.st_size);

err_out:;
    return nerrs;
}
int e3sm_io_driver_hdf5::inq_get_size (int fid, MPI_Offset *size) {
    return this->inq_put_size (fid, size);
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

int e3sm_io_driver_hdf5::def_var (
    int fid, std::string name, MPI_Datatype type, int ndim, int *dimids, int *did) {
    int nerrs = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    int i;
    hid_t h5did;
    hid_t sid    = -1;
    hid_t dcplid = -1;
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK], mdims[E3SM_IO_DRIVER_MAX_RANK];

    dcplid = H5Pcreate (H5P_DATASET_CREATE);
    CHECK_HID (dcplid)

    for (i = 0; i < ndim; i++) { dims[i] = mdims[i] = fp->dsizes[dimids[i]]; }
    if (ndim) {
        if (dims[0] == H5S_UNLIMITED) {
            dims[0] = 1;

            herr = H5Pset_chunk (dcplid, ndim, dims);
            CHECK_HERR
            dims[0] = 0;
        }
    }

    sid = H5Screate_simple (ndim, dims, mdims);
    CHECK_HID (sid);

    h5did = H5Dcreate2 (fp->id, name.c_str (), mpi_type_to_hdf5_type (type), sid, H5P_DEFAULT,
                        dcplid, H5P_DEFAULT);
    CHECK_HID (h5did)

    *did = fp->dids.size ();
    fp->dids.push_back (h5did);

err_out:;
    if (sid != -1) H5Sclose (sid);
    if (dcplid != -1) H5Pclose (dcplid);
    return nerrs;
}

int e3sm_io_driver_hdf5::inq_var (int fid, std::string name, int *did) {
    int nerrs = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    hid_t h5did;

    h5did = H5Dopen2 (fp->id, name.c_str (), H5P_DEFAULT);
    CHECK_HID (h5did)

    *did = fp->dids.size ();
    fp->dids.push_back (h5did);

err_out:;
    return nerrs;
}

int e3sm_io_driver_hdf5::inq_var_off (int fid, int vid, MPI_Offset *off) {
    throw "Function not supported";
    return 1;
}
int e3sm_io_driver_hdf5::def_dim (int fid, std::string name, MPI_Offset size, int *dimid) {
    int nerrs = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    hid_t sid     = -1;
    hid_t aid     = -1;
    hsize_t hsize;
    char aname[128];

    hsize = (hsize_t)size;
    if (hsize == NC_UNLIMITED) size = H5S_UNLIMITED;

    sid = H5Screate (H5S_SCALAR);
    CHECK_HID (sid)

    sprintf (aname, "_DIM_%s", name.c_str ());
    aid = H5Acreate2 (fp->id, aname, H5T_NATIVE_HSIZE, sid, H5P_DEFAULT, H5P_DEFAULT);
    CHECK_HID (aid)

    herr = H5Awrite (aid, H5T_NATIVE_HSIZE, &hsize);

    *dimid = fp->dsizes.size ();
    fp->dsizes.push_back (hsize);

err_out:;
    if (aid != -1) H5Aclose (aid);
    if (sid != -1) H5Sclose (sid);
    return nerrs;
}
int e3sm_io_driver_hdf5::inq_dim (int fid, std::string name, int *dimid) {
    int nerrs = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    hid_t aid     = -1;
    hsize_t hsize;
    char aname[128];

    sprintf (aname, "_DIM_%s", name.c_str ());
    aid = H5Aopen (fp->id, aname, H5P_DEFAULT);
    CHECK_HID (aid)

    herr = H5Aread (aid, H5T_NATIVE_HSIZE, &hsize);

    *dimid = fp->dsizes.size ();
    fp->dsizes.push_back (hsize);

err_out:;
    if (aid != -1) H5Aclose (aid);
    return nerrs;
}

int e3sm_io_driver_hdf5::inq_dimlen (int fid, int dimid, MPI_Offset *size) {
    hdf5_file *fp = this->files[fid];

    *size = (MPI_Offset) (fp->dsizes[dimid]);

    return 0;
}

int e3sm_io_driver_hdf5::enddef (int fid) { return 0; }
int e3sm_io_driver_hdf5::redef (int fid) { return 0; }
int e3sm_io_driver_hdf5::wait (int fid) {
    int nerrs = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];

    herr = H5Fflush (fp->id, H5F_SCOPE_GLOBAL);
    CHECK_HERR

err_out:;
    return nerrs;
}

int e3sm_io_driver_hdf5::put_att (
    int fid, int vid, std::string name, MPI_Datatype type, MPI_Offset size, void *buf) {
    int nerrs = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    hid_t asid = -1, aid = -1;
    hid_t did;
    hsize_t asize;
    htri_t exists;

    asize = (size_t)size;
    asid  = H5Screate_simple (1, &asize, &asize);
    CHECK_HID (asid)

    if (vid == E3SM_IO_GLOBAL_ATTR)
        did = fp->id;
    else
        did = fp->dids[vid];

    exists = H5Aexists (did, name.c_str ());
    if (!exists) {
        aid = H5Acreate2 (did, name.c_str (), mpi_type_to_hdf5_type (type), asid, H5P_DEFAULT,
                          H5P_DEFAULT);
    } else {
        aid = H5Aopen (did, name.c_str (), H5P_DEFAULT);
    }
    CHECK_HID (aid)

    herr = H5Awrite (aid, mpi_type_to_hdf5_type (type), buf);
    CHECK_HERR

err_out:;
    if (asid >= 0) H5Sclose (asid);
    if (aid >= 0) H5Aclose (aid);
    return nerrs;
}

int e3sm_io_driver_hdf5::get_att (int fid, int vid, std::string name, void *buf) {
    int nerrs = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    hid_t asid = -1, aid = -1;
    hid_t did;
    hid_t tid;
    hsize_t asize;
    htri_t exists;

    if (vid == E3SM_IO_GLOBAL_ATTR)
        did = fp->id;
    else
        did = fp->dids[vid];

    aid = H5Aopen (did, name.c_str (), H5P_DEFAULT);
    CHECK_HID (aid)
    tid = H5Aget_type (aid);
    CHECK_HID (tid)
    herr = H5Aread (aid, tid, buf);
    CHECK_HERR

err_out:;
    if (asid >= 0) H5Sclose (asid);
    if (aid >= 0) H5Aclose (aid);
    if (tid >= 0) H5Tclose (tid);
    return nerrs;
}

int e3sm_io_driver_hdf5::put_vara (int fid,
                                   int vid,
                                   MPI_Datatype type,
                                   MPI_Offset *start,
                                   MPI_Offset *count,
                                   void *buf,
                                   e3sm_io_op_mode mode) {
    int nerrs = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    int i;
    int ndim = -1;
    double ts, te;
    hid_t dsid = -1, msid = -1;
    hid_t did;
    hid_t dxplid;
    hsize_t hstart[E3SM_IO_DRIVER_MAX_RANK], hblock[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK];

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
    ts = MPI_Wtime ();
    if (dims[0] < hstart[0] + hblock[0]) {
        dims[0] = hstart[0] + hblock[0];
        if (fp->recsize < dims[0]) { fp->recsize = dims[0]; }

        H5Sclose (dsid);
        herr = H5Dset_extent (did, dims);
        CHECK_HERR
        dsid = H5Dget_space (did);
        CHECK_HID (dsid)
    }
    this->text = MPI_Wtime () - ts;

    // Init memory space
#ifdef ENABLE_LOGVOL
    if (this->use_logvol) {
        msid = H5S_CONTIG;
    } else
#endif
    {
        msid = H5Screate_simple (ndim, hblock, hblock);
        CHECK_HID (msid)
    }

    ts   = MPI_Wtime ();
    herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, hstart, NULL, this->one, hblock);
    CHECK_HERR
    te = MPI_Wtime ();
    this->tsel += te - ts;

    switch (mode) {
        case indep: {
            dxplid = this->dxplid_indep;
            break;
        }
        case coll: {
            dxplid = this->dxplid_coll;
            break;
        }
        case nb: {
            dxplid = this->dxplid_indep_nb;
            break;
        }
        case nbe: {
            if (this->use_logvol) {
                dxplid = this->dxplid_indep;
            } else {
                throw "bput not supported";
            }
            break;
        }
        default:
            throw "Unrecognized mode";
    }

    herr = H5Dwrite (did, mpi_type_to_hdf5_type (type), msid, dsid, dxplid, buf);
    CHECK_HERR
    this->twrite += MPI_Wtime () - te;

err_out:;
    if (dsid >= 0) H5Sclose (dsid);
#ifdef ENABLE_LOGVOL
    if (!(this->use_logvol)) {
#endif
        if (msid >= 0) H5Sclose (msid);
#ifdef ENABLE_LOGVOL
    }
#endif
    return nerrs;
}

int e3sm_io_driver_hdf5::put_vars (int fid,
                                   int vid,
                                   MPI_Datatype type,
                                   MPI_Offset *start,
                                   MPI_Offset *count,
                                   MPI_Offset *stride,
                                   void *buf,
                                   e3sm_io_op_mode mode) {
    int nerrs = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    int i;
    int ndim = -1;
    double ts, te;
    hid_t dsid = -1, msid = -1;
    hid_t did;
    hid_t dxplid;
    hsize_t hstart[E3SM_IO_DRIVER_MAX_RANK], hblock[E3SM_IO_DRIVER_MAX_RANK], hstride[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK];

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
    ts = MPI_Wtime ();
    if (dims[0] < hstart[0] + (hblock[0] - 1) * hstride[0] + 1) {
        dims[0] = hstart[0] + (hblock[0] - 1) * hstride[0] + 1;
        if (fp->recsize < dims[0]) { fp->recsize = dims[0]; }

        H5Sclose (dsid);
        herr = H5Dset_extent (did, dims);
        CHECK_HERR
        dsid = H5Dget_space (did);
        CHECK_HID (dsid)
    }
    this->text = MPI_Wtime () - ts;

    // Init memory space
#ifdef ENABLE_LOGVOL
    if (this->use_logvol) {
        msid = H5S_CONTIG;
    } else
#endif
    {
        msid = H5Screate_simple (ndim, hblock, hblock);
        CHECK_HID (msid)
    }

    ts   = MPI_Wtime ();
    herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, hstart, hstride, this->one, hblock);
    CHECK_HERR
    te = MPI_Wtime ();
    this->tsel += te - ts;

    switch (mode) {
        case indep: {
            dxplid = this->dxplid_indep;
            break;
        }
        case coll: {
            dxplid = this->dxplid_coll;
            break;
        }
        case nb: {
            dxplid = this->dxplid_indep_nb;
            break;
        }
        case nbe: {
            if (this->use_logvol) {
                dxplid = this->dxplid_indep;
            } else {
                throw "bput not supported";
            }
            break;
        }
        default:
            throw "Unrecognized mode";
    }

    herr = H5Dwrite (did, mpi_type_to_hdf5_type (type), msid, dsid, dxplid, buf);
    CHECK_HERR
    this->twrite += MPI_Wtime () - te;

err_out:;
    if (dsid >= 0) H5Sclose (dsid);
#ifdef ENABLE_LOGVOL
    if (!(this->use_logvol)) {
#endif
        if (msid >= 0) H5Sclose (msid);
#ifdef ENABLE_LOGVOL
    }
#endif
    return nerrs;
}
int e3sm_io_driver_hdf5::put_varn (int fid,
                                   int vid,
                                   MPI_Datatype type,
                                   int nreq,
                                   MPI_Offset **starts,
                                   MPI_Offset **counts,
                                   void *buf,
                                   e3sm_io_op_mode mode) {
    int err;
    int nerrs = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    int i, j;
    double ts, te;
    hsize_t esize, rsize, rsize_old = 0;
    int ndim;
    hid_t dsid = -1, msid = -1;
    hid_t mtype;
    char *bufp = (char *)buf;
    hid_t did;
    hid_t dxplid;
    hsize_t **hstarts = NULL, **hcounts = NULL;
    hsize_t start[E3SM_IO_DRIVER_MAX_RANK], block[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK], mdims[E3SM_IO_DRIVER_MAX_RANK];

    did = fp->dids[vid];

    mtype = mpi_type_to_hdf5_type (type);
    esize = (hsize_t)H5Tget_size (mtype);
    CHECK_HID (esize)

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim = H5Sget_simple_extent_dims (dsid, dims, mdims);
    CHECK_HID (ndim)

    // Extend rec dim if needed
    ts = MPI_Wtime ();
    if (ndim && mdims[0] == H5S_UNLIMITED) {
        MPI_Offset max_rec = 0;
        for (i = 0; i < nreq; i++) {
            if (max_rec < starts[i][0] + counts[i][0]) { max_rec = starts[i][0] + counts[i][0]; }
        }
        if (dims[0] < (hsize_t)max_rec) {
            dims[0] = (hsize_t)max_rec;
            if (fp->recsize < dims[0]) { fp->recsize = dims[0]; }

            H5Sclose (dsid);
            herr = H5Dset_extent (did, dims);
            CHECK_HERR
            dsid = H5Dget_space (did);
            CHECK_HID (dsid)
        }
    }
    text += MPI_Wtime () - ts;

    switch (mode) {
        case indep: {
            dxplid = this->dxplid_indep;
            break;
        }
        case coll: {
            dxplid = this->dxplid_coll;
            break;
        }
        case nb: {
            dxplid = this->dxplid_indep_nb;
            break;
        }
        case nbe: {
            if (this->use_logvol) {
                dxplid = this->dxplid_indep;
            } else {
                throw "bput not supported";
            }
            break;
        }
        default:
            throw "Unrecognized mode";
    }

#ifdef ENABLE_LOGVOL
    if (this->use_logvol) { msid = H5S_CONTIG; }
    if ((this->use_logvol) && (this->use_logvol_varn)) {
        ts = MPI_Wtime ();
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
        te = MPI_Wtime ();
        tsel += te - ts;
        herr = H5Dwrite_n (did, mtype, nreq, hstarts, hcounts, dxplid, buf);
        CHECK_HERR
        twrite += MPI_Wtime () - te;
    } else
#endif
    {
        // Call H5DWrite
        for (i = 0; i < nreq; i++) {
            rsize = esize;
            for (j = 0; j < ndim; j++) { rsize *= counts[i][j]; }

            if (rsize) {
                for (j = 0; j < ndim; j++) {
                    start[j] = (hsize_t)starts[i][j];
                    block[j] = (hsize_t)counts[i][j];
                }

#ifdef ENABLE_LOGVOL
                if (!(this->use_logvol)) {
#endif
                    // Recreate only when size mismatch
                    if (rsize != rsize_old) {
                        if (msid >= 0) H5Sclose (msid);
                        msid = H5Screate_simple (1, &rsize, &rsize);
                        CHECK_HID (msid)

                        rsize_old = rsize;
                    }
#ifdef ENABLE_LOGVOL
                }
#endif

                ts   = MPI_Wtime ();
                herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, start, NULL, one, block);
                CHECK_HERR
                te = MPI_Wtime ();
                tsel += te - ts;

                herr = H5Dwrite (did, mtype, msid, dsid, dxplid, bufp);
                CHECK_HERR

                twrite += MPI_Wtime () - te;
                bufp += rsize;
            }
        }
    }

err_out:;
    if (dsid >= 0) H5Sclose (dsid);
#ifdef ENABLE_LOGVOL
    if (!(this->use_logvol)) {
#endif
        if (msid >= 0) H5Sclose (msid);
#ifdef ENABLE_LOGVOL
    } else if (this->use_logvol_varn) {
        free (hstarts[0]);
        free (hstarts);
    }
#endif
    return nerrs;
}

int e3sm_io_driver_hdf5::put_vard (int fid,
                                   int vid,
                                   MPI_Datatype type,
                                   MPI_Offset nelem,
                                   MPI_Datatype ftype,
                                   void *buf,
                                   e3sm_io_op_mode mode) {
    throw "vard API not supported for HDF5 driver";
    return 1;
}

int e3sm_io_driver_hdf5::get_vara (int fid,
                                   int vid,
                                   MPI_Datatype type,
                                   MPI_Offset *start,
                                   MPI_Offset *count,
                                   void *buf,
                                   e3sm_io_op_mode mode) {
    int nerrs = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    int i;
    int ndim = -1;
    double ts, te;
    hid_t dsid = -1, msid = -1;
    hid_t did;
    hid_t dxplid;
    hsize_t hstart[E3SM_IO_DRIVER_MAX_RANK], hblock[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK];

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
    ts = MPI_Wtime ();
    if (dims[0] < hstart[0] + hblock[0]) {
        dims[0] = hstart[0] + hblock[0];
        if (fp->recsize < dims[0]) { fp->recsize = dims[0]; }

        H5Sclose (dsid);
        herr = H5Dset_extent (did, dims);
        CHECK_HERR
        dsid = H5Dget_space (did);
        CHECK_HID (dsid)
    }
    this->text = MPI_Wtime () - ts;

    // Init memory space
#ifdef ENABLE_LOGVOL
    if (this->use_logvol) {
        msid = H5S_CONTIG;
    } else
#endif
    {
        msid = H5Screate_simple (ndim, hblock, hblock);
        CHECK_HID (msid)
    }

    ts   = MPI_Wtime ();
    herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, hstart, NULL, this->one, hblock);
    CHECK_HERR
    te = MPI_Wtime ();
    this->tsel += te - ts;

    switch (mode) {
        case indep: {
            dxplid = this->dxplid_indep;
            break;
        }
        case coll: {
            dxplid = this->dxplid_coll;
            break;
        }
        case nb: {
            dxplid = this->dxplid_indep_nb;
            break;
        }
        case nbe: {
            throw "bput not supported for read";
            break;
        }
        default:
            throw "Unrecognized mode";
    }

    herr = H5Dread (did, mpi_type_to_hdf5_type (type), msid, dsid, dxplid, buf);
    CHECK_HERR
    this->tread += MPI_Wtime () - te;

err_out:;
    if (dsid >= 0) H5Sclose (dsid);
#ifdef ENABLE_LOGVOL
    if (!(this->use_logvol)) {
#endif
        if (msid >= 0) H5Sclose (msid);
#ifdef ENABLE_LOGVOL
    }
#endif
    return nerrs;
}
int e3sm_io_driver_hdf5::get_vars (int fid,
                                   int vid,
                                   MPI_Datatype type,
                                   MPI_Offset *start,
                                   MPI_Offset *count,
                                   MPI_Offset *stride,
                                   void *buf,
                                   e3sm_io_op_mode mode) {
    int nerrs = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    int i;
    int ndim = -1;
    double ts, te;
    hid_t dsid = -1, msid = -1;
    hid_t did;
    hid_t dxplid;
    hsize_t hstart[E3SM_IO_DRIVER_MAX_RANK], hblock[E3SM_IO_DRIVER_MAX_RANK], hstride[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK];

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
    ts = MPI_Wtime ();
    if (dims[0] < hstart[0] + (hblock[0] - 1) * hstride[0] + 1) {
        dims[0] = hstart[0] + (hblock[0] - 1) * hstride[0] + 1;
        if (fp->recsize < dims[0]) { fp->recsize = dims[0]; }

        H5Sclose (dsid);
        herr = H5Dset_extent (did, dims);
        CHECK_HERR
        dsid = H5Dget_space (did);
        CHECK_HID (dsid)
    }
    this->text = MPI_Wtime () - ts;

    // Init memory space
#ifdef ENABLE_LOGVOL
    if (this->use_logvol) {
        msid = H5S_CONTIG;
    } else
#endif
    {
        msid = H5Screate_simple (ndim, hblock, hblock);
        CHECK_HID (msid)
    }

    ts   = MPI_Wtime ();
    herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, hstart, hstride, this->one, hblock);
    CHECK_HERR
    te = MPI_Wtime ();
    this->tsel += te - ts;

    switch (mode) {
        case indep: {
            dxplid = this->dxplid_indep;
            break;
        }
        case coll: {
            dxplid = this->dxplid_coll;
            break;
        }
        case nb: {
            dxplid = this->dxplid_indep_nb;
            break;
        }
        case nbe: {
            throw "bput not supported for read";
            break;
        }
        default:
            throw "Unrecognized mode";
    }

    herr = H5Dread (did, mpi_type_to_hdf5_type (type), msid, dsid, dxplid, buf);
    CHECK_HERR
    this->tread += MPI_Wtime () - te;

err_out:;
    if (dsid >= 0) H5Sclose (dsid);
#ifdef ENABLE_LOGVOL
    if (!(this->use_logvol)) {
#endif
        if (msid >= 0) H5Sclose (msid);
#ifdef ENABLE_LOGVOL
    }
#endif
    return nerrs;
}
int e3sm_io_driver_hdf5::get_varn (int fid,
                                   int vid,
                                   MPI_Datatype type,
                                   int nreq,
                                   MPI_Offset **starts,
                                   MPI_Offset **counts,
                                   void *buf,
                                   e3sm_io_op_mode mode) {
    int err;
    int nerrs = 0;
    herr_t herr;
    hdf5_file *fp = this->files[fid];
    int i, j;
    double ts, te;
    hsize_t esize, rsize, rsize_old = 0;
    int ndim;
    hid_t dsid = -1, msid = -1;
    hid_t mtype;
    char *bufp = (char *)buf;
    hid_t did;
    hid_t dxplid;
    hsize_t **hstarts = NULL, **hcounts = NULL;
    hsize_t start[E3SM_IO_DRIVER_MAX_RANK], block[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK], mdims[E3SM_IO_DRIVER_MAX_RANK];

    did = fp->dids[vid];

    mtype = mpi_type_to_hdf5_type (type);
    esize = (hsize_t)H5Tget_size (mtype);
    CHECK_HID (esize)

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim = H5Sget_simple_extent_dims (dsid, dims, mdims);
    CHECK_HID (ndim)

    // Extend rec dim if needed
    ts = MPI_Wtime ();
    if (ndim && mdims[0] == H5S_UNLIMITED) {
        MPI_Offset max_rec = 0;
        for (i = 0; i < nreq; i++) {
            if (max_rec < starts[i][0] + counts[i][0]) { max_rec = starts[i][0] + counts[i][0]; }
        }
        if (dims[0] < (hsize_t)max_rec) {
            dims[0] = (hsize_t)max_rec;
            if (fp->recsize < dims[0]) { fp->recsize = dims[0]; }

            H5Sclose (dsid);
            herr = H5Dset_extent (did, dims);
            CHECK_HERR
            dsid = H5Dget_space (did);
            CHECK_HID (dsid)
        }
    }
    text += MPI_Wtime () - ts;

    switch (mode) {
        case indep: {
            dxplid = this->dxplid_indep;
            break;
        }
        case coll: {
            dxplid = this->dxplid_coll;
            break;
        }
        case nb: {
            dxplid = this->dxplid_indep_nb;
            break;
        }
        case nbe: {
            throw "bput not supported in read";
            break;
        }
        default:
            throw "Unrecognized mode";
    }

#ifdef ENABLE_LOGVOL
    if (this->use_logvol) { msid = H5S_CONTIG; }
    if (0 && (this->use_logvol) && (this->use_logvol_varn)) {
        ts = MPI_Wtime ();
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
        te = MPI_Wtime ();
        tsel += te - ts;
        // herr = H5Dread_n (did, mtype, nreq, hstarts, hcounts, dxplid, buf);
        // CHECK_HERR
        twrite += MPI_Wtime () - te;
    } else
#endif
    {
        // Call H5DWrite
        for (i = 0; i < nreq; i++) {
            rsize = esize;
            for (j = 0; j < ndim; j++) { rsize *= counts[i][j]; }

            if (rsize) {
                for (j = 0; j < ndim; j++) {
                    start[j] = (hsize_t)starts[i][j];
                    block[j] = (hsize_t)counts[i][j];
                }

#ifdef ENABLE_LOGVOL
                if (!(this->use_logvol)) {
#endif
                    // Recreate only when size mismatch
                    if (rsize != rsize_old) {
                        if (msid >= 0) H5Sclose (msid);
                        msid = H5Screate_simple (1, &rsize, &rsize);
                        CHECK_HID (msid)

                        rsize_old = rsize;
                    }
#ifdef ENABLE_LOGVOL
                }
#endif

                ts   = MPI_Wtime ();
                herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, start, NULL, one, block);
                CHECK_HERR
                te = MPI_Wtime ();
                tsel += te - ts;

                herr = H5Dread (did, mtype, msid, dsid, dxplid, bufp);
                CHECK_HERR

                twrite += MPI_Wtime () - te;
                bufp += rsize;
            }
        }
    }

err_out:;
    if (dsid >= 0) H5Sclose (dsid);
#ifdef ENABLE_LOGVOL
    if (!(this->use_logvol)) {
#endif
        if (msid >= 0) H5Sclose (msid);
#ifdef ENABLE_LOGVOL
    } else if (this->use_logvol_varn) {
        free (hstarts[0]);
        free (hstarts);
    }
#endif
    return nerrs;
}

int e3sm_io_driver_hdf5::get_vard (int fid,
                                   int vid,
                                   MPI_Datatype type,
                                   MPI_Offset nelem,
                                   MPI_Datatype ftype,
                                   void *buf,
                                   e3sm_io_op_mode mode) {
    throw "vard API not supported for HDF5 driver";
    return 1;
}
