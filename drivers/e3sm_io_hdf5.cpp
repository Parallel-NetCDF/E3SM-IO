#include "e3sm_io_hdf5.hpp"

#include <string.h>

// hid_t dxplid_nb  = -1;
hid_t dxplid_coll     = -1;
hid_t dxplid_indep    = -1;
hid_t dxplid_coll_nb  = -1;
hid_t dxplid_indep_nb = -1;
hsize_t one[H5S_MAX_RANK];
MPI_Offset mone[H5S_MAX_RANK];
static int f_ndim;
hsize_t f_dims[1048576];
hid_t f_dids[1048576];
static int f_nd;

#ifdef ENABLE_LOGVOL
hid_t log_vlid = -1;
#endif

hid_t nc_type_to_hdf5_type (nc_type nctype) {
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
    }

    return -1;
}

hid_t mpi_type_to_hdf5_type (MPI_Datatype mpitype) {
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
    }

    return -1;
}

double tsel, twrite, tread, text;
int hdf5_wrap_init () {
    herr_t herr = 0;
    int i;

    dxplid_coll = H5Pcreate (H5P_DATASET_XFER);
    CHECK_HID (dxplid_coll)
    herr = H5Pset_dxpl_mpio (dxplid_coll, H5FD_MPIO_COLLECTIVE);
    CHECK_HERR
    dxplid_coll_nb = H5Pcreate (H5P_DATASET_XFER);
    CHECK_HID (dxplid_coll_nb)
    herr = H5Pset_dxpl_mpio (dxplid_coll_nb, H5FD_MPIO_COLLECTIVE);
    CHECK_HERR

    dxplid_indep = H5Pcreate (H5P_DATASET_XFER);
    CHECK_HID (dxplid_indep)
    dxplid_indep_nb = H5Pcreate (H5P_DATASET_XFER);
    CHECK_HID (dxplid_indep_nb)
    // dxplid_nb = H5Pcreate (H5P_DATASET_XFER);
    // CHECK_HID (dxplid_nb)

#ifdef ENABLE_LOGVOL
    herr = H5Pset_nonblocking (dxplid_coll, H5VL_LOG_REQ_BLOCKING);
    CHECK_HERR
    herr = H5Pset_nonblocking (dxplid_indep, H5VL_LOG_REQ_BLOCKING);
    CHECK_HERR
    herr = H5Pset_nonblocking (dxplid_coll_nb, H5VL_LOG_REQ_NONBLOCKING);
    CHECK_HERR
    herr = H5Pset_nonblocking (dxplid_indep_nb, H5VL_LOG_REQ_NONBLOCKING);
    CHECK_HERR

    // Register LOG VOL plugin
    log_vlid = H5VLregister_connector (&H5VL_log_g, H5P_DEFAULT);
    CHECK_HID (log_vlid)
#endif

    for (i = 0; i < H5S_MAX_RANK; i++) {
        one[i]  = 1;
        mone[i] = 1;
    }

    f_ndim = 0;
    f_nd   = 0;

    tsel = twrite = tread = text = 0;

fn_exit:;
    return (int)herr;
}

void hdf5_wrap_finalize () {
    int rank;
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

int hdf5_put_vara (
    int vid, hid_t mtype, hid_t dxplid, MPI_Offset *mstart, MPI_Offset *mcount, void *buf) {
    herr_t herr = 0;
    int i;
    int ndim = -1;
    double ts, te;
    hid_t dsid = -1, msid = -1;
    hid_t did;
    hsize_t start[H5S_MAX_RANK], block[H5S_MAX_RANK];
    hsize_t dims[H5S_MAX_RANK];

    did = f_dids[vid];

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim = H5Sget_simple_extent_dims (dsid, dims, NULL);
    CHECK_HID (ndim)

    for (i = 0; i < ndim; i++) {
        start[i] = (hsize_t)mstart[i];
        block[i] = (hsize_t)mcount[i];
    }

    // Extend rec dim
    ts = MPI_Wtime ();
    if (dims[0] < start[0] + block[0]) {
        dims[0] = start[0] + block[0];

        H5Sclose (dsid);
        herr = H5Dset_extent (did, dims);
        CHECK_HERR
        dsid = H5Dget_space (did);
        CHECK_HID (dsid)
    }
    text = MPI_Wtime () - ts;

#ifndef ENABLE_LOGVOL
    msid = H5Screate_simple (ndim, block, block);
    CHECK_HID (msid)
#endif

    ts   = MPI_Wtime ();
    herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, start, NULL, one, block);
    CHECK_HERR
    te = MPI_Wtime ();
    tsel += te - ts;
#ifdef ENABLE_LOGVOL
    herr = H5Dwrite (did, mtype, H5S_CONTIG, dsid, dxplid, buf);
    CHECK_HERR
#else
    herr = H5Dwrite (did, mtype, msid, dsid, dxplid, buf);
    CHECK_HERR
#endif
    twrite += MPI_Wtime () - te;

fn_exit:;
    if (dsid >= 0) H5Sclose (dsid);
#ifndef ENABLE_LOGVOL
    if (msid >= 0) H5Sclose (msid);
#endif
    return (int)herr;
}

int hdf5_put_vars (int vid,
                   hid_t mtype,
                   hid_t dxplid,
                   MPI_Offset *mstart,
                   MPI_Offset *mcount,
                   MPI_Offset *mstride,
                   void *buf) {
    herr_t herr = 0;
    int i;
    int ndim = -1;
    double ts, te;
    hid_t dsid = -1, msid = -1;
    hid_t did;
    hsize_t start[H5S_MAX_RANK], block[H5S_MAX_RANK], stride[H5S_MAX_RANK];
    hsize_t dims[H5S_MAX_RANK];

    did = f_dids[vid];

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim = H5Sget_simple_extent_dims (dsid, dims, NULL);
    CHECK_HID (ndim)

    for (i = 0; i < ndim; i++) {
        start[i]  = (hsize_t)mstart[i];
        block[i]  = (hsize_t)mcount[i];
        stride[i] = (hsize_t)mstride[i];
    }

    // Extend rec dim
    ts = MPI_Wtime ();
    if (dims[0] < start[0] + block[0]) {
        dims[0] = start[0] + (block[0] - 1) * stride[0] + 1;

        H5Sclose (dsid);
        herr = H5Dset_extent (did, dims);
        CHECK_HERR
        dsid = H5Dget_space (did);
        CHECK_HID (dsid)
    }
    text = MPI_Wtime () - ts;

#ifndef ENABLE_LOGVOL
    msid = H5Screate_simple (ndim, block, block);
    CHECK_HID (msid)
#endif

    ts   = MPI_Wtime ();
    herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, start, stride, one, block);
    CHECK_HERR
    te = MPI_Wtime ();
    tsel += te - ts;
#ifdef ENABLE_LOGVOL
    herr = H5Dwrite (did, mtype, H5S_CONTIG, dsid, dxplid, buf);
    CHECK_HERR
#else
    herr = H5Dwrite (did, mtype, msid, dsid, dxplid, buf);
    CHECK_HERR
#endif
    twrite += MPI_Wtime () - te;

fn_exit:;
    if (dsid >= 0) H5Sclose (dsid);
#ifndef ENABLE_LOGVOL
    if (msid >= 0) H5Sclose (msid);
#endif
    return (int)herr;
}

int hdf5_put_var1 (int vid, hid_t mtype, hid_t dxplid, MPI_Offset *mstart, void *buf) {
    return hdf5_put_vara (vid, mtype, dxplid, mstart, mone, buf);
}

#ifdef ENABLE_LOGVOL
#define USE_LOGVOL_VARN
#endif

#ifdef USE_LOGVOL_VARN
int hdf5_put_varn (int vid,
                   MPI_Datatype mpitype,
                   hid_t dxplid,
                   int cnt,
                   MPI_Offset **mstarts,
                   MPI_Offset **mcounts,
                   void *buf) {
    int err;
    herr_t herr = 0;
    int i, j;
    double ts, te;
    hsize_t esize, rsize, rsize_old = 0;
    int ndim;
    hid_t dsid = -1, msid = -1;
    hid_t mtype;
    char *bufp = (char*)buf;
    hid_t did;
    hsize_t **starts = NULL, **counts;
    hsize_t dims[H5S_MAX_RANK], mdims[H5S_MAX_RANK];

    did = f_dids[vid];

    mtype = mpi_type_to_hdf5_type (mpitype);
    CHECK_HID (mtype)

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim = H5Sget_simple_extent_dims (dsid, dims, mdims);
    CHECK_HID (ndim)

    // Extend rec dim if needed
    ts = MPI_Wtime ();
    if (ndim && mdims[0] == H5S_UNLIMITED) {
        MPI_Offset max_rec = 0;
        for (i = 0; i < cnt; i++) {
            if (max_rec < mstarts[i][0] + mcounts[i][0]) {
                max_rec = mstarts[i][0] + mcounts[i][0];
            }
        }
        if (dims[0] < (hsize_t)max_rec) {
            dims[0] = (hsize_t)max_rec;

            H5Sclose (dsid);
            herr = H5Dset_extent (did, dims);
            CHECK_HERR
            dsid = H5Dget_space (did);
            CHECK_HID (dsid)
        }
    }
    text += MPI_Wtime () - ts;

    // Call H5Dwriten
    te = MPI_Wtime ();

    // Convert starts and counts;
    starts    = (hsize_t **)malloc (sizeof (hsize_t *) * cnt * 2);
    counts    = starts + cnt;
    starts[0] = (hsize_t *)malloc (sizeof (hsize_t ) * cnt * ndim * 2);
    counts[0] = starts[0] + cnt * ndim;
    for (i = 1; i < cnt; i++) {
        starts[i] = starts[i - 1] + ndim;
        counts[i] = counts[i - 1] + ndim;
    }
    for (i = 1; i < cnt; i++) {
        for (j = 0; j < ndim; j++) {
            starts[i][j] = (hsize_t)mstarts[i][j];
            counts[i][j] = (hsize_t)mcounts[i][j];
        }
    }
    ts = MPI_Wtime ();
    tsel += ts - te;

    herr = H5Dwrite_n (did, mtype, cnt, starts, counts, dxplid, buf);
    CHECK_HERR
    twrite += MPI_Wtime () - ts;

fn_exit:;
    if (dsid >= 0) H5Sclose (dsid);

    if (starts) {
        free (starts[0]);
        free (starts);
    }

    return (int)herr;
}
#else
int hdf5_put_varn (int vid,
                   MPI_Datatype mpitype,
                   hid_t dxplid,
                   int cnt,
                   MPI_Offset **mstarts,
                   MPI_Offset **mcounts,
                   void *buf) {
    int err;
    herr_t herr = 0;
    int i, j;
    double ts, te;
    hsize_t esize, rsize, rsize_old = 0;
    int ndim;
    hid_t dsid = -1, msid = -1;
    hid_t mtype;
    char *bufp = buf;
    hid_t did;
    hsize_t start[H5S_MAX_RANK], block[H5S_MAX_RANK];
    hsize_t dims[H5S_MAX_RANK], mdims[H5S_MAX_RANK];

    did = f_dids[vid];

    mtype = mpi_type_to_hdf5_type (mpitype);
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
        for (i = 0; i < cnt; i++) {
            if (max_rec < mstarts[i][0] + mcounts[i][0]) {
                max_rec = mstarts[i][0] + mcounts[i][0];
            }
        }
        if (dims[0] < (hsize_t)max_rec) {
            dims[0] = (hsize_t)max_rec;

            H5Sclose (dsid);
            herr = H5Dset_extent (did, dims);
            CHECK_HERR
            dsid = H5Dget_space (did);
            CHECK_HID (dsid)
        }
    }
    text += MPI_Wtime () - ts;

    // Call H5DWrite
    for (i = 0; i < cnt; i++) {
        rsize = esize;
        for (j = 0; j < ndim; j++) { rsize *= mcounts[i][j]; }

        if (rsize) {
            // err = hdf5_put_vara (vid, mtype, dxplid, mstarts[i], mcounts[i], bufp);
            // CHECK_ERR

            for (j = 0; j < ndim; j++) {
                start[j] = (hsize_t)mstarts[i][j];
                block[j] = (hsize_t)mcounts[i][j];
            }

#ifndef ENABLE_LOGVOL
            // Recreate only when size mismatch
            if (rsize != rsize_old) {
                if (msid >= 0) H5Sclose (msid);
                msid = H5Screate_simple (1, &rsize, &rsize);
                CHECK_HID (msid)

                rsize_old = rsize;
            }
#endif

            ts = MPI_Wtime ();
            herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, start, NULL, one, block);
            CHECK_HERR
            te = MPI_Wtime ();
            tsel += te - ts;
#ifdef ENABLE_LOGVOL
            herr = H5Dwrite (did, mtype, H5S_CONTIG, dsid, dxplid, bufp);
            CHECK_HERR
#else
            herr = H5Dwrite (did, mtype, msid, dsid, dxplid, bufp);
            CHECK_HERR
#endif
            twrite += MPI_Wtime () - te;
            bufp += rsize;
        }
    }

fn_exit:;
    if (dsid >= 0) H5Sclose (dsid);
#ifndef ENABLE_LOGVOL
    if (msid >= 0) H5Sclose (msid);
#endif
    return (int)herr;
}
#endif

int hdf5_get_vara (
    int vid, hid_t mtype, hid_t dxplid, MPI_Offset *mstart, MPI_Offset *mcount, void *buf) {
    herr_t herr = 0;
    int i;
    double ts, te;
    int ndim   = -1;
    hid_t dsid = -1, msid = -1;
    hid_t did;
    hsize_t start[H5S_MAX_RANK], block[H5S_MAX_RANK];

    did = f_dids[vid];

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim = H5Sget_simple_extent_ndims (dsid);
    CHECK_HID (ndim)

    for (i = 0; i < ndim; i++) {
        start[i] = (hsize_t)mstart[i];
        block[i] = (hsize_t)mcount[i];
    }

#ifndef ENABLE_LOGVOL
    msid = H5Screate_simple (ndim, block, block);
    CHECK_HID (msid)
#endif

    ts   = MPI_Wtime ();
    herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, start, NULL, one, block);
    CHECK_HERR
    te = MPI_Wtime ();
    tsel += te - ts;
#ifdef ENABLE_LOGVOL
    herr = H5Dread (did, mtype, H5S_CONTIG, dsid, dxplid, buf);
    CHECK_HERR
#else
    herr = H5Dread (did, mtype, msid, dsid, dxplid, buf);
    CHECK_HERR
#endif
    tread += MPI_Wtime () - te;

fn_exit:;
    if (dsid >= 0) H5Sclose (dsid);
#ifndef ENABLE_LOGVOL
    if (msid >= 0) H5Sclose (msid);
#endif
    return (int)herr;
}

int hdf5_get_vars (int vid,
                   hid_t mtype,
                   hid_t dxplid,
                   MPI_Offset *mstart,
                   MPI_Offset *mcount,
                   MPI_Offset *mstride,
                   void *buf) {
    herr_t herr = 0;
    int i;
    double ts, te;
    int ndim   = -1;
    hid_t dsid = -1, msid = -1;
    hid_t did;
    hsize_t start[H5S_MAX_RANK], block[H5S_MAX_RANK], stride[H5S_MAX_RANK];

    did = f_dids[vid];

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim = H5Sget_simple_extent_ndims (dsid);
    CHECK_HID (ndim)

    for (i = 0; i < ndim; i++) {
        start[i]  = (hsize_t)mstart[i];
        block[i]  = (hsize_t)mcount[i];
        stride[i] = (hsize_t)mstride[i];
    }

#ifndef ENABLE_LOGVOL
    msid = H5Screate_simple (ndim, block, block);
    CHECK_HID (msid)
#endif

    ts   = MPI_Wtime ();
    herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, start, stride, one, block);
    CHECK_HERR
    te = MPI_Wtime ();
    tsel += te - ts;
#ifdef ENABLE_LOGVOL
    herr = H5Dread (did, mtype, H5S_CONTIG, dsid, dxplid, buf);
    CHECK_HERR
#else
    herr = H5Dread (did, mtype, msid, dsid, dxplid, buf);
    CHECK_HERR
#endif
    tread += MPI_Wtime () - te;

fn_exit:;
    if (dsid >= 0) H5Sclose (dsid);
#ifndef ENABLE_LOGVOL
    if (msid >= 0) H5Sclose (msid);
#endif
    return (int)herr;
}

int hdf5_get_var1 (int vid, hid_t mtype, hid_t dxplid, MPI_Offset *mstart, void *buf) {
    return hdf5_get_vara (vid, mtype, dxplid, mstart, mone, buf);
}

int hdf5_get_varn (int vid,
                   MPI_Datatype mpitype,
                   hid_t dxplid,
                   int cnt,
                   MPI_Offset **mstarts,
                   MPI_Offset **mcounts,
                   void *buf) {
    int err;
    herr_t herr = 0;
    int i, j;
    size_t esize, rsize;
    int ndim;
    hid_t dsid = -1;
    hid_t mtype;
    hid_t did;
    char *bufp = (char*)buf;

    did = f_dids[vid];

    mtype = mpi_type_to_hdf5_type (mpitype);

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim = H5Sget_simple_extent_ndims (dsid);
    CHECK_HID (ndim)

    esize = H5Tget_size (mtype);
    CHECK_HID (esize)
    for (i = 0; i < cnt; i++) {
        err = hdf5_get_vara (vid, mtype, dxplid, mstarts[i], mcounts[i], bufp);
        CHECK_ERR
        rsize = esize;
        for (j = 0; j < ndim; j++) { rsize *= mcounts[i][j]; }
        bufp += rsize;
    }

fn_exit:;
    if (dsid >= 0) H5Sclose (dsid);

    return (int)herr;
}

int hdf5_put_att (
    hid_t fid, int vid, const char *name, hid_t atype, MPI_Offset nelems, const void *buf) {
    herr_t herr = 0;
    hid_t asid = -1, aid = -1;
    hid_t did;
    hsize_t asize;
    htri_t exists;

    asize = (size_t)nelems;
    asid  = H5Screate_simple (1, &asize, &asize);
    CHECK_HID (asid)

    if (vid == NC_GLOBAL)
        did = fid;
    else
        did = f_dids[vid];

    exists = H5Aexists (did, name);
    CHECK_HID (exists)
    if (!exists) {
        aid = H5Acreate2 (did, name, atype, asid, H5P_DEFAULT, H5P_DEFAULT);
    } else {
        aid = H5Aopen (did, name, H5P_DEFAULT);
    }
    CHECK_HID (aid)

    herr = H5Awrite (aid, atype, buf);
    CHECK_HERR

fn_exit:;
    if (asid >= 0) H5Sclose (asid);
    if (aid >= 0) H5Aclose (aid);
    return (int)herr;
}

int hdf5_get_att (hid_t fid, int vid, const char *name, hid_t atype, void *buf) {
    herr_t herr = 0;
    hid_t aid   = -1;
    hid_t did;

    if (vid == NC_GLOBAL)
        did = fid;
    else
        did = f_dids[vid];

    aid = H5Aopen (did, name, H5P_DEFAULT);
    CHECK_HID (aid)
    herr = H5Aread (aid, atype, buf);
    CHECK_HERR

fn_exit:;
    if (aid >= 0) H5Aclose (aid);
    return (int)herr;
}

int hdf5_def_var (hid_t fid, const char *name, nc_type nctype, int ndim, int *dimids, int *vid) {
    herr_t herr = 0;
    int i;
    hid_t did;
    hid_t sid    = -1;
    hid_t dcplid = -1;
    hsize_t dims[H5S_MAX_RANK], mdims[H5S_MAX_RANK];

    dcplid = H5Pcreate (H5P_DATASET_CREATE);
    CHECK_HID (dcplid)

    for (i = 0; i < ndim; i++) { dims[i] = mdims[i] = f_dims[dimids[i]]; }
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

    did = H5Dcreate2 (fid, name, nc_type_to_hdf5_type (nctype), sid, H5P_DEFAULT, dcplid,
                      H5P_DEFAULT);
    CHECK_HID (did)

    f_dids[f_nd] = did;
    *vid         = f_nd++;

fn_exit:;
    if (sid != -1) H5Sclose (sid);
    if (dcplid != -1) H5Pclose (dcplid);
    return (int)herr;
}

int hdf5_inq_varid (hid_t fid, const char *name, int *vid) {
    herr_t herr = 0;
    hid_t did;

    did = H5Dopen2 (fid, name, H5P_DEFAULT);
    CHECK_HID (did)

    f_dids[f_nd] = did;
    *vid         = f_nd++;

fn_exit:;
    return (int)herr;
}

int hdf5_def_dim (hid_t fid, const char *name, MPI_Offset msize, int *did) {
    herr_t herr = 0;
    int i;
    hid_t sid = -1;
    hid_t aid = -1;
    hsize_t size;
    char aname[128];

    size = (hsize_t)msize;
    if (size == NC_UNLIMITED) size = H5S_UNLIMITED;

    sid = H5Screate (H5S_SCALAR);
    CHECK_HID (sid)

    sprintf (aname, "_NCDIM_%s", name);
    aid = H5Acreate2 (fid, aname, H5T_NATIVE_HSIZE, sid, H5P_DEFAULT, H5P_DEFAULT);
    CHECK_HID (aid)

    herr = H5Awrite (aid, H5T_NATIVE_HSIZE, &size);

    f_dims[f_ndim] = size;
    *did           = f_ndim++;

fn_exit:;
    if (aid != -1) H5Aclose (aid);
    if (sid != -1) H5Sclose (sid);
    return (int)herr;
}

int hdf5_inq_dimid (hid_t fid, const char *name, int *did) {
    herr_t herr = 0;
    int i;
    hid_t sid = -1;
    hid_t aid;
    hsize_t size;
    char aname[128];

    sprintf (aname, "_NCDIM_%s", name);
    aid = H5Aopen (fid, aname, H5P_DEFAULT);
    CHECK_HID (aid)

    herr = H5Aread (aid, H5T_NATIVE_HSIZE, &size);
    CHECK_HERR

    f_dims[f_ndim] = size;
    *did           = f_ndim++;

fn_exit:;
    if (aid != -1) H5Aclose (aid);
    return (int)herr;
}

int hdf5_inq_dimlen (hid_t fid, int did, MPI_Offset *msize) {
    *msize = f_dims[did];
    return 0;
}

int hdf5_inq_file_info (hid_t fid, MPI_Info *info) {
    herr_t herr = 0;
    hid_t pid;

    pid = H5Fget_access_plist (fid);
    CHECK_HID (pid);
    herr = H5Pget_fapl_mpio (pid, NULL, info);
    CHECK_HERR

fn_exit:;
    if (pid != -1) H5Pclose (pid);
    return (int)herr;
}

int hdf5_inq_put_size (hid_t fid, size_t *size) {
    herr_t herr = 0;

    *size = 0;

fn_exit:;

    return (int)herr;
}

int hdf5_inq_get_size (hid_t fid, size_t *size) {
    herr_t herr = 0;

    *size = 0;

fn_exit:;

    return (int)herr;
}

int hdf5_close_vars (hid_t fid) {
    herr_t herr = 0;
    int i;

    for (i = 0; i < f_nd; i++) {
        herr = H5Dclose (f_dids[i]);
        CHECK_HERR
    }

fn_exit:;
    return (int)herr;
}
