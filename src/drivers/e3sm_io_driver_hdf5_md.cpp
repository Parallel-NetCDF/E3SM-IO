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

#include <hdf5.h>

#include <e3sm_io.h>
#include <e3sm_io_err.h>
#include <e3sm_io_driver_hdf5.hpp>
#include <e3sm_io_profile.hpp>

int e3sm_io_driver_hdf5::hdf5_file::flush_multidatasets () {
#ifndef HDF5_HAVE_MULTI_DATASET_API
    return -1;
#else
    int err = 0;
    herr_t herr=0;
    size_t i, ndsets;

    if (! driver.use_dwrite_multi) return 0;

    /* flush pending write requests */
    ndsets = wdset_ids.size();
    if (ndsets > 0) {
#if 1
        herr = H5Dwrite_multi(ndsets,
                              wdset_ids.data(),
                              wmem_type_ids.data(),
                              wmem_space_ids.data(),
                              wdset_space_ids.data(),
                              driver.dxplid_coll,
                              (const void **)(wbufs.data()));
        CHECK_HERR
#else
/* Current Issues:
 *     Cannot mix fixed-size and record-size datasets
 *     Cannot write across multiple time dimensions
 *     Cannot interleave writes with and without type conversion
 */
        if (driver.nfixVars > 0) {
            /* first 'driver.nfixVars' datasets are fixed-size */
            herr = H5Dwrite_multi(driver.nfixVars,
                                  wdset_ids.data(),
                                  wmem_type_ids.data(),
                                  wmem_space_ids.data(),
                                  wdset_space_ids.data(),
                                  driver.dxplid_coll,
                                  (const void **)(wbufs.data()));
            CHECK_HERR
        }

        /* the remaining datasets are record-size, i.e. chunked with flexible dimension */
        herr = H5Dwrite_multi(ndsets-driver.nfixVars,
                              wdset_ids.data()+driver.nfixVars,
                              wmem_type_ids.data()+driver.nfixVars,
                              wmem_space_ids.data()+driver.nfixVars,
                              wdset_space_ids.data()+driver.nfixVars,
                              driver.dxplid_coll,
                              (const void **)(wbufs.data()+driver.nfixVars));
        CHECK_HERR
#endif

        for (i=0; i<ndsets; i++) {
            H5Sclose (wmem_space_ids[i]);
            CHECK_HERR
            H5Sclose (wdset_space_ids[i]);
            CHECK_HERR
        }
        wdset_ids.clear();
        wmem_type_ids.clear();
        wmem_space_ids.clear();
        wdset_space_ids.clear();
        wbufs.clear();
    }

    /* flush pending read requests */
    ndsets = rdset_ids.size();
    if (ndsets > 0) {
        herr = H5Dread_multi(ndsets,
                             rdset_ids.data(),
                             rmem_type_ids.data(),
                             rmem_space_ids.data(),
                             rdset_space_ids.data(),
                             driver.dxplid_coll,
                             (void **)(rbufs.data()));
        CHECK_HERR

        for (i=0; i<ndsets; i++) {
            H5Sclose (rmem_space_ids[i]);
            CHECK_HERR
            H5Sclose (rdset_space_ids[i]);
            CHECK_HERR
        }
        rdset_ids.clear();
        rmem_type_ids.clear();
        rmem_space_ids.clear();
        rdset_space_ids.clear();
        rbufs.clear();
    }

err_out:
    return err;
#endif
}

/* post a varn request, which will later be flushed at wait() */
int e3sm_io_driver_hdf5::post_varn(int            fid,
                                   int            vid,
                                   MPI_Datatype   itype,
                                   int            nreqs,
                                   MPI_Offset   **starts,
                                   MPI_Offset   **counts,
                                   void          *buf,
                                   bool           isWrite)
{
    int i, j, err = 0, ndim;
    herr_t herr;
    hid_t dsid = -1, did;
    hsize_t start[E3SM_IO_DRIVER_MAX_RANK], block[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t dims[E3SM_IO_DRIVER_MAX_RANK];
    hsize_t rsize, rsize_all = 0;
    hdf5_file *fp = this->files[fid];
    H5S_seloper_t op = H5S_SELECT_SET;

    did = fp->dids[vid];

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim = H5Sget_simple_extent_dims (dsid, dims, NULL);
    CHECK_HID (ndim)

    /* set filespace selection */
    if (ndim > 0) { /* select hyperslab filespace only for non-scalar datasets */
        for (i=0; i<nreqs; i++) {
            rsize = 1;
            for (j = 0; j < ndim; j++) { rsize *= counts[i][j]; }
            if (rsize == 0) continue;
            rsize_all += rsize;

            /* type cast from MPI_Offset to hsize_t */
            for (j = 0; j < ndim; j++) {
                start[j] = (hsize_t)starts[i][j];
                block[j] = (hsize_t)counts[i][j];
            }

            /* union all nreqs hyperslabs */
            herr = H5Sselect_hyperslab (dsid, op, start, NULL, this->one, block);
            CHECK_HERR
            op = H5S_SELECT_OR;
        }
    }

    if (rsize_all > 0) {
        size_t tsize;
        hid_t mtype, msid, tid;

        mtype = e3sm_io_type_mpi2hdf5 (itype);

        /* create memory space */
        msid = H5Screate_simple (1, &rsize_all, NULL);
        CHECK_HID (msid)

        /* queue this varn request */
        if (isWrite) {
            fp->wdset_ids.push_back(did);
            fp->wmem_type_ids.push_back(mtype);
            fp->wmem_space_ids.push_back(msid);
            fp->wdset_space_ids.push_back(dsid);
            fp->wbufs.push_back(buf);
        } else {
            fp->rdset_ids.push_back(did);
            fp->rmem_type_ids.push_back(mtype);
            fp->rmem_space_ids.push_back(msid);
            fp->rdset_space_ids.push_back(dsid);
            fp->rbufs.push_back(buf);
        }

        /* retrieve external data type */
        tid   = H5Dget_type (did);
        tsize = H5Tget_size (tid);
        if (isWrite)
            this->amount_WR += tsize * rsize_all;
        else
            this->amount_RD += tsize * rsize_all;
        H5Tclose (tid);

        if (! fp->dset_isRec[vid]) this->nfixVars++;
    }
    else {
        herr = H5Sclose(dsid);
        CHECK_HERR
    }

err_out:
    return err;
}

