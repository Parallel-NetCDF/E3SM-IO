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
#ifdef ENABLE_LOGVOL
#include "H5VL_log.h"
#endif
//
#include <e3sm_io.h>
#include <e3sm_io_err.h>

#include <e3sm_io_driver_hdf5.hpp>
#include <e3sm_io_driver_hdf5_int.hpp>
#include <e3sm_io_profile.hpp>

int e3sm_io_driver_hdf5::index_order_cmp (const void *a, const void *b) {
    return (((Index_order *)a)->index - ((Index_order *)b)->index);
}

static inline herr_t count_data (int cnt, int ndim, MPI_Offset **blocks, int *index) {
    int j, k;
    hsize_t rsize;
    index[0] = 0;
    for (k = 0; k < cnt; ++k) {
        rsize = 1;
        for (j = 0; j < ndim; j++) { rsize *= blocks[k][j]; }
        if (rsize) {
            if (ndim == 1) {
                index[0]++;
            } else if (ndim == 2) {
                index[0] += blocks[k][0];
            } else if (ndim == 3) {
                index[0] += blocks[k][0] * blocks[k][1];
            } else {
                printf ("critical error, dimension is greater than 3.\n");
                return -1;
            }
        }
    }
    return 0;
}

inline herr_t e3sm_io_driver_hdf5::pack_data (Index_order *index_order,
                                              int *index,
                                              char *src,
                                              hsize_t esize,
                                              int ndim,
                                              const hsize_t *dims,
                                              const hsize_t *start,
                                              const hsize_t *block) {
    unsigned i, j;
    if (ndim == 1) {
        index_order[index[0]].index    = start[0];
        index_order[index[0]].coverage = esize * block[0];
        index_order[index[0]].data     = src;
        index[0]++;
    } else if (ndim == 2) {
        for (i = 0; i < block[0]; ++i) {
            index_order[index[0]].index    = start[1] + start[0] * dims[1];
            index_order[index[0]].coverage = esize * block[1];
            index_order[index[0]].data     = src;
            src += index_order[index[0]].coverage;
            index[0]++;
        }
    } else if (ndim == 3) {
        for (i = 0; i < block[0]; ++i) {
            for (j = 0; j < block[1]; ++j) {
                index_order[index[0]].index =
                    start[0] * dims[1] * dims[2] + start[1] * dims[2] + start[2];
                index_order[index[0]].coverage = esize * block[2];
                index_order[index[0]].data     = src;
                src += index_order[index[0]].coverage;
                index[0]++;
            }
        }
    } else {
        printf ("critical error, dimension is greater than 3.\n");
        return -1;
    }
    return 0;
}

inline herr_t e3sm_io_driver_hdf5::copy_index_buf (Index_order *index_order,
                                                   int total_blocks,
                                                   char *out_buf) {
    hsize_t displs = 0;
    int i;
    for (i = 0; i < total_blocks; ++i) {
        memcpy (out_buf + displs, index_order[i].data, index_order[i].coverage);
        displs += index_order[i].coverage;
    }
    return 0;
}

herr_t e3sm_io_driver_hdf5::hdf5_file::register_multidataset (
    void *buf, hid_t did, hid_t dsid, hid_t msid, hid_t mtype, int write) {
    H5D_rw_multi_t dset;

    dset.dset_id       = did;
    dset.mem_space_id  = msid;
    dset.dset_space_id = dsid;
    dset.mem_type_id   = mtype;
    dset.buf = buf;
    if (write) {
        wreqs.push_back (dset);
    } else {
        rreqs.push_back (dset);
    }

    return 0;
}

static inline void print_no_collective_cause (uint32_t local_no_collective_cause,
                                              uint32_t global_no_collective_cause) {
    switch (local_no_collective_cause) {
        case H5D_MPIO_COLLECTIVE: {
            // printf("MPI-IO collective successful\n");
            break;
        }
        case H5D_MPIO_SET_INDEPENDENT: {
            printf ("local flag: MPI-IO independent flag is on\n");
            break;
        }
        case H5D_MPIO_DATATYPE_CONVERSION: {
            printf ("local flag: MPI-IO datatype conversion needed\n");
            break;
        }
        case H5D_MPIO_DATA_TRANSFORMS: {
            printf ("local flag: MPI-IO H5D_MPIO_DATA_TRANSFORMS.\n");
            break;
        }
            /*
                case H5D_MPIO_SET_MPIPOSIX: {
                    printf("local flag: MPI-IO H5D_MPIO_SET_MPIPOSIX \n");
                }
            */
        case H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES: {
            printf ("local flag: MPI-IO NOT_SIMPLE_OR_SCALAR_DATASPACES\n");
            break;
        }
            /*
                case H5D_MPIO_POINT_SELECTIONS: {
                    printf("local flag: MPI-IO H5D_MPIO_POINT_SELECTIONS\n");
                }
            */
        case H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET: {
            printf ("local flag: MPI-IO H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET\n");
            break;
        }
            /*
                case H5D_MPIO_FILTERS: {
                    printf("local flag: MPI-IO H5D_MPIO_FILTERS\n");
                    break;
                }
            */
        default: {
            printf ("undefined label for collective cause\n");
            break;
        }
    }

    switch (global_no_collective_cause) {
        case H5D_MPIO_COLLECTIVE: {
            // printf("MPI-IO collective successful\n");
            break;
        }
        case H5D_MPIO_SET_INDEPENDENT: {
            printf ("global flag: MPI-IO independent flag is on\n");
            break;
        }
        case H5D_MPIO_DATATYPE_CONVERSION: {
            printf ("global flag: MPI-IO datatype conversion needed\n");
            break;
        }
        case H5D_MPIO_DATA_TRANSFORMS: {
            printf ("global flag: MPI-IO H5D_MPIO_DATA_TRANSFORMS.\n");
            break;
        }
            /*
                case H5D_MPIO_SET_MPIPOSIX: {
                    printf("global flag: MPI-IO H5D_MPIO_SET_MPIPOSIX \n");
                }
            */
        case H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES: {
            printf ("global flag: MPI-IO NOT_SIMPLE_OR_SCALAR_DATASPACES\n");
            break;
        }
            /*
                case H5D_MPIO_POINT_SELECTIONS: {
                    printf("global flag: MPI-IO H5D_MPIO_POINT_SELECTIONS\n");
                }
            */
        case H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET: {
            printf ("global flag: MPI-IO H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET\n");
            break;
        }
            /*
                case H5D_MPIO_FILTERS: {
                    printf("global flag: MPI-IO H5D_MPIO_FILTERS\n");
                    break;
                }
            */
        default: {
            printf ("undefined label for collective cause\n");
            break;
        }
    }
}

int e3sm_io_driver_hdf5::hdf5_file::flush_multidatasets () {
    herr_t herr=0;
    int err = 0, rank;
    size_t i, j;
    size_t esize;   // Element size of the memory type
    hsize_t dims[H5S_MAX_RANK], mdims[H5S_MAX_RANK];    // Memory space dims
    hid_t dsid = -1; // Dataset space ID for following dummy call
    hid_t tid = -1; // Dataset type ID for following dummy call
    size_t ndset = dids.size();    // # datasets
    std::vector<std::vector<H5D_rw_multi_t>> reqs(ndset);   // Requests organized by datasaet ID
    int *nreqs = NULL;  // # requests per dataset
    int *nreqs_all;  // max # requests per dataset across all processes

    MPI_Comm_rank (comm, &rank);

#ifdef HDF5_HAVE_DWRITE_MULTI
    if (this->driver.use_dwrite_multi) {    
        H5Dwrite_multi (this->driver.dxplid_coll, wreqs.size (), wreqs.data());
    }
    else
#endif
    {
        if (this->driver.collective_flush) {
            // Organize by did
            for(auto &r: wreqs){
                reqs[inv_dids[r.dset_id]].push_back(r);
            }

            // Count #req per dataset
            nreqs = (int*)malloc(sizeof(int) * ndset * 2);
            nreqs_all = nreqs + ndset;
            for(i = 0; i < ndset; i++){
                nreqs[i] = reqs[i].size();
            }
            err = MPI_Allreduce(nreqs, nreqs_all, ndset, MPI_INT, MPI_MAX, comm);
            CHECK_MPIERR

            for(j = 0; j < ndset; j++){
                // Dummy space for dummy call
                dsid = H5Dget_space (dids[j]);
                CHECK_HID(dsid);
                herr = H5Sselect_none (dsid);
                CHECK_HERR

                // Dummy type for dummy call
                tid = H5Dget_type (dids[j]);
                CHECK_HID(tid)

                for (i = 0; i < (size_t)nreqs_all[j]; ++i) {
                    if (i < (size_t)nreqs[j]){
                        herr = H5Dwrite (dids[j], reqs[j][i].mem_type_id,
                                        reqs[j][i].mem_space_id, reqs[j][i].dset_space_id,
                                        driver.dxplid_coll, reqs[j][i].buf);
                    }
                    else{   // Follow collective I/O with dummy call
                        herr = H5Dwrite (dids[j], tid,
                                        H5S_ALL, dsid,
                                        driver.dxplid_coll, NULL);
                    }
                    CHECK_HERR
                    /* Lagacy code to study why HDF5 does not follow collective dxpl
                    if (!rank) {
                        uint32_t local_no_collective_cause, global_no_collective_cause;
                        H5Pget_mpio_no_collective_cause (this->driver.dxplid_coll, &local_no_collective_cause,
                                                        &global_no_collective_cause);
                        print_no_collective_cause (local_no_collective_cause, global_no_collective_cause);
                    }
                    */
                }

                H5Sclose (dsid);
                dsid = -1;

                H5Tclose (tid);
                tid = -1;
            }
        }
        else{
            for(auto &r: wreqs){
                herr = H5Dwrite (r.dset_id, r.mem_type_id, r.mem_space_id, r.dset_space_id,
                                 driver.dxplid_indep, r.buf);
                CHECK_HERR
            }
        }
    }

    // Count data size
    for (i = 0; i < wreqs.size (); ++i) {
        H5Sget_simple_extent_dims (wreqs[i].mem_space_id, dims, mdims);
        esize = H5Tget_size (wreqs[i].mem_type_id);
        E3SM_IO_TIMER_ADD (E3SM_IO_TIMER_HDF5_DSIZE, (dims[0] * esize))
    }

    // Free data spaces
    for (auto &req : wreqs) {
        H5Sclose (req.dset_space_id);
        H5Sclose (req.mem_space_id);
    }
    wreqs.clear ();

err_out:;
    if (dsid >= 0) {
        H5Sclose (dsid);
    }
    if (tid >= 0) {
        H5Tclose (tid);
    }
    return err;
}

herr_t e3sm_io_driver_hdf5::hdf5_file::pull_multidatasets () {
    herr_t herr=0;
    int err = 0, rank;
    char *temp_buf   = NULL, *temp_buf_ptr;
    size_t i, j, temp_size = 0;
    size_t esize;   // Element size of the memory type
    hsize_t dims[H5S_MAX_RANK], mdims[H5S_MAX_RANK];    // Memory space dims
    hid_t dsid = -1; // Dataset space ID for following dummy call
    hid_t tid = -1; // Dataset type ID for following dummy call
    size_t ndset = dids.size();    // # datasets
    std::vector<std::vector<H5D_rw_multi_t>> reqs(ndset);   // Requests organized by datasaet ID
    int *nreqs = NULL;  // # requests per dataset
    int *nreqs_all;  // max # requests per dataset across all processes
    
    MPI_Comm_rank (comm, &rank);

    // printf("Rank %d number of datasets to be written %d\n", rank, multi_datasets.size());
#ifdef HDF5_HAVE_DWRITE_MULTI
    if (this->driver.use_dwrite_multi) {    
        H5Dread_multi (this->driver.dxplid_coll, rreqs.size (), rreqs.data());
    }
    else
#endif
    {
        if (this->driver.collective_flush) {
            // Organize by did
            for(auto &r: rreqs){
                reqs[inv_dids[r.dset_id]].push_back(r);
            }

            // Count #req per dataset
            nreqs = (int*)malloc(sizeof(int) * ndset * 2);
            nreqs_all = nreqs + ndset;
            for(i = 0; i < ndset; i++){
                nreqs[i] = reqs[i].size();
            }
            err = MPI_Allreduce(nreqs, nreqs_all, ndset, MPI_INT, MPI_MAX, comm);
            CHECK_MPIERR

            for(j = 0; j < ndset; j++){
                // Dummy space for dummy call
                dsid = H5Dget_space (dids[j]);
                CHECK_HID(dsid);
                herr = H5Sselect_none (dsid);
                CHECK_HERR

                // Dummy type for dummy call
                tid = H5Dget_type (dids[j]);
                CHECK_HID(tid)

                for (i = 0; i < (size_t)nreqs_all[j]; ++i) {
                    if (i < (size_t)nreqs[j]){
                        herr = H5Dread (dids[j], reqs[j][i].mem_type_id,
                                        reqs[j][i].mem_space_id, reqs[j][i].dset_space_id,
                                        driver.dxplid_coll, reqs[j][i].buf);
                    }
                    else{   // Follow collective I/O with dummy call
                        herr = H5Dread (dids[j], tid,
                                        H5S_ALL, dsid,
                                        driver.dxplid_coll, NULL);
                    }
                    CHECK_HERR
                    /* Lagacy code to study why HDF5 does not follow collective dxpl
                    if (!rank) {
                        uint32_t local_no_collective_cause, global_no_collective_cause;
                        H5Pget_mpio_no_collective_cause (this->driver.dxplid_coll, &local_no_collective_cause,
                                                        &global_no_collective_cause);
                        print_no_collective_cause (local_no_collective_cause, global_no_collective_cause);
                    }
                    */
                }

                H5Sclose (dsid);
                dsid = -1;

                H5Tclose (tid);
                tid = -1;
            }
        }
        else{
            for (auto &r: rreqs){
                herr = H5Dread (r.dset_id, r.mem_type_id, r.mem_space_id, r.dset_space_id,
                                 driver.dxplid_indep, r.buf);
                CHECK_HERR
            } 
        }
    }

    // Count data size
    for (i = 0; i < rreqs.size (); ++i) {
        H5Sget_simple_extent_dims (rreqs[i].mem_space_id, dims, mdims);
        esize = H5Tget_size (rreqs[i].mem_type_id);
        // data using the recorded dataset_segments.
        E3SM_IO_TIMER_ADD (E3SM_IO_TIMER_HDF5_DSIZE, (dims[0] * esize))
    }

    // printf("rank %d number of hyperslab called %d\n", rank, hyperslab_count);

    // Read is completed here, but data is out-of-order in user buffer. We need to rearrange all
    // data using the recorded dataset_segments.
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_CPY)
    for (i = 0; i < rreqs.size (); ++i) {
        // First, we make a copy of data in the current dataset.
        H5Sget_simple_extent_dims (rreqs[i].mem_space_id, dims, mdims);

        esize = H5Tget_size (rreqs[i].mem_type_id);
        if (dims[0] * esize > temp_size) {
            if (temp_size) {
                free (temp_buf);
                temp_size = 0;
            }
            temp_size = dims[0] * esize;
            temp_buf  = (char *)malloc (sizeof (char) * temp_size);
        }
        memcpy (temp_buf, rreqs[i].buf, esize * dims[0]);
        // Copy data back from temp_buf to user memory defined by data_segments. This array is
        // sorted previously to align the HDF5 memory space.
        temp_buf_ptr = temp_buf;
        for (j = 0; j < dataset_segments[i].size (); ++j) {
            memcpy (dataset_segments[i][j].data, temp_buf_ptr, dataset_segments[i][j].coverage);
            temp_buf_ptr += dataset_segments[i][j].coverage;
        }
    }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_CPY)
    if (temp_size) { free (temp_buf); }

    rreqs.clear ();

err_out:;
    if (dsid >= 0) {
        H5Sclose (dsid);
    }
    if (tid >= 0) {
        H5Tclose (tid);
    }
    return 0;
}

// TODO: Fix bug of not extending the datasets when needed
// This function was not in use for now
int e3sm_io_driver_hdf5::put_varn_merge (int fid,
                                         int vid,
                                         MPI_Datatype type,
                                         int nreq,
                                         MPI_Offset **starts,
                                         MPI_Offset **counts,
                                         void *buf,
                                         e3sm_io_op_mode mode) {
    herr_t herr   = 0;
    int err = 0;
    hdf5_file *fp = this->files[fid];
    int i, j;
    hsize_t esize, rsize, rsize_old = 0, memspace_size, total_memspace_size, hyperslab_set;
    int ndim;
    hid_t dsid = -1, msid = -1;
    hid_t mtype;
    char *bufp = (char *)buf;
    hid_t did;
    hid_t dxplid;
    hsize_t start[H5S_MAX_RANK], block[H5S_MAX_RANK];
    hsize_t dims[H5S_MAX_RANK], mdims[H5S_MAX_RANK];
    char *buf2;
    int index;
    Index_order *index_order;
    int total_blocks;
    hid_t tid = -1;
    MPI_Offset putsize;
    MPI_Offset **count;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    did = fp->dids[vid];

    mtype = mpi_type_to_hdf5_type (type);
    esize = (hsize_t)H5Tget_size (mtype);

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim = H5Sget_simple_extent_dims (dsid, dims, mdims);

    herr = count_data (nreq, ndim, counts, &total_blocks);
    CHECK_HERR
    index_order = (Index_order *)malloc (sizeof (Index_order) * total_blocks);
    index       = 0;

    // Call H5DWrite
    int rank;
    MPI_Comm_rank (fp->comm, &rank);

#if defined(DEBUG) && DEBUG == 1
    char filename[128];
    FILE *stream;
    if (rank == 0) {
        sprintf (filename, "e3sm_hdf5_reqs.txt");
        stream = fopen (filename, "r");
        if (stream) {
            fclose (stream);
            stream = fopen (filename, "a");
        } else {
            stream = fopen (filename, "w");
        }
    }
#endif
    hyperslab_set       = 0;
    total_memspace_size = 0;
    for (i = 0; i < nreq; i++) {
        rsize = esize;
        for (j = 0; j < ndim; j++) { rsize *= counts[i][j]; }
        memspace_size = rsize / esize;
        total_memspace_size += memspace_size;
        if (rsize) {
            for (j = 0; j < ndim; j++) {
                start[j] = (hsize_t)starts[i][j];
                block[j] = (hsize_t)counts[i][j];
#if defined(DEBUG) && DEBUG == 1
                if (rank == 0) { fprintf (stream, "%llu+%llu;", start[j], block[j]); }
#endif
            }
#if defined(DEBUG) && DEBUG == 1
            if (rank == 0) { fprintf (stream, ","); }
#endif
            // Recreate only when size mismatch
            if (rsize != rsize_old) { rsize_old = rsize; }

            E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)
            if (!hyperslab_set) {
                herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, start, NULL, this->one, block);
                hyperslab_set = 1;
            } else {
                herr = H5Sselect_hyperslab (dsid, H5S_SELECT_OR, start, NULL, this->one, block);
            }
            CHECK_HERR
            E3SM_IO_TIMER_SWAP (E3SM_IO_TIMER_HDF5_SEL, E3SM_IO_TIMER_HDF5_WR)
            herr = pack_data (index_order, &index, bufp, esize, ndim, dims, start, block);
            CHECK_HERR
            E3SM_IO_TIMER_ADD (E3SM_IO_TIMER_HDF5_NSLAB, 1)
            E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_WR)

            bufp += rsize;
        }
    }
#if defined(DEBUG) && DEBUG == 1
    if (rank == 0) {
        fprintf (stream, "\n");
        fclose (stream);
    }
#endif

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SORT_REQ)
    qsort (index_order, total_blocks, sizeof (Index_order), index_order_cmp);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_SORT_REQ)

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_CPY)
    buf2 = (char *)malloc (esize * total_memspace_size);
    herr = copy_index_buf (index_order, total_blocks, buf2);
    CHECK_HERR
    memcpy (buf, buf2, esize * total_memspace_size);
    free (index_order);
    free (buf2);

    msid = H5Screate_simple (1, &total_memspace_size, &total_memspace_size);
    CHECK_HID (msid)

    mode = indep;
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

    // Call H5Dread immediately if blocking or Log VOL is used
    if ((mode == indep) || (mode == coll)) {
        herr = H5Dwrite (did, mpi_type_to_hdf5_type (type), msid, dsid, dxplid, buf);
        CHECK_HERR
    } else {  // Otherwier, queue request in driver
        herr = fp->register_multidataset (buf, did, dsid, msid, mpi_type_to_hdf5_type (type), 1);
        CHECK_HERR
        // Prevent freeing of dsid and msid, they will be freed after flush
        dsid = msid = -1;
    }
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_CPY)
    /* The folowing code is to place dummy H5Dwrite for collective call.*/

    tid     = H5Dget_type (did);
    putsize = H5Tget_size (tid);
    for (count = counts; count < counts + nreq; count++) {
        for (i = 0; i < ndim; i++) { putsize *= (*count)[i]; }
    }
    fp->putsize += putsize;

err_out:
    if (tid >= 0) H5Tclose (tid);
    if (dsid >= 0) H5Sclose (dsid);
    if (msid >= 0) H5Sclose (msid);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}

// TODO: Fix bug of not extending the datasets when needed
// This function was not in use for now
int e3sm_io_driver_hdf5::get_varn_merge (int fid,
                                         int vid,
                                         MPI_Datatype type,
                                         int nreq,
                                         MPI_Offset **starts,
                                         MPI_Offset **counts,
                                         void *buf,
                                         e3sm_io_op_mode mode) {
    herr_t herr   = 0;
    int err = 0;
    hdf5_file *fp = this->files[fid];
    int i, j, index;
    size_t esize, rsize;
    hsize_t total_mem_size = 0;
    int ndim;
    char *bufp = (char *)buf;
    hid_t dsid = -1, msid = -1;
    hid_t mtype;
    hid_t did;
    hid_t dxplid;
    hsize_t dims[H5S_MAX_RANK], mdims[H5S_MAX_RANK];
    hsize_t start[H5S_MAX_RANK], block[H5S_MAX_RANK];
    int total_blocks;
    hid_t tid = -1;
    MPI_Offset getsize;
    MPI_Offset **count;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5)

    did = fp->dids[vid];

    mtype = mpi_type_to_hdf5_type (type);

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim  = H5Sget_simple_extent_dims (dsid, dims, mdims);
    esize = H5Tget_size (mtype);

    herr = count_data (nreq, ndim, counts, &total_blocks);
    CHECK_HERR
    if (fp->dataset_segments.size () <= fp->rreqs.size ()) {
        fp->dataset_segments.resize (fp->rreqs.size ());
    }
    fp->dataset_segments[fp->rreqs.size ()].resize (total_blocks);

    index = 0;
    for (i = 0; i < nreq; ++i) {
        for (j = 0; j < ndim; ++j) {
            start[j] = (hsize_t)starts[i][j];
            block[j] = (hsize_t)counts[i][j];
        }

        herr = pack_data (fp->dataset_segments[fp->rreqs.size ()].data (), &index, bufp,
                          esize, ndim, dims, start, block);
        CHECK_HERR
        E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SEL)
        if (i) {
            herr = H5Sselect_hyperslab (dsid, H5S_SELECT_OR, start, NULL, this->one, block);
        } else {
            herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, start, NULL, this->one, block);
        }
        CHECK_HERR
        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_SEL)
        rsize = 1;
        for (j = 0; j < ndim; ++j) { rsize *= block[j]; }
        total_mem_size += rsize;
        bufp += rsize * esize;
    }

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_HDF5_SORT_REQ)
    qsort (fp->dataset_segments[fp->rreqs.size ()].data (), total_blocks,
           sizeof (Index_order), index_order_cmp);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5_SORT_REQ)

    msid = H5Screate_simple (1, &total_mem_size, &total_mem_size);
    CHECK_HID (msid)

    // dataset_size is incremented in the following function.

    mode = indep;
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

    // Call H5Dread immediately if blocking or Log VOL is used
    if ((mode == indep) || (mode == coll)) {
        herr = H5Dread (did, mpi_type_to_hdf5_type (type), msid, dsid, dxplid, buf);
        CHECK_HERR
    } else {  // Otherwier, queue request in driver
        herr = fp->register_multidataset (buf, did, dsid, msid, mpi_type_to_hdf5_type (type), 0);
        CHECK_HERR
        // Prevent freeing of dsid and msid, they will be freed after flush
        dsid = msid = -1;
    }

    tid     = H5Dget_type (did);
    getsize = H5Tget_size (tid);
    for (count = counts; count < counts + nreq; count++) {
        for (i = 0; i < ndim; i++) { getsize *= (*count)[i]; }
    }
    fp->getsize += getsize;

err_out:
    if (tid >= 0) H5Tclose (tid);
    if (dsid >= 0) H5Sclose (dsid);
    if (msid >= 0) H5Sclose (msid);
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_HDF5)
    return err;
}
