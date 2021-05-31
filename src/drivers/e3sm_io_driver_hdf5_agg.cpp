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
    if (write) {
        dset.u.wbuf = buf;
    } else {
        dset.u.rbuf = buf;
    }

    multi_datasets.push_back (dset);

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
    herr_t herr;
    int nerrs = 0;
    int i;
    uint32_t local_no_collective_cause, global_no_collective_cause;
    int rank;
    size_t esize;
    hsize_t dims[H5S_MAX_RANK], mdims[H5S_MAX_RANK];

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

#ifdef HDF5_HAVE_DWRITE_MULTI
    H5Dwrite_multi (dxplid_coll, multi_datasets.size (), multi_datasets);
#else
    for (i = 0; i < multi_datasets.size (); ++i) {
        // MPI_Barrier(MPI_COMM_WORLD);
        herr = H5Dwrite (multi_datasets[i].dset_id, multi_datasets[i].mem_type_id,
                         multi_datasets[i].mem_space_id, multi_datasets[i].dset_space_id,
                         driver.dxplid_coll, multi_datasets[i].u.wbuf);
        CHECK_HERR

        if (!rank) {
            H5Pget_mpio_no_collective_cause (this->driver.dxplid_coll, &local_no_collective_cause,
                                             &global_no_collective_cause);
            print_no_collective_cause (local_no_collective_cause, global_no_collective_cause);
        }
    }
#endif

    // Count data size
    for (i = 0; i < multi_datasets.size (); ++i) {
        H5Sget_simple_extent_dims (multi_datasets[i].mem_space_id, dims, mdims);
        esize = H5Tget_size (multi_datasets[i].mem_type_id);
        driver.total_data_size += dims[0] * esize;
    }

    // Free data spaces
    for (auto &req : multi_datasets) {
        H5Sclose (req.dset_space_id);
        H5Sclose (req.mem_space_id);
    }
    multi_datasets.clear ();

err_out:;
    return herr;
}

herr_t e3sm_io_driver_hdf5::hdf5_file::pull_multidatasets () {
    int i;
    unsigned j;
    uint32_t local_no_collective_cause, global_no_collective_cause;
    int rank;
    char *temp_buf   = NULL, *temp_buf_ptr;
    size_t temp_size = 0, esize;
    double start;
    hsize_t dims[H5S_MAX_RANK], mdims[H5S_MAX_RANK];

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    // printf("Rank %d number of datasets to be written %d\n", rank, multi_datasets.size());
#ifdef HDF5_HAVE_DWRITE_MULTI
    H5Dread_multi (dxplid_coll, multi_datasets.size (), multi_datasets);
#else
    for (i = 0; i < multi_datasets.size (); ++i) {
        // MPI_Barrier(MPI_COMM_WORLD);
        /*
                hsize_t total_data_size = 1, total_mem_size = 1;
                ndim = H5Sget_simple_extent_dims (multi_datasets[i].dset_space_id, dims, mdims);
                for ( j = 0; j < ndim; ++j ) {
                    printf("dataspace: dims[%d]=%lld\n", j, (long long int) dims[j]);
                    total_data_size *= dims[j];
                }
                ndim = H5Sget_simple_extent_dims (multi_datasets[i].mem_space_id, dims, mdims);
                for ( j = 0; j < ndim; ++j ) {
                    printf("memspace: dims[%d]=%lld\n", j, (long long int) dims[j]);
                    total_mem_size *= dims[j];
                }
                printf("total_data_size = %lld, total_mem_size = %lld\n", (long long int)
           total_data_size, (long long int) total_mem_size);
        */

        H5Dread (multi_datasets[i].dset_id, multi_datasets[i].mem_type_id,
                 multi_datasets[i].mem_space_id, multi_datasets[i].dset_space_id,
                 driver.dxplid_coll, multi_datasets[i].u.rbuf);

        if (!rank) {
            H5Pget_mpio_no_collective_cause (this->driver.dxplid_coll, &local_no_collective_cause,
                                             &global_no_collective_cause);
            print_no_collective_cause (local_no_collective_cause, global_no_collective_cause);
        }
    }
#endif

    // Count data size
    for (i = 0; i < multi_datasets.size (); ++i) {
        H5Sget_simple_extent_dims (multi_datasets[i].mem_space_id, dims, mdims);
        esize = H5Tget_size (multi_datasets[i].mem_type_id);
        driver.total_data_size += dims[0] * esize;
    }

    // printf("rank %d number of hyperslab called %d\n", rank, hyperslab_count);

    // Read is completed here, but data is out-of-order in user buffer. We need to rearrange all
    // data using the recorded dataset_segments.
    start = MPI_Wtime ();
    for (i = 0; i < multi_datasets.size (); ++i) {
        // First, we make a copy of data in the current dataset.
        H5Sget_simple_extent_dims (multi_datasets[i].mem_space_id, dims, mdims);

        esize = H5Tget_size (multi_datasets[i].mem_type_id);
        if (dims[0] * esize > temp_size) {
            if (temp_size) {
                free (temp_buf);
                temp_size = 0;
            }
            temp_size = dims[0] * esize;
            temp_buf  = (char *)malloc (sizeof (char) * temp_size);
        }
        memcpy (temp_buf, multi_datasets[i].u.rbuf, esize * dims[0]);
        // Copy data back from temp_buf to user memory defined by data_segments. This array is
        // sorted previously to align the HDF5 memory space.
        temp_buf_ptr = temp_buf;
        for (j = 0; j < dataset_segments[i].size (); ++j) {
            memcpy (dataset_segments[i][j].data, temp_buf_ptr, dataset_segments[i][j].coverage);
            temp_buf_ptr += dataset_segments[i][j].coverage;
        }
    }
    driver.tcpy = MPI_Wtime () - start;
    if (temp_size) { free (temp_buf); }

    multi_datasets.clear ();

    return 0;
}

int e3sm_io_driver_hdf5::put_varn_merge (int fid,
                                         int vid,
                                         MPI_Datatype type,
                                         int nreq,
                                         MPI_Offset **starts,
                                         MPI_Offset **counts,
                                         void *buf,
                                         e3sm_io_op_mode mode) {
    herr_t herr   = 0;
    int nerrs     = 0;
    hdf5_file *fp = this->files[fid];
    int i, j;
    double ts, te;
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
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

#if DEBUG == 1
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
#if DEBUG == 1
                if (rank == 0) { fprintf (stream, "%llu+%llu;", start[j], block[j]); }
#endif
            }
#if DEBUG == 1
            if (rank == 0) { fprintf (stream, ","); }
#endif
            // Recreate only when size mismatch
            if (rsize != rsize_old) { rsize_old = rsize; }

            ts = MPI_Wtime ();
            if (!hyperslab_set) {
                herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, start, NULL, this->one, block);
                hyperslab_set = 1;
            } else {
                herr = H5Sselect_hyperslab (dsid, H5S_SELECT_OR, start, NULL, this->one, block);
            }
            CHECK_HERR
            te = MPI_Wtime ();
            this->tsel += te - ts;

            herr = pack_data (index_order, &index, bufp, esize, ndim, dims, start, block);
            CHECK_HERR
            this->hyperslab_count++;

            this->twrite += MPI_Wtime () - te;
            bufp += rsize;
        }
    }
#if DEBUG == 1
    if (rank == 0) {
        fprintf (stream, "\n");
        fclose (stream);
    }
#endif

    ts = MPI_Wtime ();
    qsort (index_order, total_blocks, sizeof (Index_order), index_order_cmp);
    this->tsort += MPI_Wtime () - ts;

    ts   = MPI_Wtime ();
    buf2 = (char *)malloc (esize * total_memspace_size);
    herr = copy_index_buf (index_order, total_blocks, buf2);
    CHECK_HERR
    memcpy (buf, buf2, esize * total_memspace_size);
    free (index_order);
    free (buf2);

    msid = H5Screate_simple (1, &total_memspace_size, &total_memspace_size);
    CHECK_HID (msid)

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
    if ((mode == indep) || (mode == coll)
#ifdef ENABLE_LOGVOL
        || this->use_logvol
#endif
    ) {
        herr = H5Dwrite (did, mpi_type_to_hdf5_type (type), msid, dsid, dxplid, buf);
        CHECK_HERR
    } else {  // Otherwier, queue request in driver
        herr = fp->register_multidataset (buf, did, dsid, msid, mpi_type_to_hdf5_type (type), 1);
        CHECK_HERR
        // Prevent freeing of dsid and msid, they will be freed after flush
        dsid = msid = -1;
    }
    this->tcpy += MPI_Wtime () - ts;
    /* The folowing code is to place dummy H5Dwrite for collective call.*/

err_out:;
    if (dsid >= 0) H5Sclose (dsid);
    if (msid >= 0) H5Sclose (msid);
    return nerrs;
}

int e3sm_io_driver_hdf5::get_varn_merge (int fid,
                                         int vid,
                                         MPI_Datatype type,
                                         int nreq,
                                         MPI_Offset **starts,
                                         MPI_Offset **counts,
                                         void *buf,
                                         e3sm_io_op_mode mode) {
    herr_t herr   = 0;
    int nerrs     = 0;
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
    double ts;

    did = fp->dids[vid];

    mtype = mpi_type_to_hdf5_type (type);

    dsid = H5Dget_space (did);
    CHECK_HID (dsid)

    ndim  = H5Sget_simple_extent_dims (dsid, dims, mdims);
    esize = H5Tget_size (mtype);

    herr = count_data (nreq, ndim, counts, &total_blocks);
    CHECK_HERR
    if (fp->dataset_segments.size () <= fp->multi_datasets.size ()) {
        fp->dataset_segments.resize (fp->multi_datasets.size ());
    }
    fp->dataset_segments[fp->multi_datasets.size ()].resize (total_blocks);

    index = 0;
    for (i = 0; i < nreq; ++i) {
        for (j = 0; j < ndim; ++j) {
            start[j] = (hsize_t)starts[i][j];
            block[j] = (hsize_t)counts[i][j];
        }

        herr = pack_data (fp->dataset_segments[fp->multi_datasets.size ()].data (), &index, bufp,
                          esize, ndim, dims, start, block);
        CHECK_HERR
        ts = MPI_Wtime ();
        if (i) {
            herr = H5Sselect_hyperslab (dsid, H5S_SELECT_OR, start, NULL, this->one, block);
        } else {
            herr = H5Sselect_hyperslab (dsid, H5S_SELECT_SET, start, NULL, this->one, block);
        }
        CHECK_HERR
        this->tsel += MPI_Wtime () - ts;
        rsize = 1;
        for (j = 0; j < ndim; ++j) { rsize *= block[j]; }
        total_mem_size += rsize;
        bufp += rsize * esize;
    }

    ts = MPI_Wtime ();
    qsort (fp->dataset_segments[fp->multi_datasets.size ()].data (), total_blocks,
           sizeof (Index_order), index_order_cmp);
    this->tsort += MPI_Wtime () - ts;

    msid = H5Screate_simple (1, &total_mem_size, &total_mem_size);
    CHECK_HID (msid)

    // dataset_size is incremented in the following function.

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
    if ((mode == indep) || (mode == coll)
#ifdef ENABLE_LOGVOL
        || this->use_logvol
#endif
    ) {
        herr = H5Dread (did, mpi_type_to_hdf5_type (type), msid, dsid, dxplid, buf);
        CHECK_HERR
    } else {  // Otherwier, queue request in driver
        herr = fp->register_multidataset (buf, did, dsid, msid, mpi_type_to_hdf5_type (type), 0);
        CHECK_HERR
        // Prevent freeing of dsid and msid, they will be freed after flush
        dsid = msid = -1;
    }

err_out:;
    if (dsid >= 0) H5Sclose (dsid);
    if (msid >= 0) H5Sclose (msid);
    return nerrs;
}