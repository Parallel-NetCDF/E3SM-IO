/* Test program for HDF5 PR # 1859
 * https://github.com/HDFGroup/hdf5/issues/1859
 */

#include "hdf5.h"

#define NREQS 2

int
main(int argc, char *argv[])
{
    herr_t err;
    hid_t   file, dataset, cparms;
    hid_t   dataspace[NREQS], ids[NREQS];
    hid_t   filespace[NREQS], mType[NREQS];
    hsize_t rec_size = 3;
    hsize_t dims[NREQS], maxdims[NREQS], chunk_dims[NREQS];
    int i, fillvalue   = 0;
    int data1[3] = {1, 1, 1};
    int data2[3] = {2, 2, 2};
    const void *buf[NREQS] = {data1, data2};

    MPI_Init(&argc, &argv);

    /* Create a new file */
    hid_t acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(acc_tpl1, MPI_COMM_WORLD, MPI_INFO_NULL);
    file = H5Fcreate("testfile.h5", H5F_ACC_TRUNC, H5P_DEFAULT, acc_tpl1);
    H5Pclose(acc_tpl1);

    /* set dataset dimensions */
    dims[0]    = 0;
    dims[1]    = rec_size;
    maxdims[0] = H5S_UNLIMITED;
    maxdims[1] = rec_size;
    dataspace[0] = H5Screate_simple(2, dims, maxdims);

    /* set chunking */
    cparms = H5Pcreate(H5P_DATASET_CREATE);
    chunk_dims[0] = 1;
    chunk_dims[1] = rec_size;
    err = H5Pset_chunk(cparms, 2, chunk_dims);
    err = H5Pset_fill_value(cparms, H5T_NATIVE_INT, &fillvalue);

    /* Create a new dataset */
    dataset = H5Dcreate2(file, "var", H5T_NATIVE_INT, dataspace[0],
                         H5P_DEFAULT, cparms, H5P_DEFAULT);
    H5Sclose(dataspace[0]);
    H5Pclose(cparms);

    /* HDF5 multidataset APIs require filespace to be updated with the latest
     * dimensions at the moment H5Dwrite_multi is called.
     */
#define FINAL_DIMS
#ifdef FINAL_DIMS
    /* Extend the dataset */
    dims[0] = NREQS;
    H5Dset_extent(dataset, dims);
#endif

    for (i=0; i<NREQS; i++) {
        hsize_t start[NREQS] = {0, 0};
        hsize_t count[NREQS] = {1, 1};
        hsize_t block[NREQS] = {1, rec_size};

#ifndef FINAL_DIMS
        /* Extend the dataset */
        dims[0] = i+1;
        H5Dset_extent(dataset, dims);
#endif

        /* Select a hyperslab */
        filespace[i] = H5Dget_space(dataset);
        start[0] = i;
        H5Sselect_hyperslab(filespace[i], H5S_SELECT_SET, start, NULL, count, block);

        /* Define memory space */
        dataspace[i] = H5Screate_simple(1, &rec_size, NULL);

        ids[i] = dataset;
        mType[i] = H5T_NATIVE_INT;
    }
#ifndef FINAL_DIMS
    for (i=0; i<NREQS; i++)
        H5Sset_extent_simple(filespace[i], 2, dims, maxdims);
#endif

    /* collective I/O */
    hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);

    err = H5Dwrite_multi(NREQS, ids, mType, dataspace, filespace, xfer_plist, buf);

    /* Close/release resources */
    H5Pclose(xfer_plist);
    for (i=0; i<NREQS; i++) {
        H5Sclose(dataspace[i]);
        H5Sclose(filespace[i]);
    }
    H5Dclose(dataset);
    H5Fclose(file);

    MPI_Finalize();
    return 0;
}
