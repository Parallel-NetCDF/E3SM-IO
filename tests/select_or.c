#include <stdlib.h>
#include <assert.h>
#include <hdf5.h>

/*
 * H5S_SELECT_OR will flatten each hyperslab into offset-length pairs and
 * merge/sort all pairs of all hyperslabs into an increasing order of offsets.
 * Thus, there is no way in HDF5 to create a filespace that describes offsets
 * in an non-increasing order. This example program creates 2 filespaces. Both
 * calls H5Sselect_hyperslab() twice. The 1st filespace sets the left column
 * block of an 2D array and union it with the right column block. The 2nd
 * filespace is create in a reverse order.
 *
 */

#define RANK   2
#define NY     3
#define NX     10

int main (int argc, char **argv)
{
    int     i, j, *buf;
    hid_t   file_id, dset1, dset2, filespace1, filespace2;
    hsize_t dims[2], count[2], start[2];
    herr_t  err;

    file_id = H5Fcreate("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert(file_id >= 0);

    dims[0] = NY;
    dims[1] = NX;
    filespace1 = H5Screate_simple(RANK, dims, NULL); 
    assert(filespace1 >= 0);

    dset1 = H5Dcreate(file_id, "var1", H5T_NATIVE_INT, filespace1,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(dset1 >= 0);
    dset2 = H5Dcreate(file_id, "var2", H5T_NATIVE_INT, filespace1,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(dset2 >= 0);
    err = H5Sclose(filespace1);
    assert(err >= 0);

    buf = (int*) malloc(sizeof(int)*dims[0]*dims[1]);
    for (i=0; i<dims[0]; i++)
        for (j=0; j<dims[1]; j++)
            buf[i*dims[1] + j] = i*dims[1] + j + 10;

    start[0] = 0;
    count[0] = dims[0];
    count[1] = dims[1]/2;

    /* create a union of 2 column-wise hyperspaces */
    filespace1 = H5Dget_space(dset1);
    assert(filespace1 >= 0);

    /* right column block first */
    start[1] = dims[1]/2;
    err = H5Sselect_hyperslab(filespace1, H5S_SELECT_SET, start, NULL, count, NULL);
    assert(err >= 0);

    /* left column block second */
    start[1] = 0;
    err = H5Sselect_hyperslab(filespace1, H5S_SELECT_OR, start, NULL, count, NULL);
    assert(err >= 0);

    err = H5Dwrite(dset1, H5T_NATIVE_INT, H5S_BLOCK, filespace1, H5P_DEFAULT, buf);
    assert(err >= 0);

    err = H5Dclose(dset1);
    assert(err >= 0);

    /* create a union of 2 column-wise hyperspaces */
    filespace2 = H5Dget_space(dset2);
    assert(filespace2 >= 0);

    /* left column block first */
    start[1] = 0;
    err = H5Sselect_hyperslab(filespace2, H5S_SELECT_SET, start, NULL, count, NULL);
    assert(err >= 0);

    /* right column block second */
    start[1] = dims[1]/2;
    err = H5Sselect_hyperslab(filespace2, H5S_SELECT_OR, start, NULL, count, NULL);
    assert(err >= 0);

    err = H5Dwrite(dset2, H5T_NATIVE_INT, H5S_BLOCK, filespace2, H5P_DEFAULT, buf);
    assert(err >= 0);

    printf("whether two file spaces are equal? %d\n",
           H5Sextent_equal(filespace1, filespace2));

    free(buf);
    err = H5Sclose(filespace1);
    assert(err >= 0);
    err = H5Sclose(filespace2);
    assert(err >= 0);
    H5Fclose(file_id);
 
    return 0;
}
