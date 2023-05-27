/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM-IO benchmark software package.
 *
 * See README.md for compile and run instructions.
 *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy() */
#include <unistd.h> /* getopt(), access() */

#include <assert.h>

#include <mpi.h>
#include <pnetcdf.h>

#ifndef MAX
#define MAX(a, b) ((a) > (b)) ? (a) : (b)
#endif
#ifndef MIN
#define MIN(a, b) ((a) < (b)) ? (a) : (b)
#endif

#define CHECK_NC_ERR { \
    if (err != NC_NOERR) { \
        fprintf(stderr, "Error at line %d: %s\n", __LINE__, \
                ncmpi_strerrno(err)); \
        goto err_out; \
    } \
}

/*----< usage() >------------------------------------------------------------*/
static void usage (char *argv0) {
    char *help =
    "Usage: %s [OPTION]... FILE\n"
    "    [-h] Print help\n"
    "    FILE decomposition map file name\n";
    fprintf (stderr, help, argv0);
}

/*----< main() >-------------------------------------------------------------*/
int main (int argc, char **argv)
{
    char name[64], filename[1024];
    int i, j, k, n, err=0, rank, nprocs, ncid, dimid, varid;
    int len, *nreqs, *wr_amount;
    MPI_Offset num_decomp, decomp_nprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* command-line arguments */
    while ((i = getopt(argc, argv, "h")) != EOF)
        switch(i) {
            case 'h': if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }

    if (optind >= argc || argv[optind] == NULL) { /* input file is mandatory */
        if (rank==0) usage(argv[0]);
        MPI_Finalize();
        return 1;
    }
    strcpy(filename, argv[optind]);

    if (filename[0] == '\0') { /* input file name is mandatory */
        if (!rank) {
            fprintf(stderr, "Error: input file is missing\n");
            usage(argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    /* open input file */
    err = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    CHECK_NC_ERR

    err = ncmpi_inq_dimid(ncid, "num_decomp", &dimid);
    CHECK_NC_ERR
    err = ncmpi_inq_dimlen(ncid, dimid, &num_decomp);
    CHECK_NC_ERR
    printf("Number of decomposition maps = %lld\n", num_decomp);

    err = ncmpi_inq_dimid(ncid, "decomp_nprocs", &dimid);
    CHECK_NC_ERR
    err = ncmpi_inq_dimlen(ncid, dimid, &decomp_nprocs);
    CHECK_NC_ERR
    printf("Number of processes = %lld\n", decomp_nprocs);

    nreqs = (int*) malloc(decomp_nprocs * sizeof(int));
    wr_amount = (int*) malloc(decomp_nprocs * sizeof(int));

    for (i=0; i<num_decomp; i++) {

        printf("Map %d\n", i);

        MPI_Offset ndims;
        int *dims;
        sprintf(name, "D%d.dims", i+1);
        err = ncmpi_inq_attlen(ncid, NC_GLOBAL, name, &ndims);
        CHECK_NC_ERR
        printf("\tD%d.ndims           = %d\n", i+1, ndims);
        dims = (int*) malloc(ndims * sizeof(int));
        err = ncmpi_get_att_int(ncid, NC_GLOBAL, name, dims);
        CHECK_NC_ERR
        printf("\tD%d.dims            = %d", i+1, dims[0]);
        len = dims[0];
        for (j=1; j<ndims; j++) {
            len *= dims[j];        
            printf(" x %d", dims[j]);
        }
        printf("\n");
        printf("\tD%d.dim size        = %d\n", i+1, len);
        free(dims);

        sprintf(name, "D%d.nreqs", i+1);
        err = ncmpi_inq_varid(ncid, name, &varid);
        CHECK_NC_ERR
        err = ncmpi_get_var_int_all(ncid, varid, nreqs);
        CHECK_NC_ERR

        MPI_Offset total_nreqs;
        sprintf(name, "D%d.total_nreqs", i+1);
        err = ncmpi_inq_dimid(ncid, name, &dimid);
        CHECK_NC_ERR
        err = ncmpi_inq_dimlen(ncid, dimid, &total_nreqs);
        CHECK_NC_ERR
        printf("\tD%d.total_nreqs     = %d\n", i+1, total_nreqs);

        int *lengths;
        lengths = (int*) malloc(total_nreqs * sizeof(int));
        sprintf(name, "D%d.lengths", i+1);
        err = ncmpi_inq_varid(ncid, name, &varid);
        CHECK_NC_ERR
        err = ncmpi_get_var_int_all(ncid, varid, lengths);
        CHECK_NC_ERR

        for (j=0; j<decomp_nprocs; j++) wr_amount[j] = 0;

        int max_req=0;
        n = 0;
        j = 0;
        for (k=0; k<total_nreqs; k++) {
            max_req = MAX(max_req, lengths[k]);
            wr_amount[j] += lengths[k];
            n++;
            if (n == nreqs[j]) {
                j++;
                n = 0;
            }
        }

        int max_amnt=wr_amount[0];
        for (j=1; j<decomp_nprocs; j++)
            max_amnt = MAX(max_amnt, wr_amount[j]);

        printf("\tD%d.lengths MAX req = %d\n", i+1, max_req);
        printf("\tD%d.lengths MAX sum = %d\n", i+1, max_amnt);

        free(lengths);
    }
    free(nreqs);
    free(wr_amount);

err_out:
    err = ncmpi_close(ncid);
    CHECK_NC_ERR

    MPI_Finalize();

    return (err != 0);
}
