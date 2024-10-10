/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * See README.md for compile and run instructions.
 *
 *********************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strdup(), strncpy(), strstr() */
#include <unistd.h> /* getopt() */

#include <mpi.h>
#include <pnetcdf.h>

#define ERR {                                                            \
    if (err < 0) {                                                       \
        printf ("Error in %s line %d function %s\n", __FILE__, __LINE__, \
                __func__);                                               \
        nerrs++;                                                         \
        goto err_out;                                                    \
    }                                                                    \
}

int verbose;

/*----< read_io_decomp() >-------------------------------------------------*/
int read_decomp(char  *map_file,
                int    dims[2],
                int   *nreqs,
                int  **offsets,
                int  **lengths)
{
    char name[128];
    int err, nerrs=0, rank, nprocs, ncid, varid;
    int i, *all_nreqs=NULL, dimid, ndims;
    MPI_Offset decomp_nprocs, total_nreqs, start, count;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* open input file that contains I/O decomposition information */
    err = ncmpi_open(MPI_COMM_WORLD, map_file, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    ERR

    /* number of processes used when the decomposition was produced */
    err = ncmpi_inq_dimid(ncid, "decomp_nprocs", &dimid);
    ERR
    err = ncmpi_inq_dimlen(ncid, dimid, &decomp_nprocs);
    ERR

    if (decomp_nprocs != nprocs) {
        printf("Error: Number of MPI processes must be %lld\n", decomp_nprocs);
        err = -1;
        goto err_out;
    }

    /* total number of noncontiguous requests of all processes */
    err = ncmpi_inq_dimid(ncid, "D3.total_nreqs", &dimid);
    ERR
    err = ncmpi_inq_dimlen(ncid, dimid, &total_nreqs);
    ERR

    /* obtain the number of dimensions of this decomposition */
    err = ncmpi_get_att_int(ncid, NC_GLOBAL, "D3.ndims", &ndims);
    ERR
    /* obtain the dimension lengths of this decomposition */
    err = ncmpi_get_att_int(ncid, NC_GLOBAL, "D3.dims", dims);
    ERR

    err = ncmpi_inq_varid(ncid, "D3.nreqs", &varid);
    ERR

    /* read all process's numbers of requests */
    all_nreqs = (int*)malloc(decomp_nprocs * sizeof (int));
    err = ncmpi_get_var_int_all(ncid, varid, all_nreqs);
    ERR

    /* calculate start index in D3.offsets for this process */
    start = 0;
    for (i=0; i<rank; i++) start += all_nreqs[i];

    /* calculate number of requests for this process */
    count = all_nreqs[rank];
    *nreqs = count;

    if (verbose)
        printf ("D3 rank %d: start=%lld count=%lld\n", rank, start, count);

    *offsets = (int*)malloc(*nreqs * sizeof(int));
    *lengths = (int*)malloc(*nreqs * sizeof (int));

    /* read starting offsets of requests into offsets[] */
    err = ncmpi_inq_varid(ncid, "D3.offsets", &varid);
    ERR
    err = ncmpi_get_vara_int_all(ncid, varid, &start, &count, *offsets);
    ERR

    /* read lengths of requests into lengths[] */
    err = ncmpi_inq_varid(ncid, "D3.lengths", &varid);
    ERR
    err = ncmpi_get_vara_int_all(ncid, varid, &start, &count, *lengths);
    ERR

    err = ncmpi_close(ncid);
    ERR

err_out:
    if (all_nreqs != NULL) free(all_nreqs);

    return err;
}

/*----< e3sm_io_write() >-------------------------------------------------*/
int e3sm_io_write(char *out_file,
                  int   nvars,
                  int   dims[2],
                  int   nreqs,
                  int  *offsets,
                  int  *lengths)
{
    int i, j, err, nerrs=0, ncid, cmode, *varids=NULL, dimids[2], amnt=0;
    float *buf=NULL, *bufptr;
    MPI_Offset start[2], count[2];

    cmode = NC_64BIT_DATA | NC_CLOBBER;
    err = ncmpi_create(MPI_COMM_WORLD, out_file, cmode, MPI_INFO_NULL, &ncid); ERR

    err = ncmpi_def_dim(ncid, "Y", dims[0], &dimids[0]); ERR
    err = ncmpi_def_dim(ncid, "X", dims[1], &dimids[1]); ERR

    varids = (int*)malloc(sizeof(int) * nvars);
    for (i=0; i<nvars; i++) {
        char name[64];
        sprintf(name, "var_%d", i);
        err = ncmpi_def_var(ncid, name, NC_FLOAT, 2, dimids, &varids[i]); ERR
    }

    err = ncmpi_enddef(ncid); ERR

    for (amnt=0, j=0; j<nreqs; j++) amnt += lengths[j];
    amnt *= nvars;

    buf = (float*)malloc(sizeof(float) * amnt);

    bufptr = buf;
    for (i=0; i<nvars; i++) {
        for (j=0; j<nreqs; j++) {
            start[0] = offsets[j] / dims[1];
            start[1] = offsets[j] % dims[1];
            count[0] = 1;
            count[1] = lengths[j];

            err = ncmpi_iput_vara_float(ncid, varids[i], start, count, bufptr, NULL); ERR
            bufptr += count[0] * count[1];
        }
    }
    err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); ERR

    err = ncmpi_close(ncid); ERR

err_out:
    if (varids != NULL) free(varids);
    if (buf != NULL) free(buf);

    return nerrs;
}

/*----< usage() >------------------------------------------------------------*/
static void usage (char *argv0) {
    char *help = "Usage: %s [OPTION] FILE\n\
       [-h] Print this help message\n\
       [-v] Verbose mode\n\
       [-n nvars] number of variables\n\
       [-o path] Output file path\n\
       FILE: Name of input file storing data decomposition maps.\n";
    fprintf (stderr, help, argv0);
}

/*----< main() >-------------------------------------------------------------*/
int main (int argc, char **argv) {
    char *out_path=NULL, *map_file=NULL;
    int i, err, nerrs=0, rank, nprocs, nvars=1;
    int nreqs, dims[2], *offsets, *lengths;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    verbose = 0;

    /* command-line arguments */
    while ((i = getopt (argc, argv, "hvn:o:")) != EOF)
        switch (i) {
            case 'v':
                verbose = 1;
                break;
            case 'n':
                nvars = atoi(optarg);
                break;
            case 'o':
                out_path = strdup(optarg);
                break;
            case 'h':
            default:
                if (rank == 0) usage(argv[0]);
                goto err_out;
        }

    if (out_path == NULL) { /* input file is mandatory */
        if (!rank) usage (argv[0]);
        printf("Error: Decomposition file is required\n");
        goto err_out;
    }
    map_file = strdup(argv[optind]);

    /* read request information from decomposition file */
    err = read_decomp(map_file, dims, &nreqs, &offsets, &lengths);
    ERR

    /* write to output file */
    err = e3sm_io_write(out_path, nvars, dims, nreqs, offsets, lengths);
    ERR

err_out:
    free(offsets);
    free(lengths);
    free(out_path);
    free(map_file);
    MPI_Finalize ();

    return (nerrs > 0);
}


