/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Copyright (C) 2024, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * E3SM-IO F case on Perlmutter:
 * 1. total write amount = 28 GB
 * 2. number of MPI processes = 21600
 * 3. number of Lustre OSTs = 128
 * 4. Lustre stripe size = 1 MB
 * 5. Two-phase I/O takes 28GB/(128*1MB)=132 rounds
 * 6. In each round, each process sends 1MB/21600=48.5B to each aggregator
 *
 * To mimic the two-phase I/O in E3SM-IO, run with command-line options:
 *    -n 132 -l 12 -r 128
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <mpi.h>

#define ERR \
    if (err != MPI_SUCCESS) { \
        int errorStringLen; \
        char errorString[MPI_MAX_ERROR_STRING]; \
        MPI_Error_string(err, errorString, &errorStringLen); \
        printf("Error at line %d: %s\n",__LINE__,errorString); \
        nerrs++; \
        goto err_out; \
    }

/*----< usage() >------------------------------------------------------------*/
static void usage (char *argv0) {
    char *help = "Usage: %s [OPTION]\n\
       [-h] Print this help message\n\
       [-v] Verbose mode (default: no)\n\
       [-n num] number of iterations (default: 1)\n\
       [-l len] individual message size (default: 12 int)\n\
       [-r ratio] number of processes to aggregators ratio (default: 1)\n";
    fprintf (stderr, help, argv0);
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv) {
    extern int optind;
    extern char *optarg;
    int i, j, rank, nprocs, err, nerrs=0, verbose, len, ntimes, ratio;
    int is_aggr, num_aggrs, *aggr_ranks, *buf;
    MPI_Request *reqs;
    MPI_Status *st;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    verbose = 0;
    len = 12;
    ntimes = 1;
    ratio = 1;

    /* command-line arguments */
    while ((i = getopt (argc, argv, "hvl:n:r:")) != EOF)
        switch (i) {
            case 'v':
                verbose = 1;
                break;
            case 'l':
                len = atoi(optarg);
                break;
            case 'n':
                ntimes = atoi(optarg);
                break;
            case 'r':
                ratio = atoi(optarg);
                break;
            case 'h':
            default:
                if (rank == 0) usage(argv[0]);
                goto err_out;
        }

    is_aggr = 0;
    num_aggrs = nprocs / ratio;

    if (verbose && rank == 0) {
        printf("nprocs    = %d\n", nprocs);
        printf("len       = %d\n", len);
        printf("ntimes    = %d\n", ntimes);
        printf("ratio     = %d\n", ratio);
        printf("num_aggrs = %d\n", num_aggrs);
    }

    aggr_ranks = (int*) malloc(sizeof(int) * num_aggrs);
    if (verbose && rank == 0) printf("aggr_ranks: ");
    for (i=0; i<num_aggrs; i++) {
        aggr_ranks[i] = i * ratio;
        if (rank == aggr_ranks[i]) is_aggr = 1;
        if (verbose && rank == 0) printf(" %d", aggr_ranks[i]);
    }
    if (verbose && rank == 0) printf("\n");

    buf = (int*) malloc(sizeof(int) * (nprocs + num_aggrs) * len);

    reqs = (MPI_Request*) calloc(nprocs + num_aggrs, sizeof(MPI_Request));
    st = (MPI_Status*)malloc(sizeof(MPI_Status) * (nprocs + num_aggrs));

    for (i=0; i<ntimes; i++) {
        int nreqs=0;
        int *ptr = buf;

        /* post recv requests */
        if (is_aggr) {
            for (j=0; j<nprocs; j++) {
                err = MPI_Irecv(ptr, len, MPI_INT, j, 0, MPI_COMM_WORLD, &reqs[nreqs++]);
                ERR
                ptr += len;
            }
        }

        /* post send requests */
        for (j=0; j<num_aggrs; j++) {
            err = MPI_Issend(ptr, len, MPI_INT, aggr_ranks[j], 0, MPI_COMM_WORLD, &reqs[nreqs++]);
            ERR
            ptr += len;
        }

        err = MPI_Waitall(nreqs, reqs, st);
        ERR
    }

    free(st);
    free(reqs);
    free(buf);
    free(aggr_ranks);

err_out:
    MPI_Finalize();
    return 0;
}


