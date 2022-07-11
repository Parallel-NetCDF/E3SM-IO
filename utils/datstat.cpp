/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program uses the E3SM I/O patterns recorded by the PIO library to
 * evaluate the performance of high-level I/O libraries.
 * The E3SM I/O patterns consist of a large number of small,
 * noncontiguous requests on each MPI process, which presents a challenge for
 * achieving a good performance.
 *
 * See README.md for compile and run instructions.
 *
 *********************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
/**/
#include <assert.h>
#include <errno.h>
#include <libgen.h> /* basename() */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>
/**/
#include <unistd.h> /* getopt() */
/**/
#include <mpi.h>
/**/
#include <e3sm_io_err.h>

#define CHECK_RET(A)                                          \
    {                                                         \
        if (ret != A) { ERR_OUT ("Unexpected return value") } \
    }

/*----< add_decomp() >------------------------------------------------------------*/
static int add_decomp (const char *infname) {
    int err = 0;
    int ret;
    int i, j;
    int ver, np, ndim;
    int rank, ndecom;
    int max_ndecom, max_ndecom_rank;
    long long max_zero;
    int max_zero_rank;
    long long total_ndecom, total_zero, total_zero_rank;
    long long off;
    std::vector<long long> dsizes;
    FILE *fd;

    printf ("=========================================================================\n");
    printf ("Deomposition file %s\n", infname);
    printf ("=========================================================================\n");

    fd = fopen (infname, "r");
    if (fd == NULL) {
        printf ("Error: fail to open file %s (%s)\n", infname, strerror (errno));
        MPI_Finalize ();
        exit (1);
    }

    ret = fscanf (fd, "version %d npes %d ndims %d", &ver, &np, &ndim);
    CHECK_RET (3)

    printf ("version: %d\nnpes: %d\nndims: %d\n", ver, np, ndim);
    CHECK_RET (3)

    dsizes.resize (ndim);
    for (i = 0; i < ndim; i++) {
        ret = fscanf (fd, "%lld", &(dsizes[i]));
        CHECK_RET (1)
    }

    printf ("Dims: ");
    for (i = 0; i < ndim; i++) { printf ("%lld, ", dsizes[i]); }
    printf ("\n");

    total_ndecom = total_zero = max_ndecom = max_zero = 0;
    max_ndecom_rank = max_zero_rank = 0;
    for (i = 0; i < np; i++) {
        ret = fscanf (fd, "%d %d", &rank, &ndecom);
        CHECK_RET (2)

        total_ndecom += ndecom;
        if (max_ndecom < ndecom) {
            max_ndecom      = ndecom;
            max_ndecom_rank = rank;
        }

        total_zero_rank = 0;
        for (j = 0; j < ndecom; j++) {
            ret = fscanf (fd, "%lld", &off);
            CHECK_RET (1)

            if (off == 0) { total_zero_rank++; }
        }

        if (total_zero_rank > max_zero) {
            max_zero      = total_zero_rank;
            max_zero_rank = rank;
        }
        total_zero += total_zero_rank;
    }

    printf ("total_cells: %lld\n", total_ndecom);
    printf ("total_zero_cells: %lld\n", total_zero);
    printf ("max_cells: %d @ rank %d\n", max_ndecom, max_ndecom_rank);
    printf ("max_zero_cells: %lld @ rank %d\n", max_zero, max_zero_rank);
    printf ("zero ratio: %lf\n", (double)total_zero / (double)total_ndecom);

err_out:
    fclose (fd);
    return err;
}

static void usage (char *argv0) {
    char *help =(char*)
    "Usage: %s [OPTION]...\n"
    "       -h             Print help\n"
    "       -d input_file  decomposition map .dat file to be analyzed\n";
    fprintf (stderr, help, argv0);
}

/*----< main() >------------------------------------------------------------*/
int main (int argc, char **argv) {
    int err, i;
    std::vector<char *> infname;

    /* get command-line arguments */
    while ((i = getopt (argc, argv, "hd:")) != EOF) switch (i) {
            case 'd':
                infname.push_back (optarg);
                break;
            case 'h':
            default:
                usage (argv[0]);
                return 1;
        }

    if (infname.size () == 0) {
        printf ("Error: input decomposition files are required.\n");
        printf ("       add -h to see all available command-line options.\n");
        return 1;
    }

    for (auto &fname : infname) { err = add_decomp (fname); }

    return (err < 0) ? 1 : 0;
}
