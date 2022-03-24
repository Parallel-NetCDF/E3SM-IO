/*********************************************************************
 *
 * Copyright (C) 2022, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 *********************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <e3sm_io_err.h>
#include <errno.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> /* getopt() */

#include <e3sm_io.h>
#include <e3sm_io_driver.hpp>
#ifdef ENABLE_PNC
#include <e3sm_io_driver_pnc.hpp>
#endif
#ifdef ENABLE_ADIOS2
#include <e3sm_io_driver_adios2.hpp>
#endif
#ifdef ENABLE_HDF5
#include <e3sm_io_driver_hdf5.hpp>
#endif
#ifdef ENABLE_NETCDF4
#include <e3sm_io_driver_nc4.hpp>
#endif

#define MAX_NFILES 6
#define LINE_SIZE  4692802

static int verbose, line_sz, raw_decom;

/*----< intcompare() >------------------------------------------------------*/
/* This subroutine is used in qsort() */
static int intcompare (const void *p1, const void *p2) {
    int i = *((int *)p1);
    int j = *((int *)p2);
    if (i > j) return (1);
    if (i < j) return (-1);
    return (0);
}

/*----< add_decomp() >-------------------------------------------------------*/
static int add_decomp (
    int ncid, const char *infname, int label, e3sm_io_driver *driver, bool &have_decom_dim) {
    char *buf, name[128], *map, *str;
    FILE *fd;
    int i, j, rank, nprocs, dimid, ndims, err = 0;
    int varid[5], *nreqs, *raw_nreqs = NULL, *off, *len, *raw_off;
    int total_nreqs, max_nreqs, min_nreqs, maxlen, minlen, total_raw_nreqs;
    MPI_Offset k, gsize, *dims, start, count, raw_start, raw_count;
    int *dims_C;

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    fd = fopen (infname, "r");
    if (fd == NULL) {
        printf ("Error: fail to open file %s(%s)\n", infname, strerror (errno));
        MPI_Finalize ();
        exit (1);
    }

    /* buffer stores one line read from input file */
    buf = (char *)malloc (line_sz);

    /* header lines:(first 2 lines of the decomposition file), for example
     *     version 2001 npes 43200 ndims 2
     *     777602 72
     * version is 2001
     * number of MPI processes used to generate this decomposition is 43200
     * number of dimensions of decomposed variable is 2
     * 1st dimension is of size 777602 and 2nd is 72(in Fortran order)
     */
    fgets (buf, LINE_SIZE, fd);
    if (verbose) printf ("header:\n\t%s", buf);
    strtok (buf, " ");
    if (strcmp (buf, "version")) {
        printf ("Error at line %d: wrong header of input file %s\n", __LINE__, infname);
        fclose (fd);
        MPI_Finalize ();
        exit (1);
    }

    strtok (NULL, " "); /* skip version number */
    strtok (NULL, " "); /* skip toke "npes" */
    nprocs = atoi (strtok (NULL, " "));
    strtok (NULL, " "); /* token "ndims" */
    ndims = atoi (strtok (NULL, " "));
    dims  = (MPI_Offset *)malloc (ndims * sizeof (MPI_Offset));

    /* get dimension sizes
     * Note the dimensions are in Fortran order in the decomposition file
     * generated from PIO library
     */
    fgets (buf, LINE_SIZE, fd);
    if (verbose) printf ("\t%s", buf);
    dims[0] = atoll (strtok (buf, " "));
    for (i = 1; i < ndims; i++) dims[i] = atoll (strtok (NULL, " "));
    /* Note dims[] is in Fortran order */
    if (verbose) {
        if (ndims == 1)
            printf ("lable D%d: dims = %lld\n", label, dims[0]);
        else if (ndims == 2)
            printf ("lable D%d: dims = %lld %lld\n", label, dims[0], dims[1]);
    }

    /* gsize is total number of elements in the global array */
    gsize = dims[0];
    for (i = 1; i < ndims; i++) gsize *= dims[i];

    /* map is used to check whether the entire array is covered by requests
     * of all processes.
     */
    map = (char *)calloc (gsize, 1);

    /* nreqs[i] is the number of elements accessed by process i */
    nreqs     = (int *)calloc (nprocs, sizeof (int));
    raw_nreqs = (int *)calloc (nprocs, sizeof (int));

    /* decomposition data format:
     *  (process.rank.ID)(number.of.requests)
     *     a list of element offsets accessed by this process(one element each)
     * Offsets are indices of array flattened into 1D and start with 1, i.e.
     * Fortran index based. Note the offsets are not sorted in an increasing
     * order and may contain 0s which should be ignored.
     */
    for (rank = 0; rank < nprocs; rank++) {
        int ncontig = 0, decomp_rank;
        char *ret;

        /* reads the first line of rank ID and no. requests */
        while (NULL != (ret = fgets (buf, LINE_SIZE, fd))) {
            if (buf[0] != '\n') break; /* non-empty line */
        }
        if (ret == NULL || strncmp (buf, "Obtained", 8) == 0) {
            /* there is no request for remaining ranks */
            for (decomp_rank = rank; decomp_rank < nprocs; decomp_rank++) nreqs[rank] = 0;
            break; /* loop of rank */
        }

        decomp_rank = atoi (strtok (buf, " ")); /* rank ID */
        while (rank < decomp_rank)              /* this rank has no request */
            nreqs[rank++] = 0;

        nreqs[rank] = atoi (strtok (NULL, " ")); /* number of requests */
        if (nreqs[rank] == 0)                    /* this rank has zero request */
            continue;                            /* loop of rank */

        /* Record number of raw offsets before it is merged */
        raw_nreqs[rank] = nreqs[rank];

        off = (int *)malloc (nreqs[rank] * sizeof (int));
        fgets (buf, LINE_SIZE + 1, fd); /* 2nd line: list of offsets */
        if (strlen (buf) >= LINE_SIZE) {
            printf ("Error: line size is larger than default %d\n", LINE_SIZE);
            printf ("       use command-line option -l to use a larger size\n");
            goto err_out;
        }

        /* construct the offset list */
        off[0] = atoi (strtok (buf, " "));
        j = 1;
        while (off[0] == 0) { /* skip leading 0 values, if there is any */
            off[0] = atoi (strtok (NULL, " "));
            j++;
        }
        off[0]--;
        k = 1;
        for (; j < nreqs[rank]; j++) {
            off[k] = atoi (strtok (NULL, " "));
            if (off[k] == 0) continue; /* skip 0 values */
            off[k]--;                  /* offset is 1 based */
            k++;
        }

        /* sort off[] into an increasing order */
        qsort ((void *)off, k, sizeof (int), intcompare);

        ncontig = 1;
        for (j = 1; j < k; j++) {
            /* break contiguity at dimension boundaries or noncontiguous */
            if (off[j] % dims[0] == 0 || off[j] > off[j - 1] + 1) ncontig++;
        }
        nreqs[rank] = ncontig;

        for (j = 0; j < k; j++) map[off[j]] = 1;
        free (off);
    }

    /* find total, max, and min nreqs amount all processes */
    total_nreqs = max_nreqs = min_nreqs = nreqs[0];
    total_raw_nreqs = raw_nreqs[0];
    for (i = 1; i < nprocs; i++) {
        total_nreqs += nreqs[i];
        total_raw_nreqs += raw_nreqs[i];
        max_nreqs = (nreqs[i] > max_nreqs) ? nreqs[i] : max_nreqs;
        min_nreqs = (nreqs[i] < min_nreqs) ? nreqs[i] : min_nreqs;
    }

    if (verbose) printf ("max_nreqs=%d min_nreqs=%d\n", max_nreqs, min_nreqs);

    /* check if the entire array is covered */
    for (k = 0, j = 0; j < gsize; j++) k += map[j];
    if (k != gsize) {
        printf ("Warning: decomposition %d does not cover the entire array\n", label);
        printf ("\tglobal %dD array size: %lld", ndims, dims[0]);
        for (i = 1; i < ndims; i++) printf (" x %lld", dims[i]);
        printf (" = %lld\n", gsize);
        printf ("\tthe decomposition map covers only %lld of them\n", k);
    }
    free (map);

    /* check if dimension decomp_nprocs has been defined in the netCDF file */
    if (have_decom_dim) {
        /* if decomp_nprocs already exist, check if value matches */
        MPI_Offset decomp_nprocs;
        err = driver->inq_dim (ncid, "decomp_nprocs", &dimid);
        CHECK_ERR
        err = driver->inq_dimlen (ncid, dimid, &decomp_nprocs);
        CHECK_ERR
        if (decomp_nprocs != nprocs) {
            printf ("Error: decomp_nprocs=%lld mismatches among input files %d\n", decomp_nprocs,
                    nprocs);
            MPI_Finalize ();
            exit (1);
        }
    } else {
        err = driver->def_dim (ncid, "decomp_nprocs", nprocs, &dimid);
        CHECK_ERR
        have_decom_dim = true;
    }

    /* define variable nreqs for this decomposition */
    sprintf (name, "D%d.nreqs", label);
    err = driver->def_var (ncid, name, NC_INT, 1, &dimid, &varid[0]);
    CHECK_ERR

    /* define attribute description for this variable */
    str = (char *)"Number of noncontiguous requests per process";
    err = driver->put_att (ncid, varid[0], "description", NC_CHAR, strlen (str), str);
    CHECK_ERR

    /* define variable raw_nreqs for this decomposition */
    if (raw_decom) {
        sprintf (name, "D%d.raw_nreqs", label);
        err = driver->def_var (ncid, name, NC_INT, 1, &dimid, &varid[3]);
        CHECK_ERR

        /* define attribute description for this variable */
        str = (char *)"Number of file offsets accessed per process before merge";
        err = driver->put_att (ncid, varid[3], "description", NC_CHAR, strlen (str), str);
        CHECK_ERR
    }

    /* define dimension total_nreqs for this decomposition */
    sprintf (name, "D%d.total_nreqs", label);
    err = driver->def_dim (ncid, name, total_nreqs, &dimid);
    CHECK_ERR

    /* define variable offsets(store starting element indices) */
    sprintf (name, "D%d.offsets", label);
    err = driver->def_var (ncid, name, NC_INT, 1, &dimid, &varid[1]);
    CHECK_ERR
    str = (char *)"Flattened starting indices of noncontiguous requests";
    err = driver->put_att (ncid, varid[1], "description", NC_CHAR, strlen (str), str);
    CHECK_ERR

    /* define variable lengths(store number of elements) */
    sprintf (name, "D%d.lengths", label);
    err = driver->def_var (ncid, name, NC_INT, 1, &dimid, &varid[2]);
    CHECK_ERR
    str = (char *)"Lengths of noncontiguous requests";
    err = driver->put_att (ncid, varid[2], "description", NC_CHAR, strlen (str), str);
    CHECK_ERR

    maxlen = minlen = 0;
    /*
    err = driver->put_att (ncid, varid[2], "max", NC_INT, 1, &maxlen);
    CHECK_ERR
    err = driver->put_att (ncid, varid[2], "min", NC_INT, 1, &minlen);
    CHECK_ERR
    */

    if (raw_decom) {
        /* define dimension total_raw_nreqs for this decomposition */
        sprintf (name, "D%d.total_raw_nreqs", label);
        err = driver->def_dim (ncid, name, total_raw_nreqs, &dimid);
        CHECK_ERR

        /* define variable raw_offsets */
        sprintf (name, "D%d.raw_offsets", label);
        err = driver->def_var (ncid, name, NC_INT, 1, &dimid, &varid[4]);
        CHECK_ERR
        str = (char *)"File offsets accessed before merge";
        err = driver->put_att (ncid, varid[4], "description", NC_CHAR, strlen (str), str);
        CHECK_ERR
    }

    /* add attribute to describe dimensionality */
    sprintf (name, "D%d.ndims", label);
    err = driver->put_att (ncid, NC_GLOBAL, name, NC_INT, 1, &ndims);
    CHECK_ERR

    /* swap dims in Fortran order to dims_C in C order */
    dims_C = (int *)malloc (ndims * sizeof (int));
    for (i = 0; i < ndims; i++) dims_C[i] = (int)(dims[ndims - i - 1]);
    sprintf (name, "D%d.dims", label);
    err = driver->put_att (ncid, NC_GLOBAL, name, NC_INT, ndims, dims_C);
    CHECK_ERR
    free (dims_C);

    sprintf (name, "D%d.max_nreqs", label);
    err = driver->put_att (ncid, NC_GLOBAL, name, NC_INT, 1, &max_nreqs);
    CHECK_ERR
    sprintf (name, "D%d.min_nreqs", label);
    err = driver->put_att (ncid, NC_GLOBAL, name, NC_INT, 1, &min_nreqs);
    CHECK_ERR

    /* exit define mode */
    err = driver->enddef (ncid);
    CHECK_ERR

    /* write variable containing number of requests for each process */
    err = driver->put_vara (ncid, varid[0], MPI_INT, NULL, NULL, nreqs, coll);
    CHECK_ERR

    /* write variable containing number of requests before merge for each process */
    if (raw_decom) {
        err = driver->put_vara (ncid, varid[3], MPI_INT, NULL, NULL, raw_nreqs, coll);
        CHECK_ERR
    }

    /* read the offsets again into allocated array off */
    start = 0;
    raw_start = 0;
    rewind (fd);
    fgets (buf, LINE_SIZE, fd);
    fgets (buf, LINE_SIZE, fd);
    for (rank = 0; rank < nprocs; rank++) {
        int prev, ncontig = 0, decomp_rank;
        char *ret;

        /* reads the first line of rank ID and no. requests */
        while (NULL != (ret = fgets (buf, LINE_SIZE, fd))) {
            if (buf[0] != '\n') break; /* non-empty line */
        }
        if (ret == NULL || strncmp (buf, "Obtained", 8) == 0) {
            /* there is no request for remaining ranks */
            for (decomp_rank = rank; decomp_rank < nprocs; decomp_rank++) nreqs[rank] = 0;
            break; /* loop of rank */
        }

        decomp_rank = atoi (strtok (buf, " ")); /* rank ID */
        while (rank < decomp_rank)                /* this rank has no request */
            rank++;

        nreqs[rank] = atoi (strtok (NULL, " ")); /* number of requests */
        if (nreqs[rank] == 0)                     /* this rank has zero request */
            continue;                             /* loop of rank */

        off = (int *)malloc (nreqs[rank] * sizeof (int));
        len = (int *)malloc (nreqs[rank] * sizeof (int));
        raw_off = (int *)malloc (raw_nreqs[rank] * sizeof (int));
        fgets (buf, LINE_SIZE, fd);
        i = 0;
        off[0] = raw_off[i++] = atoi (strtok (buf, " "));
        j = 1;
        while (off[0] == 0) {
            off[0] = raw_off[i++] = atoi (strtok (NULL, " "));
            j++;
        }
        off[0]--;
        k = 1;
        for (; j < nreqs[rank]; j++) {
            off[k] = raw_off[i++] = atoi (strtok (NULL, " "));
            if (off[k] == 0) continue; /* skip 0 values */
            off[k]--;                   /* offset is 1 based */
            k++;
        }

        /* sort off[] into an increasing order */
        qsort ((void *)off, k, sizeof (int), intcompare);

        ncontig = 1;
        prev    = 0;
        len[0]  = 1;
        for (j = 1; j < k; j++) {
            /* break contiguity at dimension boundaries or noncontiguous */
            if (off[j] % dims[0] == 0 || off[j] > off[j - 1] + 1) ncontig++;

            if (off[j] % dims[0] == 0 || off[j] > off[prev] + len[prev]) {
                prev++;
                if (prev < j) off[prev] = off[j];
                len[prev] = 1;
            } else
                len[prev]++;
        }
        assert (prev + 1 == ncontig);

        if (rank == 0) maxlen = minlen = len[0];
        for (j = 0; j < ncontig; j++) {
            maxlen = (len[j] > maxlen) ? len[j] : maxlen;
            minlen = (len[j] < minlen) ? len[j] : minlen;
        }

        /* write/append to variables offsets and lengths */
        count = ncontig;
        err = driver->put_vara (ncid, varid[1], MPI_INT, &start, &count, off, coll);
        CHECK_ERR
        err = driver->put_vara (ncid, varid[2], MPI_INT, &start, &count, len, coll);
        CHECK_ERR
        start += ncontig;

        /* write/append raw file offsets before merge */
        if (raw_decom) {
            raw_count = raw_nreqs[rank];
            err = driver->put_vara (ncid, varid[4], MPI_INT, &raw_start, &raw_count, raw_off, coll);
            CHECK_ERR
            raw_start += raw_count;
        }

        free (off);
        free (len);
        free (raw_off);
    }
    free (dims);

    err = driver->put_att (ncid, varid[2], "max", NC_INT, 1, &maxlen);
    CHECK_ERR
    err = driver->put_att (ncid, varid[2], "min", NC_INT, 1, &minlen);
    CHECK_ERR

err_out:
    fclose (fd);
    free (nreqs);
    free (raw_nreqs);
    free (buf);

    return err;
}

/*----< extract_file_names() >-----------------------------------------------*/
static void extract_file_names (const char *inList, int *num_decomp, char **infname) {
    FILE *fd;
    int i, num_files;
    char line[1024];

    for (i = 0; i < MAX_NFILES; i++) infname[i] = NULL;

    fd = fopen (inList, "r");
    if (fd == NULL) {
        printf ("Error: open fails on file %s(%s)\n", inList, strerror (errno));
        exit (1);
    }

    /* count number of input files */
    num_files = 0;
    while (fgets (line, 1024, fd)) {
        if (strlen (line) == 0) continue; /* skip empty lines */
        if (line[0] == '#') continue;      /* skip comment line(start with #) */
        num_files++;
    }

    /* read input file names */
    rewind (fd);
    i = 0;
    while (fgets (line, 1024, fd)) {
        char *tail;
        if (strlen (line) == 0) continue; /* skip empty lines */
        if (line[0] == '#') continue;      /* skip comment line(start with #) */
        /* remove blanks at tail. Note fgets stores newline to the buffer */
        tail = line + strlen (line) - 1;
        while (*tail == ' ' || *tail == '\t' || *tail == '\n') tail--;
        tail[1] = '\0';
        /* save file name to in_list */
        infname[i] = strdup (line);
        i++;
    }
    assert (i == num_files);
    fclose (fd);

    *num_decomp = num_files;
}

static void usage (char *argv0) {
    const char *help =
    "Usage: ./dat2decom [-h|-v|-r|-l num] -a fmt -i input_file -o out_file\n"
    "  [-h]           Print help\n"
    "  [-v]           Verbose mode\n"
    "  [-r]           Include original decomposition map\n"
    "  [-l num]       max number of characters per line in input file\n"
    "  -a fmt         output file format, fmt is one of the followings\n"
    "     cdf5:       NetCDF classic 64-bit data format\n"
    "     netcdf4:    NetCDF-4 (HDF5-based) format\n"
    "     hdf5:       HDF5 format\n"
    "     bp:         ADIOS2 BP format\n"
    "  -i input_file  input file containing a list of decomposition .dat file\n"
    "  -o out_file    output file name\n";
    fprintf (stderr, help, argv0);
}

/*----< main() >------------------------------------------------------------*/
int main (int argc, char **argv) {
    char *inList = NULL, *infname[MAX_NFILES], *outfname = NULL, cmd_line[4096];
    int i, rank, ncid, num_decomp = 0, dimid, err = 0;
    bool have_decom_dim = false;
    MPI_Info info;
    e3sm_io_api api = pnetcdf;
    e3sm_io_driver *driver = NULL;
    e3sm_io_config cfg;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    cmd_line[0] = '\0';
    for (i = 0; i < argc; i++) {
        strcat (cmd_line, argv[i]);
        strcat (cmd_line, " ");
    }

    for (i = 0; i < MAX_NFILES; i++) infname[i] = NULL;
    line_sz   = LINE_SIZE;
    verbose   = 0;
    raw_decom = 0;

    /* get command-line arguments */
    while ((i = getopt (argc, argv, "hvro:l:i:a:")) != EOF) {
        switch (i) {
            case 'a':;
                if (std::string (optarg) == "cdf5") {
                    api = pnetcdf;
                } else if (std::string (optarg) == "hdf5") {
                    api = hdf5;
                } else if (std::string (optarg) == "bp") {
                    api = adios;
                } else if (std::string (optarg) == "netcdf4") {
                    api = netcdf4;
                } else {
                    api = undef_api;
                }
                break;
            case 'v':
                verbose = 1;
                break;
            case 'o':
                outfname = strdup (optarg);
                break;
            case 'l':
                line_sz = atoi (optarg);
                break;
            case 'r':
                raw_decom = 1;
                break;
            case 'i':
                inList = strdup (optarg);
                break;
            case 'h':
            default:
                if (rank == 0) usage (argv[0]);
                err = -1;
                goto err_out;
        }
    }

    if (inList == NULL) { /* input file name is mandatory */
        if (rank == 0) {
            printf ("input file is missing\n");
            usage (argv[0]);
        }
        MPI_Finalize ();
        return 1;
    }
    if (outfname == NULL) { /* output file name is mandatory */
        if (rank == 0) {
            printf ("output file is missing\n");
            usage (argv[0]);
        }
        MPI_Finalize ();
        return 1;
    }

    extract_file_names (inList, &num_decomp, infname);

    if (verbose && rank == 0) {
        printf ("input  file: %s\n", inList);
        printf ("output file: %s\n", outfname);
        printf ("Number of decomposition files: %d\n", num_decomp);
        for (i = 0; i < num_decomp; i++) printf ("decomposition file %d: %s\n", i, infname[i]);
    }

    if (num_decomp == 0) {
        if (rank == 0) printf ("Error: number of input decomposition files is zero.\n");
        MPI_Finalize ();
        return 1;
    }

    if (num_decomp != 3 && num_decomp != 6 && num_decomp != 5) {
        if (rank == 0) {
            printf ("Error: e3sm_io currently supports 3 case studies:\n");
            printf ("       F case(3 decomposition files),\n");
            printf ("       G case(6 decomposition files), and\n");
            printf ("       I case(5 decomposition files).\n");
        }
        MPI_Finalize ();
        return 1;
    }

    // Set up dummy config for the driver
    cfg.io_comm        = MPI_COMM_WORLD;
    cfg.info           = MPI_INFO_NULL;
    cfg.num_iotasks    = cfg.np;
    cfg.num_group      = 1;
    cfg.out_path[0]    = '\0';
    cfg.in_path[0]     = '\0';
    cfg.cfg_path[0]    = '\0';
    cfg.hx             = -1;
    cfg.wr             = 0;
    cfg.rd             = 0;
    cfg.nvars          = 0;
    cfg.strategy       = undef_io;
    cfg.api            = api;
    cfg.chunksize      = 0;
    cfg.filter         = none;
    cfg.verbose        = 0;
    cfg.keep_outfile   = 0;
    cfg.profiling      = 0;
    cfg.two_buf        = 0;
    cfg.non_contig_buf = 0;
    cfg.io_stride      = 1;
    cfg.sub_comm       = MPI_COMM_NULL;
    cfg.rank           = rank;

    // Initialize driver
    driver = e3sm_io_get_driver (NULL, &cfg);
    CHECK_PTR (driver)

    err = MPI_Info_create (&info);
    CHECK_MPIERR
    err = MPI_Info_set (info, "nc_header_align_size", "8192");
    CHECK_MPIERR

    /* create a new NC file */
    err = driver->create (outfname, MPI_COMM_WORLD, info, &ncid);
    CHECK_ERR
    err = MPI_Info_free (&info);
    CHECK_MPIERR

    /* add the number of decompositions */
    err = driver->def_dim (ncid, "num_decomp", num_decomp, &dimid);
    CHECK_ERR

    /* add command line used */
    err = driver->put_att (ncid, NC_GLOBAL, "command_line", NC_CHAR, strlen (cmd_line), cmd_line);
    CHECK_ERR

    for (i = 0; i < num_decomp; i++) {
        err = add_decomp (ncid, infname[i], i + 1, driver, have_decom_dim);
        CHECK_ERR;
        err = driver->redef (ncid);
        CHECK_ERR;
    }

    err = driver->close (ncid);
    CHECK_ERR

err_out:
    if (driver) { delete driver; }
    for (i = 0; i < MAX_NFILES; i++)
        if (infname[i] != NULL) free (infname[i]);
    free (outfname);
    free (inList);

    MPI_Finalize ();
    return err;
}
