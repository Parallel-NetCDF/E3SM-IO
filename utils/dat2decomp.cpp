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
#include <string.h> /* strdup() */
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

// #define REDUCED_RATIO 100

static int verbose, line_sz, raw_decom, fill_gaps, sort_off;

typedef struct {
    int off;
    int len;
} off_len;

/*----< off_len_compare() >--------------------------------------------------*/
/* This subroutine is used in qsort() */
static
int off_len_compare(const void *p1, const void *p2) {
    int off1 = ((off_len*)p1)->off;
    int off2 = ((off_len*)p2)->off;
    if (off1 > off2) return (1);
    if (off1 < off2) return (-1);
    return (0);
}

/*----< intcompare() >------------------------------------------------------*/
/* This subroutine is used in qsort() */
static
int intcompare(const void *p1, const void *p2) {
    int i = *((int *)p1);
    int j = *((int *)p2);
    if (i > j) return (1);
    if (i < j) return (-1);
    return (0);
}

#define TRUE "true"
#define FALSE "false"

/*----< add_decomp() >-------------------------------------------------------*/
static
int add_decomp(int             ncid,
               const char     *infname,
               int             label,
               e3sm_io_driver *driver,
               bool           &have_decom_dim)
{
    char *buf, name[128], *map, *str;
    FILE *fd;
    int i, j, rank, nprocs, dimid, ndims, dimX, ngaps, cur, err=0;
    int varid[6], *nreqs, **off, **len, *fill_starts, dim_nprocs;
    int total_nreqs, max_nreqs, min_nreqs, maxlen, minlen;
    MPI_Offset k, gsize, *dims, start, count;
    int *raw_nreqs, **raw_off, total_raw_nreqs;
    MPI_Offset raw_start;
    int *dims_C;

    fd = fopen(infname, "r");
    if (fd == NULL) {
        printf("Error: fail to open file %s(%s)\n", infname, strerror(errno));
        MPI_Finalize();
        exit(1);
    }

    /* buffer stores one line read from input file */
    buf = (char*) malloc(line_sz);

    /* header lines:(first 2 lines of the decomposition file), for example
     *     version 2001 npes 43200 ndims 2
     *     777602 72
     * version is 2001
     * number of MPI processes used to generate this decomposition is 43200
     * number of dimensions of decomposed variable is 2
     * 1st dimension is of size 777602 and 2nd is 72(in Fortran order)
     */
    fgets(buf, LINE_SIZE, fd);
    if (verbose) printf("header:\n\t%s", buf);
    strtok(buf, " ");
    if (strcmp(buf, "version")) {
        printf("Error at line %d: wrong header of input file %s\n",
               __LINE__, infname);
        fclose(fd);
        MPI_Finalize();
        exit(1);
    }

    strtok(NULL, " "); /* skip version number */
    strtok(NULL, " "); /* skip toke "npes" */
    nprocs = atoi(strtok(NULL, " "));
    strtok(NULL, " "); /* token "ndims" */
    ndims = atoi(strtok(NULL, " "));
    dims = (MPI_Offset*) malloc(ndims * sizeof(MPI_Offset));

    /* get dimension sizes
     * Note the dimensions are in Fortran order in the decomposition file
     * generated from PIO library
     */
    fgets(buf, LINE_SIZE, fd);
    if (verbose) printf("\t%s", buf);
    dims[0] = atoll(strtok(buf, " "));
    for (i=1; i<ndims; i++) dims[i] = atoll(strtok(NULL, " "));

#ifdef REDUCED_RATIO
    printf("label D%d: dim[%d] is reduced from %lld to %lld\n",
           label, ndims-1, dims[ndims-1], dims[ndims-1]/REDUCED_RATIO);
    dims[ndims-1] /= REDUCED_RATIO;
#endif

    /* Note dims[] read from the PIO decomposition file are in Fortran order */
    if (verbose) {
        if (ndims == 1)
            printf("label D%d: dims = %lld\n", label, dims[0]);
        else if (ndims == 2)
            printf("label D%d: dims = %lld x %lld (in C order)\n", label, dims[1], dims[0]);
        else if (ndims == 3)
            printf("label D%d: dims = %lld x %lld x %lld (in C order)\n", label, dims[2], dims[1], dims[0]);
    }
    dimX = dims[0]; /* the least significant dimension */

    /* gsize is total number of elements in the global array */
    gsize = dims[0];
    for (i=1; i<ndims; i++) gsize *= dims[i];

    /* map is used to check whether the entire array is covered by requests
     * of all processes.
     */
    map = (char*) calloc(gsize, 1);

    /* nreqs[i] is the number of elements accessed by process i */
    nreqs = (int*)  malloc(nprocs * sizeof(int));
    off   = (int**) malloc(nprocs * sizeof(int*));
    len   = (int**) malloc(nprocs * sizeof(int*));

    /* start index in nreqs[] for requests to be written with fill values */
    fill_starts = (int*) malloc(nprocs * sizeof(int));

    /* decomposition data format:
     *  (process.rank.ID)(number.of.requests)
     *     a list of element offsets accessed by this process(one element each)
     * Offsets are indices of array flattened into 1D and start with 1, i.e.
     * Fortran index based. Note the offsets are not sorted in an increasing
     * order and may contain 0s which should be ignored.
     */
    for (rank=0; rank<nprocs; rank++) {
        int prev, decomp_rank;
        char *ret;

        /* reads the first line of rank ID and no. requests */
        while (NULL != (ret = fgets(buf, LINE_SIZE, fd))) {
            if (buf[0] != '\n') break; /* non-empty line */
        }
        if (ret == NULL || strncmp(buf, "Obtained", 8) == 0) {
            /* there is no request for remaining ranks */
            for (decomp_rank=rank; decomp_rank<nprocs; decomp_rank++)
                nreqs[rank] = 0;
            break; /* loop of rank */
        }

        decomp_rank = atoi(strtok(buf, " ")); /* rank ID */
        while (rank < decomp_rank)            /* this rank has no request */
            nreqs[rank++] = 0;

        nreqs[rank] = atoi(strtok(NULL, " ")); /* number of requests */
        if (nreqs[rank] == 0) /* this rank has zero request */
            continue;         /* loop of rank */

        /* read and store offsets */
        off[rank] = (int *)malloc(nreqs[rank] * sizeof(int));
        len[rank] = (int *)malloc(nreqs[rank] * sizeof(int));

        fgets(buf, LINE_SIZE + 1, fd); /* 2nd line: list of offsets */
        if (strlen(buf) >= LINE_SIZE) {
            printf("Error: line size is larger than default %d\n", LINE_SIZE);
            printf("       use command-line option -l to use a larger size\n");
            goto err_out;
        }

        /* construct the offset list:
         * The offsets stored in the dat files are in Fortran 1-based index.
         * There are 0s that need to be removed.
         */
        off[rank][0] = atoi(strtok(buf, " "));
        j = 1;
        while (off[rank][0] == 0 || off[rank][0] > gsize) {
            /* skip leading 0 values and values > gsize, if there is any */
            off[rank][0] = atoi(strtok(NULL, " "));
            j++;
        }
        off[rank][0]--;
        k = 1;
        for (; j<nreqs[rank]; j++) {
            off[rank][k] = atoi(strtok(NULL, " "));
            if (off[rank][k] == 0 || off[rank][k] > gsize)
                continue;   /* skip 0 values and > gsize */
            off[rank][k]--; /* offset is 1 based */
            k++;
        }
        /* k now is the number of non-zero offsets */

        if (sort_off)
            /* sort off[rank][] into an increasing order */
            qsort((void *)off[rank], k, sizeof(int), intcompare);

        /* build a map for checking if the decompositions cover all elements */
        for (j=0; j<k; j++) map[off[rank][j]] = 1;

        /* coalescing contiguous offsets (must break boundary at dimension X) */
        prev = 0;
        len[rank][0] = 1;
        for (j=1; j<k; j++) {
            if (off[rank][j] != off[rank][prev] + len[rank][prev] || /* not contiguous from previous offset-length pair */
                off[rank][j] % dimX == 0) { /* break at dimension boundaries */
                prev++;
                if (prev < j) off[rank][prev] = off[rank][j];
                len[rank][prev] = 1;
            } else
                len[rank][prev]++;
        }
        /* set nreqs[] to the number of offset-length pairs */
        nreqs[rank] = prev+1;

        /* find max and min contiguous length among all offset-length pairs and
         * among all processes
         */
        if (rank == 0) maxlen = minlen = len[rank][0];
        for (j=0; j<nreqs[rank]; j++) {
            maxlen = (len[rank][j] > maxlen) ? len[rank][j] : maxlen;
            minlen = (len[rank][j] < minlen) ? len[rank][j] : minlen;
        }
    }

    /* start index for requests to be written with fill values */
    for (rank=0; rank<nprocs; rank++) fill_starts[rank] = nreqs[rank];

    /* calculate number of gaps (unconverted consecutive blocks of elements) */
    ngaps = 0;
    cur = 1;
    for (k=0, j=0; j<gsize; j++) {
        k += map[j];
        if (map[j] == 0) {
            if (cur == 1) {
                ngaps++;
                cur = 0;
            }
            /* Note dims[] is in Fortran order */
            else if (j % dimX == 0) /* end of dim X */
                ngaps++;
        }
        else cur = 1;
    }

    /* check if the entire array is covered */
    if (k != gsize) {
        printf("Warning: decomposition %d does not cover the entire array\n",
               label);
        printf("\tglobal %dD array size: %lld", ndims, dims[ndims-1]);
        for (i=ndims-2; i>=0; i--) printf(" x %lld", dims[i]);
        printf(" = %lld\n", gsize);
        printf("\tthe decomposition map covers only %lld of them\n", k);
        printf("\tnumber of missing elements are:   %lld\n", gsize-k);
        printf("\tnumber of gaps:                   %d\n", ngaps);
        if (fill_gaps) printf("\tOption to fill in gaps is enabled\n");
    }
    if (fill_gaps && k != gsize) {
        /* assign gaps to all processes */
        int *fill_nreqs;
        fill_nreqs = (int*) malloc(nprocs * sizeof(int));

        /* divide number of gaps evenly among all processes */
        int rem = ngaps / nprocs;
        for (rank=0; rank<nprocs; rank++) {
            int nelems = (rank < ngaps % nprocs) ? rem + 1 : rem;
            fill_nreqs[rank] = nelems;
            /* extend requests */
            off[rank] = (int*) realloc(off[rank], (nreqs[rank] + nelems) * sizeof(int));
            len[rank] = (int*) realloc(len[rank], (nreqs[rank] + nelems) * sizeof(int));
        }

        /* assign contiguous gaps to the same process */
        rank = 0;
        int prev = 1;
        int ncontig = nreqs[rank] - 1; /* add at the end */
        for (j=0; j<gsize; j++) {
            if (map[j] == 0) {
                /* Note dims[] is in Fortran order */
                if (prev == 1 || j % dimX == 0) { /* end of dim X */
                    ncontig++;
                    if (ncontig == nreqs[rank] + fill_nreqs[rank]) {
                        rank++;
                        ncontig = nreqs[rank]; /* add at the end */
                    }
                    assert(rank <= nprocs);
                    off[rank][ncontig] = j;
                    len[rank][ncontig] = 1;
                    prev = 0;
                }
                else
                    len[rank][ncontig]++;
            }
            else prev = 1;
        }

        for (rank=0; rank<nprocs; rank++) {
            nreqs[rank] += fill_nreqs[rank];

            /* sort off-len pairs into an increasing order */
            if (sort_off) {
                off_len *pairs = (off_len*) malloc(nreqs[rank] * sizeof(off_len));
                for (i=0; i<nreqs[rank]; i++) {
                    pairs[i].off = off[rank][i];
                    pairs[i].len = len[rank][i];
                }
                qsort((void*)pairs, nreqs[rank], sizeof(off_len), off_len_compare);
                for (i=0; i<nreqs[rank]; i++) {
                    off[rank][i] = pairs[i].off;
                    len[rank][i] = pairs[i].len;
                }
                free(pairs);
            }
        }
        free(fill_nreqs);
    }
    free(map);

    if (raw_decom) {
        /* populate raw_off */
        raw_nreqs = (int*)  malloc(nprocs * sizeof(int));
        raw_off   = (int**) malloc(nprocs * sizeof(int*));
        for (rank=0; rank<nprocs; rank++) {
            raw_nreqs[rank] = 0;
            for (i=0; i<fill_starts[rank]; i++)
                raw_nreqs[rank] += len[rank][i];
            raw_off[rank] = (int*) malloc(raw_nreqs[rank] * sizeof(int));
            k = 0;
            for (i=0; i<fill_starts[rank]; i++)
                for (j=0; j<len[rank][i]; j++)
                    /* raw offsets are Fortran based, starting from 1 */
                    raw_off[rank][k++] = off[rank][i] + j + 1;
        }
        total_raw_nreqs = raw_nreqs[0];
        for (rank=1; rank<nprocs; rank++)
            total_raw_nreqs += raw_nreqs[rank];
    }

    /* find total, max, and min nreqs amount all processes */
    total_nreqs = max_nreqs = min_nreqs = nreqs[0];
    for (i=1; i<nprocs; i++) {
        total_nreqs += nreqs[i];
        max_nreqs = (nreqs[i] > max_nreqs) ? nreqs[i] : max_nreqs;
        min_nreqs = (nreqs[i] < min_nreqs) ? nreqs[i] : min_nreqs;
    }

    if (verbose) printf("total_nreqs=%d max_nreqs=%d min_nreqs=%d\n", total_nreqs, max_nreqs, min_nreqs);

    /* check if dimension decomp_nprocs has been defined in the file */
    if (have_decom_dim) {
        /* if decomp_nprocs already exist, check if value matches */
        MPI_Offset decomp_nprocs;
        err = driver->inq_dim(ncid, "decomp_nprocs", &dim_nprocs);
        CHECK_ERR
        err = driver->inq_dimlen(ncid, dim_nprocs, &decomp_nprocs);
        CHECK_ERR
        if (decomp_nprocs != nprocs) {
            printf("Error: decomp_nprocs=%lld mismatches among input files %d\n",
                   decomp_nprocs, nprocs);
            MPI_Finalize();
            exit(1);
        }
    } else {
        err = driver->def_dim(ncid, "decomp_nprocs", nprocs, &dim_nprocs);
        CHECK_ERR
        have_decom_dim = true;
    }

    /* define variable nreqs for this decomposition */
    sprintf(name, "D%d.nreqs", label);
    err = driver->def_var(ncid, name, NC_INT, 1, &dim_nprocs, &varid[0]);
    CHECK_ERR

    /* add an attribute to describe this variable */
    str = (char *)"Number of noncontiguous requests per process";
    err = driver->put_att(ncid, varid[0], "description", NC_CHAR, strlen(str), str);
    CHECK_ERR

    /* define variable fill_starts for this decomposition */
    sprintf(name, "D%d.fill_starts", label);
    err = driver->def_var(ncid, name, NC_INT, 1, &dim_nprocs, &varid[1]);
    CHECK_ERR

    /* add an attribute to describe this variable */
    str = (char*)"Start index in offsets[] and lengths[] for requests to be written with fill values";
    err = driver->put_att(ncid, varid[1], "description", NC_CHAR, strlen(str), str);
    CHECK_ERR

    /* define dimension total_nreqs for this decomposition */
    sprintf(name, "D%d.total_nreqs", label);
    err = driver->def_dim(ncid, name, total_nreqs, &dimid);
    CHECK_ERR

    /* define variable offsets(store starting element indices) */
    sprintf(name, "D%d.offsets", label);
    err = driver->def_var(ncid, name, NC_INT, 1, &dimid, &varid[2]);
    CHECK_ERR
    str = (char *)"Flattened starting indices of noncontiguous requests";
    err = driver->put_att(ncid, varid[2], "description", NC_CHAR, strlen(str), str);
    CHECK_ERR

    /* define variable lengths(store number of elements) */
    sprintf(name, "D%d.lengths", label);
    err = driver->def_var(ncid, name, NC_INT, 1, &dimid, &varid[3]);
    CHECK_ERR
    str = (char *)"Lengths of noncontiguous requests";
    err = driver->put_att(ncid, varid[3], "description", NC_CHAR, strlen(str), str);
    CHECK_ERR

    err = driver->put_att(ncid, varid[3], "max", NC_INT, 1, &maxlen);
    CHECK_ERR
    err = driver->put_att(ncid, varid[3], "min", NC_INT, 1, &minlen);
    CHECK_ERR

    if (raw_decom) {
        /* define variable raw_nreqs for this decomposition */
        sprintf(name, "D%d.raw_nreqs", label);
        err = driver->def_var(ncid, name, NC_INT, 1, &dim_nprocs, &varid[4]);
        CHECK_ERR

        /* define attribute description for this variable */
        str = (char *)"Number of file offsets accessed per process before merge";
        err = driver->put_att(ncid, varid[4], "description", NC_CHAR, strlen(str), str);
        CHECK_ERR

        /* define dimension total_raw_nreqs for this decomposition */
        sprintf(name, "D%d.total_raw_nreqs", label);
        err = driver->def_dim(ncid, name, total_raw_nreqs, &dimid);
        CHECK_ERR

        /* define variable raw_offsets */
        sprintf(name, "D%d.raw_offsets", label);
        err = driver->def_var(ncid, name, NC_INT, 1, &dimid, &varid[5]);
        CHECK_ERR
        str = (char *)"File offsets accessed before merge";
        err = driver->put_att(ncid, varid[5], "description", NC_CHAR, strlen(str), str);
        CHECK_ERR
    }

    /* add attribute to describe dimensionality */
    sprintf(name, "D%d.ndims", label);
    err = driver->put_att(ncid, NC_GLOBAL, name, NC_INT, 1, &ndims);
    CHECK_ERR

    /* swap dims in Fortran order to dims_C in C order */
    dims_C = (int*) malloc(ndims * sizeof(int));
    for (i=0; i<ndims; i++) dims_C[i] = (int)(dims[ndims - i - 1]);
    sprintf(name, "D%d.dims", label);
    err = driver->put_att(ncid, NC_GLOBAL, name, NC_INT, ndims, dims_C);
    CHECK_ERR
    free(dims_C);
    free(dims);

    sprintf(name, "D%d.max_nreqs", label);
    err = driver->put_att(ncid, NC_GLOBAL, name, NC_INT, 1, &max_nreqs);
    CHECK_ERR
    sprintf(name, "D%d.min_nreqs", label);
    err = driver->put_att(ncid, NC_GLOBAL, name, NC_INT, 1, &min_nreqs);
    CHECK_ERR

    sprintf(name, "D%d.fully_covered", label);
    if (ngaps)
        err = driver->put_att(ncid, NC_GLOBAL, name, NC_CHAR, strlen(FALSE), FALSE);
    else
        err = driver->put_att(ncid, NC_GLOBAL, name, NC_CHAR, strlen(TRUE), TRUE);
    CHECK_ERR

    /* exit define mode */
    err = driver->enddef(ncid);
    CHECK_ERR

    /* write variable containing number of requests for each process */
    err = driver->put_vara(ncid, varid[0], MPI_INT, NULL, NULL, nreqs, coll);
    CHECK_ERR

    /* write variable containing number of requests for each process */
    err = driver->put_vara(ncid, varid[1], MPI_INT, NULL, NULL, fill_starts, coll);
    CHECK_ERR

    /* write variable containing number of requests before coalescing */
    if (raw_decom) {
        err = driver->put_vara(ncid, varid[4], MPI_INT, NULL, NULL, raw_nreqs, coll);
        CHECK_ERR
    }

    /* read the offsets again into allocated array off */
    start = 0;
    raw_start = 0;
    for (rank=0; rank<nprocs; rank++) {
        /* write/append to variables offsets and lengths */
        count = nreqs[rank];
        err = driver->put_vara(ncid, varid[2], MPI_INT, &start, &count, off[rank], coll);
        CHECK_ERR
        err = driver->put_vara(ncid, varid[3], MPI_INT, &start, &count, len[rank], coll);
        CHECK_ERR
        start += count;

        /* write/append raw file offsets before merge */
        if (raw_decom) {
            count = raw_nreqs[rank];
            err = driver->put_vara(ncid, varid[5], MPI_INT, &raw_start, &count, raw_off[rank], coll);
            CHECK_ERR
            raw_start += count;
        }
    }

err_out:
    fclose(fd);

    free(nreqs);
    free(fill_starts);
    for (rank=0; rank<nprocs; rank++) {
        free(off[rank]);
        free(len[rank]);
    }
    free(off);
    free(len);
    if (raw_decom) {
        for (rank=0; rank<nprocs; rank++)
            free(raw_off[rank]);
        free(raw_off);
        free(raw_nreqs);
    }
    free(buf);

    return err;
}

/*----< extract_file_names() >-----------------------------------------------*/
static
void extract_file_names(const char  *inList,
                        int         *num_decomp,
                        char       **infname)
{
    FILE *fd;
    int i, num_files;
    char line[1024];

    for (i=0; i<MAX_NFILES; i++) infname[i] = NULL;

    fd = fopen(inList, "r");
    if (fd == NULL) {
        printf("Error: open fails on file %s(%s)\n", inList, strerror(errno));
        exit(1);
    }

    /* count number of input files */
    num_files = 0;
    while (fgets(line, 1024, fd)) {
        char *ptr = line;
        while (*ptr == ' ' || *ptr == '\t') ptr++; /* skip blank space */
        if (*ptr == '\n') continue; /* skip empty lines */
        if (*ptr == '#')  continue; /* skip comment line(start with #) */
        /* remove blanks at tail. Note fgets stores newline to the buffer */
        while (ptr[strlen(ptr)-1] == '\n' || ptr[strlen(ptr)-1] == ' ')
            ptr[strlen(ptr)-1] = '\0';
        /* save file name to in_list */
        infname[num_files++] = strdup(ptr);
    }
    fclose(fd);

    *num_decomp = num_files;
}

static void usage(char *argv0) {
    const char *help =
    "Usage: %s [-h|-v|-r|-f|-s|-l num] -a fmt -i input_file -o out_file\n"
    "  -h             Print help\n"
    "  -v             Verbose mode\n"
    "  -r             Include original decomposition maps\n"
    "  -f             Fill in unassigned elements in decomposition maps\n"
    "  -s             Sort offsets of decomposition maps increasingly\n"
    "  -l num         max number of characters per line in input file\n"
    "  -a fmt         output file format, fmt is one of the followings\n"
    "     cdf5:       NetCDF classic 64-bit data format\n"
    "     netcdf4:    NetCDF-4 (HDF5-based) format\n"
    "     hdf5:       HDF5 format\n"
    "     bp:         ADIOS2 BP format\n"
    "  -i input_file  a text file containing a list of decomposition\n"
    "                 map .dat file names\n"
    "  -o out_file    output file name\n";
    fprintf(stderr, help, argv0);
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv) {
    char *inList = NULL, *infname[MAX_NFILES], *outfname = NULL, cmd_line[4096];
    int i, rank, nprocs, ncid, num_decomp=0, dimid, err=0;
    bool have_decom_dim = false;
    MPI_Info info;
    e3sm_io_api api = pnetcdf;
    e3sm_io_driver *driver = NULL;
    e3sm_io_config cfg;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    for (i=0; i<MAX_NFILES; i++) infname[i] = NULL;

    if (nprocs > 1) {
        nprocs = 1;
        comm = MPI_COMM_SELF;
        if (rank == 0)
            printf("Warning: %s is for sequential run. Run on 1 process now.\n",
                   argv[0]);
        else
            goto err_out;
    }

    cmd_line[0] = '\0';
    for (i=0; i<argc; i++) {
        strcat(cmd_line, argv[i]);
        strcat(cmd_line, " ");
    }

    line_sz   = LINE_SIZE;
    verbose   = 0;
    raw_decom = 0;
    fill_gaps = 0;
    sort_off  = 0;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hvrfso:l:i:a:")) != EOF) {
        switch (i) {
            case 'a':;
                if (std::string(optarg) == "cdf5") {
                    api = pnetcdf;
                } else if (std::string(optarg) == "hdf5") {
                    api = hdf5;
                } else if (std::string(optarg) == "bp") {
                    api = adios;
                } else if (std::string(optarg) == "netcdf4") {
                    api = netcdf4;
                } else {
                    api = undef_api;
                }
                break;
            case 'v': verbose = 1;
                      break;
            case 'o': outfname = strdup(optarg);
                      break;
            case 'l': line_sz = atoi(optarg);
                      break;
            case 'f': fill_gaps = 1;
                      break;
            case 'r': raw_decom = 1;
                      break;
            case 's': sort_off = 1;
                      break;
            case 'i': inList = strdup(optarg);
                      break;
            case 'h':
            default:
                if (rank == 0) usage(argv[0]);
                MPI_Finalize();
                return 1;
        }
    }

    if (inList == NULL) { /* input file name is mandatory */
        if (rank == 0) {
            printf("input file is missing\n");
            usage(argv[0]);
        }
        MPI_Finalize();
        return 1;
    }
    if (outfname == NULL) { /* output file name is mandatory */
        if (rank == 0) {
            printf("output file is missing\n");
            usage(argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    extract_file_names(inList, &num_decomp, infname);

    if (verbose && rank == 0) {
        printf("input  file: %s\n", inList);
        printf("output file: %s\n", outfname);
        printf("Number of decomposition files: %d\n", num_decomp);
        for (i=0; i<num_decomp; i++)
            printf("decomposition file %d: %s\n", i, infname[i]);
    }

    if (num_decomp == 0) {
        if (rank == 0)
            printf("Error: number of input decomposition files is zero.\n");
        MPI_Finalize();
        return 1;
    }

    if (num_decomp != 3 && num_decomp != 6 && num_decomp != 5) {
        if (rank == 0) {
            printf("Error: e3sm_io currently supports 3 case studies:\n");
            printf("       F case(3 decomposition files),\n");
            printf("       G case(6 decomposition files), and\n");
            printf("       I case(5 decomposition files).\n");
        }
        MPI_Finalize();
        return 1;
    }

    // Set up dummy config for the driver
    cfg.io_comm        = comm;
    cfg.info           = MPI_INFO_NULL;
    cfg.num_iotasks    = 1;
    cfg.num_subfiles   = 0;
    cfg.out_path[0]    = '\0';
    cfg.in_path[0]     = '\0';
    cfg.decomp_path[0] = '\0';
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
    driver = e3sm_io_get_driver(NULL, &cfg);
    CHECK_PTR(driver)

    err = MPI_Info_create(&info);
    CHECK_MPIERR
    err = MPI_Info_set(info, "nc_header_align_size", "8192");
    CHECK_MPIERR

    /* create a new NC file */
    err = driver->create(outfname, comm, info, &ncid);
    CHECK_ERR
    err = MPI_Info_free(&info);
    CHECK_MPIERR

    /* add the number of decompositions */
    err = driver->def_dim(ncid, "num_decomp", num_decomp, &dimid);
    CHECK_ERR

    /* add command line used */
    err = driver->put_att(ncid, NC_GLOBAL, "command_line", NC_CHAR, strlen(cmd_line), cmd_line);
    CHECK_ERR

    for (i=0; i<num_decomp; i++) {
        err = add_decomp(ncid, infname[i], i+1, driver, have_decom_dim);
        CHECK_ERR
        err = driver->redef(ncid);
        CHECK_ERR
    }

    err = driver->close(ncid);
    CHECK_ERR

err_out:
    if (driver) { delete driver; }
    for (i=0; i<MAX_NFILES; i++)
        if (infname[i] != NULL) free(infname[i]);
    if (outfname != NULL) free(outfname);
    if (inList != NULL) free(inList);

    MPI_Finalize();
    return err;
}
