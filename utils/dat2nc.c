/*********************************************************************
 *
 * Copyright(C) 2018, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program converts a PIO decomposition file(text) to a NetCDF file,
 * to be used by e3sm_io.c to benchmark PnetCDF performance.
 */

#include <assert.h>
#include <errno.h>
#include <libgen.h> /* basename() */
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strdup() */
#include <unistd.h> /* getopt() */

#include <mpi.h>
#include <pnetcdf.h>

#define MAX_NFILES 6
#define LINE_SIZE  4692802

// #define REDUCED_RATIO 100

#define ERR { \
    if (err != NC_NOERR) { \
        printf("Error at line %d in %s: %s\n", __LINE__, __FILE__, \
               ncmpi_strerrno(err)); \
        goto fn_exit; \
    } \
}

#define CHECK_MPIERR {                                                       \
    if (err != MPI_SUCCESS) {                                                \
        int el = 256;                                                        \
        char errstr[256];                                                    \
        MPI_Error_string (err, errstr, &el);                                 \
        printf ("Error in %s line %d function %s: %s\n", __FILE__, __LINE__, \
                __func__, errstr);                                           \
        err = -1;                                                            \
        goto fn_exit;                                                        \
    }                                                                        \
}

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
int add_decomp(int ncid, const char *infname, int label) {
    char *buf, name[128], *str;
    FILE *fd;
    int i, j, rank, nprocs, dimid, ndims, dimX, ngaps, cur, err=NC_NOERR;
    int varid[6], *nreqs, **off, **len, *fill_starts, dim_nprocs;
    int *map, total_nreqs, max_nreqs, min_nreqs, maxlen, minlen;
    MPI_Offset k, gsize, *dims, *dims_C, start, count;
    int *raw_nreqs, **raw_off, total_raw_nreqs;
    MPI_Offset raw_start;

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
    map = (int*) malloc(gsize * sizeof(int));
    for (i=0; i<gsize; i++) map[i] = -1;

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
            goto fn_exit;
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
        for (j=0; j<k; j++) map[off[rank][j]] = rank;

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
        if (verbose) printf("%2d: decomposition %d nreqs=%d\n", rank, label, nreqs[rank]);

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
        k += (map[j] >= 0) ? 1 : 0;
        if (map[j] == -1) {
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
            if (map[j] == -1) {
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
            if (verbose) printf("%2d: after fille: decomposition %d nreqs=%d\n", rank, label, nreqs[rank]);

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

    /* check if dimension decomp_nprocs has been defined in the netCDF file */
    err = ncmpi_inq_dimid(ncid, "decomp_nprocs", &dim_nprocs);
    if (err == NC_EBADDIM) { /* not defined */
        err = ncmpi_def_dim(ncid, "decomp_nprocs", nprocs, &dim_nprocs);
        ERR
    } else {
        /* if decomp_nprocs already exist, check if value matches */
        MPI_Offset decomp_nprocs;
        err = ncmpi_inq_dimlen(ncid, dim_nprocs, &decomp_nprocs);
        ERR
        if (decomp_nprocs != nprocs) {
            printf("Error: decomp_nprocs=%lld mismatches among input files %d\n",
                   decomp_nprocs, nprocs);
            MPI_Finalize();
            exit(1);
        }
    }

    /* define variable nreqs for this decomposition */
    sprintf(name, "D%d.nreqs", label);
    err = ncmpi_def_var(ncid, name, NC_INT, 1, &dim_nprocs, &varid[0]);
    ERR

    /* add an attribute to describe this variable */
    str = "Number of noncontiguous requests per process";
    err = ncmpi_put_att_text(ncid, varid[0], "description", strlen(str), str);
    ERR

    /* define variable fill_starts for this decomposition */
    sprintf(name, "D%d.fill_starts", label);
    err = ncmpi_def_var(ncid, name, NC_INT, 1, &dim_nprocs, &varid[1]);
    ERR

    /* add an attribute to describe this variable */
    str = "Start index in offsets[] and lengths[] for requests to be written with fill values";
    err = ncmpi_put_att_text(ncid, varid[1], "description", strlen(str), str);
    ERR

    /* define dimension total_nreqs for this decomposition */
    sprintf(name, "D%d.total_nreqs", label);
    err = ncmpi_def_dim(ncid, name, total_nreqs, &dimid);
    ERR

    /* define variable offsets(store starting element indices) */
    sprintf(name, "D%d.offsets", label);
    err = ncmpi_def_var(ncid, name, NC_INT, 1, &dimid, &varid[2]);
    ERR
    str = "Flattened starting indices of noncontiguous requests";
    err = ncmpi_put_att_text(ncid, varid[2], "description", strlen(str), str);
    ERR

    /* define variable lengths(store number of elements) */
    sprintf(name, "D%d.lengths", label);
    err = ncmpi_def_var(ncid, name, NC_INT, 1, &dimid, &varid[3]);
    ERR
    str = "Lengths of noncontiguous requests";
    err = ncmpi_put_att_text(ncid, varid[3], "description", strlen(str), str);
    ERR

    err = ncmpi_put_att_int(ncid, varid[3], "max", NC_INT, 1, &maxlen);
    ERR
    err = ncmpi_put_att_int(ncid, varid[3], "min", NC_INT, 1, &minlen);
    ERR

    if (raw_decom) {
        /* define variable raw_nreqs for this decomposition */
        sprintf(name, "D%d.raw_nreqs", label);
        err = ncmpi_def_var(ncid, name, NC_INT, 1, &dim_nprocs, &varid[4]);
        ERR

        /* define attribute description for this variable */
        str = "Number of file offsets accessed per process before merge";
        err = ncmpi_put_att_text(ncid, varid[4], "description", strlen(str), str);
        ERR

        /* define dimension total_raw_nreqs for this decomposition */
        sprintf(name, "D%d.total_raw_nreqs", label);
        err = ncmpi_def_dim(ncid, name, total_raw_nreqs, &dimid);
        ERR

        /* define variable raw_offsets */
        sprintf(name, "D%d.raw_offsets", label);
        err = ncmpi_def_var(ncid, name, NC_INT, 1, &dimid, &varid[5]);
        ERR
        str = "File offsets accessed before merge";
        err = ncmpi_put_att_text(ncid, varid[5], "description", strlen(str), str);
        ERR
    }

    /* add attribute to describe dimensionality */
    sprintf(name, "D%d.ndims", label);
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, name, NC_INT, 1, &ndims);
    ERR

    /* swap dims in Fortran order to dims_C in C order */
    dims_C = (MPI_Offset *)malloc(ndims * sizeof(MPI_Offset));
    for (i=0; i<ndims; i++) dims_C[i] = dims[ndims - i - 1];
    sprintf(name, "D%d.dims", label);
    err = ncmpi_put_att_longlong(ncid, NC_GLOBAL, name, NC_INT, ndims, dims_C);
    ERR
    free(dims_C);
    free(dims);

    sprintf(name, "D%d.max_nreqs", label);
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, name, NC_INT, 1, &max_nreqs);
    ERR
    sprintf(name, "D%d.min_nreqs", label);
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, name, NC_INT, 1, &min_nreqs);
    ERR

    sprintf(name, "D%d.fully_covered", label);
    if (ngaps)
        err = ncmpi_put_att_text(ncid, NC_GLOBAL, name, strlen(FALSE), FALSE);
    else
        err = ncmpi_put_att_text(ncid, NC_GLOBAL, name, strlen(TRUE), TRUE);
    ERR

    /* exit define mode */
    err = ncmpi_enddef(ncid);
    ERR

    /* write variable containing number of requests for each process */
    err = ncmpi_put_var_int_all(ncid, varid[0], nreqs);
    ERR

    /* write variable containing number of requests for each process */
    err = ncmpi_put_var_int_all(ncid, varid[1], fill_starts);
    ERR

    /* write variable containing number of requests before coalescing */
    if (raw_decom) {
        err = ncmpi_put_var_int_all(ncid, varid[4], raw_nreqs);
        ERR
    }

    start = 0;
    raw_start = 0;
    for (rank=0; rank<nprocs; rank++) {
        /* write/append to variables offsets and lengths */
        count = nreqs[rank];
        err = ncmpi_put_vara_int_all(ncid, varid[2], &start, &count, off[rank]);
        ERR
        err = ncmpi_put_vara_int_all(ncid, varid[3], &start, &count, len[rank]);
        ERR
        start += count;

        /* write/append raw file offsets before merge */
        if (raw_decom) {
            count = raw_nreqs[rank];
            err = ncmpi_put_vara_int_all(ncid, varid[5], &raw_start, &count, raw_off[rank]);
            ERR
            raw_start += count;
        }
    }

fn_exit:
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
    char *help =
    "Usage: %s [-h|-v|-r|-f|-s|-l num] -i input_file -o out_file\n"
    "  -h             Print help\n"
    "  -v             Verbose mode\n"
    "  -r             Include original decomposition maps\n"
    "  -f             Fill in unassigned elements in decomposition maps\n"
    "  -s             Sort offsets of decomposition maps increasingly\n"
    "  -l num         max number of characters per line in input file\n"
    "  -i input_file  a text file containing a list of decomposition\n"
    "                 map .dat file names\n"
    "  -o out_file    name of output netCDF file\n";
    fprintf(stderr, help, argv0);
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv) {
    char *inList=NULL, *infname[MAX_NFILES], *outfname=NULL, cmd_line[4096];
    int i, rank, ncid, num_decomp=0, dimid, err=NC_NOERR;
    MPI_Info info;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    cmd_line[0] = '\0';
    for (i=0; i<argc; i++) {
        strcat(cmd_line, argv[i]);
        strcat(cmd_line, " ");
    }

    for (i=0; i<MAX_NFILES; i++) infname[i] = NULL;
    line_sz   = LINE_SIZE;
    verbose   = 0;
    raw_decom = 0;
    fill_gaps = 0;
    sort_off  = 0;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hvrfso:l:i:")) != EOF)
        switch (i) {
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

    err = MPI_Info_create(&info);
    CHECK_MPIERR
    err = MPI_Info_set(info, "nc_header_align_size", "8192");
    CHECK_MPIERR

    /* create a new NC file */
    err = ncmpi_create(MPI_COMM_WORLD, outfname, NC_NOCLOBBER|NC_64BIT_DATA,
                       info, &ncid);
    if (err != NC_NOERR) {
        printf("Error at line %d in %s: %s\n", __LINE__, __FILE__,
               ncmpi_strerrno(err));
        goto fn_exit;
    }
    err = MPI_Info_free(&info);
    CHECK_MPIERR

    /* add the number of decompositions */
    err = ncmpi_def_dim(ncid, "num_decomp", num_decomp, &dimid);
    ERR

    /* add command line used */
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "command_line",
                             strlen(cmd_line), cmd_line);
    ERR

    for (i=0; i<num_decomp; i++) {
        err = add_decomp(ncid, infname[i], i+1); ERR
        err = ncmpi_redef(ncid); ERR
    }

    err = ncmpi_close(ncid);
    ERR

fn_exit:
    for (i=0; i<MAX_NFILES; i++)
        if (infname[i] != NULL) free(infname[i]);
    free(outfname);
    free(inList);

    MPI_Finalize();
    return(err != NC_NOERR);
}
