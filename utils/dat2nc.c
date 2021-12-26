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
#include <string.h>
#include <unistd.h> /* getopt() */

#include <mpi.h>
#include <pnetcdf.h>

#define MAX_NFILES 6
#define LINE_SIZE 4692802

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

static int verbose, line_sz, raw_decom;

/*----< intcompare() >------------------------------------------------------*/
/* This subroutine is used in qsort() */
static
int intcompare(const void *p1, const void *p2) {
    int i = *((int *)p1);
    int j = *((int *)p2);
    if (i > j) return(1);
    if (i < j) return(-1);
    return(0);
}

/*----< add_decomp() >-------------------------------------------------------*/
static
int add_decomp(int ncid, const char *infname, int label) {
    char *buf, name[128], *map, *str;
    FILE *fd;
    int i, j, rank, nprocs, dimid, ndims, err=NC_NOERR;
    int varid[5], *nreqs, *raw_nreqs=NULL, *off, *len, *raw_off;
    int total_nreqs, max_nreqs, min_nreqs, maxlen, minlen, total_raw_nreqs;
    MPI_Offset k, gsize, *dims, *dims_C, start, count, raw_start, raw_count;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    fd = fopen(infname, "r");
    if (fd == NULL) {
        printf("Error: fail to open file %s(%s)\n", infname, strerror(errno));
        MPI_Finalize();
        exit(1);
    }

    /* buffer stores one line read from input file */
    buf = (char *)malloc(line_sz);

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
    dims = (MPI_Offset *)malloc(ndims * sizeof(MPI_Offset));

    /* get dimension sizes
     * Note the dimensions are in Fortran order in the decomposition file
     * generated from PIO library
     */
    fgets(buf, LINE_SIZE, fd);
    if (verbose) printf("\t%s", buf);
    dims[0] = atoll(strtok(buf, " "));
    for (i = 1; i < ndims; i++) dims[i] = atoll(strtok(NULL, " "));
    /* Note dims[] is in Fortran order */
    if (verbose) {
        if (ndims == 1)
            printf("lable D%d: dims = %lld\n", label, dims[0]);
        else if (ndims == 2)
            printf("lable D%d: dims = %lld %lld\n", label, dims[0], dims[1]);
    }

    /* gsize is total number of elements in the global array */
    gsize = dims[0];
    for (i=1; i<ndims; i++) gsize *= dims[i];

    /* map is used to check whether the entire array is covered by requests
     * of all processes.
     */
    map = (char *)calloc(gsize, 1);

    /* nreqs[i] is the number of elements accessed by process i */
    nreqs = (int *)calloc(nprocs, sizeof(int));
    raw_nreqs = (int *)calloc(nprocs, sizeof(int));

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
        while (NULL != (ret = fgets(buf, LINE_SIZE, fd))) {
            if (buf[0] != '\n') break; /* non-empty line */
        }
        if (ret == NULL || strncmp(buf, "Obtained", 8) == 0) {
            /* there is no request for remaining ranks */
            for (decomp_rank = rank; decomp_rank < nprocs; decomp_rank++)
                nreqs[rank] = 0;
            break; /* loop of rank */
        }

        decomp_rank = atoi(strtok(buf, " ")); /* rank ID */
        while (rank < decomp_rank)            /* this rank has no request */
            nreqs[rank++] = 0;

        nreqs[rank] = atoi(strtok(NULL, " ")); /* number of requests */
        if (nreqs[rank] == 0) /* this rank has zero request */
            continue;         /* loop of rank */

        /* Record number of raw offsets before it is merged */
        raw_nreqs[rank] = nreqs[rank];

        off = (int *)malloc(nreqs[rank] * sizeof(int));
        fgets(buf, LINE_SIZE + 1, fd); /* 2nd line: list of offsets */
        if (strlen(buf) >= LINE_SIZE) {
            printf("Error: line size is larger than default %d\n", LINE_SIZE);
            printf("       use command-line option -l to use a larger size\n");
            goto fn_exit;
        }

        /* construct the offset list */
        off[0] = atoi(strtok(buf, " "));
        j = 1;
        while (off[0] == 0) { /* skip leading 0 values, if there is any */
            off[0] = atoi(strtok(NULL, " "));
            j++;
        }
        off[0]--;
        k = 1;
        for (; j < nreqs[rank]; j++) {
            off[k] = atoi(strtok(NULL, " "));
            if (off[k] == 0) continue; /* skip 0 values */
            off[k]--;                  /* offset is 1 based */
            k++;
        }

        /* sort off[] into an increasing order */
        qsort((void *)off, k, sizeof(int), intcompare);

        ncontig = 1;
        for (j=1; j<k; j++) {
            /* break contiguity at dimension boundaries or noncontiguous */
            if (off[j] % dims[0] == 0 || off[j] > off[j - 1] + 1)
               ncontig++;
        }
        nreqs[rank] = ncontig;

        for (j = 0; j < k; j++) map[off[j]] = 1;
        free(off);
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

    if (verbose) printf("max_nreqs=%d min_nreqs=%d\n", max_nreqs, min_nreqs);

    /* check if the entire array is covered */
    for (k = 0, j = 0; j < gsize; j++) k += map[j];
    if (k != gsize) {
        printf("Warning: decomposition %d does not cover the entire array\n",
               label);
        printf("\tglobal %dD array size: %lld", ndims, dims[0]);
        for (i=1; i<ndims; i++) printf(" x %lld", dims[i]);
        printf(" = %lld\n", gsize);
        printf("\tthe decomposition map covers only %lld of them\n", k);
    }
    free(map);

    /* check if dimension decomp_nprocs has been defined in the netCDF file */
    err = ncmpi_inq_dimid(ncid, "decomp_nprocs", &dimid);
    if (err == NC_EBADDIM) { /* not defined */
        err = ncmpi_def_dim(ncid, "decomp_nprocs", nprocs, &dimid);
        ERR
    } else {
        /* if decomp_nprocs already exist, check if value matches */
        MPI_Offset decomp_nprocs;
        err = ncmpi_inq_dimlen(ncid, dimid, &decomp_nprocs);
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
    err = ncmpi_def_var(ncid, name, NC_INT, 1, &dimid, &varid[0]);
    ERR

    /* define attribute description for this variable */
    str = "Number of noncontiguous requests per process";
    err = ncmpi_put_att_text(ncid, varid[0], "description", strlen(str), str);
    ERR

    /* define variable raw_nreqs for this decomposition */
    if (raw_decom) {
        sprintf(name, "D%d.raw_nreqs", label);
        err = ncmpi_def_var(ncid, name, NC_INT, 1, &dimid, &varid[3]);
        ERR

        /* define attribute description for this variable */
        str = "Number of file offsets accessed per process before merge";
        err = ncmpi_put_att_text(ncid, varid[3], "description", strlen(str), str);
        ERR
    }

    /* define dimension total_nreqs for this decomposition */
    sprintf(name, "D%d.total_nreqs", label);
    err = ncmpi_def_dim(ncid, name, total_nreqs, &dimid);
    ERR

    /* define variable offsets(store starting element indices) */
    sprintf(name, "D%d.offsets", label);
    err = ncmpi_def_var(ncid, name, NC_INT, 1, &dimid, &varid[1]);
    ERR
    str = "Flattened starting indices of noncontiguous requests";
    err = ncmpi_put_att_text(ncid, varid[1], "description", strlen(str), str);
    ERR

    /* define variable lengths(store number of elements) */
    sprintf(name, "D%d.lengths", label);
    err = ncmpi_def_var(ncid, name, NC_INT, 1, &dimid, &varid[2]);
    ERR
    str = "Lengths of noncontiguous requests";
    err = ncmpi_put_att_text(ncid, varid[2], "description", strlen(str), str);
    ERR

    maxlen = minlen = 0;
    err             = ncmpi_put_att_int(ncid, varid[2], "max", NC_INT, 1, &maxlen);
    ERR
    err = ncmpi_put_att_int(ncid, varid[2], "min", NC_INT, 1, &minlen);
    ERR

    if (raw_decom) {
        /* define dimension total_raw_nreqs for this decomposition */
        sprintf(name, "D%d.total_raw_nreqs", label);
        err = ncmpi_def_dim(ncid, name, total_raw_nreqs, &dimid);
        ERR

        /* define variable raw_offsets */
        sprintf(name, "D%d.raw_offsets", label);
        err = ncmpi_def_var(ncid, name, NC_INT, 1, &dimid, &varid[4]);
        ERR
        str = "File offsets accessed before merge";
        err = ncmpi_put_att_text(ncid, varid[4], "description", strlen(str), str);
        ERR
    }

    /* add attribute to describe dimensionality */
    sprintf(name, "D%d.ndims", label);
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, name, NC_INT, 1, &ndims);
    ERR

    /* swap dims in Fortran order to dims_C in C order */
    dims_C = (MPI_Offset *)malloc(ndims * sizeof(MPI_Offset));
    for (i=0; i < ndims; i++) dims_C[i] = dims[ndims - i - 1];
    sprintf(name, "D%d.dims", label);
    err = ncmpi_put_att_longlong(ncid, NC_GLOBAL, name, NC_INT, ndims, dims_C);
    ERR
    free(dims_C);

    sprintf(name, "D%d.max_nreqs", label);
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, name, NC_INT, 1, &max_nreqs);
    ERR
    sprintf(name, "D%d.min_nreqs", label);
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, name, NC_INT, 1, &min_nreqs);
    ERR

    /* exit define mode */
    err = ncmpi_enddef(ncid);
    ERR

    /* write variable containing number of requests for each process */
    err = ncmpi_put_var_int_all(ncid, varid[0], nreqs);
    ERR

    /* write variable containing number of requests before merge for each process */
    if (raw_decom){
        err = ncmpi_put_var_int_all(ncid, varid[3], raw_nreqs);
        ERR
    }

    /* read the offsets again into allocated array off */
    start = 0;
    raw_start = 0;
    rewind(fd);
    fgets(buf, LINE_SIZE, fd);
    fgets(buf, LINE_SIZE, fd);
    for (rank = 0; rank < nprocs; rank++) {
        int prev, ncontig = 0, decomp_rank;
        char *ret;

        /* reads the first line of rank ID and no. requests */
        while (NULL != (ret = fgets(buf, LINE_SIZE, fd))) {
            if (buf[0] != '\n') break; /* non-empty line */
        }
        if (ret == NULL || strncmp(buf, "Obtained", 8) == 0) {
            /* there is no request for remaining ranks */
            for (decomp_rank = rank; decomp_rank < nprocs; decomp_rank++) nreqs[rank] = 0;
            break; /* loop of rank */
        }

        decomp_rank = atoi(strtok(buf, " ")); /* rank ID */
        while (rank < decomp_rank)            /* this rank has no request */
            rank++;

        nreqs[rank] = atoi(strtok(NULL, " ")); /* number of requests */
        if (nreqs[rank] == 0)                  /* this rank has zero request */
            continue;                          /* loop of rank */

        off = (int *)malloc(nreqs[rank] * sizeof(int));
        len = (int *)malloc(nreqs[rank] * sizeof(int));
        raw_off = (int *)malloc(raw_nreqs[rank] * sizeof(int));
        fgets(buf, LINE_SIZE, fd);
        i=0;
        off[0] = raw_off[i++] = atoi(strtok(buf, " "));
        j      = 1;
        while (off[0] == 0) {
            off[0] = raw_off[i++] = atoi(strtok(NULL, " "));
            j++;
        }
        off[0]--;
        k = 1;
        for (; j < nreqs[rank]; j++) {
            off[k] = raw_off[i++] = atoi(strtok(NULL, " "));
            if (off[k] == 0) continue; /* skip 0 values */
            off[k]--;                  /* offset is 1 based */
            k++;
        }

        /* sort off[] into an increasing order */
        qsort((void *)off, k, sizeof(int), intcompare);

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
        assert(prev + 1 == ncontig);

        if (rank == 0) maxlen = minlen = len[0];
        for (j = 0; j < ncontig; j++) {
            maxlen = (len[j] > maxlen) ? len[j] : maxlen;
            minlen = (len[j] < minlen) ? len[j] : minlen;
        }

        /* write/append to variables offsets and lengths */
        count = ncontig;
        err   = ncmpi_put_vara_int_all(ncid, varid[1], &start, &count, off);
        ERR
        err = ncmpi_put_vara_int_all(ncid, varid[2], &start, &count, len);
        ERR
        start += ncontig;

        /* write/append raw file offsets before merge */
        if (raw_decom){
            raw_count = raw_nreqs[rank];
            err = ncmpi_put_vara_int_all(ncid, varid[4], &raw_start, &raw_count, raw_off);
            ERR
            raw_start += raw_count; 
        }

        free(off);
        free(len);
        free(raw_off);
    }
    free(dims);

    err = ncmpi_put_att_int(ncid, varid[2], "max", NC_INT, 1, &maxlen);
    ERR
    err = ncmpi_put_att_int(ncid, varid[2], "min", NC_INT, 1, &minlen);
    ERR

fn_exit:
    fclose(fd);
    free(nreqs);
    free(raw_nreqs);
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
        printf("Error: open fails on file %s(%s)\n",inList,strerror(errno));
        exit(1);
    }

    /* count number of input files */
    num_files = 0;
    while (fgets(line, 1024, fd)) {
        if (strlen(line) == 0)
            continue; /* skip empty lines */
        if (line[0] == '#')
            continue; /* skip comment line(start with #) */
        num_files++;
    }

    /* read input file names */
    rewind(fd);
    i=0;
    while (fgets(line, 1024, fd)) {
        char *tail;
        if (strlen(line) == 0)
            continue; /* skip empty lines */
        if (line[0] == '#')
            continue; /* skip comment line(start with #) */
        /* remove blanks at tail. Note fgets stores newline to the buffer */
        tail = line + strlen(line) - 1;
        while (*tail == ' ' || *tail == '\t' || *tail == '\n') tail--;
        tail[1] = '\0';
        /* save file name to in_list */
        infname[i] = strdup(line);
        i++;
    }
    assert(i == num_files);
    fclose(fd);

    *num_decomp = num_files;
}

static void usage(char *argv0) {
    char *help =
        "Usage: %s [-h|-v|-r|-l] -i input_file -o out_file\n"
        "   -h               Print help\n"
        "   -v               Verbose mode\n"
        "   -r               Include original decomposition map\n"
        "   -l num           max number of characters per line in input file\n"
        "   -i input_file    list of decomposition file names\n"
        "   -o out_file      name of output netCDF file\n";
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
    for (i=0; i < argc; i++) {
        strcat(cmd_line, argv[i]);
        strcat(cmd_line, " ");
    }

    for (i=0; i<MAX_NFILES; i++) infname[i] = NULL;
    line_sz   = LINE_SIZE;
    verbose   = 0;
    raw_decom = 0;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hvro:l:i:")) != EOF)
        switch (i) {
            case 'v': verbose = 1;
                      break;
            case 'o': outfname = strdup(optarg);
                      break;
            case 'l': line_sz = atoi(optarg);
                      break;
            case 'r': raw_decom = 1;
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
        err = add_decomp(ncid, infname[i], i + 1); ERR;
        err = ncmpi_redef(ncid); ERR;
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
