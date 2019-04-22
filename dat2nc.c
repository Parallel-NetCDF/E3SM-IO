/*********************************************************************
 *
 * Copyright (C) 2018, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program converts a PIO decomposition file (text) to a NetCDF file,
 * to be used by e3sm_io.c to benchmark PnetCDF performance.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <unistd.h> /* getopt() */

#include <assert.h>
#include <errno.h>
#include <mpi.h>
#include <pnetcdf.h>

#define LINE_SIZE 4692802

#define ERR {if(err!=NC_NOERR){printf("Error at line %d in %s: %s\n", __LINE__,__FILE__, ncmpi_strerrno(err));nerrs++;}}

static int verbose, line_sz;

/*----< intcompare() >------------------------------------------------------*/
/* This subroutine is used in qsort() */
static int intcompare(const void *p1, const void *p2)
{
    int i = *((int *)p1);
    int j = *((int *)p2);
    if (i > j) return (1);
    if (i < j) return (-1);
    return (0);
}

/*----< add_decomp() >------------------------------------------------------------*/
static int
add_decomp(int         ncid,
           const char *infname,
           int         label)
{
    char  *buf, name[128], *map, *str;
    FILE *fd;
    int   rank, nprocs, ndims;
    int   i, j, dimid, varid[3], *nreqs, err, nerrs=0, *off, *len;
    int   total_nreqs, max_nreqs, min_nreqs;
    MPI_Offset k, gsize, *dims, *dims_C, start, count;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    fd = fopen(infname, "r");
    if (fd == NULL) {
        printf("Error: fail to open file %s (%s)\n",infname,strerror(errno));
        MPI_Finalize();
	exit(1);
    }

    /* buffer stores one line read from input file */
    buf = (char*) malloc(line_sz);

    /* header lines: (first 2 lines of the decomposition file), for example
     *     version 2001 npes 43200 ndims 2
     *     777602 72
     * version is 2001
     * number of MPI processes used to generate this decomposition is 43200
     * number of dimensions of decomposed variable is 2
     * 1st dimension is of size 777602 and 2nd is 72 (in Fortran order)
     */
    fgets(buf, LINE_SIZE, fd);
    if (verbose) printf("header:\n\t%s",buf);
    strtok(buf, " ");
    if (strcmp(buf, "version")) {
        printf("Error at line %d: wrong header of input file %s\n",
               __LINE__,infname);
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
    if (verbose) printf("\t%s",buf);
    dims[0] = atoll(strtok(buf, " "));
    for (i=1; i<ndims; i++)
        dims[i] = atoll(strtok(NULL, " "));
    /* Note dims[] is in Fortran order */
    if (verbose) {
        if (ndims == 1) printf("lable D%d: dims = %lld\n",label, dims[0]);
        if (ndims == 2) printf("lable D%d: dims = %lld %lld\n",label, dims[0], dims[1]);
    }

    /* gsize is total number of elements in the global array */
    gsize=dims[0];
    for (i=1; i<ndims; i++) gsize *= dims[i];

    /* map is used to check whether the entire array is covered by requests
     * of all processes.
     */
    map = (char*) calloc(gsize, 1);

    /* nreqs[i] is the number of elements accessed by process i */
    nreqs = (int*) calloc(nprocs, sizeof(int));

    /* decomposition data format:
     *     (process.rank.ID) (number.of.requests)
     *     a list of element offsets accessed by this process (one element each)
     * Offsets are indices of array flattened into 1D and start with 1, i.e.
     * Fortran index based. Note the offsets are not sorted in an increasing
     * order and may contain 0s which should be ignored.
     */
    for (rank=0; rank<nprocs; rank++) {
        int ncontig=0, decomp_rank;
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
        while (rank < decomp_rank) /* this rank has no request */
            nreqs[rank++] = 0;

        nreqs[rank] = atoi(strtok(NULL, " "));  /* number of requests */
        if (nreqs[rank] == 0) /* this rank has zero request */
            continue; /* loop of rank */

        off = (int*) malloc(nreqs[rank] * sizeof(int));
        fgets(buf, LINE_SIZE+1, fd); /* 2nd line: list of offsets */
        if (strlen(buf) >= LINE_SIZE) {
            printf("Error: line size is larger than default %d\n",LINE_SIZE);
            printf("       use command-line option -l to use a larger size\n");
            goto fn_exit;
        }

        /* construct the offset list */
        off[0] = atoi(strtok(buf, " "));
        j=1;
        while (off[0] == 0) { /* skip leading 0 values, if there is any */
            off[0] = atoi(strtok(NULL, " "));
            j++;
        }
        off[0]--;
        k = 1;
        for (; j<nreqs[rank]; j++) {
            off[k] = atoi(strtok(NULL, " "));
            if (off[k] == 0) continue; /* skip 0 values */
            off[k]--;  /* offset is 1 based */
            k++;
        }

        /* sort off[] into an increasing order */
        qsort((void*)off, k, sizeof(int), intcompare);

        ncontig = 1;
        for (j=1; j<k; j++) {
            /* break contiguity at dimension boundaries or noncontiguous */
            if (off[j] % dims[0] == 0 || off[j] > off[j-1] + 1)
                ncontig++;
        }
        nreqs[rank] = ncontig;

        for (j=0; j<k; j++) map[off[j]]=1;
        free(off);
    }

    /* find total, max, and min nreqs amount all processes */
    total_nreqs = max_nreqs = min_nreqs = nreqs[0];
    for (i=1; i<nprocs; i++) {
        total_nreqs += nreqs[i];
        max_nreqs = (nreqs[i] > max_nreqs) ? nreqs[i] : max_nreqs;
        min_nreqs = (nreqs[i] < min_nreqs) ? nreqs[i] : min_nreqs;
    }

    if (verbose)
        printf("max_nreqs=%d min_nreqs=%d\n",max_nreqs, min_nreqs);

    /* check if the entire array is covered */
    for (k=0,j=0; j<gsize; j++) k+=map[j];
    if (k != gsize)
        printf("Warning: global array size = %lld but only %lld written\n",gsize,k);
    free(map);

    /* check if dimension decomp_nprocs has been defined in the netCDF file */
    err = ncmpi_inq_dimid(ncid, "decomp_nprocs", &dimid);
    if (err == NC_EBADDIM) {  /* not defined */
        err = ncmpi_def_dim(ncid, "decomp_nprocs", nprocs, &dimid); ERR
    }
    else {
        /* if decomp_nprocs already exist, check if value matches */
        MPI_Offset decomp_nprocs;
        err = ncmpi_inq_dimlen(ncid, dimid, &decomp_nprocs); ERR
        if (decomp_nprocs != nprocs) {
            printf("Error: decomp_nprocs=%lld mismatches among input files %d\n", decomp_nprocs, nprocs);
            MPI_Finalize();
            exit(1);
        }
    }

    /* define variable nreqs for this decomposition */
    sprintf(name, "D%d.nreqs", label);
    err = ncmpi_def_var(ncid, name, NC_INT, 1, &dimid, &varid[0]); ERR

    /* define attribute description for this variable */
    str = "Number of noncontiguous requests per process";
    err = ncmpi_put_att_text(ncid, varid[0], "description", strlen(str), str); ERR

    /* define dimension total_nreqs for this decomposition */
    sprintf(name, "D%d.total_nreqs", label);
    err = ncmpi_def_dim(ncid, name, total_nreqs, &dimid); ERR

    /* define variable offsets (store starting element indices) */
    sprintf(name, "D%d.offsets", label);
    err = ncmpi_def_var(ncid, name, NC_INT, 1, &dimid, &varid[1]); ERR
    str = "Flattened starting indices of noncontiguous requests";
    err = ncmpi_put_att_text(ncid, varid[1], "description", strlen(str), str); ERR

    /* define variable lengths (store number of elements) */
    sprintf(name, "D%d.lengths", label);
    err = ncmpi_def_var(ncid, name, NC_INT, 1, &dimid, &varid[2]); ERR
    str = "Lengths of noncontiguous requests";
    err = ncmpi_put_att_text(ncid, varid[2], "description", strlen(str), str); ERR

    /* add attribute to describe dimensionality */
    sprintf(name, "D%d.ndims", label);
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, name, NC_INT, 1, &ndims); ERR

    /* swap dims in Fortran order to dims_C in C order */
    dims_C = (MPI_Offset*) malloc(ndims * sizeof(MPI_Offset));
    for (i=0; i<ndims; i++) dims_C[i] = dims[ndims-i-1];
    sprintf(name, "D%d.dims", label);
    err = ncmpi_put_att_longlong(ncid, NC_GLOBAL, name, NC_INT, ndims, dims_C);
    ERR
    free(dims_C);

    sprintf(name, "D%d.max_nreqs", label);
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, name, NC_INT, 1, &max_nreqs); ERR
    sprintf(name, "D%d.min_nreqs", label);
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, name, NC_INT, 1, &min_nreqs); ERR

    /* exit define mode */
    err = ncmpi_enddef(ncid); ERR

    /* write variable containing number of requests for each process */
    err = ncmpi_put_var_int_all(ncid, varid[0], nreqs); ERR

    /* read the offsets again into allocated array off */
    start = 0;
    rewind(fd);
    fgets(buf, LINE_SIZE, fd);
    fgets(buf, LINE_SIZE, fd);
    for (rank=0; rank<nprocs; rank++) {
        int prev, ncontig=0, decomp_rank;
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
        while (rank < decomp_rank) /* this rank has no request */
            rank++;

        nreqs[rank] = atoi(strtok(NULL, " "));  /* number of requests */
        if (nreqs[rank] == 0) /* this rank has zero request */
            continue; /* loop of rank */

        off = (int*) malloc(nreqs[rank] * sizeof(int));
        len = (int*) malloc(nreqs[rank] * sizeof(int));
        fgets(buf, LINE_SIZE, fd);
        off[0] = atoi(strtok(buf, " "));
        j=1;
        while (off[0] == 0) {
            off[0] = atoi(strtok(NULL, " "));
            j++;
        }
        off[0]--;
        k = 1;
        for (; j<nreqs[rank]; j++) {
            off[k] = atoi(strtok(NULL, " "));
            if (off[k] == 0) continue; /* skip 0 values */
            off[k]--;  /* offset is 1 based */
            k++;
        }

        /* sort off[] into an increasing order */
        qsort((void*)off, k, sizeof(int), intcompare);

        ncontig = 1;
        prev = 0;
        len[0] = 1;
        for (j=1; j<k; j++) {
            /* break contiguity at dimension boundaries or noncontiguous */
            if (off[j] % dims[0] == 0 || off[j] > off[j-1] + 1) ncontig++;

            if (off[j] % dims[0] == 0 || off[j] > off[prev] + len[prev]) {
                prev++;
                if (prev < j) off[prev] = off[j];
                len[prev] = 1;
            }
            else len[prev]++;
        }
        assert(prev+1 == ncontig);

        /* write/append to variables offsets and lengths */
        count = ncontig;
        err = ncmpi_put_vara_int_all(ncid, varid[1], &start, &count, off); ERR
        err = ncmpi_put_vara_int_all(ncid, varid[2], &start, &count, len); ERR
        start += ncontig;

        free(off);
        free(len);
    }
    free(dims);

fn_exit:
    fclose(fd);
    free(nreqs);
    free(buf);

    return nerrs;
}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [OPTION]...\n"
    "       -h               Print help\n"
    "       -v               Verbose mode\n"
    "       -l num           max number of characters per line in input file\n"
    "       -o out_file      name of output netCDF file\n"
    "       -1 input_file    name of 1st decomposition file\n"
    "       -2 input_file    name of 2nd decomposition file\n"
    "       -3 input_file    name of 3rd decomposition file\n"
    "       -4 input_file    name of 4th decomposition file\n"
    "       -5 input_file    name of 5th decomposition file\n"
    "       -6 input_file    name of 6th decomposition file\n";
    fprintf(stderr, help, argv0);
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv) {
    char *infname[6], outfname[1024], cmd_line[4096];
    int  i, rank, ncid, num_decomp=0, dimid, err, nerrs=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    cmd_line[0] = '\0';
    for (i=0; i<argc; i++) {
        strcat(cmd_line, argv[i]);
        strcat(cmd_line, " ");
    }

    for (i=0; i<6; i++) infname[i] = NULL;
    outfname[0] = '\0';
    line_sz = LINE_SIZE;
    verbose = 0;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hvo:l:1:2:3:4:5:6:")) != EOF)
        switch(i) {
            case 'v': verbose = 1;
                      break;
            case 'o': strcpy(outfname, optarg);
                      break;
            case 'l': line_sz = atoi(optarg);
                      break;
            case '1': infname[0] = optarg;
                      num_decomp++;
                      break;
            case '2': infname[1] = optarg;
                      num_decomp++;
                      break;
            case '3': infname[2] = optarg;
                      num_decomp++;
                      break;
            case '4': infname[3] = optarg;
                      num_decomp++;
                      break;
            case '5': infname[4] = optarg;
                      num_decomp++;
                      break;
            case '6': infname[5] = optarg;
                      num_decomp++;
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }

    if (num_decomp == 0) {
        printf("Error: input decomposition files are required.\n");
        printf("       add -h to see all available command-line options.\n");
        MPI_Finalize();
        return 1;
    }

    if (num_decomp != 3 && num_decomp != 6) {
        printf("Error: e3sm_io currently supports F case (3 decomposition files)\n");
        printf("       and G case (6 decomposition files) only.\n");
        MPI_Finalize();
        return 1;
    }

    if (verbose) {
        for (i=0; i<6; i++) {
            if (infname[i] == NULL) continue;
            printf("input file %d = %s\n",i, infname[i]);
        }
    }

    if (outfname[0] == '\0') {
        strcpy(outfname, infname[0]);
        char *ptr = strstr(infname[0], ".dat");
        if (ptr == NULL) {
            if (rank == 0) {
                printf("Error at line %d: expected input file name extension is \".dat\"\n",__LINE__);
                printf("Abort ...\n");
            }
            MPI_Finalize();
            return 1;
        }
        sprintf(strstr(outfname, ".dat"),".nc");
    }
    if (verbose) printf("output file name =%s\n",outfname);

    if ((infname[0] != NULL && strcmp(infname[0], outfname) == 0) ||
        (infname[1] != NULL && strcmp(infname[1], outfname) == 0) ||
        (infname[2] != NULL && strcmp(infname[2], outfname) == 0) ||
        (infname[3] != NULL && strcmp(infname[3], outfname) == 0) ||
        (infname[4] != NULL && strcmp(infname[4], outfname) == 0) ||
        (infname[5] != NULL && strcmp(infname[5], outfname) == 0)) {
        if (rank == 0) {
            printf("Error: input and output file names conflict\n");
            printf("Abort ...\n");
        }
        MPI_Finalize();
        return 1;
    }

    /* create a new NC file */
    err = ncmpi_create(MPI_COMM_WORLD, outfname, NC_NOCLOBBER, MPI_INFO_NULL, &ncid);
    if (err != NC_NOERR) {
        printf("Error at line %d in %s: %s\n", __LINE__,__FILE__, ncmpi_strerrno(err));
        printf("Abort\n");
        nerrs++;
        goto fn_exit;
    }

    /* add the number of decompositions */
    err = ncmpi_def_dim(ncid, "num_decomp", num_decomp, &dimid); ERR

    /* add command line used */
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "command_line", strlen(cmd_line), cmd_line); ERR

    for (i=0; i<6; i++) {
        if (infname[i] == NULL) continue;
        nerrs += add_decomp(ncid, infname[i], i+1);
        err = ncmpi_redef(ncid); ERR;
    }

    err = ncmpi_close(ncid); ERR

fn_exit:
    MPI_Finalize();
    return (nerrs > 0);
}
