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

#define LINE_SIZE 2097152

#define ERR {if(err!=NC_NOERR){printf("Error at line %d in %s: %s\n", __LINE__,__FILE__, ncmpi_strerrno(err));nerrs++;}}

static int verbose, line_sz;

/*----< add_decomp() >------------------------------------------------------------*/
int add_decomp(int   ncid,
               char *infname,
               char *label)
{
    char  *buf, name[128], *map, *str;
    FILE *fd;
    int   rank, nprocs, ndims;
    int   i, j, dimids[2], varid[2], *nreqs, err, nerrs=0, *off;
    int max_nreqs, min_nreqs;
    MPI_Offset k, gsize, *dims, start[2], count[2], dim_len;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    fd = fopen(infname, "r");
    if (fd == NULL) {
        printf("Error: fail to open file %s (%s)\n",infname,strerror(errno));
        MPI_Finalize();
	exit(1);
    }

    buf = (char*) malloc(line_sz);

    /* process header lines:
     * version 2001 npes 43200 ndims 2 
     * 777602 72 
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
    dims = (MPI_Offset*) malloc(3 * ndims * sizeof(MPI_Offset));

    /* get dimension sizes
     * Note the dimensions are in Fortran order in the decomposition file
     * generated from PIO library
     */
    fgets(buf, LINE_SIZE, fd);
    if (verbose) printf("\t%s",buf);
    dims[ndims-1] = atoll(strtok(buf, " "));
    for (i=ndims-2; i>=0; i--)
        dims[i] = atoll(strtok(NULL, " "));

    gsize=dims[0];
    for (i=1; i<ndims; i++) gsize *= dims[i];
    
    map = (char*) calloc(gsize, 1);
    nreqs = (int*) malloc(nprocs * sizeof(int));

    /* calculate max_nreqs, min_nreqs and nreqs[] */
    for (rank=0; rank<nprocs; rank++) {
        fgets(buf, LINE_SIZE, fd);
        assert(rank == atoi(strtok(buf, " ")));
        nreqs[rank] = atoi(strtok(NULL, " "));
        off = (int*) malloc(nreqs[rank] * sizeof(int));
        fgets(buf, LINE_SIZE+1, fd);
        if (strlen(buf) >= LINE_SIZE) {
            printf("Error: line size is larger than default %d\n",LINE_SIZE);
            printf("       use command-line option -l to use a larger size\n");
            goto fn_exit;
        }

        /* read offset list */
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
        if (rank == 0) {
            max_nreqs = min_nreqs = k;
        }
        else {
            max_nreqs = (k > max_nreqs) ? k : max_nreqs;
            min_nreqs = (k < min_nreqs) ? k : min_nreqs;
        }
        nreqs[rank]=k;

        for (j=0; j<k; j++) map[off[j]]=1;
        free(off);
    }
    if (verbose) printf("max_nreqs=%d min_nreqs=%d\n",max_nreqs, min_nreqs);
    for (k=0,j=0; j<gsize; j++) k+=map[j];
    if (k != gsize)
        printf("Warning: global array size = %lld but only %lld written\n",gsize,k);
    free(map);

    /* check if dimension num_procs has been defined */
    err = ncmpi_inq_dimid(ncid, "num_procs", &dimids[0]);
    if (err == NC_EBADDIM) {  /* not defined */
        err = ncmpi_def_dim(ncid, "num_procs", nprocs, &dimids[0]); ERR
    }
    else {
        MPI_Offset num_procs;
        err = ncmpi_inq_dimlen(ncid, dimids[0], &num_procs); ERR
        if (num_procs != nprocs) {
            printf("Error: num_procs mismatches among input files\n");
            MPI_Finalize();
            exit(1);
        }
    }

    sprintf(name, "%s.max_nreqs", label);
    err = ncmpi_def_dim(ncid, name, max_nreqs, &dimids[1]); ERR
    sprintf(name, "%s.nreqs", label);
    err = ncmpi_def_var(ncid, name, NC_INT, 1, dimids, &varid[0]); ERR
    str = "Number of noncontiguous subarray requests by each MPI process";
    err = ncmpi_put_att_text(ncid, varid[0], "description", strlen(str), str); ERR

    sprintf(name, "%s.offsets", label);
    err = ncmpi_def_var(ncid, name, NC_INT, 2, dimids, &varid[1]); ERR
    str = "Flattened starting indices of noncontiguous requests. Each row corresponds to requests by an MPI process.";
    err = ncmpi_put_att_text(ncid, varid[1], "description", strlen(str), str); ERR

    if (ndims > 1) {
        err = ncmpi_get_att_longlong(ncid, NC_GLOBAL, "dim_len_Y", &dim_len);
        if (err == NC_ENOTATT) {
            err = ncmpi_put_att_longlong(ncid, NC_GLOBAL, "dim_len_Y", NC_INT, 1, &dims[0]);
            ERR
        }
        else {
            if (dim_len != dims[0]) {
                printf("Error: dim_len_Y mismatches among input files\n");
                MPI_Finalize();
                exit(1);
            }
        }
    }
    err = ncmpi_get_att_longlong(ncid, NC_GLOBAL, "dim_len_X", &dim_len);
    if (err == NC_ENOTATT) {
        dim_len = (ndims == 1) ? dims[0] : dims[1];
        err = ncmpi_put_att_longlong(ncid, NC_GLOBAL, "dim_len_X", NC_INT, 1, &dim_len);
        ERR
    }
    else {
        if ((ndims == 1 && dim_len != dims[0]) || (ndims == 2 && dim_len != dims[1])) {
            printf("Error: dim_len_X mismatches among input files\n");
            MPI_Finalize();
            exit(1);
        }
    }
    free(dims);

    sprintf(name, "%s.max_nreqs", label);
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, name, NC_INT, 1, &max_nreqs); ERR
    sprintf(name, "%s.min_nreqs", label);
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, name, NC_INT, 1, &min_nreqs); ERR
    err = ncmpi_enddef(ncid); ERR

    start[0] = 0;
    count[0] = nprocs;
    err = ncmpi_put_var_int_all(ncid, varid[0], nreqs); ERR

    off = (int*) malloc(max_nreqs * sizeof(int));

    rewind(fd);
    fgets(buf, LINE_SIZE, fd);
    fgets(buf, LINE_SIZE, fd);
    for (rank=0; rank<nprocs; rank++) {
        fgets(buf, LINE_SIZE, fd);
        assert(rank == atoi(strtok(buf, " ")));
        nreqs[rank] = atoi(strtok(NULL, " "));
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

        start[0] = rank;
        start[1] = 0;
        count[0] = 1;
        count[1] = k;
        err = ncmpi_put_vara_int_all(ncid, varid[1], start, count, off); ERR
    }
    free(off);

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
    "Usage: %s [OPTION]... [FILE]...\n"
    "       -h               Print help\n"
    "       -q               Quiet mode (reports when fail)\n"
    "       -l num           max number of characters per line in input file\n"
    "       -o out_file      name of output netCDF file\n"
    "       -1 input_file    name of 1st 1D decomposition file\n"
    "       -2 input_file    name of 2nd 1D decomposition file\n"
    "       -3 input_file    name of     2D decomposition file\n";
    fprintf(stderr, help, argv0);
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv) {
    char *infname_D1=NULL, *infname_D2=NULL, *infname_D3=NULL, outfname[1024];
    int  i, rank, ncid, err, nerrs=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    outfname[0] = '\0';
    line_sz = LINE_SIZE;
    verbose = 1;

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hqo:l:1:2:3:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'o': strcpy(outfname, optarg);
                      break;
            case 'l': line_sz = atoi(optarg);
                      break;
            case '1': infname_D1 = optarg;
                      break;
            case '2': infname_D2 = optarg;
                      break;
            case '3': infname_D3 = optarg;
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }

    if (infname_D1 == NULL || infname_D2 == NULL || infname_D3 == NULL) {
        /* 3 input data decomposition files are mandatory */
        usage(argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (verbose) printf("input files =%s %s %s\n",infname_D1,infname_D2,infname_D3);

    if (outfname[0] == '\0') {
        strcpy(outfname, infname_D1);
        char *ptr = strstr(infname_D1, ".dat");
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

    if (strcmp(infname_D1, outfname) == 0) {
        if (rank == 0) {
            printf("Error: input and output file names are the same\n");
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

    nerrs += add_decomp(ncid, infname_D3, "D3");
    err = ncmpi_redef(ncid); ERR;
    nerrs += add_decomp(ncid, infname_D2, "D2");
    err = ncmpi_redef(ncid); ERR;
    nerrs += add_decomp(ncid, infname_D1, "D1");

    err = ncmpi_close(ncid); ERR

fn_exit:
    MPI_Finalize();
    return (nerrs > 0);
}
