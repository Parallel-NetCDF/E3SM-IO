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

#define LINE_SIZE 1048576

#define ERR {if(err!=NC_NOERR){printf("Error at line %d in %s: %s\n", __LINE__,__FILE__, ncmpi_strerrno(err));nerrs++;}}

static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] [-q] [-o out_file] input_file\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [-o out_file] output netCDF file name\n"
    "       in_file input decomposition file\n";
    fprintf(stderr, help, argv0);
}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv) {
    char  buf[LINE_SIZE+1], *infname=NULL, outfname[1024], *map, *str;
    FILE *fd;
    int   rank, nprocs, verbose=1, ndims;
    int   i, j, ncid, dimids[2], varid[2], *nreqs, err, nerrs=0, *off;
    int max_nreqs, min_nreqs;
    MPI_Offset k, gsize, *dims, start[2], count[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    outfname[0] = '\0';

    /* get command-line arguments */
    while ((i = getopt(argc, argv, "hqo:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'o': strcpy(outfname, optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }

    if (argv[optind] == NULL) { /* input file is mandatory */
        usage(argv[0]);
        MPI_Finalize();
        return 1;
    }
    infname = argv[optind];
    if (verbose) printf("input file name =%s\n",infname);

    if (outfname[0] == '\0') {
        char *ptr = strstr(infname, ".dat");
        if (ptr == NULL) {
            if (rank == 0)
                printf("Error at line %d: expected input file name extension is \".dat\"\n",__LINE__);
            MPI_Finalize();
            return 1;
        }
        strcpy(outfname, infname);
        sprintf(strstr(outfname, ".dat"),".nc");
    }
    if (verbose) printf("output file name =%s\n",outfname);

    if (strcmp(infname, outfname) == 0) {
        if (rank == 0)
            printf("Error: input and output file names are the same\n");
        MPI_Finalize();
        return 1;
    }

    fd = fopen(infname, "r");
    if (fd == NULL) {
        printf("Error: fail to open file %s (%s)\n",infname,strerror(errno));
        MPI_Finalize();
	exit(1);
    }

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
    atoi(strtok(NULL, " ")); /* token "ndims" */
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
        assert(strlen(buf) < LINE_SIZE);

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

    err = ncmpi_create(MPI_COMM_WORLD, outfname, NC_NOCLOBBER, MPI_INFO_NULL, &ncid);
    if (err != NC_NOERR) {
        printf("Error at line %d in %s: %s\n", __LINE__,__FILE__,
        ncmpi_strerrno(err));
        nerrs++;
        goto fn_exit;
    }
    err = ncmpi_def_dim(ncid, "num_procs", nprocs, &dimids[0]); ERR
    err = ncmpi_def_dim(ncid, "max_nreqs", max_nreqs, &dimids[1]); ERR
    err = ncmpi_def_var(ncid, "nreqs",   NC_INT, 1, dimids, &varid[0]); ERR
    str = "Number of noncontiguous subarray requests by each MPI process";
    err = ncmpi_put_att_text(ncid, varid[0], "description", strlen(str), str); ERR

    err = ncmpi_def_var(ncid, "offsets", NC_INT, 2, dimids, &varid[1]); ERR
    str = "Flattened starting indices of noncontiguous requests. Each row corresponds to requests by an MPI process.";
    err = ncmpi_put_att_text(ncid, varid[1], "description", strlen(str), str); ERR

    err = ncmpi_put_att_int(ncid, NC_GLOBAL, "var_ndims", NC_INT, 1, &ndims); ERR
    for (i=0; i<ndims; i++) {
        char dim_name[32];
        sprintf(dim_name, "dim_len_%d", i);
        err = ncmpi_put_att_longlong(ncid, NC_GLOBAL, dim_name, NC_INT, 1, &dims[i]); ERR
    }
    free(dims);
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, "max_nreqs", NC_INT, 1, &max_nreqs); ERR
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, "min_nreqs", NC_INT, 1, &min_nreqs); ERR
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
    err = ncmpi_close(ncid); ERR

fn_exit:
    free(nreqs);
    fclose(fd);

    MPI_Finalize();
    return 0;
}
