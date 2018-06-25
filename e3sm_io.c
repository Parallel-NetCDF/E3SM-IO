/*********************************************************************
 *
 * Copyright (C) 2018, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program uses the E3SM I/O patterns recorded by the PIO library
 * to evaluate the performance of three PnetCDF APIs: ncmpi_vard_all(),
 * ncmpi_iput_varn(), and ncmpi_iput_vara(). The E3SM I/O patterns consist
 * of a large number of small, noncontiguous requests on each MPI process,
 * which presents a challenge for achieving a good performance.
 *
 * To compile:
 *    mpicc e3sm_io.c -o e3sm_io -I/path/PnetCDF/include -L/path/PnetCDF/lib -lpnetcdf
 *
 * To run, for example:
 *    mpiexec -n 8 ./e3sm_io -q datasets/48602x72_64p.nc -n 10
 *
 *    ---- benchmarking vard API -----------------------
 *    -----------------------------------------------------------
 *    MAX heap memory allocated by PnetCDF internally is 0.27 MiB
 *    Total number of variables          = 10
 *    Total write amount                 = 133.49 MiB = 0.13 GiB
 *    Max number of requests             = 22144
 *    Max Time of open + metadata define = 0.0005 sec
 *    Max Time of I/O preparing          = 0.0354 sec
 *    Max Time of ncmpi_put_vard         = 0.4521 sec
 *    Max Time of close                  = 0.0037 sec
 *    Max Time of TOTAL                  = 0.4916 sec
 *    I/O bandwidth                      = 271.5230 MiB/sec
 *
 *    ---- benchmarking varn API -----------------------
 *    -----------------------------------------------------------
 *    MAX heap memory allocated by PnetCDF internally is 20.47 MiB
 *    Total number of variables          = 10
 *    Total write amount                 = 133.49 MiB = 0.13 GiB
 *    Max number of requests             = 22144
 *    Max Time of open + metadata define = 0.0004 sec
 *    Max Time of I/O preparing          = 0.0211 sec
 *    Max Time of ncmpi_iput_varn        = 0.0497 sec
 *    Max Time of ncmpi_wait_all         = 0.4593 sec
 *    Max Time of close                  = 0.0079 sec
 *    Max Time of TOTAL                  = 0.5385 sec
 *    I/O bandwidth                      = 247.9046 MiB/sec
 *
 *    ---- benchmarking vara API -----------------------
 *    -----------------------------------------------------------
 *    MAX heap memory allocated by PnetCDF internally is 45.78 MiB
 *    Total number of variables          = 10
 *    Total write amount                 = 133.49 MiB = 0.13 GiB
 *    Max number of requests             = 22144
 *    Max Time of open + metadata define = 0.0004 sec
 *    Max Time of I/O preparing          = 0.0210 sec
 *    Max Time of ncmpi_iput_vara        = 0.9343 sec
 *    Max Time of ncmpi_wait_all         = 0.6547 sec
 *    Max Time of close                  = 0.0141 sec
 *    Max Time of TOTAL                  = 1.6247 sec
 *    I/O bandwidth                      = 82.1640 MiB/sec
 *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h> /* strtoll() */
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() unlink() */

#include <mpi.h>
#include <pnetcdf.h>

#define USE_DOUBLE
#ifdef USE_DOUBLE
#define itype double
#define xtype NC_DOUBLE
#define mpitype MPI_DOUBLE
#else
#define itype float
#define xtype NC_FLOAT
#define mpitype MPI_FLOAT
#endif

#define MAX(a,b) ((a) > (b)) ? (a) : (b)
#define MIN(a,b) ((a) < (b)) ? (a) : (b)

#define ERR { \
    if (err != NC_NOERR) { \
        printf("Error in %s line %d: %s\n", __FILE__, __LINE__, \
               ncmpi_strerror(err)); \
        nerrs++; \
        goto fn_exit; \
    } \
}

static int verbose; /* verbose mode to print additional messages on screen */

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

/*----< read_io_pattern() >-------------------------------------------------*/
/* read I/O pattern file, infname. The contents of the file are, for example,
 *      netcdf 48602x72_64p {
 *      // file format: CDF-1
 *      dimensions:
 *              num_procs = 64 ;
 *              max_nreqs = 55296 ;
 *      variables:
 *              int nreqs(num_procs) ;
 *              int offsets(num_procs, max_nreqs) ;
 *
 *      // global attributes:
 *              :var_ndims = 2 ;
 *              :dim_len_0 = 48602 ;
 *              :dim_len_1 = 72 ;
 *              :max_nreqs = 55296 ;
 *              :min_nreqs = 51840 ;
 *      }
 */
static int
read_io_pattern(int          verbose,
                char        *infname,
                int         *ndims,       /* number of dimensions of variable */
                MPI_Offset **dims,        /* dimension lengths */
                int         *contig_nreqs,/* num of contiguous requests */
                int        **disps,
                int        **blocklens)
{
    int err, nerrs=0, rank, nprocs, ncid, varid, proc_start, proc_numb;
    int i, j, k, nreqs, dimids[2];
    MPI_Offset num_procs, max_nreqs, start[2], count[2];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* read start-count from input file */
    err = ncmpi_open(MPI_COMM_WORLD, infname, NC_NOWRITE, MPI_INFO_NULL,
                     &ncid); ERR
    err = ncmpi_inq_dimid(ncid, "num_procs", &dimids[0]); ERR
    err = ncmpi_inq_dimid(ncid, "max_nreqs", &dimids[1]); ERR
    err = ncmpi_inq_dimlen(ncid, dimids[0], &num_procs); ERR
    err = ncmpi_inq_dimlen(ncid, dimids[1], &max_nreqs); ERR

    err = ncmpi_get_att_int(ncid, NC_GLOBAL, "var_ndims", ndims); ERR
    *dims = (MPI_Offset*) malloc(*ndims * sizeof(MPI_Offset));
    for (i=0; i<*ndims; i++) {
        char dim_name[32];
        sprintf(dim_name, "dim_len_%d", i);
        err = ncmpi_get_att_longlong(ncid, NC_GLOBAL, dim_name, (*dims)+i); ERR
    }

    /* num_procs is the number of processes used to generate the input I/O
     * pattern file. nprocs is the number of processes running this benchmark.
     * This benchmark allows the two to be different. When nprocs is smaller
     * than num_procs, some of nprocs processes will carry out the requests
     * from more than one of num_procs processes. The requests responsible by
     * this process starts from proc_start with the number proc_numb. When
     * nprocs is bigger than num_procs, then those processes with rank
     * ID >= num_procs will have no data to write and they will just
     * participate the collective subroutines.
     */
    proc_numb = num_procs / nprocs;
    proc_start = rank * proc_numb;
    if (rank < num_procs % nprocs) {
        proc_start += rank;
        proc_numb++;
    }
    else
        proc_start += num_procs % nprocs;

    /* read number of requests for this process */
    err = ncmpi_inq_varid(ncid, "nreqs", &varid); ERR
    start[0] = proc_start;
    count[0] = proc_numb;
    int *num_reqs = (int*) malloc(proc_numb * sizeof(int));
    err = ncmpi_get_vara_int_all(ncid, varid, start, count, num_reqs); ERR

    /* calculate the total number of requests responsible by this process */
    nreqs = 0;
    for (i=0; i<proc_numb; i++) nreqs += num_reqs[i];
    if (verbose) printf("rank %d: nreqs=%d\n",rank,nreqs);

    /* read the starting offsets of all requests into disps[] */
    *disps = (int*) malloc(proc_numb * max_nreqs * sizeof(int));
    err = ncmpi_inq_varid(ncid, "offsets", &varid); ERR
    start[0] = proc_start;
    start[1] = 0;
    count[0] = proc_numb;
    count[1] = max_nreqs;
    err = ncmpi_get_vara_int_all(ncid, varid, start, count, *disps); ERR
    err = ncmpi_close(ncid); ERR

    /* coalesce the disps array into contiguous requests */
    k = 0;
    for (i=0; i<proc_numb; i++)
        for (j=0; j<num_reqs[i]; j++)
            (*disps)[k++] = (*disps)[i*max_nreqs + j];
    free(num_reqs);

    /* sort disps[] in increasing order (this is to satisfy the MPI fileview
     * or monotonically nondecreasing file offset requirement) */
    qsort((void*)(*disps), nreqs, sizeof(int), intcompare);

    /* define MPI derived datatype (each start-count pair is 1D contiguous) */
    *blocklens = (int*)malloc(nreqs*sizeof(int));
    (*blocklens)[0] = 1;
    for (j=0, i=1; i<nreqs; i++) {
        if ((*disps)[i] == (*disps)[i-1]+1 && (*disps)[i]%((*dims)[*ndims-1]))
            /* contiguous disp && each blocklens[i] is no longer than the last
             * dimension length */
            (*blocklens)[j]++;
        else {
            j++;
            (*disps)[j] = (*disps)[i];
            (*blocklens)[j] = 1;
        }
    }
    *contig_nreqs = j+1; /* the true number of contiguous requests */

    if (verbose) {
        int min_blocklen = (*blocklens)[0];
        int max_blocklen = (*blocklens)[0];
        for (i=1; i<*contig_nreqs; i++) {
            max_blocklen = MAX((*blocklens)[i], max_blocklen);
            min_blocklen = MIN((*blocklens)[i], min_blocklen);
        }
        printf("%3d nreqs=%d contig nreqs=%4d max_blocklen=%d min_blocklen=%d\n",
               rank, nreqs, *contig_nreqs, max_blocklen, min_blocklen);
    }

fn_exit:
    if (nerrs && *disps != NULL) {
        free(*disps);
        *disps = NULL;
    }
    if (nerrs && *blocklens != NULL) {
        free(*blocklens);
        *blocklens = NULL;
    }
    return nerrs;
}

/*----< print_info() >------------------------------------------------------*/
static
void print_info(MPI_Info *info_used)
{
    int  i, nkeys;

    MPI_Info_get_nkeys(*info_used, &nkeys);
    printf("MPI File Info: nkeys = %d\n",nkeys);
    for (i=0; i<nkeys; i++) {
        char key[MPI_MAX_INFO_KEY], value[MPI_MAX_INFO_VAL];
        int  valuelen, flag;

        MPI_Info_get_nthkey(*info_used, i, key);
        MPI_Info_get_valuelen(*info_used, key, &valuelen, &flag);
        MPI_Info_get(*info_used, key, valuelen+1, value, &flag);
        printf("MPI File Info: [%2d] key = %25s, value = %s\n",i,key,value);
    }
}

/*----< run_vard() >--------------------------------------------------------*/
static
int run_vard(char       *out_dir,   /* output folder name */
             MPI_Info    info,
             int         ndims,     /* number of dimensions of variables */
             MPI_Offset *dims,      /* dimension lengths */
             int         nvars,     /* number of variables */
             int         nreqs,     /* number of request per variable */
             int        *disps,     /* [nreqs] request's displacements */
             int        *blocklens) /* [nreqs] request's block lengths */
{
    char outfname[512];
    int i, err, nerrs=0, rank, cmode, ncid, *varids, *dimids;
    int *var_blocklens, max_nreqs;
    size_t buflen=0;
    itype *buffer;
    double pre_timing, open_timing, io_timing, close_timing;
    double timing, total_timing,  max_timing;
    MPI_Aint *var_disps;
    MPI_Offset *var_offsets, put_size, total_size;
    MPI_Datatype vartype, filetype;
    MPI_Info info_used=MPI_INFO_NULL;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Barrier(comm); /*-----------------------------------------*/
    total_timing = open_timing = MPI_Wtime();

    /* set output file name */
    sprintf(outfname, "%s/testfile_vard.nc",out_dir);

    /* create a new CDF-5 file for writing */
    cmode = NC_CLOBBER|NC_64BIT_DATA;
    err = ncmpi_create(comm, outfname, cmode, info, &ncid); ERR

    /* define dimensions */
    dimids = (int*) malloc(ndims * sizeof(int));
    for (i=0; i<ndims; i++) {
        char dim_name[32];
        sprintf(dim_name, "dim_len_%d", i);
        err = ncmpi_def_dim(ncid, dim_name, dims[i], &dimids[i]); ERR
    }

    /* define variables */
    varids = (int*) malloc(nvars * sizeof(int));
    for (i=0; i<nvars; i++) {
        char varname[128];
        sprintf(varname, "var%04d",i);
        err = ncmpi_def_var(ncid, varname, xtype, ndims, dimids, &varids[i]);
        ERR
    }
    free(dimids);
    err = ncmpi_enddef(ncid); ERR
    err = ncmpi_inq_file_info(ncid, &info_used); ERR
    open_timing = MPI_Wtime() - open_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    pre_timing = MPI_Wtime();

    /* define MPI datatype for one variable */
    MPI_Type_indexed(nreqs, blocklens, disps, mpitype, &vartype);
    MPI_Type_commit(&vartype);

    /* concatenate multiple vartype into filetype */
    var_blocklens = (int*) malloc(nvars * sizeof(int));
    var_disps = (MPI_Aint*) malloc(nvars * sizeof(MPI_Aint));
    var_offsets = (MPI_Offset*) malloc(nvars * sizeof(MPI_Offset));
    var_blocklens[0] = 1;
    var_disps[0] = 0;
    err = ncmpi_inq_varoffset(ncid, varids[0], &var_offsets[0]); ERR

    /* calculate distance between consecutive variables */
    for (i=1; i<nvars; i++) {
        var_blocklens[i] = 1;
        err = ncmpi_inq_varoffset(ncid, varids[i], &var_offsets[i]); ERR
        var_disps[i] = var_offsets[i] - var_offsets[0];
    }

    /* construct filetype */
    MPI_Type_create_hindexed(nvars, var_blocklens, var_disps, vartype,
                             &filetype);
    MPI_Type_commit(&filetype);
    MPI_Type_free(&vartype);
    free(var_blocklens);
    free(var_disps);
    free(var_offsets);

    MPI_Comm_rank(comm, &rank);

    /* allocate and initialize write buffer */
    for (i=0; i<nreqs; i++) buflen += blocklens[i];
    buffer = (itype*) malloc(buflen * nvars * sizeof(itype));
    for (i=0; i<buflen*nvars; i++) buffer[i] = rank;

    pre_timing = MPI_Wtime() - pre_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    io_timing = MPI_Wtime();

    /* write all variables in one vard call */
    err = ncmpi_put_vard_all(ncid, varids[0], filetype, buffer, buflen*nvars,
                             mpitype); ERR
    io_timing = MPI_Wtime() - io_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    close_timing = MPI_Wtime();

    MPI_Type_free(&filetype);
    err = ncmpi_inq_put_size(ncid, &put_size); ERR
    err = ncmpi_close(ncid); ERR

    free(varids);
    free(buffer);
    timing = MPI_Wtime();
    close_timing = timing - close_timing;
    total_timing = timing - total_timing;

    MPI_Reduce(&nreqs, &max_nreqs, 1, MPI_INT, MPI_MAX, 0, comm);
    MPI_Reduce(&put_size, &total_size, 1, MPI_OFFSET, MPI_SUM, 0, comm);

    MPI_Reduce(&open_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    open_timing = max_timing;
    MPI_Reduce(&pre_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    pre_timing = max_timing;
    MPI_Reduce(&io_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    io_timing = max_timing;
    MPI_Reduce(&close_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    close_timing = max_timing;
    MPI_Reduce(&total_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    total_timing = max_timing;

    /* check if there is any PnetCDF internal malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, comm);
        if (rank == 0 && sum_size > 0) {
            printf("-----------------------------------------------------------\n");
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
        }
    }
    MPI_Offset m_alloc=0, max_alloc;
    ncmpi_inq_malloc_max_size(&m_alloc);
    MPI_Reduce(&m_alloc, &max_alloc, 1, MPI_OFFSET, MPI_MAX, 0, comm);
    if (rank == 0) {
        printf("-----------------------------------------------------------\n");
        printf("MAX heap memory allocated by PnetCDF internally is %.2f MiB\n",
               (float)max_alloc/1048576);
        printf("Total number of variables          = %d\n",nvars);
        printf("Total write amount                 = %.2f MiB = %.2f GiB\n",
               (double)total_size/1048576,(double)total_size/1073741824);
        printf("Max number of requests             = %d\n",max_nreqs);
        printf("Max Time of open + metadata define = %.4f sec\n",open_timing);
        printf("Max Time of I/O preparing          = %.4f sec\n",pre_timing);
        printf("Max Time of ncmpi_put_vard         = %.4f sec\n",io_timing);
        printf("Max Time of close                  = %.4f sec\n",close_timing);
        printf("Max Time of TOTAL                  = %.4f sec\n",total_timing);
        printf("I/O bandwidth                      = %.4f MiB/sec\n",
               (double)total_size/1048576.0/total_timing);
        if (verbose) print_info(&info_used);
    }
fn_exit:
    if (info_used != MPI_INFO_NULL) MPI_Info_free(&info_used);
    unlink(outfname);
    MPI_Barrier(comm);
    return nerrs;
}

/*----< run_varn() >--------------------------------------------------------*/
static
int run_varn(char       *out_dir,   /* output folder name */
             MPI_Info    info,
             int         ndims,     /* number of dimensions of variables */
             MPI_Offset *dims,      /* dimension lengths */
             int         nvars,     /* number of variables */
             int         nreqs,     /* number of request per variable */
             int        *disps,     /* [nreqs] request's displacements */
             int        *blocklens) /* [nreqs] request's block lengths */
{
    char outfname[512];
    int i, j, err, nerrs=0, rank, cmode, ncid, *varids, *dimids;
    int max_nreqs;
    size_t buflen=0;
    itype *buffer, *buf_ptr;
    double pre_timing, open_timing, post_timing, wait_timing, close_timing;
    double timing, total_timing,  max_timing;
    MPI_Offset **starts, **counts, put_size, total_size;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Barrier(comm); /*-----------------------------------------*/
    total_timing = open_timing = MPI_Wtime();

    /* set output file name */
    sprintf(outfname, "%s/testfile_varn.nc",out_dir);

    /* create a new CDF-5 file for writing */
    cmode = NC_CLOBBER|NC_64BIT_DATA;
    err = ncmpi_create(comm, outfname, cmode, info, &ncid); ERR

    /* define dimensions */
    dimids = (int*) malloc(ndims * sizeof(int));
    for (i=0; i<ndims; i++) {
        char dim_name[32];
        sprintf(dim_name, "dim_len_%d", i);
        err = ncmpi_def_dim(ncid, dim_name, dims[i], &dimids[i]); ERR
    }

    /* define variables */
    varids = (int*) malloc(nvars * sizeof(int));
    for (i=0; i<nvars; i++) {
        char varname[128];
        sprintf(varname, "var%04d",i);
        err = ncmpi_def_var(ncid, varname, xtype, ndims, dimids, &varids[i]);
        ERR
    }
    free(dimids);
    err = ncmpi_enddef(ncid); ERR
    open_timing = MPI_Wtime() - open_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    pre_timing = MPI_Wtime();

    /* construct varn API arguments starts[][] and counts[][] */
    starts = (MPI_Offset**) malloc(2 * nreqs * sizeof(MPI_Offset*));
    counts = starts + nreqs;
    starts[0] = (MPI_Offset*) malloc(2 * nreqs * ndims * sizeof(MPI_Offset));
    counts[0] = starts[0] + nreqs * ndims;
    for (i=1; i<nreqs; i++) {
        starts[i] = starts[i-1] + ndims;
        counts[i] = counts[i-1] + ndims;
    }

    for (i=0; i<nreqs; i++) {
        MPI_Offset disp=disps[i];
        for (j=ndims-1; j>=0; j--) {
            starts[i][j] = disp % dims[j];
            disp /= dims[j];
            counts[i][j] = 1;
        }
        counts[i][ndims-1] = blocklens[i]; /* each blocklens[i] is no bigger than dims[ndims-1] */
    }

    MPI_Comm_rank(comm, &rank);

    /* allocate and initialize write buffer */
    for (i=0; i<nreqs; i++) buflen += blocklens[i];
    buffer = (itype*) malloc(buflen * nvars * sizeof(itype));
    for (i=0; i<buflen*nvars; i++) buffer[i] = rank;

    pre_timing = MPI_Wtime() - pre_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    post_timing = MPI_Wtime();

    buf_ptr = buffer;
    for (i=0; i<nvars; i++) {
        err = ncmpi_iput_varn(ncid, varids[i], nreqs, starts, counts,
                              buf_ptr, -1, mpitype, NULL); ERR
        buf_ptr += buflen;
    }
    post_timing = MPI_Wtime() - post_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    wait_timing = MPI_Wtime();

    err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); ERR
    wait_timing = MPI_Wtime() - wait_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    close_timing = MPI_Wtime();

    err = ncmpi_inq_put_size(ncid, &put_size); ERR
    err = ncmpi_close(ncid); ERR

    free(varids);
    free(buffer);
    free(starts[0]); free(starts);

    timing = MPI_Wtime();
    close_timing = timing - close_timing;
    total_timing = timing - total_timing;

    MPI_Reduce(&nreqs, &max_nreqs, 1, MPI_INT, MPI_MAX, 0, comm);
    MPI_Reduce(&put_size, &total_size, 1, MPI_OFFSET, MPI_SUM, 0, comm);

    MPI_Reduce(&open_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    open_timing = max_timing;
    MPI_Reduce(&pre_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    pre_timing = max_timing;
    MPI_Reduce(&post_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    post_timing = max_timing;
    MPI_Reduce(&wait_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    wait_timing = max_timing;
    MPI_Reduce(&close_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    close_timing = max_timing;
    MPI_Reduce(&total_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    total_timing = max_timing;

    /* check if there is any PnetCDF internal malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, comm);
        if (rank == 0 && sum_size > 0) {
            printf("-----------------------------------------------------------\n");
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
        }
    }
    MPI_Offset m_alloc=0, max_alloc;
    ncmpi_inq_malloc_max_size(&m_alloc);
    MPI_Reduce(&m_alloc, &max_alloc, 1, MPI_OFFSET, MPI_MAX, 0, comm);
    if (rank == 0) {
        printf("-----------------------------------------------------------\n");
        printf("MAX heap memory allocated by PnetCDF internally is %.2f MiB\n",
               (float)max_alloc/1048576);
        printf("Total number of variables          = %d\n",nvars);
        printf("Total write amount                 = %.2f MiB = %.2f GiB\n",
               (double)total_size/1048576,(double)total_size/1073741824);
        printf("Max number of requests             = %d\n",max_nreqs);
        printf("Max Time of open + metadata define = %.4f sec\n",open_timing);
        printf("Max Time of I/O preparing          = %.4f sec\n",pre_timing);
        printf("Max Time of ncmpi_iput_varn        = %.4f sec\n",post_timing);
        printf("Max Time of ncmpi_wait_all         = %.4f sec\n",wait_timing);
        printf("Max Time of close                  = %.4f sec\n",close_timing);
        printf("Max Time of TOTAL                  = %.4f sec\n",total_timing);
        printf("I/O bandwidth                      = %.4f MiB/sec\n",
               (double)total_size/1048576.0/total_timing);
    }
fn_exit:
    unlink(outfname);
    MPI_Barrier(comm);
    return nerrs;
}

/*----< run_vara() >--------------------------------------------------------*/
static
int run_vara(char       *out_dir,   /* output folder name */
             MPI_Info    info,
             int         ndims,     /* number of dimensions of variables */
             MPI_Offset *dims,      /* dimension lengths */
             int         nvars,     /* number of variables */
             int         nreqs,     /* number of request per variable */
             int        *disps,     /* [nreqs] request's displacements */
             int        *blocklens) /* [nreqs] request's block lengths */
{
    char outfname[512];
    int i, j, err, nerrs=0, rank, cmode, ncid, *varids, *dimids;
    int max_nreqs;
    size_t buflen=0;
    itype *buffer, *buf_ptr;
    double pre_timing, open_timing, post_timing, wait_timing, close_timing;
    double timing, total_timing,  max_timing;
    MPI_Offset **starts, **counts, put_size, total_size;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Barrier(comm); /*-----------------------------------------*/
    total_timing = open_timing = MPI_Wtime();

    sprintf(outfname, "%s/testfile_vara.nc",out_dir);

    /* create a new CDF-5 file for writing */
    cmode = NC_CLOBBER|NC_64BIT_DATA;
    err = ncmpi_create(comm, outfname, cmode, info, &ncid); ERR

    /* define dimensions */
    dimids = (int*) malloc(ndims * sizeof(int));
    for (i=0; i<ndims; i++) {
        char dim_name[32];
        sprintf(dim_name, "dim_len_%d", i);
        err = ncmpi_def_dim(ncid, dim_name, dims[i], &dimids[i]); ERR
    }

    /* define variables */
    varids = (int*) malloc(nvars * sizeof(int));
    for (i=0; i<nvars; i++) {
        char varname[128];
        sprintf(varname, "var%04d",i);
        err = ncmpi_def_var(ncid, varname, xtype, ndims, dimids, &varids[i]);
        ERR
    }
    free(dimids);
    err = ncmpi_enddef(ncid); ERR
    open_timing = MPI_Wtime() - open_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    pre_timing = MPI_Wtime();

    /* construct starts[][] and counts[][] */
    starts = (MPI_Offset**) malloc(2 * nreqs * sizeof(MPI_Offset*));
    counts = starts + nreqs;
    starts[0] = (MPI_Offset*) malloc(2 * nreqs * ndims * sizeof(MPI_Offset));
    counts[0] = starts[0] + nreqs * ndims;
    for (i=1; i<nreqs; i++) {
        starts[i] = starts[i-1] + ndims;
        counts[i] = counts[i-1] + ndims;
    }

    for (i=0; i<nreqs; i++) {
        MPI_Offset disp=disps[i];
        for (j=ndims-1; j>=0; j--) {
            starts[i][j] = disp % dims[j];
            disp /= dims[j];
            counts[i][j] = 1;
        }
        counts[i][ndims-1] = blocklens[i]; /* each blocklens[i] is no bigger than dims[ndims-1] */
    }

    MPI_Comm_rank(comm, &rank);

    /* allocate write buffer */
    for (i=0; i<nreqs; i++) buflen += blocklens[i];
    buffer = (itype*) malloc(buflen * nvars * sizeof(itype));
    for (i=0; i<buflen*nvars; i++) buffer[i] = rank;

    pre_timing = MPI_Wtime() - pre_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    post_timing = MPI_Wtime();

    buf_ptr = buffer;
    for (i=0; i<nvars; i++) {
        for (j=0; j<nreqs; j++) {
            err = ncmpi_iput_vara(ncid, varids[i], starts[j], counts[j],
                                  buf_ptr, -1, mpitype, NULL); ERR
            buf_ptr += blocklens[j];
        }
    }
    post_timing = MPI_Wtime() - post_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    wait_timing = MPI_Wtime();

    err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); ERR
    wait_timing = MPI_Wtime() - wait_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    close_timing = MPI_Wtime();

    err = ncmpi_inq_put_size(ncid, &put_size); ERR
    err = ncmpi_close(ncid); ERR

    free(varids);
    free(buffer);
    free(starts[0]); free(starts);

    timing = MPI_Wtime();
    close_timing = timing - close_timing;
    total_timing = timing - total_timing;

    MPI_Reduce(&nreqs, &max_nreqs, 1, MPI_INT, MPI_MAX, 0, comm);
    MPI_Reduce(&put_size, &total_size, 1, MPI_OFFSET, MPI_SUM, 0, comm);

    MPI_Reduce(&open_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    open_timing = max_timing;
    MPI_Reduce(&pre_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    pre_timing = max_timing;
    MPI_Reduce(&post_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    post_timing = max_timing;
    MPI_Reduce(&wait_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    wait_timing = max_timing;
    MPI_Reduce(&close_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    close_timing = max_timing;
    MPI_Reduce(&total_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    total_timing = max_timing;

    /* check if there is any PnetCDF internal malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, comm);
        if (rank == 0 && sum_size > 0) {
            printf("-----------------------------------------------------------\n");
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
        }
    }
    MPI_Offset m_alloc=0, max_alloc;
    ncmpi_inq_malloc_max_size(&m_alloc);
    MPI_Reduce(&m_alloc, &max_alloc, 1, MPI_OFFSET, MPI_MAX, 0, comm);
    if (rank == 0) {
        printf("-----------------------------------------------------------\n");
        printf("MAX heap memory allocated by PnetCDF internally is %.2f MiB\n",
               (float)max_alloc/1048576);
        printf("Total number of variables          = %d\n",nvars);
        printf("Total write amount                 = %.2f MiB = %.2f GiB\n",
               (double)total_size/1048576,(double)total_size/1073741824);
        printf("Max number of requests             = %d\n",max_nreqs);
        printf("Max Time of open + metadata define = %.4f sec\n",open_timing);
        printf("Max Time of I/O preparing          = %.4f sec\n",pre_timing);
        printf("Max Time of ncmpi_iput_vara        = %.4f sec\n",post_timing);
        printf("Max Time of ncmpi_wait_all         = %.4f sec\n",wait_timing);
        printf("Max Time of close                  = %.4f sec\n",close_timing);
        printf("Max Time of TOTAL                  = %.4f sec\n",total_timing);
        printf("I/O bandwidth                      = %.4f MiB/sec\n",
               (double)total_size/1048576.0/total_timing);
    }
fn_exit:
    unlink(outfname);
    MPI_Barrier(comm);
    return nerrs;
}

/*----< usage() >------------------------------------------------------------*/
static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [-h] [-q] [-n nvars] [-o output_file] input_file\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode\n"
    "       [-n nvars]: number of variables (default 1)\n"
    "       [-o output_file]: output file name (default ./testfile.nc)\n"
    "       input_file: intput netCDF file name\n";
    fprintf(stderr, help, argv0);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char** argv)
{
    extern int optind;
    char *infname, out_dir[1024];
    int i, ndims, rank, nprocs, err, nerrs=0, nvars=0;
    int contig_nreqs, *disps=NULL, *blocklens=NULL;
    MPI_Offset *dims;
    MPI_Info info=MPI_INFO_NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    out_dir[0] = '\0';
    verbose = 1;

    /* command-line arguments */
    while ((i = getopt(argc, argv, "hqn:o:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'n': nvars = atoi(optarg);
                      break;
            case 'o': strcpy(out_dir, optarg);
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

    /* input file conatins number of write requests and their file access
     * offsets (per array element) */
    infname = argv[optind];
    if (verbose && rank==0) printf("input file name =%s\n",infname);

    /* set the output folder name */
    if (out_dir[0] == '\0') {
        strcpy(out_dir, ".");
    }
    if (verbose && rank==0) printf("output folder name =%s\n",out_dir);

    /* set the number of variables (same size, type, pattern) */
    nvars = (nvars <= 0) ? 1 : nvars;
    if (verbose && rank==0) printf("number variables = %d\n",nvars);

    /* read I/O pattern from input file */
    err = read_io_pattern(verbose, infname, &ndims, &dims, &contig_nreqs,
                          &disps, &blocklens);
    if (err) goto fn_exit;

    /* set MPI-IO hints */
    MPI_Info_create(&info);
    MPI_Info_set(info, "romio_ds_write", "disable"); /* MPI-IO data sieving */
    MPI_Info_set(info, "romio_cb_write", "enable");  /* collective write */
    MPI_Info_set(info, "romio_no_indep_rw", "true"); /* no independent MPI-IO */

    /* set PnetCDF I/O hints */
    MPI_Info_set(info, "nc_var_align_size", "1"); /* no gap between variables */
    MPI_Info_set(info, "nc_in_place_swap", "enable"); /* in-place byte swap */

    if (!rank) printf("\n---- benchmarking vard API -----------------------\n");
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    nerrs += run_vard(out_dir, info, ndims, dims, nvars, contig_nreqs,
                      disps, blocklens);

    if (!rank) printf("\n---- benchmarking varn API -----------------------\n");
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    nerrs += run_varn(out_dir, info, ndims, dims, nvars, contig_nreqs,
                      disps, blocklens);

    if (!rank) printf("\n---- benchmarking vara API -----------------------\n");
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    nerrs += run_vara(out_dir, info, ndims, dims, nvars, contig_nreqs,
                      disps, blocklens);

fn_exit:
    if (info != MPI_INFO_NULL) MPI_Info_free(&info);

    if (dims != NULL) free(dims);
    if (disps != NULL) free(disps);
    if (blocklens != NULL) free(blocklens);
    MPI_Finalize();
    return (nerrs > 0);
}

