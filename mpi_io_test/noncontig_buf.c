/*
 * Copyright (C) 2018, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program evaluates the performance of MPI collective write given a
 * fileview containing a long list of noncontiguous, small layout in file.
 * Under the same fileview, it compares two cases. One uses a contiguous
 * allocated user buffer and the other noncontiguous buffer.
 */
#include <stdio.h>
#include <stdlib.h> /* strtoll() */
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() unlink() */
#include <limits.h>  /* INT_MAX */
#include <assert.h>

#include <mpi.h>
#include <pnetcdf.h>

static int verbose;

#ifndef MAX
#define MAX(a,b) ((a) > (b)) ? (a) : (b)
#endif
#ifndef MIN
#define MIN(a,b) ((a) < (b)) ? (a) : (b)
#endif

#define ERR { \
    if (err != NC_NOERR) { \
        printf("Error in %s:%d: %s\n", __FILE__, __LINE__, \
               ncmpi_strerror(err)); \
        nerrs++; \
        goto fn_exit; \
    } \
}

/*----< read_decomp() >------------------------------------------------------*/
static int
read_decomp(const char   *infname,
            MPI_Offset   *array_size,  /* OUT */
            int          *nreqs,       /* OUT */
            int         **blocklens,   /* OUT */
            MPI_Aint    **disps)       /* OUT */
{
    char name[128];
    int err, nerrs=0, rank, ncid, varid;
    int i, dimid, *offs=NULL, *lens=NULL;
    MPI_Offset num_procs, tmp;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* read from input file */
    err = ncmpi_open(MPI_COMM_WORLD, infname, NC_NOWRITE, MPI_INFO_NULL,
                     &ncid); ERR

    err = ncmpi_begin_indep_data(ncid); ERR

    /* read num_procs */
    err = ncmpi_inq_dimid(ncid, "num_procs", &dimid); ERR
    err = ncmpi_inq_dimlen(ncid, dimid, &num_procs); ERR

    /* read array_size */
    err = ncmpi_inq_dimid(ncid, "array_size", &dimid); ERR
    err = ncmpi_inq_dimlen(ncid, dimid, array_size); ERR

    if (rank >= num_procs) {
        nreqs = 0;
        blocklens = NULL;
        disps = NULL;
        err = ncmpi_close(ncid); ERR
        return 0;
    }

    /* read number of noncontiguous requests */
    sprintf(name, "nreq_%05d", rank);
    err = ncmpi_inq_dimid(ncid, name, &dimid); ERR
    err = ncmpi_inq_dimlen(ncid, dimid, &tmp); ERR
    *nreqs = (int)tmp;
    offs = (int*) malloc(*nreqs * sizeof(int));
    lens = (int*) malloc(*nreqs * sizeof(int));

    /* get offsets */
    sprintf(name, "off_%05d", rank);
    err = ncmpi_inq_varid(ncid, name, &varid); ERR
    err = ncmpi_get_var_int(ncid, varid, offs); ERR

    /* get lengths */
    sprintf(name, "len_%05d", rank);
    err = ncmpi_inq_varid(ncid, name, &varid); ERR
    err = ncmpi_get_var_int(ncid, varid, lens); ERR

    err = ncmpi_close(ncid); ERR

    /* cast int to MPI_Aint */
    *disps = (MPI_Aint*) malloc(*nreqs * sizeof(MPI_Aint));
    for (i=0; i<*nreqs; i++) (*disps)[i] = offs[i];

fn_exit:
    if (offs != NULL) free(offs);
    *blocklens = lens;
    return nerrs;
}

/*----< err_handler() >------------------------------------------------------*/
void err_handler(int err, char *err_msg)
{
    int rank, errorStringLen;
    char errorString[MPI_MAX_ERROR_STRING];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Error_string(err, errorString, &errorStringLen);
    if (err_msg == NULL) err_msg = "";
    printf("rank %d: MPI error (%s) : %s\n", rank, err_msg, errorString);
}

/*----< construct_filetype() >-----------------------------------------------*/
static int
construct_filetype(int             nvars,
                   MPI_Offset      array_size,
                   int             nreqs,
                   const int      *blocklens,
                   const MPI_Aint *disps,
                   MPI_Datatype   *filetype)  /* OUT */
{
    int i, j, err, *nlens, *lens;
    MPI_Aint *ndisps=NULL, *disp=NULL, var_off, extent;
    MPI_Datatype ftype;

    if (nvars == 0 || nreqs == 0) {
        *filetype = MPI_BYTE;
        return 0;
    }

    /* PnetCDF wait API checks if each varn is contiguous and if yes, it uses
     * one MPI_Type_create_hindexed() to construct fileviews for all varn,
     * across all variables.
     */
    ndisps = (MPI_Aint*) malloc(nvars * nreqs * sizeof(MPI_Aint));
    nlens = (int*) malloc(nvars * nreqs * sizeof(int));
    for (j=0; j<nreqs; j++) {
        ndisps[j] = disps[j] * sizeof(float);
        nlens[j]  = blocklens[j];
    }

    /* propagate into nvars */
    array_size *= sizeof(float);
    var_off = array_size;
    disp = ndisps + nreqs;
    lens = nlens + nreqs;
    for (i=1; i<nvars; i++) {
        memcpy(disp, ndisps, nreqs * sizeof(MPI_Aint));
        memcpy(lens, nlens,  nreqs * sizeof(int));
        for (j=0; j<nreqs; j++)
            disp[j] += var_off;
        disp += nreqs;
        lens += nreqs;
        var_off += array_size;
    }

    /* create filetype */
    err = MPI_Type_create_hindexed(nvars*nreqs, nlens, ndisps, MPI_FLOAT,
                                   &ftype);
    if (err != MPI_SUCCESS) {
        err_handler(err, "MPI_Type_create_hindexed()");
        return 1;
    }
    MPI_Type_commit(&ftype);

    extent = array_size * nvars;
    err = MPI_Type_create_resized(ftype, 0, extent, filetype);
    if (err != MPI_SUCCESS) {
        err_handler(err, "MPI_Type_create_resized()");
        return 1;
    }
    MPI_Type_commit(filetype);

    if (verbose) {
        int type_size;
        MPI_Type_size(*filetype, &type_size);
        printf("filetype size = %d\n", type_size);
    }

    MPI_Type_free(&ftype);
    free(nlens);
    free(ndisps);
    return 0;
}

/*----< construct_buftype() >------------------------------------------------*/
static int
construct_buftype(int           nvars,
                  int           nelems,
                  int           gap,
                  MPI_Datatype *buftype)  /* OUT */
{
    int i, err, *nlens;
    MPI_Aint *ndisps=NULL;

    if (nvars == 0 || nelems == 0) {
        *buftype = MPI_BYTE;
        return 0;
    }

    /* construct a noncontiguous buffer datatype */
    nlens = (int*) malloc(nvars * sizeof(int));
    ndisps = (MPI_Aint*) malloc(nvars * sizeof(MPI_Aint));
    nlens[0] = nelems;
    ndisps[0] = 0;
    for (i=1; i<nvars; i++) {
        nlens[i] = nelems;
        ndisps[i] = ndisps[i-1] + (nelems+gap) * sizeof(float);
    }

    err = MPI_Type_create_hindexed(nvars, nlens, ndisps, MPI_FLOAT, buftype);
    if (err != MPI_SUCCESS) {
        err_handler(err, "MPI_Type_create_hindexed()");
        return 1;
    }
    MPI_Type_commit(buftype);

    if (verbose) {
        int type_size;
        MPI_Type_size(*buftype, &type_size);
        printf("buftype size = %d\n", type_size);
    }

    free(nlens);
    free(ndisps);
    return 0;
}

/*----< test_mpi_io() >------------------------------------------------------*/
static int
test_mpi_io(const char  *outfname,
            int          nvars,
            MPI_Offset   array_size,
            int          nreqs,
            int         *blocklens,
            MPI_Aint    *disps)
{
    size_t buf_len;
    int i, err, rank, nprocs, tmp, pos, nelems, max_nelems, min_nelems;
    int gap, max_nreqs, min_nreqs, max_blocklen, min_blocklen, count;
    float *buf, *pack_buf;
    double t_pack, t_contig, t_noncontig, max_t;

    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Info info=MPI_INFO_NULL;
    MPI_File fh;
    MPI_Datatype filetype, buftype;
    MPI_Status status;
    MPI_Offset w_len, w_sum;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    /* construct filetype */
    err = construct_filetype(nvars, array_size, nreqs, blocklens, disps,
                             &filetype);
    if (err) return err;

    for (nelems=0, i=0; i<nreqs; i++) nelems += blocklens[i];
    w_len = nelems;

    /* construct buftype of a noncontiguous layout */
    gap = 10;
    err = construct_buftype(nvars, nelems, gap, &buftype);
    if (err) return err;

    /* allocate a contiguous buffer space */
    buf_len = (nelems + gap) * nvars;
    buf = (float*) malloc(buf_len * sizeof(float));
    for (i=0; i<buf_len; i++) buf[i] = rank;

    buf_len = nvars * nelems * sizeof(float);
    pack_buf = (float*)malloc(buf_len);

    MPI_Barrier(comm); /*-------------------------------------------*/
    t_pack = MPI_Wtime();
    pos = 0;
    if (nreqs > 0)
        MPI_Pack(buf, 1, buftype, pack_buf, buf_len, &pos, MPI_COMM_SELF);
    t_pack = MPI_Wtime() - t_pack;
    free(pack_buf);

    /* set MPI-IO hints */
    MPI_Info_create(&info);
    MPI_Info_set(info, "romio_ds_write", "disable"); /* MPI-IO data sieving */
    MPI_Info_set(info, "romio_cb_write", "enable");  /* collective write */
    MPI_Info_set(info, "romio_no_indep_rw", "true"); /* no independent MPI-IO */

    /* open output file */
    err = MPI_File_open(comm, outfname, MPI_MODE_CREATE | MPI_MODE_RDWR,
                        info, &fh);
    if (err != MPI_SUCCESS) {
        err_handler(err, "MPI_File_open()");
        return 1;
    }

    /* set fileview */
    err = MPI_File_set_view(fh, 0, MPI_BYTE, filetype, "native", info);
    if (err != MPI_SUCCESS) {
        err_handler(err, "MPI_File_set_view()");
        return 1;
    }
    MPI_Info_free(&info);

    /* collective write to file using a contiguous buffer ----------*/
    MPI_Barrier(comm); /*-------------------------------------------*/
    t_contig = MPI_Wtime();
    err = MPI_File_write_all(fh, buf, nvars * nelems, MPI_FLOAT, &status);
    if (err != MPI_SUCCESS) {
        err_handler(err, "MPI_File_write_all()");
        return 1;
    }
    t_contig = MPI_Wtime() - t_contig;
    MPI_Get_count(&status, MPI_FLOAT, &count);
    assert(count == nvars * nelems);

    /* collective write to file using a noncontiguous buffer -------*/
    MPI_Barrier(comm); /*-------------------------------------------*/
    t_noncontig = MPI_Wtime();
    err = MPI_File_write_all(fh, buf, 1, buftype, &status);
    if (err != MPI_SUCCESS) {
        err_handler(err, "MPI_File_write_all()");
        return 1;
    }
    t_noncontig = MPI_Wtime() - t_noncontig;
    MPI_Get_count(&status, MPI_FLOAT, &count);
    assert(count == nvars * nelems);

    /* close file */
    err = MPI_File_close(&fh);
    if (err != MPI_SUCCESS) {
        err_handler(err, "MPI_File_close()");
        return 1;
    }

    free(buf);
    if (buftype != MPI_BYTE) MPI_Type_free(&buftype);
    if (buftype != MPI_BYTE) MPI_Type_free(&filetype);

    /* output timing results ----------------------------------------*/
    max_blocklen = -1;
    min_blocklen = INT_MAX;
    for (i=0; i<nreqs; i++) {
        max_blocklen = MAX(blocklens[i], max_blocklen);
        min_blocklen = MIN(blocklens[i], min_blocklen);
    }
    MPI_Reduce(&nreqs, &max_nreqs, 1, MPI_INT, MPI_MAX, 0, comm);
    if (nreqs == 0) nreqs = INT_MAX;
    MPI_Reduce(&nreqs, &min_nreqs, 1, MPI_INT, MPI_MIN, 0, comm);
    MPI_Reduce(&max_blocklen, &tmp, 1, MPI_INT, MPI_MAX, 0, comm);
    max_blocklen = tmp;
    if (min_blocklen == 0) min_blocklen = INT_MAX;
    MPI_Reduce(&min_blocklen, &tmp, 1, MPI_INT, MPI_MIN, 0, comm);
    min_blocklen = tmp;
    MPI_Reduce(&nelems, &max_nelems, 1, MPI_INT, MPI_MAX, 0, comm);
    if (nelems == 0) nelems = INT_MAX;
    MPI_Reduce(&nelems, &min_nelems, 1, MPI_INT, MPI_MIN, 0, comm);
    MPI_Reduce(&t_contig, &max_t, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    t_contig = max_t;
    MPI_Reduce(&t_noncontig, &max_t, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    t_noncontig = max_t;
    MPI_Reduce(&t_pack, &max_t, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    t_pack = max_t;
    MPI_Reduce(&w_len, &w_sum, 1, MPI_OFFSET, MPI_SUM, 0, comm);

    if (rank == 0) {
        float mea;
        double total_size = nvars * w_sum * sizeof(float);
        printf("-----------------------------------------------------------\n");
        printf("Total number of MPI processes        = %d\n",nprocs);
        printf("Total number of variables            = %d\n",nvars);
        printf("Total write amount                   = %.2f MiB = %.2f GiB\n",
               total_size/1048576,total_size/1073741824);
        printf("Max no. noncontig requests per var   = %d\n",max_nreqs);
        printf("Min no. noncontig requests per var   = %d\n",min_nreqs);
        max_blocklen *= sizeof(float);
        printf("Max length of contig request         = %d bytes\n",max_blocklen);
        min_blocklen *= sizeof(float);
        printf("Min length of contig request         = %d bytes\n",min_blocklen);
        mea = (float)max_nelems * sizeof(float) / 1048576.0;
        printf("Max write amount per variable        = %.2f MiB\n",mea);
        mea = (float)min_nelems * sizeof(float) / 1048576.0;
        printf("Min write amount per variable        = %.2f MiB\n",mea);
        printf("Max write time when buf is contig    = %.4f sec\n",t_contig);
        printf("Max write time when buf is noncontig = %.4f sec\n",t_noncontig);
        printf("Max time of MPI_Pack()               = %.4f sec\n",t_pack);
        printf("-----------------------------------------------------------\n");
    }

    return 0;
}

/*----< usage() >------------------------------------------------------------*/
static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [OPTION]... [FILE]...\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode\n"
    "       [-n] number of variables\n"
    "       [-o] output file name (default \"./testfile\")\n"
    "       input_file: name of input netCDF file describing data decompositions\n";
    fprintf(stderr, help, argv0);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char** argv)
{
    extern int optind;
    char *infname, outfname[1024];
    int i, nvars=1, rank, nprocs, err, nreqs, *blocklens;
    MPI_Aint *disps;
    MPI_Offset array_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    verbose = 1;

    strcpy(outfname, "testfile");

    /* command-line arguments */
    while ((i = getopt(argc, argv, "hqn:o:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'n': nvars = atoi(optarg);
                      break;
            case 'o': strcpy(outfname, optarg);
                      break;
            case 'h':
            default:  if (rank==0) usage(argv[0]);
                      MPI_Finalize();
                      return 1;
        }

    if (argv[optind] == NULL) { /* input file is mandatory */
        if (rank == 0) usage(argv[0]);
        MPI_Finalize();
        return 1;
    }

    /* input file contains number of write requests and their file access
     * offsets (per array element) */
    infname = argv[optind];
    if (rank==0) {
        printf("input  file name = %s\n",infname);
        printf("output file name = %s\n",outfname);
    }

    /* read I/O decomposition from input file */
    err = read_decomp(infname, &array_size, &nreqs, &blocklens, &disps);
    if (err) goto all_exit;

    /* test MPI-IO */
    err = test_mpi_io(outfname, nvars, array_size, nreqs, blocklens, disps);
    
    if (blocklens != NULL) free(blocklens);
    if (disps != NULL) free(disps);

all_exit:
    MPI_Finalize();
    return (err > 0);
}

