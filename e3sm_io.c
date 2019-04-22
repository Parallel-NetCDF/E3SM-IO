/*********************************************************************
 *
 * Copyright (C) 2018, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program uses the E3SM I/O patterns recorded by the PIO library to
 * evaluate the performance of two PnetCDF APIs: ncmpi_vard_all(), and
 * ncmpi_iput_varn(). The E3SM I/O patterns consist of a large number of small,
 * noncontiguous requests on each MPI process, which presents a challenge for
 * achieving a good performance.
 *
 * See README.md for compile and run instructions.
 *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h> /* strtoll() */
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() unlink() */

#include <e3sm_io.h>

/* minimum PnetCDF version required is 1.10.0 */
#if (PNETCDF_VERSION_MAJOR*1000000 + PNETCDF_VERSION_MINOR*1000 + PNETCDF_VERSION_SUB < 1010000)
#error "PnetCDF 1.10.0 and later is required to build this program"
#endif

static int verbose; /* verbose mode to print additional messages on screen */
static int keep_outfile; /* whether to keep the output files when exits */

static void print_info(MPI_Info *info_used);

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

/*----< read_io_decomp() >-------------------------------------------------*/
/* read I/O decomposition file, infname. The contents of the file are, for
 * example, 866x72_16p.nc.
 *   % ncmpidump -h 866x72_16p.nc
 *   netcdf 866x72_16p {
 *   // file format: CDF-1
 *   dimensions:
 *           num_procs = 16 ;
 *           D3.max_nreqs = 4032 ;
 *           D2.max_nreqs = 56 ;
 *           D1.max_nreqs = 70 ;
 *   variables:
 *           int D3.nreqs(num_procs) ;
 *           int D3.offsets(num_procs, D3.max_nreqs) ;
 *           int D2.nreqs(num_procs) ;
 *           int D2.offsets(num_procs, D2.max_nreqs) ;
 *           int D1.nreqs(num_procs) ;
 *           int D1.offsets(num_procs, D1.max_nreqs) ;
 *
 *   // global attributes:
 *                   :dim_len_X = 866 ;
 *                   :dim_len_Y = 72 ;
 *                   :D3.max_nreqs = 4032 ;
 *                   :D3.min_nreqs = 3744 ;
 *                   :D2.max_nreqs = 56 ;
 *                   :D2.min_nreqs = 52 ;
 *                   :D1.max_nreqs = 70 ;
 *                   :D1.min_nreqs = 40 ;
 *   }
 */
static int
read_io_decomp(const char  *infname,
               const char  *label,       /* name label */
               int          ndims,       /* number of dimensions of decomposition */
               MPI_Offset  *dims,        /* [2] dimension lengths */
               int         *contig_nreqs,/* num of contiguous requests */
               int        **disps,       /* displacements */
               int        **blocklens)   /* lengths of contiguous request */
{
    char name[128];
    int err, nerrs=0, rank, nprocs, ncid, varid, proc_start, proc_numb;
    int i, j, k, nreqs, dimids[2];
    MPI_Offset num_procs, max_nreqs, start[2], count[2];

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* read start-count from input file */
    err = ncmpi_open(MPI_COMM_WORLD, infname, NC_NOWRITE, MPI_INFO_NULL,
                     &ncid); ERR

    err = ncmpi_inq_dimid(ncid, "num_procs", &dimids[0]); ERR
    sprintf(name, "%s.max_nreqs", label);
    err = ncmpi_inq_dimid(ncid, name, &dimids[1]); ERR
    err = ncmpi_inq_dimlen(ncid, dimids[0], &num_procs); ERR
    err = ncmpi_inq_dimlen(ncid, dimids[1], &max_nreqs); ERR

    /* ndims is either 1 or 2, decomposition dimensions, not variable dimensions */
    if (ndims > 1) {
        err = ncmpi_get_att_longlong(ncid, NC_GLOBAL, "dim_len_Y", &dims[0]); ERR
        err = ncmpi_get_att_longlong(ncid, NC_GLOBAL, "dim_len_X", &dims[1]); ERR
    }
    else {
        err = ncmpi_get_att_longlong(ncid, NC_GLOBAL, "dim_len_X", &dims[0]); ERR
    }

    /* num_procs is the number of processes used to generate the E3SM data
     * decomposition files. nprocs is the number of processes running this
     * benchmark. This benchmark allows the two to be different. When nprocs is
     * less than num_procs, some of nprocs processes will carry out the
     * requests from more than one of num_procs processes. The requests
     * responsible by this process starts from proc_start with the number
     * proc_numb. When nprocs is bigger than num_procs, then those processes
     * with rank ID >= num_procs will have no data to write and will simply
     * participate the collective I/O subroutines.
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
    sprintf(name, "%s.nreqs", label);
    err = ncmpi_inq_varid(ncid, name, &varid); ERR
    start[0] = proc_start;
    count[0] = proc_numb;
    int *num_reqs = (int*) malloc(proc_numb * sizeof(int));
    err = ncmpi_get_vara_int_all(ncid, varid, start, count, num_reqs); ERR

    /* calculate the total number of requests responsible by this process */
    nreqs = 0;
    for (i=0; i<proc_numb; i++) nreqs += num_reqs[i];

    /* read the starting offsets of all requests into disps[] */
    *disps = (int*) malloc(proc_numb * max_nreqs * sizeof(int));
    sprintf(name, "%s.offsets", label);
    err = ncmpi_inq_varid(ncid, name, &varid); ERR
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
        if ((*disps)[i] == (*disps)[i-1]+1 && (*disps)[i]%dims[ndims-1])
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

    /* no request for this process */
    if (nreqs == 0) *contig_nreqs = 0;

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

/*----< write_small_vars() >-------------------------------------------------*/
static int
write_small_vars(int          ncid,
                 int          vid,    /* starting variable ID */
                 int         *varids,
                 int          rec_no,
                 int          gap,
                 MPI_Offset   lev,
                 MPI_Offset   ilev,
                 MPI_Offset   nbnd,
                 MPI_Offset   nchars,
                 int        **int_buf,
                 char       **txt_buf,
                 double     **dbl_buf)
{
    int i, err, nerrs=0;
    MPI_Offset start[2], count[2];

    /* scalar and small variables are written by rank 0 only */
    i = vid;

    if (rec_no == 0) {
        /* lev */
        err = ncmpi_iput_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += lev + gap;
        /* hyam */
        err = ncmpi_iput_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += lev + gap;
        /* hybm */
        err = ncmpi_iput_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += lev + gap;
        /* P0 */
        err = ncmpi_iput_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += 1 + gap;
        /* ilev */
        err = ncmpi_iput_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += ilev + gap;
        /* hyai */
        err = ncmpi_iput_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += ilev + gap;
        /* hybi */
        err = ncmpi_iput_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += ilev + gap;
    }
    else
        i += 7;

    /* time */
    start[0] = rec_no;
    err = ncmpi_iput_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* date */
    start[0] = rec_no;
    err = ncmpi_iput_var1_int(ncid, varids[i++], start, *int_buf, NULL); ERR
    *int_buf += 1;
    /* datesec */
    start[0] = rec_no;
    err = ncmpi_iput_var1_int(ncid, varids[i++], start, *int_buf, NULL); ERR
    *int_buf += 1;
    /* time_bnds */
    start[0] = rec_no; start[1] = 0;
    count[0] = 1;      count[1] = nbnd;
    err = ncmpi_iput_vara_double(ncid, varids[i++], start, count, *dbl_buf, NULL); ERR
    *dbl_buf += nbnd + gap;
    /* date_written */
    start[0] = rec_no; start[1] = 0;
    count[0] = 1;      count[1] = nchars;
    err = ncmpi_iput_vara_text(ncid, varids[i++], start, count, *txt_buf, NULL); ERR
    *txt_buf += nchars;
    /* time_written */
    start[0] = rec_no; start[1] = 0;
    count[0] = 1;      count[1] = nchars;
    err = ncmpi_iput_vara_text(ncid, varids[i++], start, count, *txt_buf, NULL); ERR
    *txt_buf += nchars;

    if (rec_no == 0) {
        /* ndbase */
        err = ncmpi_iput_var_int(ncid, varids[i++], *int_buf, NULL); ERR
        *int_buf += 1;
        /* nsbase */
        err = ncmpi_iput_var_int(ncid, varids[i++], *int_buf, NULL); ERR
        *int_buf += 1;
        /* nbdate */
        err = ncmpi_iput_var_int(ncid, varids[i++], *int_buf, NULL); ERR
        *int_buf += 1;
        /* nbsec */
        err = ncmpi_iput_var_int(ncid, varids[i++], *int_buf, NULL); ERR
        *int_buf += 1;
        /* mdt */
        err = ncmpi_iput_var_int(ncid, varids[i++], *int_buf, NULL); ERR
        *int_buf += 1;
    }
    else
        i += 5;

    /* ndcur */
    start[0] = rec_no;
    err = ncmpi_iput_var1_int(ncid, varids[i++], start, *int_buf, NULL); ERR
    *int_buf += 1;
    /* nscur */
    start[0] = rec_no;
    err = ncmpi_iput_var1_int(ncid, varids[i++], start, *int_buf, NULL); ERR
    *int_buf += 1;
    /* co2vmr */
    start[0] = rec_no;
    err = ncmpi_iput_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* ch4vmr */
    start[0] = rec_no;
    err = ncmpi_iput_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* n2ovmr */
    start[0] = rec_no;
    err = ncmpi_iput_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* f11vmr */
    start[0] = rec_no;
    err = ncmpi_iput_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* f12vmr */
    start[0] = rec_no;
    err = ncmpi_iput_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* sol_tsi */
    start[0] = rec_no;
    err = ncmpi_iput_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* nsteph */
    start[0] = rec_no;
    err = ncmpi_iput_var1_int(ncid, varids[i++], start, *int_buf, NULL); ERR
    *int_buf += 1;
fn_exit:
    return err;
}

/*----< write_small_vars() >-------------------------------------------------*/
static int
read_small_vars(int          ncid,
                 int          vid,    /* starting variable ID */
                 int         *varids,
                 int          rec_no,
                 int          gap,
                 MPI_Offset   lev,
                 MPI_Offset   ilev,
                 MPI_Offset   nbnd,
                 MPI_Offset   nchars,
                 int        **int_buf,
                 char       **txt_buf,
                 double     **dbl_buf)
{
    int i, err, nerrs=0;
    MPI_Offset start[2], count[2];

    /* scalar and small variables are written by rank 0 only */
    i = vid;

    if (rec_no == 0) {
        /* lev */
        err = ncmpi_iget_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += lev + gap;
        /* hyam */
        err = ncmpi_iget_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += lev + gap;
        /* hybm */
        err = ncmpi_iget_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += lev + gap;
        /* P0 */
        err = ncmpi_iget_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += 1 + gap;
        /* ilev */
        err = ncmpi_iget_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += ilev + gap;
        /* hyai */
        err = ncmpi_iget_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += ilev + gap;
        /* hybi */
        err = ncmpi_iget_var_double(ncid, varids[i++], *dbl_buf, NULL); ERR
        *dbl_buf += ilev + gap;
    }
    else
        i += 7;

    /* time */
    start[0] = rec_no;
    err = ncmpi_iget_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* date */
    start[0] = rec_no;
    err = ncmpi_iget_var1_int(ncid, varids[i++], start, *int_buf, NULL); ERR
    *int_buf += 1;
    /* datesec */
    start[0] = rec_no;
    err = ncmpi_iget_var1_int(ncid, varids[i++], start, *int_buf, NULL); ERR
    *int_buf += 1;
    /* time_bnds */
    start[0] = rec_no; start[1] = 0;
    count[0] = 1;      count[1] = nbnd;
    err = ncmpi_iget_vara_double(ncid, varids[i++], start, count, *dbl_buf, NULL); ERR
    *dbl_buf += nbnd + gap;
    /* date_written */
    start[0] = rec_no; start[1] = 0;
    count[0] = 1;      count[1] = nchars;
    err = ncmpi_iget_vara_text(ncid, varids[i++], start, count, *txt_buf, NULL); ERR
    *txt_buf += nchars;
    /* time_written */
    start[0] = rec_no; start[1] = 0;
    count[0] = 1;      count[1] = nchars;
    err = ncmpi_iget_vara_text(ncid, varids[i++], start, count, *txt_buf, NULL); ERR
    *txt_buf += nchars;

    if (rec_no == 0) {
        /* ndbase */
        err = ncmpi_iget_var_int(ncid, varids[i++], *int_buf, NULL); ERR
        *int_buf += 1;
        /* nsbase */
        err = ncmpi_iget_var_int(ncid, varids[i++], *int_buf, NULL); ERR
        *int_buf += 1;
        /* nbdate */
        err = ncmpi_iget_var_int(ncid, varids[i++], *int_buf, NULL); ERR
        *int_buf += 1;
        /* nbsec */
        err = ncmpi_iget_var_int(ncid, varids[i++], *int_buf, NULL); ERR
        *int_buf += 1;
        /* mdt */
        err = ncmpi_iget_var_int(ncid, varids[i++], *int_buf, NULL); ERR
        *int_buf += 1;
    }
    else
        i += 5;

    /* ndcur */
    start[0] = rec_no;
    err = ncmpi_iget_var1_int(ncid, varids[i++], start, *int_buf, NULL); ERR
    *int_buf += 1;
    /* nscur */
    start[0] = rec_no;
    err = ncmpi_iget_var1_int(ncid, varids[i++], start, *int_buf, NULL); ERR
    *int_buf += 1;
    /* co2vmr */
    start[0] = rec_no;
    err = ncmpi_iget_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* ch4vmr */
    start[0] = rec_no;
    err = ncmpi_iget_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* n2ovmr */
    start[0] = rec_no;
    err = ncmpi_iget_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* f11vmr */
    start[0] = rec_no;
    err = ncmpi_iget_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* f12vmr */
    start[0] = rec_no;
    err = ncmpi_iget_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* sol_tsi */
    start[0] = rec_no;
    err = ncmpi_iget_var1_double(ncid, varids[i++], start, *dbl_buf, NULL); ERR
    *dbl_buf += 1 + gap;
    /* nsteph */
    start[0] = rec_no;
    err = ncmpi_iget_var1_int(ncid, varids[i++], start, *int_buf, NULL); ERR
    *int_buf += 1;
fn_exit:
    return err;
}

#define SET_TYPE(kind) { \
    var_types[i] = type[kind]; \
    err = ncmpi_inq_varoffset(ncid, varids[i], &var_offset); ERR \
    var_disps[i] = var_offset - offset_rec; \
    if (kind == 2) { \
        nreqs += nreqs_D2; \
        if (i < nvars-1) buf_disps[i+1] = buf_disps[i] + (nelems_D2+gap) * sizeof(dtype); \
        buf_blocklens[i] = nelems_D2; \
    } else { /* kind == 3 */ \
        nreqs += nreqs_D3; \
        if (i < nvars-1) buf_disps[i+1] = buf_disps[i] + (nelems_D3+gap) * sizeof(dtype); \
        buf_blocklens[i] = nelems_D3; \
    } \
    i++; \
}

#define SET_TYPES(kind, num) \
    for (j=0; j<num; j++) { \
        var_types[i] = type[kind]; \
        err = ncmpi_inq_varoffset(ncid, varids[i], &var_offset); ERR \
        var_disps[i] = var_offset - offset_rec; \
        if (kind == 2) { \
            nreqs += nreqs_D2; \
            if (i < nvars-1) buf_disps[i+1] = buf_disps[i] + (nelems_D2+gap) * sizeof(dtype); \
            buf_blocklens[i] = nelems_D2; \
        } else { /* kind == 3 */ \
            nreqs += nreqs_D3; \
            if (i < nvars-1) buf_disps[i+1] = buf_disps[i] + (nelems_D3+gap) * sizeof(dtype); \
            buf_blocklens[i] = nelems_D3; \
        } \
        i++; \
    }

/*----< run_vard() >--------------------------------------------------------*/
static
int run_vard(char       *out_dir,      /* output folder name */
             char       *outfile,      /* output file name */
             int         nvars,        /* number of variables 408 or 51 */
             int         num_recs,     /* number of records */
             int         noncontig_buf,/* whether to us noncontiguous buffer */
             MPI_Info    info,
             MPI_Offset *dims,         /* [2] dimension lengths */
             int         nreqs_D1,     /* no. request in decomposition 1 */
             int        *disps_D1,     /* [nreqs_D1] request's displacements */
             int        *blocklens_D1, /* [nreqs_D1] request's block lengths */
             int         nreqs_D2,     /* no. request in decomposition 2 */
             int        *disps_D2,     /* [nreqs_D2] request's displacements */
             int        *blocklens_D2, /* [nreqs_D2] request's block lengths */
             int         nreqs_D3,     /* no. request in decomposition 2 */
             int        *disps_D3,     /* [nreqs_D3] request's displacements */
             int        *blocklens_D3) /* [nreqs_D3] request's block lengths */
{
    char outfname[512], txt_buf[16], *txt_buf_ptr;
    int i, j, k, err, nerrs=0, rank, ncid, cmode, *varids;
    int *var_blocklens, *buf_blocklens, nreqs, max_nreqs, rec_no, gap=0;
    int int_buf[10], *int_buf_ptr;
    size_t fix_buflen, dbl_buflen, rec_buflen;
    size_t nelems_D1, nelems_D2, nelems_D3;
    dtype *rec_buf;
    double *dbl_buf, *dbl_buf_ptr;
    double pre_timing, open_timing, io_timing, close_timing;
    double timing, total_timing,  max_timing;
    MPI_Aint *var_disps, *buf_disps;
    MPI_Offset tmp, metadata_size, rec_size, put_size, total_size;
    MPI_Offset offset_fix, offset_rec, var_offset;
    MPI_Datatype *var_types, type[4], *filetype_rec, filetype_dbl;
    MPI_Datatype buftype_rec, buftype_dbl;
    MPI_Info info_used=MPI_INFO_NULL;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Barrier(comm); /*-----------------------------------------*/
    total_timing = pre_timing = MPI_Wtime();

    MPI_Comm_rank(comm, &rank);

    if (noncontig_buf) gap = 10;

    varids = (int*) malloc(nvars * sizeof(int));

    /* allocate arrays for constructing fileview and buffer type */
    var_types = (MPI_Datatype*) malloc(nvars * sizeof(MPI_Datatype));
    var_blocklens = (int*) malloc(nvars * 2 * sizeof(int));
    buf_blocklens = var_blocklens + nvars;
    var_disps = (MPI_Aint*) malloc(nvars * 2 * sizeof(MPI_Aint));
    buf_disps = var_disps + nvars;

    /* define MPI datatypes for 4 kinds from 3 decompositions */
    MPI_Type_indexed(nreqs_D1, blocklens_D1, disps_D1, MPI_DOUBLE, &type[0]);
    MPI_Type_commit(&type[0]);
    MPI_Type_indexed(nreqs_D2, blocklens_D2, disps_D2, MPI_DOUBLE, &type[1]);
    MPI_Type_commit(&type[1]);
    MPI_Type_indexed(nreqs_D2, blocklens_D2, disps_D2, REC_DTYPE, &type[2]);
    MPI_Type_commit(&type[2]);
    MPI_Type_indexed(nreqs_D3, blocklens_D3, disps_D3, REC_DTYPE, &type[3]);
    MPI_Type_commit(&type[3]);

    /* number of variable elements from 3 decompositions */
    for (j=0; j<nvars; j++) var_blocklens[j] = 1;
    nelems_D1 = nelems_D2 = nelems_D3 = 0;
    for (k=0; k<nreqs_D1; k++) nelems_D1 += blocklens_D1[k];
    for (k=0; k<nreqs_D2; k++) nelems_D2 += blocklens_D2[k];
    for (k=0; k<nreqs_D3; k++) nelems_D3 += blocklens_D3[k];

    if (verbose && rank == 0)
        printf("nelems_D1=%zd nelems_D2=%zd nelems_D3=%zd\n",
               nelems_D1,nelems_D2,nelems_D3);

    /* allocate and initialize write buffer for small variables */
    dbl_buflen = nelems_D2*2 + nelems_D1 + gap*3
               + 3 * dims[0] + 3 * (dims[0]+1) + 8 + 2 + 20 * gap;

    dbl_buf = (double*) malloc(dbl_buflen * sizeof(double));
    for (j=0; j<dbl_buflen; j++) dbl_buf[j] = rank;

    for (j=0; j<10; j++) int_buf[j] = rank;
    for (j=0; j<16; j++) txt_buf[j] = 'a' + rank;

    /* allocate and initialize write buffer for large variables */
    if (nvars == 408)
        rec_buflen = nelems_D2 * 315 + nelems_D3 * 63 + (315+63) * gap;
    else
        rec_buflen = nelems_D2 * 20 + nelems_D3 + (20+1) * gap;

    rec_buf = (dtype*) malloc(rec_buflen * sizeof(dtype));
    for (i=0; i<rec_buflen; i++) rec_buf[i] = rank;

    pre_timing = MPI_Wtime() - pre_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    open_timing = MPI_Wtime();

    /* set output file name */
    sprintf(outfname, "%s/%s",out_dir, outfile);

    /* create a new CDF-5 file for writing */
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(comm, outfname, cmode, info, &ncid); ERR

    /* define dimensions, variables, and attributes */
    if (nvars == 408) {
        /* for h0 file */
        err = e3sm_io_header(ncid, dims, nvars, varids); ERR
    }
    else {
        /* for h1 file */
        err = e3sm_io_header1(ncid, dims, nvars, varids); ERR
    }

    /* exit define mode and enter data mode */
    err = ncmpi_enddef(ncid); ERR

    /* I/O amount so far */
    err = ncmpi_inq_put_size(ncid, &metadata_size); ERR
    err = ncmpi_inq_file_info(ncid, &info_used); ERR
    open_timing = MPI_Wtime() - open_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    /* the first 3 variables are of type NC_DOUBLE -------------------*/
    i = 0;
    nreqs = 0;

    /* lat */
    var_types[i] = type[1];
    err = ncmpi_inq_varoffset(ncid, varids[i], &offset_fix); ERR
    var_disps[i] = 0;
    buf_disps[0] = 0;
    buf_blocklens[0] = nelems_D2;
    i++;
    nreqs += nreqs_D2;

    /* lon */
    var_types[i] = type[1];
    err = ncmpi_inq_varoffset(ncid, varids[i], &var_offset); ERR
    var_disps[i] = var_offset - offset_fix;
    buf_disps[i] = buf_disps[i-1] + (nelems_D2 + gap) * sizeof (double);
    buf_blocklens[1] = nelems_D2;
    i++;
    nreqs += nreqs_D2;

    /* area */
    var_types[i] = type[0];
    err = ncmpi_inq_varoffset(ncid, varids[i], &var_offset); ERR
    var_disps[i] = var_offset - offset_fix;
    buf_disps[i] = buf_disps[i-1] + (nelems_D2 + gap) * sizeof (double);
    buf_blocklens[2] = nelems_D1;
    i++;
    nreqs += nreqs_D1;
    fix_buflen = nelems_D2*2 + nelems_D1 + gap*3;

    /* skip next 27 small variables */
    i += 27;

    /* concatenate 3 var_types[] into filetype_dbl */
    MPI_Type_create_struct(3, var_blocklens, var_disps, var_types,
                           &filetype_dbl);
    MPI_Type_commit(&filetype_dbl);

    if (noncontig_buf) {
        /* construct buffer type for 3 variables */
        MPI_Type_create_hindexed(3, buf_blocklens, buf_disps, MPI_DOUBLE,
                                 &buftype_dbl);
        MPI_Type_commit(&buftype_dbl);
    }
    else {
        /* buffer type is contiguous */
        buftype_dbl = MPI_DOUBLE;
    }

    err = ncmpi_inq_varoffset(ncid, varids[i], &offset_rec); ERR
    err = ncmpi_inq_recsize(ncid, &rec_size); ERR
    buf_disps[i] = 0;

    if (nvars == 408) {
        SET_TYPE(2)       /* AEROD_v */
        SET_TYPES(3, 2)   /* ANRAIN and ANSNOW */
        SET_TYPES(2, 18)  /* AODABS ... ANSNOW */
        SET_TYPES(3, 2)   /* AQRAIN and AQSNOW */
        SET_TYPES(2, 6)   /* AQ_DMS ... AQ_SOAG */
        SET_TYPES(3, 5)   /* AREI ... CCN3 */
        SET_TYPES(2, 2)   /* CDNUMC and CLDHGH */
        SET_TYPES(3, 2)   /* CLDICE and CLDLIQ */
        SET_TYPES(2, 3)   /* CLDLOW ... CLDTOT */
        SET_TYPES(3, 4)   /* CLOUD ... DCQ */
        SET_TYPES(2, 11)  /* DF_DMS ... DSTSFMBL */
        SET_TYPE(3)       /* DTCOND */
        SET_TYPES(2, 2)   /* DTENDTH and DTENDTQ */
        SET_TYPES(3, 2)   /* EXTINCT and FICE */
        SET_TYPES(2, 7)   /* FLDS ... FLUTC */
        SET_TYPES(3, 4)   /* FREQI ... FREQS */
        SET_TYPES(2, 15)  /* FSDS ... ICEFRAC */
        SET_TYPES(3, 3)   /* ICIMR ... IWC */
        SET_TYPES(2, 2)   /* LANDFRAC and LHFLX */
        SET_TYPES(3, 5)   /* LINOZ_DO3 ... LINOZ_SSO3 */
        SET_TYPES(2, 3)   /* LINOZ_SZA ... LWCF */
        SET_TYPES(3, 12)  /* Mass_bc ... O3 */
        SET_TYPES(2, 2)   /* O3_SRF and OCNFRAC */
        SET_TYPE(3)       /* OMEGA */
        SET_TYPE(2)       /* OMEGA500 */
        SET_TYPE(3)       /* OMEGAT */
        SET_TYPES(2, 8)   /* PBLH ... PSL */
        SET_TYPE(3)       /* Q */
        SET_TYPES(2, 2)   /* QFLX and QREFHT */
        SET_TYPES(3, 3)   /* QRL ... RAINQM */
        SET_TYPE(2)       /* RAM1 */
        SET_TYPE(3)       /* RELHUM */
        SET_TYPES(2, 37)  /* SFDMS ... SNOWHLND */
        SET_TYPES(3, 2)   /* SNOWQM and SO2 */
        SET_TYPES(2, 10)  /* SO2_CLXF ... SWCF */
        SET_TYPE(3)       /* T */
        SET_TYPES(2, 19)  /* TAUGWX ... TVQ */
        SET_TYPE(3)       /* U */
        SET_TYPE(2)       /* U10 */
        SET_TYPES(3, 6)   /* UU ... VV */
        SET_TYPES(2, 3)   /* WD_H2O2 ... WD_SO2 */
        SET_TYPES(3, 3)   /* WSUB ... aero_water */
        SET_TYPES(2, 32)  /* airFV ... dst_c3SFWET */
        SET_TYPE(3)       /* hstobie_linoz */
        SET_TYPES(2, 129) /* mlip ... soa_c3SFWET */
    }
    else {
        SET_TYPES(2, 13)  /* CLDHGH ... T5 */
        SET_TYPE(3)       /* U */
        SET_TYPES(2, 7)   /* U250 ... Z500 */
    }

    if (noncontig_buf) {
        /* construct buffer type for record variables */
        MPI_Type_create_hindexed(nvars-30, buf_blocklens+30, buf_disps+30,
                                 REC_DTYPE, &buftype_rec);
        MPI_Type_commit(&buftype_rec);
    }
    else {
        /* all record variables are in a single contiguous buffer */
        buftype_rec = REC_DTYPE;
    }

    filetype_rec = (MPI_Datatype*)malloc(num_recs * sizeof(MPI_Datatype));
    for (j=0; j<num_recs; j++) {
        if (j > 0) {
            for (k=30; k<nvars; k++)
                var_disps[k] += rec_size;
        }
        /* concatenate nvars-30 var_types[] into filetype_rec[j] */
        MPI_Type_create_struct(nvars-30, var_blocklens+30, var_disps+30,
                               var_types+30, filetype_rec+j);
        MPI_Type_commit(filetype_rec+j);
    }

    for (j=0; j<4; j++) MPI_Type_free(&type[j]);
    free(var_types);
    free(var_disps);
    free(var_blocklens);

    pre_timing += MPI_Wtime() - timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    io_timing = MPI_Wtime();

    if (noncontig_buf) fix_buflen = rec_buflen = 1;
    else {
        if (nvars == 408)
            rec_buflen = nelems_D2 * 315 + nelems_D3 * 63;
        else
            rec_buflen = nelems_D2 * 20 + nelems_D3;
    }

    /* write first 3 NC_DOUBLE fixed-size variables in one vard call */
    err = ncmpi_put_vard_all(ncid, varids[0], filetype_dbl, dbl_buf,
                             fix_buflen, buftype_dbl); ERR

    for (rec_no=0; rec_no<num_recs; rec_no++) {
        i=3;
        dbl_buf_ptr = dbl_buf + nelems_D2*2 + nelems_D1 + gap*3;
        int_buf_ptr = int_buf;
        txt_buf_ptr = txt_buf;

        /* next 27 small variables are written by rank 0 only */
        if (rank == 0) {
            /* post nonblocking requests using ncmpi_iput_varn() */
            err = write_small_vars(ncid, i, varids, rec_no, gap, dims[0],
                                   dims[0]+1, 2, 8, &int_buf_ptr, &txt_buf_ptr,
                                   &dbl_buf_ptr); ERR
            nreqs += 27;
        }
        i += 27;

        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); ERR

        /* write remaining record variables in one vard call */
        err = ncmpi_put_vard_all(ncid, varids[30], filetype_rec[rec_no],
                                 rec_buf, rec_buflen, buftype_rec); ERR
    }
    io_timing = MPI_Wtime() - io_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    close_timing = MPI_Wtime();

    err = ncmpi_inq_put_size(ncid, &total_size); ERR
    put_size = total_size - metadata_size;
    err = ncmpi_close(ncid); ERR
    close_timing = MPI_Wtime() - close_timing;

    for (j=0; j<num_recs; j++) MPI_Type_free(filetype_rec+j);
    free(filetype_rec);
    MPI_Type_free(&filetype_dbl);

    if (noncontig_buf) {
        MPI_Type_free(&buftype_rec);
        MPI_Type_free(&buftype_dbl);
    }

    free(rec_buf);
    free(dbl_buf);
    free(varids);

    total_timing = MPI_Wtime() - total_timing;

    MPI_Reduce(&nreqs,         &max_nreqs,  1, MPI_INT,    MPI_MAX, 0, comm);
    MPI_Reduce(&put_size,      &tmp,        1, MPI_OFFSET, MPI_SUM, 0, comm);
    put_size = tmp;
    MPI_Reduce(&total_size,    &tmp,        1, MPI_OFFSET, MPI_SUM, 0, comm);
    total_size = tmp;
    MPI_Reduce(&open_timing,   &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    open_timing = max_timing;
    MPI_Reduce(&pre_timing,    &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    pre_timing = max_timing;
    MPI_Reduce(&io_timing,     &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
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
        printf("History output file                = %s\n", outfile);
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
        printf("I/O bandwidth (open-to-close)      = %.4f MiB/sec\n",
               (double)total_size/1048576.0/total_timing);
        printf("I/O bandwidth (write-only)         = %.4f MiB/sec\n",
               (double)put_size/1048576.0/io_timing);
        if (verbose) print_info(&info_used);
        printf("-----------------------------------------------------------\n");
    }
fn_exit:
    if (info_used != MPI_INFO_NULL) MPI_Info_free(&info_used);
    if (!keep_outfile) unlink(outfname);
    fflush(stdout);
    MPI_Barrier(comm);
    return nerrs;
}

#define FIX_STARTS_COUNTS(starts, counts, nreqs, disps, blocklens) { \
    starts = (MPI_Offset**) malloc(2 * nreqs * sizeof(MPI_Offset*)); \
    counts = starts + nreqs; \
    starts[0] = (MPI_Offset*) malloc(2 * nreqs * sizeof(MPI_Offset)); \
    counts[0] = starts[0] + nreqs; \
    for (i=1; i<nreqs; i++) { \
        starts[i] = starts[i-1] + 1; \
        counts[i] = counts[i-1] + 1; \
    } \
    for (i=0; i<nreqs; i++) { \
        starts[i][0] = disps[i]; \
        counts[i][0] = blocklens[i]; \
    } \
}

#define REC_STARTS_COUNTS(rec, starts, counts, nreqs, disps, blocklens) { \
    starts = (MPI_Offset**) malloc(2 * nreqs * sizeof(MPI_Offset*)); \
    counts = starts + nreqs; \
    starts[0] = (MPI_Offset*) malloc(2 * nreqs * ndims * sizeof(MPI_Offset)); \
    counts[0] = starts[0] + nreqs * ndims; \
    for (i=1; i<nreqs; i++) { \
        starts[i] = starts[i-1] + ndims; \
        counts[i] = counts[i-1] + ndims; \
    } \
    for (i=0; i<nreqs; i++) { \
        MPI_Offset disp=disps[i]; \
        for (k=1,j=ndims-1; j>0; j--,k--) { \
            starts[i][j] = disp % dims[k]; /* dims always 2D */ \
            disp /= dims[k]; \
            counts[i][j] = 1; \
        } \
        starts[i][0] = rec; /* record ID */ \
        counts[i][0] = 1;   /* one record only */ \
        /* each blocklens[i] is no bigger than dims[ndims-1] */ \
        counts[i][ndims-1] = blocklens[i]; \
    } \
}

#define POST_VARN(k, num, vid) \
    for (j=0; j<num; j++) { \
        err = ncmpi_iput_varn(ncid, vid+j, nreqs_D##k, starts_D##k, \
                              counts_D##k, rec_buf_ptr, -1, REC_DTYPE, NULL); \
        ERR \
        rec_buf_ptr += nelems_D##k + gap; \
        nreqs += nreqs_D##k; \
    }

#define POST_VARN_RD(k, num, vid) \
    for (j=0; j<num; j++) { \
        err = ncmpi_iget_varn(ncid, vid+j, nreqs_D##k, starts_D##k, \
                              counts_D##k, rec_buf_ptr, -1, REC_DTYPE, NULL); \
        ERR \
        rec_buf_ptr += nelems_D##k + gap; \
        nreqs += nreqs_D##k; \
    }

static int two_buf;

/*----< run_varn() >--------------------------------------------------------*/
static
int run_varn(char       *out_dir,      /* output folder name */
             char       *outfile,      /* output file name */
             int         nvars,        /* number of variables 408 or 51 */
             int         num_recs,     /* number of records */
             int         noncontig_buf,/* whether to us noncontiguous buffer */
             MPI_Info    info,
             MPI_Offset *dims,         /* [2] dimension lengths */
             int         nreqs_D1,     /* no. request in decomposition 1 */
             int        *disps_D1,     /* [nreqs_D1] request's displacements */
             int        *blocklens_D1, /* [nreqs_D1] request's block lengths */
             int         nreqs_D2,     /* no. request in decomposition 2 */
             int        *disps_D2,     /* [nreqs_D2] request's displacements */
             int        *blocklens_D2, /* [nreqs_D2] request's block lengths */
             int         nreqs_D3,     /* no. request in decomposition 3 */
             int        *disps_D3,     /* [nreqs_D3] request's displacements */
             int        *blocklens_D3, /* [nreqs_D3] request's block lengths */
             double     *dbl_bufp,
             dtype      *rec_bufp,
             char       *txt_buf,
             int        *int_buf)
{
    char outfname[512], *txt_buf_ptr;
    int i, j, k, err, nerrs=0, rank, ndims, ncid, cmode, *varids;
    int rec_no, gap=0, nreqs, max_nreqs, *int_buf_ptr;
    size_t dbl_buflen, rec_buflen;
    size_t nelems_D1, nelems_D2, nelems_D3;
    dtype *rec_buf, *rec_buf_ptr;
    double *dbl_buf, *dbl_buf_ptr;
    double pre_timing, open_timing, post_timing, wait_timing, close_timing;
    double timing, total_timing,  max_timing;
    MPI_Offset tmp, metadata_size, put_size, total_size;
    MPI_Offset **fix_starts_D1, **fix_counts_D1;
    MPI_Offset **fix_starts_D2, **fix_counts_D2;
    MPI_Offset **starts_D2, **counts_D2;
    MPI_Offset **starts_D3, **counts_D3;
    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Info info_used=MPI_INFO_NULL;

    MPI_Barrier(comm); /*-----------------------------------------*/
    total_timing = pre_timing = MPI_Wtime();

    open_timing = 0.0;
    post_timing = 0.0;
    wait_timing = 0.0;
    close_timing = 0.0;

    MPI_Comm_rank(comm, &rank);

    if (noncontig_buf) gap = 10;

    /* number of variable elements from 3 decompositions */
    nreqs = max_nreqs = 0;
    nelems_D1 = nelems_D2 = nelems_D3 = 0;
    for (k=0; k<nreqs_D1; k++) nelems_D1 += blocklens_D1[k];
    for (k=0; k<nreqs_D2; k++) nelems_D2 += blocklens_D2[k];
    for (k=0; k<nreqs_D3; k++) nelems_D3 += blocklens_D3[k];

    if (verbose && rank == 0)
        printf("nelems_D1=%zd nelems_D2=%zd nelems_D3=%zd\n",
               nelems_D1,nelems_D2,nelems_D3);

    /* construct varn API arguments starts[][] and counts[][] */
    ndims = 1;
    FIX_STARTS_COUNTS(fix_starts_D1, fix_counts_D1, nreqs_D1, disps_D1, blocklens_D1)
    FIX_STARTS_COUNTS(fix_starts_D2, fix_counts_D2, nreqs_D2, disps_D2, blocklens_D2)

    ndims = 2;
    REC_STARTS_COUNTS(0, starts_D2, counts_D2, nreqs_D2, disps_D2, blocklens_D2)
    ndims = 3;
    REC_STARTS_COUNTS(0, starts_D3, counts_D3, nreqs_D3, disps_D3, blocklens_D3)

    /* allocate and initialize write buffer for small variables */
    dbl_buflen = nelems_D2 * 2 + nelems_D1
               + 3 * dims[0] + 3 * (dims[0]+1) + 8 + 2
               + 20 * gap;
    if (dbl_bufp != NULL){
        dbl_buf = dbl_bufp;
    }
    else{
        dbl_buf = (double*) malloc(dbl_buflen * sizeof(double));
        for (i=0; i<dbl_buflen; i++) dbl_buf[i] = rank;
    }

    /* allocate and initialize write buffer for large variables */
    if (nvars == 408)
        rec_buflen = nelems_D2 * 315 + nelems_D3 * 63 + (315+63) * gap;
    else
        rec_buflen = nelems_D2 * 20 + nelems_D3 + (20+1) * gap;

    if (rec_bufp != NULL){
        rec_buf = rec_bufp;
    }
    else{
        rec_buf = (dtype*) malloc(rec_buflen * sizeof(dtype));
        for (i=0; i<rec_buflen; i++) rec_buf[i] = rank;

        for (i=0; i<10; i++) int_buf[i] = rank;

        for (i=0; i<16; i++) txt_buf[i] = 'a' + rank;
    }

    varids = (int*) malloc(nvars * sizeof(int));

    pre_timing = MPI_Wtime() - pre_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    /* set output file name */
    sprintf(outfname, "%s/%s",out_dir, outfile);

    /* create a new CDF-5 file for writing */
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(comm, outfname, cmode, info, &ncid); ERR

    /* define dimensions, variables, and attributes */
    if (nvars == 408) {
        /* for h0 file */
        err = e3sm_io_header(ncid, dims, nvars, varids); ERR
    }
    else {
        /* for h1 file */
        err = e3sm_io_header1(ncid, dims, nvars, varids); ERR
    }

    /* exit define mode and enter data mode */
    err = ncmpi_enddef(ncid); ERR

    /* I/O amount so far */
    err = ncmpi_inq_put_size(ncid, &metadata_size); ERR
    err = ncmpi_inq_file_info(ncid, &info_used); ERR
    open_timing += MPI_Wtime() - timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    i = 0;
    dbl_buf_ptr = dbl_buf;

    /* lat */
    err = ncmpi_iput_varn(ncid, varids[i++], nreqs_D2, fix_starts_D2, fix_counts_D2,
                          dbl_buf_ptr, nelems_D2, MPI_DOUBLE, NULL); ERR
    dbl_buf_ptr += nelems_D2 + gap;
    nreqs += nreqs_D2;

    /* lon */
    err = ncmpi_iput_varn(ncid, varids[i++], nreqs_D2, fix_starts_D2, fix_counts_D2,
                          dbl_buf_ptr, nelems_D2, MPI_DOUBLE, NULL); ERR
    dbl_buf_ptr += nelems_D2 + gap;
    nreqs += nreqs_D2;

    /* area */
    err = ncmpi_iput_varn(ncid, varids[i++], nreqs_D1, fix_starts_D1, fix_counts_D1,
                          dbl_buf_ptr, nelems_D1, MPI_DOUBLE, NULL); ERR
    dbl_buf_ptr += nelems_D1 + gap;
    nreqs += nreqs_D1;

    post_timing += MPI_Wtime() - timing;

    for (rec_no=0; rec_no<num_recs; rec_no++) {
        MPI_Barrier(comm); /*-----------------------------------------*/
        timing = MPI_Wtime();

        i=3;
        dbl_buf_ptr = dbl_buf + nelems_D2*2 + nelems_D1 + gap*3;
        int_buf_ptr = int_buf;
        txt_buf_ptr = txt_buf;

        /* next 27 small variables are written by rank 0 only */
        if (rank == 0) {
            nreqs += 27;
            /* post nonblocking requests using ncmpi_iput_varn() */
            err = write_small_vars(ncid, i, varids, rec_no, gap, dims[0],
                                   dims[0]+1, 2, 8, &int_buf_ptr, &txt_buf_ptr,
                                   &dbl_buf_ptr); ERR
        }
        i += 27;

        post_timing += MPI_Wtime() - timing;

        MPI_Barrier(comm); /*-----------------------------------------*/
        timing = MPI_Wtime();

        /* flush fixed-size and small variables */
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); ERR

        wait_timing += MPI_Wtime() - timing;

        MPI_Barrier(comm); /*-----------------------------------------*/
        timing = MPI_Wtime();

        /* high water mark of number of noncontiguous requests */
        if (nreqs > max_nreqs) max_nreqs = nreqs;
        nreqs = 0;

        rec_buf_ptr = rec_buf;

        for (j=0; j<nreqs_D2; j++) starts_D2[j][0] = rec_no;
        for (j=0; j<nreqs_D3; j++) starts_D3[j][0] = rec_no;

        if (nvars == 408) {
            if (two_buf) {
                /* write 2D variables */
                POST_VARN(2,   1,  30)   /* AEROD_v */
                POST_VARN(2,  18,  33)   /* AODABS ... ANSNOW */
                POST_VARN(2,   6,  53)   /* AQ_DMS ... AQ_SOAG */
                POST_VARN(2,   2,  64)   /* CDNUMC and CLDHGH */
                POST_VARN(2,   3,  68)   /* CLDLOW ... CLDTOT */
                POST_VARN(2,  11,  75)   /* DF_DMS ... DSTSFMBL */
                POST_VARN(2,   2,  87)   /* DTENDTH and DTENDTQ */
                POST_VARN(2,   7,  91)   /* FLDS ... FLUTC */
                POST_VARN(2,  15, 102)   /* FSDS ... ICEFRAC */
                POST_VARN(2,   2, 120)   /* LANDFRAC and LHFLX */
                POST_VARN(2,   3, 127)   /* LINOZ_SZA ... LWCF */
                POST_VARN(2,   2, 142)   /* O3_SRF and OCNFRAC */
                POST_VARN(2,   1, 145)   /* OMEGA500 */
                POST_VARN(2,   8, 147)   /* PBLH ... PSL */
                POST_VARN(2,   2, 156)   /* QFLX and QREFHT */
                POST_VARN(2,   1, 161)   /* RAM1 */
                POST_VARN(2,  37, 163)   /* SFDMS ... SNOWHLND */
                POST_VARN(2,  10, 202)   /* SO2_CLXF ... SWCF */
                POST_VARN(2,  19, 213)   /* TAUGWX ... TVQ */
                POST_VARN(2,   1, 233)   /* U10 */
                POST_VARN(2,   3, 240)   /* WD_H2O2 ... WD_SO2 */
                POST_VARN(2,  32, 246)   /* airFV ... dst_c3SFWET */
                POST_VARN(2, 129, 279)   /* mlip ... soa_c3SFWET */
                /* write 3D variables */
                POST_VARN(3,   2,  31)   /* ANRAIN and ANSNOW */
                POST_VARN(3,   2,  51)   /* AQRAIN and AQSNOW */
                POST_VARN(3,   5,  59)   /* AREI ... CCN3 */
                POST_VARN(3,   2,  66)   /* CLDICE and CLDLIQ */
                POST_VARN(3,   4,  71)   /* CLOUD ... DCQ */
                POST_VARN(3,   1,  86)   /* DTCOND */
                POST_VARN(3,   2,  89)   /* EXTINCT and FICE */
                POST_VARN(3,   4,  98)   /* FREQI ... FREQS */
                POST_VARN(3,   3, 117)   /* ICIMR ... IWC */
                POST_VARN(3,   5, 122)   /* LINOZ_DO3 ... LINOZ_SSO3 */
                POST_VARN(3,  12, 130)   /* Mass_bc ... O3 */
                POST_VARN(3,   1, 144)   /* OMEGA */
                POST_VARN(3,   1, 146)   /* OMEGAT */
                POST_VARN(3,   1, 155)   /* Q */
                POST_VARN(3,   3, 158)   /* QRL ... RAINQM */
                POST_VARN(3,   1, 162)   /* RELHUM */
                POST_VARN(3,   2, 200)   /* SNOWQM and SO2 */
                POST_VARN(3,   1, 212)   /* T */
                POST_VARN(3,   1, 232)   /* U */
                POST_VARN(3,   6, 234)   /* UU ... VV */
                POST_VARN(3,   3, 243)   /* WSUB ... aero_water */
                POST_VARN(3,   1, 278)   /* hstobie_linoz */
            } else {
                /* write variables in the same order as they defined */
                POST_VARN(2,   1,  30)   /* AEROD_v */
                POST_VARN(3,   2,  31)   /* ANRAIN and ANSNOW */
                POST_VARN(2,  18,  33)   /* AODABS ... ANSNOW */
                POST_VARN(3,   2,  51)   /* AQRAIN and AQSNOW */
                POST_VARN(2,   6,  53)   /* AQ_DMS ... AQ_SOAG */
                POST_VARN(3,   5,  59)   /* AREI ... CCN3 */
                POST_VARN(2,   2,  64)   /* CDNUMC and CLDHGH */
                POST_VARN(3,   2,  66)   /* CLDICE and CLDLIQ */
                POST_VARN(2,   3,  68)   /* CLDLOW ... CLDTOT */
                POST_VARN(3,   4,  71)   /* CLOUD ... DCQ */
                POST_VARN(2,  11,  75)   /* DF_DMS ... DSTSFMBL */
                POST_VARN(3,   1,  86)   /* DTCOND */
                POST_VARN(2,   2,  87)   /* DTENDTH and DTENDTQ */
                POST_VARN(3,   2,  89)   /* EXTINCT and FICE */
                POST_VARN(2,   7,  91)   /* FLDS ... FLUTC */
                POST_VARN(3,   4,  98)   /* FREQI ... FREQS */
                POST_VARN(2,  15, 102)   /* FSDS ... ICEFRAC */
                POST_VARN(3,   3, 117)   /* ICIMR ... IWC */
                POST_VARN(2,   2, 120)   /* LANDFRAC and LHFLX */
                POST_VARN(3,   5, 122)   /* LINOZ_DO3 ... LINOZ_SSO3 */
                POST_VARN(2,   3, 127)   /* LINOZ_SZA ... LWCF */
                POST_VARN(3,  12, 130)   /* Mass_bc ... O3 */
                POST_VARN(2,   2, 142)   /* O3_SRF and OCNFRAC */
                POST_VARN(3,   1, 144)   /* OMEGA */
                POST_VARN(2,   1, 145)   /* OMEGA500 */
                POST_VARN(3,   1, 146)   /* OMEGAT */
                POST_VARN(2,   8, 147)   /* PBLH ... PSL */
                POST_VARN(3,   1, 155)   /* Q */
                POST_VARN(2,   2, 156)   /* QFLX and QREFHT */
                POST_VARN(3,   3, 158)   /* QRL ... RAINQM */
                POST_VARN(2,   1, 161)   /* RAM1 */
                POST_VARN(3,   1, 162)   /* RELHUM */
                POST_VARN(2,  37, 163)   /* SFDMS ... SNOWHLND */
                POST_VARN(3,   2, 200)   /* SNOWQM and SO2 */
                POST_VARN(2,  10, 202)   /* SO2_CLXF ... SWCF */
                POST_VARN(3,   1, 212)   /* T */
                POST_VARN(2,  19, 213)   /* TAUGWX ... TVQ */
                POST_VARN(3,   1, 232)   /* U */
                POST_VARN(2,   1, 233)   /* U10 */
                POST_VARN(3,   6, 234)   /* UU ... VV */
                POST_VARN(2,   3, 240)   /* WD_H2O2 ... WD_SO2 */
                POST_VARN(3,   3, 243)   /* WSUB ... aero_water */
                POST_VARN(2,  32, 246)   /* airFV ... dst_c3SFWET */
                POST_VARN(3,   1, 278)   /* hstobie_linoz */
                POST_VARN(2, 129, 279)   /* mlip ... soa_c3SFWET */
            }
        }
        else {
            if (two_buf) {
                /* write 2D variables followed by 3D variables */
                POST_VARN(2, 13, 30)   /* CLDHGH ... T5 */
                POST_VARN(2,  7, 44)   /* U250 ... Z500 */
                POST_VARN(3,  1, 43)   /* U */
            } else {
                /* write variables in the same order as they defined */
                POST_VARN(2, 13, 30)   /* CLDHGH ... T5 */
                POST_VARN(3,  1, 43)   /* U */
                POST_VARN(2,  7, 44)   /* U250 ... Z500 */
            }
        }

        post_timing += MPI_Wtime() - timing;

        MPI_Barrier(comm); /*-----------------------------------------*/
        timing = MPI_Wtime();

        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); ERR

        /* high water mark of number of noncontiguous requests */
        if (nreqs > max_nreqs) max_nreqs = nreqs;
        nreqs = 0;

        wait_timing += MPI_Wtime() - timing;
    }

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    err = ncmpi_inq_put_size(ncid, &total_size); ERR
    put_size = total_size - metadata_size;
    err = ncmpi_close(ncid); ERR
    close_timing += MPI_Wtime() - timing;

    free(starts_D3[0]); free(starts_D3);
    free(starts_D2[0]); free(starts_D2);
    free(fix_starts_D2[0]); free(fix_starts_D2);
    free(fix_starts_D1[0]); free(fix_starts_D1);
    if (rec_bufp == NULL){
        free(rec_buf);
    }
    if (dbl_bufp == NULL){
        free(dbl_buf);
    }
    free(varids);

    total_timing = MPI_Wtime() - total_timing;

    nreqs = max_nreqs;
    MPI_Reduce(&nreqs,         &max_nreqs,  1, MPI_INT,    MPI_MAX, 0, comm);
    MPI_Reduce(&put_size,      &tmp,        1, MPI_OFFSET, MPI_SUM, 0, comm);
    put_size = tmp;
    MPI_Reduce(&total_size,    &tmp,        1, MPI_OFFSET, MPI_SUM, 0, comm);
    total_size = tmp;
    MPI_Reduce(&open_timing,   &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    open_timing = max_timing;
    MPI_Reduce(&pre_timing,    &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    pre_timing = max_timing;
    MPI_Reduce(&post_timing,   &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    post_timing = max_timing;
    MPI_Reduce(&wait_timing,   &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
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
        printf("History output file                = %s\n", outfile);
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
        printf("I/O bandwidth (open-to-close)      = %.4f MiB/sec\n",
               (double)total_size/1048576.0/total_timing);
        printf("I/O bandwidth (write-only)         = %.4f MiB/sec\n",
               (double)put_size/1048576.0/wait_timing);
        if (verbose) print_info(&info_used);
        printf("-----------------------------------------------------------\n");
    }
fn_exit:
    if (info_used != MPI_INFO_NULL) MPI_Info_free(&info_used);
    if (!keep_outfile) unlink(outfname);
    fflush(stdout);
    MPI_Barrier(comm);
    return nerrs;
}

/*----< run_varn() >--------------------------------------------------------*/
static
int run_varn_rd(char       *out_dir,      /* output folder name */
             char       *outfile,      /* output file name */
             int         nvars,        /* number of variables 408 or 51 */
             int         num_recs,     /* number of records */
             int         noncontig_buf,/* whether to us noncontiguous buffer */
             MPI_Info    info,
             MPI_Offset *dims,         /* [2] dimension lengths */
             int         nreqs_D1,     /* no. request in decomposition 1 */
             int        *disps_D1,     /* [nreqs_D1] request's displacements */
             int        *blocklens_D1, /* [nreqs_D1] request's block lengths */
             int         nreqs_D2,     /* no. request in decomposition 2 */
             int        *disps_D2,     /* [nreqs_D2] request's displacements */
             int        *blocklens_D2, /* [nreqs_D2] request's block lengths */
             int         nreqs_D3,     /* no. request in decomposition 3 */
             int        *disps_D3,     /* [nreqs_D3] request's displacements */
             int        *blocklens_D3, /* [nreqs_D3] request's block lengths */
             double     **dbl_bufp,
             dtype      **rec_bufp,
             char       *txt_buf,
             int        *int_buf)
{
    char outfname[512], *txt_buf_ptr;
    int i, j, k, err, nerrs=0, rank, ndims, ncid, cmode, *varids;
    int rec_no, gap=0, nreqs, max_nreqs, *int_buf_ptr;
    size_t dbl_buflen, rec_buflen;
    size_t nelems_D1, nelems_D2, nelems_D3;
    dtype *rec_buf, *rec_buf_ptr;
    double *dbl_buf, *dbl_buf_ptr;
    double pre_timing, open_timing, post_timing, wait_timing, close_timing;
    double timing, total_timing,  max_timing;
    MPI_Offset tmp, metadata_size, put_size, total_size;
    MPI_Offset **fix_starts_D1, **fix_counts_D1;
    MPI_Offset **fix_starts_D2, **fix_counts_D2;
    MPI_Offset **starts_D2, **counts_D2;
    MPI_Offset **starts_D3, **counts_D3;
    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Info info_used=MPI_INFO_NULL;

    MPI_Barrier(comm); /*-----------------------------------------*/
    total_timing = pre_timing = MPI_Wtime();

    open_timing = 0.0;
    post_timing = 0.0;
    wait_timing = 0.0;
    close_timing = 0.0;

    MPI_Comm_rank(comm, &rank);

    if (noncontig_buf) gap = 10;

    /* number of variable elements from 3 decompositions */
    nreqs = max_nreqs = 0;
    nelems_D1 = nelems_D2 = nelems_D3 = 0;
    for (k=0; k<nreqs_D1; k++) nelems_D1 += blocklens_D1[k];
    for (k=0; k<nreqs_D2; k++) nelems_D2 += blocklens_D2[k];
    for (k=0; k<nreqs_D3; k++) nelems_D3 += blocklens_D3[k];

    if (verbose && rank == 0)
        printf("nelems_D1=%zd nelems_D2=%zd nelems_D3=%zd\n",
               nelems_D1,nelems_D2,nelems_D3);

    /* construct varn API arguments starts[][] and counts[][] */
    ndims = 1;
    FIX_STARTS_COUNTS(fix_starts_D1, fix_counts_D1, nreqs_D1, disps_D1, blocklens_D1)
    FIX_STARTS_COUNTS(fix_starts_D2, fix_counts_D2, nreqs_D2, disps_D2, blocklens_D2)

    ndims = 2;
    REC_STARTS_COUNTS(0, starts_D2, counts_D2, nreqs_D2, disps_D2, blocklens_D2)
    ndims = 3;
    REC_STARTS_COUNTS(0, starts_D3, counts_D3, nreqs_D3, disps_D3, blocklens_D3)

    /* allocate and initialize write buffer for small variables */
    dbl_buflen = nelems_D2 * 2 + nelems_D1
               + 3 * dims[0] + 3 * (dims[0]+1) + 8 + 2
               + 20 * gap;

    dbl_buf = (double*) malloc(dbl_buflen * sizeof(double));
    *dbl_bufp = dbl_buf;

    /* allocate and initialize write buffer for large variables */
    if (nvars == 408)
        rec_buflen = nelems_D2 * 315 + nelems_D3 * 63 + (315+63) * gap;
    else
        rec_buflen = nelems_D2 * 20 + nelems_D3 + (20+1) * gap;

    rec_buf = (dtype*) malloc(rec_buflen * sizeof(dtype));
    *rec_bufp = rec_buf;

    varids = (int*) malloc(nvars * sizeof(int));

    pre_timing = MPI_Wtime() - pre_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    /* set output file name */
    sprintf(outfname, "%s/%s",out_dir, outfile);

    /* create a new CDF-5 file for writing */
    cmode = NC_64BIT_DATA;
    err = ncmpi_open(comm, outfname, cmode, info, &ncid); ERR

    /* define dimensions, variables, and attributes */
    if (nvars == 408) {
        /* for h0 file */
        err = e3sm_io_header_rd(ncid, dims, nvars, varids); ERR
    }
    else {
        /* for h1 file */
        err = e3sm_io_header1_rd(ncid, dims, nvars, varids); ERR
    }

    /* exit define mode and enter data mode */
    //err = ncmpi_enddef(ncid); ERR

    /* I/O amount so far */
    err = ncmpi_inq_get_size(ncid, &metadata_size); ERR
    err = ncmpi_inq_file_info(ncid, &info_used); ERR
    open_timing += MPI_Wtime() - timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    i = 0;
    dbl_buf_ptr = dbl_buf;

    /* lat */
    err = ncmpi_iget_varn(ncid, varids[i++], nreqs_D2, fix_starts_D2, fix_counts_D2,
                          dbl_buf_ptr, nelems_D2, MPI_DOUBLE, NULL); ERR
    dbl_buf_ptr += nelems_D2 + gap;
    nreqs += nreqs_D2;

    /* lon */
    err = ncmpi_iget_varn(ncid, varids[i++], nreqs_D2, fix_starts_D2, fix_counts_D2,
                          dbl_buf_ptr, nelems_D2, MPI_DOUBLE, NULL); ERR
    dbl_buf_ptr += nelems_D2 + gap;
    nreqs += nreqs_D2;

    /* area */
    err = ncmpi_iget_varn(ncid, varids[i++], nreqs_D1, fix_starts_D1, fix_counts_D1,
                          dbl_buf_ptr, nelems_D1, MPI_DOUBLE, NULL); ERR
    dbl_buf_ptr += nelems_D1 + gap;
    nreqs += nreqs_D1;

    post_timing += MPI_Wtime() - timing;

    for (rec_no=0; rec_no<num_recs; rec_no++) {
        MPI_Barrier(comm); /*-----------------------------------------*/
        timing = MPI_Wtime();

        i=3;
        dbl_buf_ptr = dbl_buf + nelems_D2*2 + nelems_D1 + gap*3;
        int_buf_ptr = int_buf;
        txt_buf_ptr = txt_buf;

        /* next 27 small variables are written by rank 0 only */
        if (rank == 0) {
            nreqs += 27;
            /* post nonblocking requests using ncmpi_iget_varn() */
            err = read_small_vars(ncid, i, varids, rec_no, gap, dims[0],
                                   dims[0]+1, 2, 8, &int_buf_ptr, &txt_buf_ptr,
                                   &dbl_buf_ptr); ERR
        }
        i += 27;

        post_timing += MPI_Wtime() - timing;

        MPI_Barrier(comm); /*-----------------------------------------*/
        timing = MPI_Wtime();

        /* flush fixed-size and small variables */
        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); ERR

        wait_timing += MPI_Wtime() - timing;

        MPI_Barrier(comm); /*-----------------------------------------*/
        timing = MPI_Wtime();

        /* high water mark of number of noncontiguous requests */
        if (nreqs > max_nreqs) max_nreqs = nreqs;
        nreqs = 0;

        rec_buf_ptr = rec_buf;

        for (j=0; j<nreqs_D2; j++) starts_D2[j][0] = rec_no;
        for (j=0; j<nreqs_D3; j++) starts_D3[j][0] = rec_no;

        if (nvars == 408) {
            if (two_buf) {
                /* write 2D variables */
                POST_VARN_RD(2,   1,  30)   /* AEROD_v */
                POST_VARN_RD(2,  18,  33)   /* AODABS ... ANSNOW */
                POST_VARN_RD(2,   6,  53)   /* AQ_DMS ... AQ_SOAG */
                POST_VARN_RD(2,   2,  64)   /* CDNUMC and CLDHGH */
                POST_VARN_RD(2,   3,  68)   /* CLDLOW ... CLDTOT */
                POST_VARN_RD(2,  11,  75)   /* DF_DMS ... DSTSFMBL */
                POST_VARN_RD(2,   2,  87)   /* DTENDTH and DTENDTQ */
                POST_VARN_RD(2,   7,  91)   /* FLDS ... FLUTC */
                POST_VARN_RD(2,  15, 102)   /* FSDS ... ICEFRAC */
                POST_VARN_RD(2,   2, 120)   /* LANDFRAC and LHFLX */
                POST_VARN_RD(2,   3, 127)   /* LINOZ_SZA ... LWCF */
                POST_VARN_RD(2,   2, 142)   /* O3_SRF and OCNFRAC */
                POST_VARN_RD(2,   1, 145)   /* OMEGA500 */
                POST_VARN_RD(2,   8, 147)   /* PBLH ... PSL */
                POST_VARN_RD(2,   2, 156)   /* QFLX and QREFHT */
                POST_VARN_RD(2,   1, 161)   /* RAM1 */
                POST_VARN_RD(2,  37, 163)   /* SFDMS ... SNOWHLND */
                POST_VARN_RD(2,  10, 202)   /* SO2_CLXF ... SWCF */
                POST_VARN_RD(2,  19, 213)   /* TAUGWX ... TVQ */
                POST_VARN_RD(2,   1, 233)   /* U10 */
                POST_VARN_RD(2,   3, 240)   /* WD_H2O2 ... WD_SO2 */
                POST_VARN_RD(2,  32, 246)   /* airFV ... dst_c3SFWET */
                POST_VARN_RD(2, 129, 279)   /* mlip ... soa_c3SFWET */
                /* write 3D variables */
                POST_VARN_RD(3,   2,  31)   /* ANRAIN and ANSNOW */
                POST_VARN_RD(3,   2,  51)   /* AQRAIN and AQSNOW */
                POST_VARN_RD(3,   5,  59)   /* AREI ... CCN3 */
                POST_VARN_RD(3,   2,  66)   /* CLDICE and CLDLIQ */
                POST_VARN_RD(3,   4,  71)   /* CLOUD ... DCQ */
                POST_VARN_RD(3,   1,  86)   /* DTCOND */
                POST_VARN_RD(3,   2,  89)   /* EXTINCT and FICE */
                POST_VARN_RD(3,   4,  98)   /* FREQI ... FREQS */
                POST_VARN_RD(3,   3, 117)   /* ICIMR ... IWC */
                POST_VARN_RD(3,   5, 122)   /* LINOZ_DO3 ... LINOZ_SSO3 */
                POST_VARN_RD(3,  12, 130)   /* Mass_bc ... O3 */
                POST_VARN_RD(3,   1, 144)   /* OMEGA */
                POST_VARN_RD(3,   1, 146)   /* OMEGAT */
                POST_VARN_RD(3,   1, 155)   /* Q */
                POST_VARN_RD(3,   3, 158)   /* QRL ... RAINQM */
                POST_VARN_RD(3,   1, 162)   /* RELHUM */
                POST_VARN_RD(3,   2, 200)   /* SNOWQM and SO2 */
                POST_VARN_RD(3,   1, 212)   /* T */
                POST_VARN_RD(3,   1, 232)   /* U */
                POST_VARN_RD(3,   6, 234)   /* UU ... VV */
                POST_VARN_RD(3,   3, 243)   /* WSUB ... aero_water */
                POST_VARN_RD(3,   1, 278)   /* hstobie_linoz */
            } else {
                /* write variables in the same order as they defined */
                POST_VARN_RD(2,   1,  30)   /* AEROD_v */
                POST_VARN_RD(3,   2,  31)   /* ANRAIN and ANSNOW */
                POST_VARN_RD(2,  18,  33)   /* AODABS ... ANSNOW */
                POST_VARN_RD(3,   2,  51)   /* AQRAIN and AQSNOW */
                POST_VARN_RD(2,   6,  53)   /* AQ_DMS ... AQ_SOAG */
                POST_VARN_RD(3,   5,  59)   /* AREI ... CCN3 */
                POST_VARN_RD(2,   2,  64)   /* CDNUMC and CLDHGH */
                POST_VARN_RD(3,   2,  66)   /* CLDICE and CLDLIQ */
                POST_VARN_RD(2,   3,  68)   /* CLDLOW ... CLDTOT */
                POST_VARN_RD(3,   4,  71)   /* CLOUD ... DCQ */
                POST_VARN_RD(2,  11,  75)   /* DF_DMS ... DSTSFMBL */
                POST_VARN_RD(3,   1,  86)   /* DTCOND */
                POST_VARN_RD(2,   2,  87)   /* DTENDTH and DTENDTQ */
                POST_VARN_RD(3,   2,  89)   /* EXTINCT and FICE */
                POST_VARN_RD(2,   7,  91)   /* FLDS ... FLUTC */
                POST_VARN_RD(3,   4,  98)   /* FREQI ... FREQS */
                POST_VARN_RD(2,  15, 102)   /* FSDS ... ICEFRAC */
                POST_VARN_RD(3,   3, 117)   /* ICIMR ... IWC */
                POST_VARN_RD(2,   2, 120)   /* LANDFRAC and LHFLX */
                POST_VARN_RD(3,   5, 122)   /* LINOZ_DO3 ... LINOZ_SSO3 */
                POST_VARN_RD(2,   3, 127)   /* LINOZ_SZA ... LWCF */
                POST_VARN_RD(3,  12, 130)   /* Mass_bc ... O3 */
                POST_VARN_RD(2,   2, 142)   /* O3_SRF and OCNFRAC */
                POST_VARN_RD(3,   1, 144)   /* OMEGA */
                POST_VARN_RD(2,   1, 145)   /* OMEGA500 */
                POST_VARN_RD(3,   1, 146)   /* OMEGAT */
                POST_VARN_RD(2,   8, 147)   /* PBLH ... PSL */
                POST_VARN_RD(3,   1, 155)   /* Q */
                POST_VARN_RD(2,   2, 156)   /* QFLX and QREFHT */
                POST_VARN_RD(3,   3, 158)   /* QRL ... RAINQM */
                POST_VARN_RD(2,   1, 161)   /* RAM1 */
                POST_VARN_RD(3,   1, 162)   /* RELHUM */
                POST_VARN_RD(2,  37, 163)   /* SFDMS ... SNOWHLND */
                POST_VARN_RD(3,   2, 200)   /* SNOWQM and SO2 */
                POST_VARN_RD(2,  10, 202)   /* SO2_CLXF ... SWCF */
                POST_VARN_RD(3,   1, 212)   /* T */
                POST_VARN_RD(2,  19, 213)   /* TAUGWX ... TVQ */
                POST_VARN_RD(3,   1, 232)   /* U */
                POST_VARN_RD(2,   1, 233)   /* U10 */
                POST_VARN_RD(3,   6, 234)   /* UU ... VV */
                POST_VARN_RD(2,   3, 240)   /* WD_H2O2 ... WD_SO2 */
                POST_VARN_RD(3,   3, 243)   /* WSUB ... aero_water */
                POST_VARN_RD(2,  32, 246)   /* airFV ... dst_c3SFWET */
                POST_VARN_RD(3,   1, 278)   /* hstobie_linoz */
                POST_VARN_RD(2, 129, 279)   /* mlip ... soa_c3SFWET */
            }
        }
        else {
            if (two_buf) {
                /* write 2D variables followed by 3D variables */
                POST_VARN_RD(2, 13, 30)   /* CLDHGH ... T5 */
                POST_VARN_RD(2,  7, 44)   /* U250 ... Z500 */
                POST_VARN_RD(3,  1, 43)   /* U */
            } else {
                /* write variables in the same order as they defined */
                POST_VARN_RD(2, 13, 30)   /* CLDHGH ... T5 */
                POST_VARN_RD(3,  1, 43)   /* U */
                POST_VARN_RD(2,  7, 44)   /* U250 ... Z500 */
            }
        }

        post_timing += MPI_Wtime() - timing;

        MPI_Barrier(comm); /*-----------------------------------------*/
        timing = MPI_Wtime();

        err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); ERR

        /* high water mark of number of noncontiguous requests */
        if (nreqs > max_nreqs) max_nreqs = nreqs;
        nreqs = 0;

        wait_timing += MPI_Wtime() - timing;
    }

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    err = ncmpi_inq_get_size(ncid, &total_size); ERR
    put_size = total_size - metadata_size;
    err = ncmpi_close(ncid); ERR
    close_timing += MPI_Wtime() - timing;

    free(starts_D3[0]); free(starts_D3);
    free(starts_D2[0]); free(starts_D2);
    free(fix_starts_D2[0]); free(fix_starts_D2);
    free(fix_starts_D1[0]); free(fix_starts_D1);
    free(varids);

    total_timing = MPI_Wtime() - total_timing;

    nreqs = max_nreqs;
    MPI_Reduce(&nreqs,         &max_nreqs,  1, MPI_INT,    MPI_MAX, 0, comm);
    MPI_Reduce(&put_size,      &tmp,        1, MPI_OFFSET, MPI_SUM, 0, comm);
    put_size = tmp;
    MPI_Reduce(&total_size,    &tmp,        1, MPI_OFFSET, MPI_SUM, 0, comm);
    total_size = tmp;
    MPI_Reduce(&open_timing,   &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    open_timing = max_timing;
    MPI_Reduce(&pre_timing,    &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    pre_timing = max_timing;
    MPI_Reduce(&post_timing,   &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    post_timing = max_timing;
    MPI_Reduce(&wait_timing,   &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
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
        printf("History input file                = %s\n", outfile);
        printf("MAX heap memory allocated by PnetCDF internally is %.2f MiB\n",
               (float)max_alloc/1048576);
        printf("Total number of variables          = %d\n",nvars);
        printf("Total read amount                 = %.2f MiB = %.2f GiB\n",
               (double)total_size/1048576,(double)total_size/1073741824);
        printf("Max number of requests             = %d\n",max_nreqs);
        printf("Max Time of open + metadata define = %.4f sec\n",open_timing);
        printf("Max Time of I/O preparing          = %.4f sec\n",pre_timing);
        printf("Max Time of ncmpi_iget_varn        = %.4f sec\n",post_timing);
        printf("Max Time of ncmpi_wait_all         = %.4f sec\n",wait_timing);
        printf("Max Time of close                  = %.4f sec\n",close_timing);
        printf("Max Time of TOTAL                  = %.4f sec\n",total_timing);
        printf("I/O bandwidth (open-to-close)      = %.4f MiB/sec\n",
               (double)total_size/1048576.0/total_timing);
        printf("I/O bandwidth (read-only)         = %.4f MiB/sec\n",
               (double)put_size/1048576.0/wait_timing);
        if (verbose) print_info(&info_used);
        printf("-----------------------------------------------------------\n");
    }
fn_exit:
    if (info_used != MPI_INFO_NULL) MPI_Info_free(&info_used);
    if (!keep_outfile) unlink(outfname);
    fflush(stdout);
    MPI_Barrier(comm);
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

/*----< usage() >------------------------------------------------------------*/
static void
usage(char *argv0)
{
    char *help =
    "Usage: %s [OPTION]... FILE\n"
    "       [-h] Print help\n"
    "       [-v] Verbose mode\n"
    "       [-k] Keep the output files when program exits\n"
    "       [-d] Run test that uses PnetCDF vard API\n"
    "       [-n] Run test that uses PnetCDF varn API\n"
    "       [-m] Run test using noncontiguous write buffer\n"
    "       [-t] Write 2D variables followed by 3D variables\n"
    "       [-r num] Number of records (default 1)\n"
    "       [-o output_dir] Output directory name (default ./)\n"
    "       FILE: Name of input netCDF file describing data decompositions\n";
    fprintf(stderr, help, argv0);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char** argv)
{
    extern int optind;
    char *infname, out_dir[1024], in_dir[1024], *outfname;
    int i, rank, nprocs, err, nerrs=0, tst_vard=0, tst_varn=0, tst_wr=0, tst_rd=0, noncontig_buf=0;
    int nvars, num_recs;
    int contig_nreqs_D1, *disps_D1=NULL, *blocklens_D1=NULL;
    int contig_nreqs_D2, *disps_D2=NULL, *blocklens_D2=NULL;
    int contig_nreqs_D3, *disps_D3=NULL, *blocklens_D3=NULL;
    MPI_Offset dims_D1[2], dims_D2[2], dims_D3[2];
    MPI_Info info=MPI_INFO_NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    out_dir[0] = '\0';
    in_dir[0] = '\0';
    verbose = 0;
    keep_outfile = 0;
    num_recs = 1;
    two_buf = 0;

    /* command-line arguments */
    while ((i = getopt(argc, argv, "hkvdnmto:r:WR:")) != EOF)
        switch(i) {
            case 'v': verbose = 1;
                      break;
            case 'k': keep_outfile = 1;
                      break;
            case 'r': num_recs = atoi(optarg);
                      break;
            case 'R': strcpy(in_dir, optarg);  
                      tst_rd = 1;
                      break;
            case 'W': tst_wr = 1;
                      break;
            case 'd': tst_vard = 1;
                      break;
            case 'n': tst_varn = 1;
                      break;
            case 'm': noncontig_buf = 1;
                      break;
            case 't': two_buf = 1;
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
    if (tst_vard == 0 && tst_varn == 0)
        /* neither command-line option -d or -n is used, run both */
        tst_vard = tst_varn = 1;

    if (tst_wr == 0 && tst_rd == 0)
        /* neither command-line option -R or -W is used, run write */
        tst_wr = 1;

    /* input file contains number of write requests and their file access
     * offsets (per array element) */
    infname = argv[optind];
    if (verbose && rank==0) printf("input file name =%s\n",infname);

    /* set the output folder name */
    if (out_dir[0] == '\0') {
        strcpy(out_dir, ".");
    }
    if (verbose && rank==0) printf("output folder name =%s\n",out_dir);

    /* set the input folder name */
    if (tst_rd){
        if (in_dir[0] == '\0') {
            strcpy(in_dir, ".");
        }
        if (verbose && rank==0) printf("input folder name =%s\n",in_dir);
    }

    /* read I/O decomposition from input file */
    err = read_io_decomp(infname, "D1", 1, dims_D1, &contig_nreqs_D1,
                         &disps_D1, &blocklens_D1);
    if (err) goto fn_exit;
    err = read_io_decomp(infname, "D2", 1, dims_D2, &contig_nreqs_D2,
                         &disps_D2, &blocklens_D2);
    if (err) goto fn_exit;
    err = read_io_decomp(infname, "D3", 2, dims_D3, &contig_nreqs_D3,
                         &disps_D3, &blocklens_D3);
    if (err) goto fn_exit;

    if (verbose && rank==0) {
        printf("number of noncontiguous requests for D1 = %d\n",contig_nreqs_D1);
        printf("number of noncontiguous requests for D2 = %d\n",contig_nreqs_D2);
        printf("number of noncontiguous requests for D3 = %d\n",contig_nreqs_D3);
    }

    /* set MPI-IO hints */
    MPI_Info_create(&info);
    MPI_Info_set(info, "romio_ds_write", "disable"); /* MPI-IO data sieving */
    MPI_Info_set(info, "romio_cb_write", "enable");  /* collective write */
    MPI_Info_set(info, "romio_no_indep_rw", "true"); /* no independent MPI-IO */

    /* set PnetCDF I/O hints */
    MPI_Info_set(info, "nc_var_align_size", "1"); /* no gap between variables */
    MPI_Info_set(info, "nc_in_place_swap", "enable"); /* in-place byte swap */

    if (!rank) {
        printf("Total number of MPI processes      = %d\n",nprocs);
        printf("Input decomposition file           = %s\n",infname);
        printf("Output file directory              = %s\n",out_dir);
        printf("Variable dimensions (C order)      = %lld x %lld\n",dims_D3[0],dims_D3[1]);
        printf("Write number of records (time dim) = %d\n",num_recs);
        printf("Using noncontiguous write buffer   = %s\n",noncontig_buf?"yes":"no");
    }

    /* vard APIs require internal data type matches external one */
    if (tst_vard) {
#if REC_XTYPE != NC_FLOAT
        if (!rank)
            printf("PnetCDF vard API requires internal and external data types match, skip\n");
#else
        if (!rank) {
            printf("\n==== benchmarking vard API ================================\n");
            printf("Variable written order: same as variables are defined\n\n");
        }
        fflush(stdout);
        if (tst_rd){
            if (!rank)
            printf("Reading not supported for vard\n");
        }
        if (tst_wr){
            MPI_Barrier(MPI_COMM_WORLD);

            nvars = 408;
            outfname = "testfile_h0_vard.nc";
            nerrs += run_vard(out_dir, outfname, nvars, num_recs,
                            noncontig_buf, info, dims_D3,
                            contig_nreqs_D1, disps_D1, blocklens_D1,
                            contig_nreqs_D2, disps_D2, blocklens_D2,
                            contig_nreqs_D3, disps_D3, blocklens_D3);

            MPI_Barrier(MPI_COMM_WORLD);

            nvars = 51;
            outfname = "testfile_h1_vard.nc";
            nerrs += run_vard(out_dir, outfname, nvars, num_recs,
                            noncontig_buf, info, dims_D3,
                            contig_nreqs_D1, disps_D1, blocklens_D1,
                            contig_nreqs_D2, disps_D2, blocklens_D2,
                            contig_nreqs_D3, disps_D3, blocklens_D3);
        }
#endif
    }

    if (tst_varn) {
        double *dbl_buf_h0 = NULL, *dbl_buf_h1 = NULL;
        dtype *rec_buf_h0 = NULL, *rec_buf_h1 = NULL;
        char txt_buf[2][16];
        int int_buf[2][10];

        if (tst_rd){        
            if (!rank) {
                printf("\n==== benchmarking varn API read ================================\n");
                printf("Variable written order: ");
                if (two_buf)
                    printf("2D variables then 3D variables\n\n");
                else
                    printf("same as variables are defined\n\n");
            }
            fflush(stdout);
            MPI_Barrier(MPI_COMM_WORLD);

            /* There are two kinds of outputs for history variables.
            * Output 1st kind history variables.
            */
            nvars = 408;
            outfname = "testfile_h0_varn.nc";
            nerrs += run_varn_rd(in_dir, outfname, nvars, num_recs,
                            noncontig_buf, info, dims_D3,
                            contig_nreqs_D1, disps_D1, blocklens_D1,
                            contig_nreqs_D2, disps_D2, blocklens_D2,
                            contig_nreqs_D3, disps_D3, blocklens_D3, &dbl_buf_h0, &rec_buf_h0, txt_buf[0], int_buf[0]);

            MPI_Barrier(MPI_COMM_WORLD);

            /* Output 2nd kind history variables. */
            nvars = 51;
            outfname = "testfile_h1_varn.nc";
            nerrs += run_varn_rd(in_dir, outfname, nvars, num_recs,
                            noncontig_buf, info, dims_D3,
                            contig_nreqs_D1, disps_D1, blocklens_D1,
                            contig_nreqs_D2, disps_D2, blocklens_D2,
                            contig_nreqs_D3, disps_D3, blocklens_D3, &dbl_buf_h1, &rec_buf_h1, txt_buf[1], int_buf[1]);
        }
        if (tst_wr){
            if (!rank) {
                printf("\n==== benchmarking varn API write ================================\n");
                printf("Variable written order: ");
                if (two_buf)
                    printf("2D variables then 3D variables\n\n");
                else
                    printf("same as variables are defined\n\n");
            }
            fflush(stdout);

            /* There are two kinds of outputs for history variables.
            * Output 1st kind history variables.
            */
            nvars = 408;
            outfname = "testfile_h0_varn.nc";
            nerrs += run_varn(out_dir, outfname, nvars, num_recs,
                            noncontig_buf, info, dims_D3,
                            contig_nreqs_D1, disps_D1, blocklens_D1,
                            contig_nreqs_D2, disps_D2, blocklens_D2,
                            contig_nreqs_D3, disps_D3, blocklens_D3, dbl_buf_h0, rec_buf_h0, txt_buf[0], int_buf[0]);

            MPI_Barrier(MPI_COMM_WORLD);

            /* Output 2nd kind history variables. */
            nvars = 51;
            outfname = "testfile_h1_varn.nc";
            nerrs += run_varn(out_dir, outfname, nvars, num_recs,
                            noncontig_buf, info, dims_D3,
                            contig_nreqs_D1, disps_D1, blocklens_D1,
                            contig_nreqs_D2, disps_D2, blocklens_D2,
                            contig_nreqs_D3, disps_D3, blocklens_D3, dbl_buf_h1, rec_buf_h1, txt_buf[1], int_buf[1]);
        }

        if (dbl_buf_h0 != NULL){
            free(dbl_buf_h0);
        }
        if (dbl_buf_h1 != NULL){
            free(dbl_buf_h1);
        }
        if (rec_buf_h0 != NULL){
            free(rec_buf_h0);
        }
        if (rec_buf_h1 != NULL){
            free(rec_buf_h1);
        }
    }

fn_exit:
    if (info != MPI_INFO_NULL) MPI_Info_free(&info);

    if (disps_D1     != NULL) free(disps_D1);
    if (blocklens_D1 != NULL) free(blocklens_D1);
    if (disps_D2     != NULL) free(disps_D2);
    if (blocklens_D2 != NULL) free(blocklens_D2);
    if (disps_D3     != NULL) free(disps_D3);
    if (blocklens_D3 != NULL) free(blocklens_D3);

    MPI_Finalize();
    return (nerrs > 0);
}

