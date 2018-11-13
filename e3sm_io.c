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
 * To compile:
 *    mpicc e3sm_io.c -o e3sm_io -I/path/PnetCDF/include -L/path/PnetCDF/lib -lpnetcdf
 *
 * To run, for example:
 *    mpiexec -n 8 ./e3sm_io -q datasets/866x72_16p.nc
 *
 *    ---- benchmarking vard API -----------------------
 *    -----------------------------------------------------------
 *    MAX heap memory allocated by PnetCDF internally is 2.16 MiB
 *    Total number of variables          = 408
 *    Total write amount                 = 16.13 MiB = 0.02 GiB
 *    Max number of requests             = 325153
 *    Max Time of open + metadata define = 0.0676 sec
 *    Max Time of I/O preparing          = 0.0260 sec
 *    Max Time of ncmpi_put_vard         = 0.6789 sec
 *    Max Time of close                  = 0.0103 sec
 *    Max Time of TOTAL                  = 0.7828 sec
 *    I/O bandwidth                      = 20.6098 MiB/sec
 *
 *    ---- benchmarking varn API -----------------------
 *    -----------------------------------------------------------
 *    MAX heap memory allocated by PnetCDF internally is 36.59 MiB
 *    Total number of variables          = 408
 *    Total write amount                 = 16.13 MiB = 0.02 GiB
 *    Max number of requests             = 325153
 *    Max Time of open + metadata define = 0.0537 sec
 *    Max Time of I/O preparing          = 0.0018 sec
 *    Max Time of ncmpi_iput_varn        = 0.0748 sec
 *    Max Time of ncmpi_wait_all         = 0.7772 sec
 *    Max Time of close                  = 0.0156 sec
 *    Max Time of TOTAL                  = 0.9231 sec
 *    I/O bandwidth                      = 17.4770 MiB/sec
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
    if (verbose) printf("rank %d: nreqs=%d\n",rank,nreqs);

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

#define SET_TYPE(kind) { \
    var_types[i] = type[kind]; \
    err = ncmpi_inq_varoffset(ncid, varids[i], &var_offset); ERR \
    var_disps[i] = var_offset - offset_flt; \
    if (kind == 2) { \
        flt_buflen += nelems_D2 + gap; \
        nreqs += nreqs_D2; \
        if (i < nvars-1) buf_disps[i+1] = buf_disps[i] + (nelems_D2+gap) * sizeof(float); \
        buf_blocklens[i] = nelems_D2; \
    } else { /* kind == 3 */ \
        flt_buflen += nelems_D3 + gap; \
        nreqs += nreqs_D3; \
        if (i < nvars-1) buf_disps[i+1] = buf_disps[i] + (nelems_D3+gap) * sizeof(float); \
        buf_blocklens[i] = nelems_D3; \
    } \
    i++; \
}

#define SET_TYPES(kind, num) \
    for (j=0; j<num; j++) { \
        var_types[i] = type[kind]; \
        err = ncmpi_inq_varoffset(ncid, varids[i], &var_offset); ERR \
        var_disps[i] = var_offset - offset_flt; \
        if (kind == 2) { \
            flt_buflen += nelems_D2 + gap; \
            nreqs += nreqs_D2; \
            if (i < nvars-1) buf_disps[i+1] = buf_disps[i] + (nelems_D2+gap) * sizeof(float); \
            buf_blocklens[i] = nelems_D2; \
        } else { /* kind == 3 */ \
            flt_buflen += nelems_D3 + gap; \
            nreqs += nreqs_D3; \
            if (i < nvars-1) buf_disps[i+1] = buf_disps[i] + (nelems_D3+gap) * sizeof(float); \
            buf_blocklens[i] = nelems_D3; \
        } \
        i++; \
    }

/*----< run_vard() >--------------------------------------------------------*/
static
int run_vard(char       *out_dir,      /* output folder name */
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
    char outfname[512];
    int i, j, k, err, nerrs=0, rank, ncid, cmode, nvars, *varids;
    int *var_blocklens, *buf_blocklens, nreqs, max_nreqs, gap=0;
    size_t dbl_buflen, flt_buflen;
    size_t nelems_D1, nelems_D2, nelems_D3;
    float *flt_buf;
    double *dbl_buf;
    double pre_timing, open_timing, io_timing, close_timing;
    double timing, total_timing,  max_timing;
    MPI_Aint *var_disps, *buf_disps;
    MPI_Offset tmp, metadata_size, put_size, total_size;
    MPI_Offset offset_dbl, offset_flt, var_offset;
    MPI_Datatype *var_types, type[4], filetype_flt, filetype_dbl;
    MPI_Datatype buftype_flt, buftype_dbl;
    MPI_Info info_used=MPI_INFO_NULL;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Barrier(comm); /*-----------------------------------------*/
    total_timing = open_timing = MPI_Wtime();

    MPI_Comm_rank(comm, &rank);

    /* set output file name */
    sprintf(outfname, "%s/testfile_vard.nc",out_dir);

    nvars = 408;
    varids = (int*) malloc(nvars * sizeof(int));

    /* create a new CDF-5 file for writing */
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(comm, outfname, cmode, info, &ncid); ERR

    /* define dimensions, variables, and attributes */
    err = e3sm_io_header(ncid, dims, nvars, varids); ERR

    /* exit define mode and enter data mode */
    err = ncmpi_enddef(ncid); ERR

    err = ncmpi_inq_put_size(ncid, &metadata_size); ERR
    err = ncmpi_inq_file_info(ncid, &info_used); ERR
    open_timing = MPI_Wtime() - open_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    pre_timing = MPI_Wtime();

    var_types = (MPI_Datatype*) malloc(nvars * sizeof(MPI_Datatype));
    var_blocklens = (int*) malloc(nvars * 2 * sizeof(int));
    buf_blocklens = var_blocklens + nvars;
    var_disps = (MPI_Aint*) malloc(nvars * 2 * sizeof(MPI_Aint));
    buf_disps = var_disps + nvars;

    /* define MPI datatypes for 4 kinds */
    MPI_Type_indexed(nreqs_D1, blocklens_D1, disps_D1, MPI_DOUBLE, &type[0]);
    MPI_Type_commit(&type[0]);
    MPI_Type_indexed(nreqs_D2, blocklens_D2, disps_D2, MPI_DOUBLE, &type[1]);
    MPI_Type_commit(&type[1]);
    MPI_Type_indexed(nreqs_D2, blocklens_D2, disps_D2, MPI_FLOAT, &type[2]);
    MPI_Type_commit(&type[2]);
    MPI_Type_indexed(nreqs_D3, blocklens_D3, disps_D3, MPI_FLOAT, &type[3]);
    MPI_Type_commit(&type[3]);

    nreqs = 0;
    for (j=0; j<nvars; j++) var_blocklens[j] = 1;
    nelems_D1 = nelems_D2 = nelems_D3 = 0;
    for (k=0; k<nreqs_D1; k++) nelems_D1 += blocklens_D1[k];
    for (k=0; k<nreqs_D2; k++) nelems_D2 += blocklens_D2[k];
    for (k=0; k<nreqs_D3; k++) nelems_D3 += blocklens_D3[k];

    if (verbose && rank == 0)
        printf("nelems_D1=%zd nelems_D2=%zd nelems_D3=%zd\n",nelems_D1,nelems_D2,nelems_D3);

    /* the first 3 variables are of type NC_DOUBLE -------------------*/
    i = 0;
    dbl_buflen = 0;
    if (noncontig_buf) gap = 10;

    /* lat */
    var_types[i] = type[1];
    err = ncmpi_inq_varoffset(ncid, varids[i], &offset_dbl); ERR
    var_disps[i] = 0;
    buf_disps[0] = 0;
    buf_blocklens[0] = nelems_D2;
    i++;
    dbl_buflen += nelems_D2 + gap;
    nreqs += nreqs_D2;

    /* lon */
    var_types[i] = type[1];
    err = ncmpi_inq_varoffset(ncid, varids[i], &var_offset); ERR
    var_disps[i] = var_offset - offset_dbl;
    buf_disps[i] = buf_disps[i-1] + (nelems_D2 + gap) * sizeof (double);
    buf_blocklens[1] = nelems_D2;
    i++;
    dbl_buflen += nelems_D2 + gap;
    nreqs += nreqs_D2;

    /* area */
    var_types[i] = type[0];
    err = ncmpi_inq_varoffset(ncid, varids[i], &var_offset); ERR
    var_disps[i] = var_offset - offset_dbl;
    buf_disps[i] = buf_disps[i-1] + (nelems_D2 + gap) * sizeof (double);
    buf_blocklens[2] = nelems_D1;
    i++;
    dbl_buflen += nelems_D1 + gap;
    nreqs += nreqs_D1;

    /* concatenate 3 var_types[] into filetype_dbl */
    MPI_Type_create_struct(3, var_blocklens, var_disps, var_types,
                           &filetype_dbl);
    MPI_Type_commit(&filetype_dbl);

    if (noncontig_buf) {
        MPI_Type_create_hindexed(3, buf_blocklens, buf_disps, MPI_DOUBLE,
                                 &buftype_dbl);
        MPI_Type_commit(&buftype_dbl);
    }
    else {
        buftype_dbl = MPI_DOUBLE;
    }

    /* allocate and initialize write buffer */
    dbl_buf = (double*) malloc(dbl_buflen * sizeof(double));
    for (j=0; j<dbl_buflen; j++) dbl_buf[j] = rank;

    /* skip next 27 variables that have no decomposition file */
    i += 27;

    /* the remaining variables are all of type NC_FLOAT --------------*/
    flt_buflen = 0;

/* TODO: currently the implementation below only handles one record.
 * Need to revise for writing more than one.
 */

    /* AEROD_v */
    var_types[i] = type[2];
    err = ncmpi_inq_varoffset(ncid, varids[i], &offset_flt); ERR
    var_disps[i] = 0;
    buf_disps[i] = 0;
    buf_disps[i+1] = nelems_D2 * sizeof(float);
    buf_blocklens[i] = nelems_D2;
    i++;
    flt_buflen += nelems_D2 + gap;
    nreqs += nreqs_D2;

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

    /* concatenate nvars-30 var_types[] into filetype_flt */
    MPI_Type_create_struct(nvars-30, var_blocklens+30, var_disps+30,
                           var_types+30, &filetype_flt);
    MPI_Type_commit(&filetype_flt);

    for (j=0; j<4; j++) MPI_Type_free(&type[j]);
    free(var_types);

    /* allocate and initialize write buffer */
    flt_buf = (float*) malloc(flt_buflen * sizeof(float));
    for (j=0; j<flt_buflen; j++) flt_buf[j] = rank;

    if (noncontig_buf) {
        MPI_Type_create_hindexed(nvars-30, buf_blocklens+30, buf_disps+30,
                                 MPI_FLOAT, &buftype_flt);
        MPI_Type_commit(&buftype_flt);
    }
    else {
        buftype_flt = MPI_FLOAT;
    }
    free(var_disps);
    free(var_blocklens);

    pre_timing = MPI_Wtime() - pre_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    io_timing = MPI_Wtime();

    if (noncontig_buf) dbl_buflen = flt_buflen = 1;

    /* write all NC_DOUBLE variables in one vard call */
    err = ncmpi_put_vard_all(ncid, varids[0], filetype_dbl, dbl_buf,
                             dbl_buflen, buftype_dbl); ERR

    /* write all NC_FLOAT variables in one vard call */
    err = ncmpi_put_vard_all(ncid, varids[30], filetype_flt, flt_buf,
                             flt_buflen, buftype_flt); ERR

    io_timing = MPI_Wtime() - io_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    close_timing = MPI_Wtime();

    MPI_Type_free(&filetype_flt);
    MPI_Type_free(&filetype_dbl);

    if (noncontig_buf) {
        MPI_Type_free(&buftype_flt);
        MPI_Type_free(&buftype_dbl);
    }

    err = ncmpi_inq_put_size(ncid, &total_size); ERR
    put_size = total_size - metadata_size;
    err = ncmpi_close(ncid); ERR

    free(flt_buf);
    free(dbl_buf);
    free(varids);

    timing = MPI_Wtime();
    close_timing = timing - close_timing;
    total_timing = timing - total_timing;

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
        printf("I/O bandwidth (open-to-close)      = %.4f MiB/sec\n",
               (double)total_size/1048576.0/total_timing);
        printf("I/O bandwidth (write-only)         = %.4f MiB/sec\n",
               (double)put_size/1048576.0/io_timing);
        if (verbose) print_info(&info_used);
    }
fn_exit:
    if (info_used != MPI_INFO_NULL) MPI_Info_free(&info_used);
    if (!keep_outfile) unlink(outfname);
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

#define POST_VARN(k, num) \
    for (j=0; j<num; j++) { \
        err = ncmpi_iput_varn(ncid, varids[i++], nreqs_D##k, starts_D##k, \
                              counts_D##k, flt_buf_ptr, -1, MPI_FLOAT, NULL); \
        ERR \
        flt_buf_ptr += nelems_D##k + gap; \
        nreqs += nreqs_D##k; \
    }

/*----< run_varn() >--------------------------------------------------------*/
static
int run_varn(char       *out_dir,      /* output folder name */
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
             int        *blocklens_D3) /* [nreqs_D3] request's block lengths */
{
    char outfname[512];
    int i, j, k, err, nerrs=0, rank, ndims, ncid, cmode, nvars, *varids;
    int gap=0, nreqs, max_nreqs;
    size_t dbl_buflen, flt_buflen;
    size_t nelems_D1, nelems_D2, nelems_D3;
    float *flt_buf, *flt_buf_ptr;
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
    total_timing = open_timing = MPI_Wtime();

    MPI_Comm_rank(comm, &rank);

    /* set output file name */
    sprintf(outfname, "%s/testfile_varn.nc",out_dir);

    nvars = 408;
    varids = (int*) malloc(nvars * sizeof(int));

    /* create a new CDF-5 file for writing */
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(comm, outfname, cmode, info, &ncid); ERR

    /* define dimensions, variables, and attributes */
    err = e3sm_io_header(ncid, dims, nvars, varids); ERR

    /* exit define mode and enter data mode */
    err = ncmpi_enddef(ncid); ERR

    err = ncmpi_inq_put_size(ncid, &metadata_size); ERR
    err = ncmpi_inq_file_info(ncid, &info_used); ERR
    open_timing = MPI_Wtime() - open_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    pre_timing = MPI_Wtime();

    nreqs = 0;
    nelems_D1 = nelems_D2 = nelems_D3 = 0;
    for (k=0; k<nreqs_D1; k++) nelems_D1 += blocklens_D1[k];
    for (k=0; k<nreqs_D2; k++) nelems_D2 += blocklens_D2[k];
    for (k=0; k<nreqs_D3; k++) nelems_D3 += blocklens_D3[k];

    if (verbose && rank == 0)
        printf("nelems_D1=%zd nelems_D2=%zd nelems_D3=%zd\n",nelems_D1,nelems_D2,nelems_D3);

    /* construct varn API arguments starts[][] and counts[][] */
    ndims = 1;
    FIX_STARTS_COUNTS(fix_starts_D1, fix_counts_D1, nreqs_D1, disps_D1, blocklens_D1)
    FIX_STARTS_COUNTS(fix_starts_D2, fix_counts_D2, nreqs_D2, disps_D2, blocklens_D2)

    ndims = 2;
    REC_STARTS_COUNTS(0, starts_D2, counts_D2, nreqs_D2, disps_D2, blocklens_D2)
    ndims = 3;
    REC_STARTS_COUNTS(0, starts_D3, counts_D3, nreqs_D3, disps_D3, blocklens_D3)

    /* allocate and initialize write buffer */
    if (noncontig_buf) gap = 10;
    dbl_buflen = nelems_D2 * 2 + nelems_D1 + 3 * gap;
    dbl_buf = (double*) malloc(dbl_buflen * sizeof(double));
    for (i=0; i<dbl_buflen; i++) dbl_buf[i] = rank;

    flt_buflen = nelems_D2 * 315 + nelems_D3 * 63 + (315+63) * gap;
    flt_buf = (float*) malloc(flt_buflen * sizeof(float));
    for (i=0; i<flt_buflen; i++) flt_buf[i] = rank;

    pre_timing = MPI_Wtime() - pre_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    post_timing = MPI_Wtime();

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
    nreqs += nreqs_D1;

    post_timing = MPI_Wtime() - post_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    wait_timing = MPI_Wtime();

    /* flush fixed-size variables */
    err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); ERR

    wait_timing = MPI_Wtime() - wait_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    /* skip next 27 variables that have no decomposition file */
    i += 27;

    flt_buf_ptr = flt_buf;

    POST_VARN(2, 1)   /* AEROD_v */
    POST_VARN(3, 2)   /* ANRAIN and ANSNOW */
    POST_VARN(2, 18)  /* AODABS ... ANSNOW */
    POST_VARN(3, 2)   /* AQRAIN and AQSNOW */
    POST_VARN(2, 6)   /* AQ_DMS ... AQ_SOAG */
    POST_VARN(3, 5)   /* AREI ... CCN3 */
    POST_VARN(2, 2)   /* CDNUMC and CLDHGH */
    POST_VARN(3, 2)   /* CLDICE and CLDLIQ */
    POST_VARN(2, 3)   /* CLDLOW ... CLDTOT */
    POST_VARN(3, 4)   /* CLOUD ... DCQ */
    POST_VARN(2, 11)  /* DF_DMS ... DSTSFMBL */
    POST_VARN(3, 1)   /* DTCOND */
    POST_VARN(2, 2)   /* DTENDTH and DTENDTQ */
    POST_VARN(3, 2)   /* EXTINCT and FICE */
    POST_VARN(2, 7)   /* FLDS ... FLUTC */
    POST_VARN(3, 4)   /* FREQI ... FREQS */
    POST_VARN(2, 15)  /* FSDS ... ICEFRAC */
    POST_VARN(3, 3)   /* ICIMR ... IWC */
    POST_VARN(2, 2)   /* LANDFRAC and LHFLX */
    POST_VARN(3, 5)   /* LINOZ_DO3 ... LINOZ_SSO3 */
    POST_VARN(2, 3)   /* LINOZ_SZA ... LWCF */
    POST_VARN(3, 12)  /* Mass_bc ... O3 */
    POST_VARN(2, 2)   /* O3_SRF and OCNFRAC */
    POST_VARN(3, 1)   /* OMEGA */
    POST_VARN(2, 1)   /* OMEGA500 */
    POST_VARN(3, 1)   /* OMEGAT */
    POST_VARN(2, 8)   /* PBLH ... PSL */
    POST_VARN(3, 1)   /* Q */
    POST_VARN(2, 2)   /* QFLX and QREFHT */
    POST_VARN(3, 3)   /* QRL ... RAINQM */
    POST_VARN(2, 1)   /* RAM1 */
    POST_VARN(3, 1)   /* RELHUM */
    POST_VARN(2, 37)  /* SFDMS ... SNOWHLND */
    POST_VARN(3, 2)   /* SNOWQM and SO2 */
    POST_VARN(2, 10)  /* SO2_CLXF ... SWCF */
    POST_VARN(3, 1)   /* T */
    POST_VARN(2, 19)  /* TAUGWX ... TVQ */
    POST_VARN(3, 1)   /* U */
    POST_VARN(2, 1)   /* U10 */
    POST_VARN(3, 6)   /* UU ... VV */
    POST_VARN(2, 3)   /* WD_H2O2 ... WD_SO2 */
    POST_VARN(3, 3)   /* WSUB ... aero_water */
    POST_VARN(2, 32)  /* airFV ... dst_c3SFWET */
    POST_VARN(3, 1)   /* hstobie_linoz */
    POST_VARN(2, 129) /* mlip ... soa_c3SFWET */

    post_timing += MPI_Wtime() - timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); ERR
    wait_timing += MPI_Wtime() - timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    close_timing = MPI_Wtime();

    err = ncmpi_inq_put_size(ncid, &total_size); ERR
    put_size = total_size - metadata_size;
    err = ncmpi_close(ncid); ERR

    free(starts_D3[0]); free(starts_D3);
    free(starts_D2[0]); free(starts_D2);
    free(fix_starts_D2[0]); free(fix_starts_D2);
    free(fix_starts_D1[0]); free(fix_starts_D1);
    free(flt_buf);
    free(dbl_buf);
    free(varids);

    timing = MPI_Wtime();
    close_timing = timing - close_timing;
    total_timing = timing - total_timing;

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
        printf("I/O bandwidth (open-to-close)      = %.4f MiB/sec\n",
               (double)total_size/1048576.0/total_timing);
        printf("I/O bandwidth (write-only)         = %.4f MiB/sec\n",
               (double)put_size/1048576.0/wait_timing);
        if (verbose) print_info(&info_used);
    }
fn_exit:
    if (info_used != MPI_INFO_NULL) MPI_Info_free(&info_used);
    if (!keep_outfile) unlink(outfname);
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
    "Usage: %s [OPTION]... [FILE]...\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode\n"
    "       [-k] Keep the output files when program exits\n"
    "       [-d] run test that uses PnetCDF vard API\n"
    "       [-n] run test that uses PnetCDF varn API\n"
    "       [-m] run test using noncontiguous write buffer\n"
    "       [-o output_dir]: output directory name (default ./)\n"
    "       input_file: name of input netCDF file describing data decompositions\n";
    fprintf(stderr, help, argv0);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char** argv)
{
    extern int optind;
    char *infname, out_dir[1024];
    int i, rank, nprocs, err, nerrs=0, tst_vard=0, tst_varn=0, noncontig_buf=0;
    int contig_nreqs_D1, *disps_D1=NULL, *blocklens_D1=NULL;
    int contig_nreqs_D2, *disps_D2=NULL, *blocklens_D2=NULL;
    int contig_nreqs_D3, *disps_D3=NULL, *blocklens_D3=NULL;
    MPI_Offset dims_D1[2], dims_D2[2], dims_D3[2];
    MPI_Info info=MPI_INFO_NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    out_dir[0] = '\0';
    verbose = 1;
    keep_outfile = 0;

    /* command-line arguments */
    while ((i = getopt(argc, argv, "hkqdnmo:")) != EOF)
        switch(i) {
            case 'q': verbose = 0;
                      break;
            case 'k': keep_outfile = 1;
                      break;
            case 'd': tst_vard = 1;
                      break;
            case 'n': tst_varn = 1;
                      break;
            case 'm': noncontig_buf = 1;
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

    /* input file contains number of write requests and their file access
     * offsets (per array element) */
    infname = argv[optind];
    if (verbose && rank==0) printf("input file name =%s\n",infname);

    /* set the output folder name */
    if (out_dir[0] == '\0') {
        strcpy(out_dir, ".");
    }
    if (verbose && rank==0) printf("output folder name =%s\n",out_dir);

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

    /* set MPI-IO hints */
    MPI_Info_create(&info);
    MPI_Info_set(info, "romio_ds_write", "disable"); /* MPI-IO data sieving */
    MPI_Info_set(info, "romio_cb_write", "enable");  /* collective write */
    MPI_Info_set(info, "romio_no_indep_rw", "true"); /* no independent MPI-IO */

    /* set PnetCDF I/O hints */
    MPI_Info_set(info, "nc_var_align_size", "1"); /* no gap between variables */
    MPI_Info_set(info, "nc_in_place_swap", "enable"); /* in-place byte swap */
    // MPI_Info_set(info, "cb_config_list", "*:*");  /* all aggregators */

    if (!rank) {
        printf("Total number of MPI processes      = %d\n",nprocs);
        printf("Input decomposition file           = %s\n",infname);
        printf("Output file directory              = %s\n",out_dir);
        printf("Variable dimensions (C order)      = %lld x %lld\n",dims_D3[0],dims_D3[1]);
        printf("Using noncontiguous write buffer   = %s\n",noncontig_buf?"yes":"no");
    }

    if (tst_vard) {
        if (!rank) printf("\n---- benchmarking vard API -----------------------\n");
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
        nerrs += run_vard(out_dir, noncontig_buf, info, dims_D3,
                          contig_nreqs_D1, disps_D1, blocklens_D1,
                          contig_nreqs_D2, disps_D2, blocklens_D2,
                          contig_nreqs_D3, disps_D3, blocklens_D3);
    }

    if (tst_varn) {
        if (!rank) printf("\n---- benchmarking varn API -----------------------\n");
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
        nerrs += run_varn(out_dir, noncontig_buf, info, dims_D3,
                          contig_nreqs_D1, disps_D1, blocklens_D1,
                          contig_nreqs_D2, disps_D2, blocklens_D2,
                          contig_nreqs_D3, disps_D3, blocklens_D3);
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

