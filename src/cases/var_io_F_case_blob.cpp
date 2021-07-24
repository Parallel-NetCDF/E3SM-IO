/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 *********************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h> /* unlink() */

#include <mpi.h>

#include <e3sm_io.h>
#include <e3sm_io_err.h>
#include <e3sm_io_case.hpp>
#include <e3sm_io_driver.hpp>

#ifdef USE_PNETCDF_DIRECTLY
#define CHECK_VAR_ERR(varid) {                                      \
    if (err != NC_NOERR) {                                          \
        char name[64];                                              \
        ncmpi_inq_varname(ncid, varid, name);                       \
        printf("Error in %s:%d: var %s - %s\n", __FILE__, __LINE__, \
               name, ncmpi_strerror(err));                          \
        goto err_out;                                               \
    }                                                               \
}
#define FILE_CREATE(filename) { \
    err = ncmpi_create(cfg.sub_comm, filename, NC_CLOBBER|NC_64BIT_DATA, cfg.info, &ncid); \
    if (err != NC_NOERR) {                                 \
        printf("Error in %s:%d: %s\n", __FILE__, __LINE__, \
               ncmpi_strerror(err));                       \
        goto err_out;                                      \
    }                                                      \
}
#define FILE_CLOSE {                                       \
    err = ncmpi_close(ncid);                               \
    if (err != NC_NOERR) {                                 \
        printf("Error in %s:%d: %s\n", __FILE__, __LINE__, \
               ncmpi_strerror(err));                       \
        goto err_out;                                      \
    }                                                      \
}
#define ENDDEF {                                           \
    err = ncmpi_enddef(ncid);                              \
    if (err != NC_NOERR) {                                 \
        printf("Error in %s:%d: %s\n", __FILE__, __LINE__, \
               ncmpi_strerror(err));                       \
        goto err_out;                                      \
    }                                                      \
}
#define INQ_PUT_SIZE(size) {                               \
    err = ncmpi_inq_put_size(ncid, &size);                 \
    if (err != NC_NOERR) {                                 \
        printf("Error in %s:%d: %s\n", __FILE__, __LINE__, \
               ncmpi_strerror(err));                       \
        goto err_out;                                      \
    }                                                      \
}
#define INQ_FILE_INFO(info) {                              \
    err = ncmpi_inq_file_info(ncid, &info);                \
    if (err != NC_NOERR) {                                 \
        printf("Error in %s:%d: %s\n", __FILE__, __LINE__, \
               ncmpi_strerror(err));                       \
        goto err_out;                                      \
    }                                                      \
}
#define IPUT_VAR_DBL(buf, adv) { \
    err = ncmpi_iput_var_double(ncid, varids[i], buf, NULL); \
    CHECK_VAR_ERR(varids[i]) \
    buf += (adv); \
    my_nreqs++; \
    i++; \
}
#define IPUT_VAR_INT(buf, adv) { \
    err = ncmpi_iput_var_int(ncid, varids[i], buf, NULL); \
    CHECK_VAR_ERR(varids[i]) \
    buf += (adv); \
    my_nreqs++; \
    i++; \
}
#define IPUT_VAR1_DBL(buf, adv) { \
    err = ncmpi_iput_var1_double(ncid, varids[i], start, buf, NULL); \
    CHECK_VAR_ERR(varids[i]) \
    buf += (adv); \
    my_nreqs++; \
    i++; \
}
#define IPUT_VAR1_INT(buf, adv) { \
    err = ncmpi_iput_var1_int(ncid, varids[i], start, buf, NULL); \
    CHECK_VAR_ERR(varids[i]) \
    buf += (adv); \
    my_nreqs++; \
    i++; \
}
#define IPUT_VARA_INT(buf, adv) { \
    err = ncmpi_iput_vara_int(ncid, varids[i], start, count, buf, NULL); \
    CHECK_VAR_ERR(varids[i]) \
    buf += (len); \
    my_nreqs++; \
    i++; \
}
#define IPUT_VARA_INT_NOADV(buf) { \
    err = ncmpi_iput_vara_int(ncid, varids[i], start, count, buf, NULL); \
    CHECK_VAR_ERR(varids[i]) \
    my_nreqs++; \
    i++; \
}
#define IPUT_VARA_INT64_NOADV(buf) { \
    err = ncmpi_iput_vara_longlong(ncid, varids[i], start, count, buf, NULL); \
    CHECK_VAR_ERR(varids[i]) \
    my_nreqs++; \
    i++; \
}
#define IPUT_VARA_DBL(buf, adv) { \
    err = ncmpi_iput_vara_double(ncid, varids[i], start, count, buf, NULL); \
    CHECK_VAR_ERR(varids[i]) \
    buf += (adv); \
    my_nreqs++; \
    i++; \
}
#define IPUT_VARA_TXT(buf, adv) { \
    err = ncmpi_iput_vara_text(ncid, varids[i], start, count, buf, NULL); \
    CHECK_VAR_ERR(varids[i]) \
    buf += (adv); \
    my_nreqs++; \
    i++; \
}
#ifdef _DOUBLE_TYPE_
#define IPUT_VARA(varid, start, count, buf) { \
    err = ncmpi_iput_vara_double(ncid, varid, start, count, buf, NULL); \
    CHECK_VAR_ERR(varid) \
    my_nreqs++; \
}
#else
#define IPUT_VARA(varid, start, count, buf) { \
    err = ncmpi_iput_vara_float(ncid, varid, start, count, buf, NULL); \
    CHECK_VAR_ERR(varid) \
    my_nreqs++; \
}
#endif
#define WAIT_ALL_REQS {                                    \
    err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);    \
    if (err != NC_NOERR) {                                 \
        printf("Error in %s:%d: %s\n", __FILE__, __LINE__, \
               ncmpi_strerror(err));                       \
        goto err_out;                                      \
    }                                                      \
    nflushes++;                                            \
}
#else
#define CHECK_VAR_ERR(varid) {                                \
    if (err != 0) {                                           \
        char var_name[64];                                    \
        driver.inq_var_name(ncid, varid, var_name);           \
        printf("Error in %s:%d: %s() var %s\n"     ,          \
               __FILE__, __LINE__, __func__, var_name);       \
        goto err_out;                                         \
    }                                                         \
}
#define FILE_CREATE(filename) { \
    err = driver.create(filename, cfg.sub_comm, cfg.info, &ncid); \
    CHECK_ERR \
}
#define FILE_CLOSE {            \
    err = driver.close(ncid);   \
    CHECK_ERR                   \
}
#define ENDDEF {                \
    err = driver.enddef(ncid);  \
    CHECK_ERR                   \
}
#define INQ_PUT_SIZE(size) {                 \
    err = driver.inq_put_size(ncid, &size);  \
    CHECK_ERR                                \
}
#define INQ_FILE_INFO(info) {                \
    err = driver.inq_file_info(ncid, &info); \
    CHECK_ERR                                \
}
#define IPUT_VAR_DBL(buf, adv) { \
    err = driver.put_vara(ncid, varids[i], MPI_DOUBLE, NULL, NULL, buf, nb); \
    CHECK_VAR_ERR(varids[i]) \
    buf += (adv); \
    my_nreqs++; \
    i++; \
}
#define IPUT_VAR_INT(buf, adv) { \
    err = driver.put_vara(ncid, varids[i], MPI_INT, NULL, NULL, buf, nb); \
    CHECK_VAR_ERR(varids[i]) \
    buf += (adv); \
    my_nreqs++; \
    i++; \
}
#define IPUT_VAR1_DBL(buf, adv) { \
    err = driver.put_vara(ncid, varids[i], MPI_DOUBLE, start, NULL, buf, nb); \
    CHECK_VAR_ERR(varids[i]) \
    buf += (adv); \
    my_nreqs++; \
    i++; \
}
#define IPUT_VAR1_INT(buf, adv) { \
    err = driver.put_vara(ncid, varids[i], MPI_INT, start, NULL, buf, nb); \
    CHECK_VAR_ERR(varids[i]) \
    buf += (adv); \
    my_nreqs++; \
    i++; \
}
#define IPUT_VARA_INT(buf, adv) { \
    err = driver.put_vara(ncid, varids[i], MPI_INT, start, count, buf, nb); \
    CHECK_VAR_ERR(varids[i]) \
    buf += (adv); \
    my_nreqs++; \
    i++; \
}
#define IPUT_VARA_INT_NOADV(buf) { \
    err = driver.put_vara(ncid, varids[i], MPI_INT, start, count, buf, nb); \
    CHECK_VAR_ERR(varids[i]) \
    my_nreqs++; \
    i++; \
}
#define IPUT_VARA_INT64_NOADV(buf) { \
    err = driver.put_vara(ncid, varids[i], MPI_LONG_LONG, start, count, buf, nb); \
    CHECK_VAR_ERR(varids[i]) \
    my_nreqs++; \
    i++; \
}
#define IPUT_VARA_DBL(buf, adv) { \
    err = driver.put_vara(ncid, varids[i], MPI_DOUBLE, start, count, buf, nb); \
    CHECK_VAR_ERR(varids[i]) \
    buf += (adv); \
    my_nreqs++; \
    i++; \
}
#define IPUT_VARA_TXT(buf, adv) { \
    err = driver.put_vara(ncid, varids[i], MPI_CHAR, start, count, buf, nb); \
    CHECK_VAR_ERR(varids[i]) \
    buf += (adv); \
    my_nreqs++; \
    i++; \
}
#define IPUT_VARA(varid, start, count, buf) { \
    err = driver.put_vara(ncid, varid, REC_ITYPE, start, count, buf, nb); \
    CHECK_VAR_ERR(varid) \
    my_nreqs++; \
}
#define WAIT_ALL_REQS { \
    err = driver.wait(ncid); \
    CHECK_ERR \
    nflushes++; \
}
#endif

/*----< write_small_fix_vars_F_case() >--------------------------------------*/
static
int write_small_fix_vars_F_case(e3sm_io_driver  &driver,
                                int              ncid,
                                int              start_varid,
                                int              gap,
                                MPI_Offset       lev,
                                MPI_Offset       ilev,
                                int            **int_buf,
                                double         **dbl_buf,
                                int             *nreqs)
{
    /* small fixed-size variable IDs */
    int varids[12] = {3, 4, 5, 6, 7, 8, 9,
                      16, 17, 18, 19, 20};

    int i, err, my_nreqs=0;

    for (i=0; i<15; i++) varids[i] += start_varid;

    i = 0;

    /* 7 of type double */
    IPUT_VAR_DBL(*dbl_buf,  lev + gap) /* lev */
    IPUT_VAR_DBL(*dbl_buf,  lev + gap) /* hyam */
    IPUT_VAR_DBL(*dbl_buf,  lev + gap) /* hybm */
    IPUT_VAR_DBL(*dbl_buf,    1 + gap) /* P0 */
    IPUT_VAR_DBL(*dbl_buf, ilev + gap) /* ilev */
    IPUT_VAR_DBL(*dbl_buf, ilev + gap) /* hyai */
    IPUT_VAR_DBL(*dbl_buf, ilev + gap) /* hybi */

    /* 5 of type int */
    IPUT_VAR_INT(*int_buf, 1) /* ndbase */
    IPUT_VAR_INT(*int_buf, 1) /* nsbase */
    IPUT_VAR_INT(*int_buf, 1) /* nbdate */
    IPUT_VAR_INT(*int_buf, 1) /* nbsec */
    IPUT_VAR_INT(*int_buf, 1) /* mdt */

    if (nreqs != NULL) *nreqs = my_nreqs;

err_out:
    return err;
}

/*----< write_small_rec_vars_F_case() >--------------------------------------*/
static
int write_small_rec_vars_F_case(e3sm_io_driver  &driver,
                                int              ncid,
                                int              start_varid,
                                int              rec_no,
                                int              gap,
                                MPI_Offset       nbnd,
                                MPI_Offset       nchars,
                                int            **int_buf,
                                char           **txt_buf,
                                double         **dbl_buf,
                                int             *nreqs)
{
    /* small record variable IDs */
    int varids[15] = {10, 11, 12, 13, 14, 15,
                      21, 22, 23, 24, 25, 26, 27, 28, 29};

    int i, err, my_nreqs=0;
    MPI_Offset start[2], count[2];

    for (i=0; i<15; i++) varids[i] += start_varid;

    start[0] = rec_no;
    start[1] = 0;
    count[0] = 1;
    i = 0;

    IPUT_VAR1_DBL(*dbl_buf,       1 + gap) /* time */
    IPUT_VAR1_INT(*int_buf,       1 + gap) /* date */
    IPUT_VAR1_INT(*int_buf,       1 + gap) /* datesec */
    count[1] = nbnd;
    IPUT_VARA_DBL(*dbl_buf,    nbnd + gap) /* time_bnds */
    count[1] = nchars;
    IPUT_VARA_TXT(*txt_buf, nchars + gap) /* date_written */
    count[1] = nchars;
    IPUT_VARA_TXT(*txt_buf, nchars + gap) /* time_written */

    IPUT_VAR1_INT(*int_buf, 1 + gap) /* ndcur */
    IPUT_VAR1_INT(*int_buf, 1 + gap) /* nscur */
    IPUT_VAR1_DBL(*dbl_buf, 1 + gap) /* co2vmr */
    IPUT_VAR1_DBL(*dbl_buf, 1 + gap) /* ch4vmr */
    IPUT_VAR1_DBL(*dbl_buf, 1 + gap) /* n2ovmr */
    IPUT_VAR1_DBL(*dbl_buf, 1 + gap) /* f11vmr */
    IPUT_VAR1_DBL(*dbl_buf, 1 + gap) /* f12vmr */
    IPUT_VAR1_DBL(*dbl_buf, 1 + gap) /* sol_tsi */
    IPUT_VAR1_INT(*int_buf, 1 + gap) /* nsteph */

    if (nreqs != NULL) *nreqs = my_nreqs;

err_out:
    return err;
}

#define POST_VARA(decomp, numVars, vid) {                                 \
    int varid = vid + num_decomp_vars;                                    \
    for (j=0; j<numVars; j++) {                                           \
        IPUT_VARA(varid+j, start_D##decomp, count_D##decomp, rec_buf_ptr) \
        rec_buf_ptr += count_D##decomp[1] + gap;                          \
        if (rec_no == 0) nvars_D[decomp - 1]++;                           \
    }                                                                     \
}

/*----< blob_F_case() >-------------------------------------------------------*/
int blob_F_case(e3sm_io_config &cfg,
                e3sm_io_decom  &decom,
                e3sm_io_driver &driver)
{
    char outfile[1040], base_name[1024], *ext;
    const char *hist;
    int i, j, err, sub_rank, global_rank, ncid=-1, nflushes=0, *varids=NULL;
    int rec_no, gap = 0, my_nreqs, num_decomp_vars, one_flush;
    int contig_nreqs[MAX_NUM_DECOMP], *nvars_D=cfg.nvars_D;
    double timing;
    MPI_Offset previous_size, metadata_size, total_size;
    MPI_Offset start[2], count[2];
    MPI_Offset start_D2[2], count_D2[2], start_D3[2], count_D3[2];
    MPI_Offset blob_start[MAX_NUM_DECOMP], blob_count[MAX_NUM_DECOMP];
    MPI_Info info_used = MPI_INFO_NULL;

    size_t ii, rec_buflen;
    size_t                 fix_int_buflen, fix_dbl_buflen;
    size_t rec_txt_buflen, rec_int_buflen, rec_dbl_buflen;
    int    *fix_int_buf=NULL, *fix_int_buf_ptr;
    double *fix_dbl_buf=NULL, *fix_dbl_buf_ptr;
    int    *rec_int_buf=NULL, *rec_int_buf_ptr;
    char   *rec_txt_buf=NULL, *rec_txt_buf_ptr;
    double *rec_dbl_buf=NULL, *rec_dbl_buf_ptr;
    itype  *rec_buf=NULL, *rec_buf_ptr;

    MPI_Barrier(cfg.io_comm); /*-----------------------------------------*/
    cfg.end2end_time = cfg.pre_time = MPI_Wtime();

    cfg.post_time  = 0.0;
    cfg.flush_time = 0.0;

    if (cfg.api == hdf5) /* I/O amount from previous I/O */
        previous_size = cfg.amount_WR;
    else
        previous_size = 0;

    MPI_Comm_rank(cfg.io_comm,  &global_rank);
    MPI_Comm_rank(cfg.sub_comm, &sub_rank);

#define FLUSH_ALL_RECORDS_AT_ONCE
#ifdef FLUSH_ALL_RECORDS_AT_ONCE
    one_flush = 1;
#else
    one_flush = 0;
    if (cfg.strategy == blob && cfg.api == hdf5) {
        printf("Error in %s:%d: %s() FLUSH_ALL_RECORDS_AT_ONCE must be enabled when using HDF5 blob I/O",
               __FILE__, __LINE__, __func__);
        return -1;
    }
#endif

    if (cfg.non_contig_buf) gap = 10;

    /* allocate write buffer for fixed-size variables */
    fix_dbl_buflen = decom.count[1] * 2          /* lat[ncol] and lon[ncol] */
                   + decom.count[0]              /* area[ncol] */
                   + 3 * decom.dims[2][0]        /* [lev]: lev, hyam, hybm */
                   + 3 * (decom.dims[2][0] + 1)  /* [ilev]: ilev, hyai, hybi */
                   + 1                           /* P0 */
                   + 5                           /* ndbase ... mdt */
                   + 15 * gap;

    fix_int_buflen = 5         /* ndbase ... mdt */
                   + 5 * gap;

    /* allocate write buffer for record variables */
    rec_txt_buflen = 2 * 8     /* [nchars]: date_written, time_written */
                   + 2 * gap;

    rec_int_buflen = 2         /* date, datesec */
                   + 2         /* ndcur, nscur */
                   + 1         /* nsteph */
                   + 5 * gap;

    rec_dbl_buflen = 1         /* time */
                   + 2         /* [nbnd]: time_bnds */
                   + 6         /* co2vmr ... sol_tsi */
                   + 8 * gap;

    if (cfg.nvars == 414)
        rec_buflen = decom.count[1] * 323
                   + decom.count[2] * 63
                   + (323 + 63) * gap;
    else
        rec_buflen = 13 * decom.count[1]   /* CLDHGH ... TS */
                   +  1 * decom.count[2]   /* U */
                   +  7 * decom.count[1]   /* U250 ... Z500 */
                   + 21 * gap;

    if (one_flush && cfg.api == pnetcdf) {
        /* write buffers should not be touched when using PnetCDF iput before
         * ncmpi_wait_all is called. For HDF5 and ADIOS blob I/O, write data
         * will be copied and cached into internally allocated buffers and user
         * buffers can be reused after put call returned.
         */
        rec_dbl_buflen *= cfg.nrecs;
        rec_txt_buflen *= cfg.nrecs;
        rec_int_buflen *= cfg.nrecs;
        rec_buflen     *= cfg.nrecs;
    }

    /* allocate and initialize write buffer */
    fix_dbl_buf = (double*) malloc(fix_dbl_buflen * sizeof(double));
    fix_int_buf = (int*)    malloc(fix_int_buflen * sizeof(int));
    rec_dbl_buf = (double*) malloc(rec_dbl_buflen * sizeof(double));
    rec_txt_buf = (char*)   malloc(rec_txt_buflen * sizeof(char));
    rec_int_buf = (int*)    malloc(rec_int_buflen * sizeof(int));
    rec_buf     = (itype*)  malloc(rec_buflen     * sizeof(itype));

    for (ii=0; ii<fix_dbl_buflen; ii++) fix_dbl_buf[ii] = global_rank;
    for (ii=0; ii<fix_int_buflen; ii++) fix_int_buf[ii] = global_rank;
    for (ii=0; ii<rec_dbl_buflen; ii++) rec_dbl_buf[ii] = global_rank;
    for (ii=0; ii<rec_txt_buflen; ii++) rec_txt_buf[ii] = 'a' + global_rank;
    for (ii=0; ii<rec_int_buflen; ii++) rec_int_buf[ii] = global_rank;
    for (ii=0; ii<rec_buflen;     ii++) rec_buf[ii]     = global_rank;

    /* there are num_decomp_vars number of decomposition variables */
    num_decomp_vars = decom.num_decomp * NVARS_DECOMP;

    /* allocate space for all variable IDs */
    varids = (int*) malloc((cfg.nvars + num_decomp_vars) * sizeof(int));

    cfg.pre_time = MPI_Wtime() - cfg.pre_time;

    MPI_Barrier(cfg.sub_comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    /* construct h0/h1 file name */
    hist = (cfg.nvars == 414) ? "_h0" :  "_h1";
    ext = strrchr(cfg.out_path, '.');
    if (ext == NULL || (strcmp(ext, ".nc") && strcmp(ext, ".h5"))) {
        sprintf(base_name, "%s%s", cfg.out_path, hist);
        sprintf(outfile, "%s%s.%04d", cfg.out_path, hist, cfg.subfile_ID);
    }
    else {
        sprintf(outfile, "%s", cfg.out_path);
        sprintf(outfile + (ext - cfg.out_path), "%s%s.%04d", hist, ext, cfg.subfile_ID);
        sprintf(base_name, "%s", cfg.out_path);
        sprintf(base_name + (ext - cfg.out_path), "%s%s", hist, ext);
    }

    /* set output subfile name */
    if (cfg.verbose && sub_rank == 0)
        printf("global_rank=%d sub_rank=%d outfile=%s\n",global_rank,sub_rank,outfile);

    /* create the output file */
    FILE_CREATE(outfile)

    cfg.open_time = MPI_Wtime() - timing;

    MPI_Barrier(cfg.sub_comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    /* define dimensions, variables, and attributes */
    if (cfg.nvars == 414) { /* for h0 file */
        err = def_F_case_h0(cfg, decom, driver, ncid, varids);
        CHECK_ERR
    }
    else { /* for h1 file */
        err = def_F_case_h1(cfg, decom, driver, ncid, varids);
        CHECK_ERR
    }

    /* exit define mode and enter data mode */
    ENDDEF

    /* I/O amount so far */
    INQ_PUT_SIZE(metadata_size)
    metadata_size -= previous_size;

    INQ_FILE_INFO(info_used)

    cfg.def_time = MPI_Wtime() - timing;

    MPI_Barrier(cfg.sub_comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    i = 0;
    my_nreqs = 0;

    /* First, write decomposition variables */
    for (j=0; j<decom.num_decomp; j++) {
        nvars_D[j] = 0; /* number of variables using decomposition j */
        start[0] = sub_rank;
        count[0] = 1;
        start[1] = 0;
        count[1] = decom.contig_nreqs[j];

        /* write to D*.nreqs, 1D array */
        contig_nreqs[j] = decom.contig_nreqs[j];
        IPUT_VARA_INT_NOADV(contig_nreqs+j)

        /* write to D*.blob_start, 1D array */
        blob_start[j] = decom.start[j];
        IPUT_VARA_INT64_NOADV(blob_start+j)

        /* write to D*.blob_count, 1D array */
        blob_count[j] = decom.count[j];
        IPUT_VARA_INT64_NOADV(blob_count+j);

        /* write to D*.offsets, 2D array */
        IPUT_VARA_INT_NOADV(decom.disps[j])

        /* write to D*.lengths, 2D array */
        IPUT_VARA_INT_NOADV(decom.blocklens[j])
    }

    /* Now, write climate variables */
    fix_dbl_buf_ptr = fix_dbl_buf;

    /* lat */
    start[0] = decom.start[1];
    count[0] = decom.count[1];
    IPUT_VARA_DBL(fix_dbl_buf_ptr, decom.count[1])
    nvars_D[1]++;

    /* lon */
    IPUT_VARA_DBL(fix_dbl_buf_ptr, decom.count[1])
    nvars_D[1]++;

    /* area */
    start[0] = decom.start[0];
    count[0] = decom.count[0];
    IPUT_VARA_DBL(fix_dbl_buf_ptr, decom.count[0])
    nvars_D[0]++;

    /* post iput for the remaining fixed-size variables */
    fix_int_buf_ptr = fix_int_buf;
    err = write_small_fix_vars_F_case(driver, ncid, num_decomp_vars, gap,
                                      decom.dims[2][0],
                                      decom.dims[2][0] + 1,
                                      &fix_int_buf_ptr, &fix_dbl_buf_ptr,
                                      &my_nreqs);
    CHECK_ERR

    rec_int_buf_ptr = rec_int_buf;
    rec_txt_buf_ptr = rec_txt_buf;
    rec_dbl_buf_ptr = rec_dbl_buf;
    rec_buf_ptr     = rec_buf;

    /* now write record variables */
    for (rec_no=0; rec_no<cfg.nrecs; rec_no++) {
        if (!one_flush || (cfg.strategy == blob && cfg.api == hdf5)) {
            /* reset the pointers to the beginning of the buffers */
            rec_int_buf_ptr = rec_int_buf;
            rec_txt_buf_ptr = rec_txt_buf;
            rec_dbl_buf_ptr = rec_dbl_buf;
            rec_buf_ptr     = rec_buf;
        }

        /* next 27 small variables are written by sub_rank 0 only */
        if (sub_rank == 0) {
            /* write 15 small record variables by rank 0 only */
            err = write_small_rec_vars_F_case(driver, ncid, num_decomp_vars,
                                              rec_no, gap, 2, 8,
                                              &rec_int_buf_ptr,
                                              &rec_txt_buf_ptr,
                                              &rec_dbl_buf_ptr,
                                              &my_nreqs);
            CHECK_ERR
        }

        start_D2[0] = rec_no; start_D2[1] = decom.start[1];
        count_D2[0] = 1;      count_D2[1] = decom.count[1];
        start_D3[0] = rec_no; start_D3[1] = decom.start[2];
        count_D3[0] = 1;      count_D3[1] = decom.count[2];

        if (cfg.nvars == 414) { /* h0 file */
            if (cfg.two_buf) {
                /* write 2D variables */
                POST_VARA(2,   1,  30)  /* AEROD_v */
                POST_VARA(2,  19,  33)  /* AODABS ... AODVIS */
                POST_VARA(2,   6,  54)  /* AQ_DMS ... AQ_SOAG */
                POST_VARA(2,   4,  64)  /* BURDEN1 ... BURDEN4 */
                POST_VARA(2,   2,  69)  /* CDNUMC and CLDHGH */
                POST_VARA(2,   3,  73)  /* CLDLOW ... CLDTOT */
                POST_VARA(2,  11,  80)  /* DF_DMS ... DSTSFMBL */
                POST_VARA(2,   2,  92)  /* DTENDTH and DTENDTQ */
                POST_VARA(2,   7,  96)  /* FLDS ... FLUTC */
                POST_VARA(2,  15, 107)  /* FSDS ... ICEFRAC */
                POST_VARA(2,   2, 125)  /* LANDFRAC and LHFLX */
                POST_VARA(2,   1, 131)  /* LINOZ_SFCSINK */
                POST_VARA(2,   3, 133)  /* LINOZ_SZA ... LWCF */
                POST_VARA(2,   2, 148)  /* O3_SRF and OCNFRAC */
                POST_VARA(2,   1, 151)  /* OMEGA500 */
                POST_VARA(2,   8, 153)  /* PBLH ... PSL */
                POST_VARA(2,   2, 162)  /* QFLX and QREFHT */
                POST_VARA(2,   1, 167)  /* RAM1 */
                POST_VARA(2,  37, 169)  /* SFDMS ... SNOWHLND */
                POST_VARA(2,  10, 208)  /* SO2_CLXF ... SWCF */
                POST_VARA(2,  19, 219)  /* TAUGWX ... TVQ */
                POST_VARA(2,   1, 239)  /* U10 */
                POST_VARA(2,   3, 246)  /* WD_H2O2 ... WD_SO2 */
                POST_VARA(2,  32, 252)  /* airFV ... dst_c3SFWET */
                POST_VARA(2, 129, 285)  /* mlip ... soa_c3SFWET */
                /* write 3D variables */
                POST_VARA(3,   2,  31)  /* ANRAIN and ANSNOW */
                POST_VARA(3,   2,  52)  /* AQRAIN and AQSNOW */
                POST_VARA(3,   4,  60)  /* AREI ... AWNI */
                POST_VARA(3,   1,  68)  /* CCN3 */
                POST_VARA(3,   2,  71)  /* CLDICE and CLDLIQ */
                POST_VARA(3,   4,  76)  /* CLOUD ... DCQ */
                POST_VARA(3,   1,  91)  /* DTCOND */
                POST_VARA(3,   2,  94)  /* EXTINCT and FICE */
                POST_VARA(3,   4, 103)  /* FREQI ... FREQS */
                POST_VARA(3,   3, 122)  /* ICIMR ... IWC */
                POST_VARA(3,   4, 127)  /* LINOZ_DO3 ... LINOZ_O3COL */
                POST_VARA(3,   1, 132)  /* LINOZ_SSO3 */
                POST_VARA(3,  12, 136)  /* Mass_bc ... O3 */
                POST_VARA(3,   1, 150)  /* OMEGA */
                POST_VARA(3,   1, 152)  /* OMEGAT */
                POST_VARA(3,   1, 161)  /* Q */
                POST_VARA(3,   3, 164)  /* QRL ... RAINQM */
                POST_VARA(3,   1, 168)  /* RELHUM */
                POST_VARA(3,   2, 206)  /* SNOWQM and SO2 */
                POST_VARA(3,   1, 218)  /* T */
                POST_VARA(3,   1, 238)  /* U */
                POST_VARA(3,   6, 240)  /* UU ... VV */
                POST_VARA(3,   3, 249)  /* WSUB ... aero_water */
                POST_VARA(3,   1, 284)  /* hstobie_linoz */
            } else {
                /* write variables in the same order as they defined */
                POST_VARA(2,   1,  30)  /* AEROD_v */
                POST_VARA(3,   2,  31)  /* ANRAIN and ANSNOW */
                POST_VARA(2,  19,  33)  /* AODABS ... AODVIS */
                POST_VARA(3,   2,  52)  /* AQRAIN and AQSNOW */
                POST_VARA(2,   6,  54)  /* AQ_DMS ... AQ_SOAG */
                POST_VARA(3,   4,  60)  /* AREI ... AWNI */
                POST_VARA(2,   4,  64)  /* BURDEN1 ... BURDEN4 */
                POST_VARA(3,   1,  68)  /* CCN3 */
                POST_VARA(2,   2,  69)  /* CDNUMC and CLDHGH */
                POST_VARA(3,   2,  71)  /* CLDICE and CLDLIQ */
                POST_VARA(2,   3,  73)  /* CLDLOW ... CLDTOT */
                POST_VARA(3,   4,  76)  /* CLOUD ... DCQ */
                POST_VARA(2,  11,  80)  /* DF_DMS ... DSTSFMBL */
                POST_VARA(3,   1,  91)  /* DTCOND */
                POST_VARA(2,   2,  92)  /* DTENDTH and DTENDTQ */
                POST_VARA(3,   2,  94)  /* EXTINCT and FICE */
                POST_VARA(2,   7,  96)  /* FLDS ... FLUTC */
                POST_VARA(3,   4, 103)  /* FREQI ... FREQS */
                POST_VARA(2,  15, 107)  /* FSDS ... ICEFRAC */
                POST_VARA(3,   3, 122)  /* ICIMR ... IWC */
                POST_VARA(2,   2, 125)  /* LANDFRAC and LHFLX */
                POST_VARA(3,   4, 127)  /* LINOZ_DO3 ... LINOZ_O3COL */
                POST_VARA(2,   1, 131)  /* LINOZ_SFCSINK */
                POST_VARA(3,   1, 132)  /* LINOZ_SSO3 */
                POST_VARA(2,   3, 133)  /* LINOZ_SZA ... LWCF */
                POST_VARA(3,  12, 136)  /* Mass_bc ... O3 */
                POST_VARA(2,   2, 148)  /* O3_SRF and OCNFRAC */
                POST_VARA(3,   1, 150)  /* OMEGA */
                POST_VARA(2,   1, 151)  /* OMEGA500 */
                POST_VARA(3,   1, 152)  /* OMEGAT */
                POST_VARA(2,   8, 153)  /* PBLH ... PSL */
                POST_VARA(3,   1, 161)  /* Q */
                POST_VARA(2,   2, 162)  /* QFLX and QREFHT */
                POST_VARA(3,   3, 164)  /* QRL ... RAINQM */
                POST_VARA(2,   1, 167)  /* RAM1 */
                POST_VARA(3,   1, 168)  /* RELHUM */
                POST_VARA(2,  37, 169)  /* SFDMS ... SNOWHLND */
                POST_VARA(3,   2, 206)  /* SNOWQM and SO2 */
                POST_VARA(2,  10, 208)  /* SO2_CLXF ... SWCF */
                POST_VARA(3,   1, 218)  /* T */
                POST_VARA(2,  19, 219)  /* TAUGWX ... TVQ */
                POST_VARA(3,   1, 238)  /* U */
                POST_VARA(2,   1, 239)  /* U10 */
                POST_VARA(3,   6, 240)  /* UU ... VV */
                POST_VARA(2,   3, 246)  /* WD_H2O2 ... WD_SO2 */
                POST_VARA(3,   3, 249)  /* WSUB ... aero_water */
                POST_VARA(2,  32, 252)  /* airFV ... dst_c3SFWET */
                POST_VARA(3,   1, 284)  /* hstobie_linoz */
                POST_VARA(2, 129, 285)  /* mlip ... soa_c3SFWET */
            }
        } else { /* h1 file */
            if (cfg.two_buf) {
                /* write 2D variables followed by 3D variables */
                POST_VARA(2, 13, 30)  /* CLDHGH ... T5 */
                POST_VARA(2,  7, 44)  /* U250 ... Z500 */
                POST_VARA(3,  1, 43)  /* U */
            } else {
                /* write variables in the same order as they defined */
                POST_VARA(2, 13, 30)  /* CLDHGH ... T5 */
                POST_VARA(3,  1, 43)  /* U */
                POST_VARA(2,  7, 44)  /* U250 ... Z500 */
            }
        }

        if (!one_flush) {
            cfg.post_time += MPI_Wtime() - timing;

            MPI_Barrier(cfg.sub_comm); /*-----------------------------------*/
            timing = MPI_Wtime();

            /* flush once per time record */
            WAIT_ALL_REQS

            cfg.flush_time += MPI_Wtime() - timing;

            timing = MPI_Wtime();
        }
    }

    if (one_flush) {
        cfg.post_time += MPI_Wtime() - timing;

        MPI_Barrier(cfg.sub_comm); /*---------------------------------------*/
        timing = MPI_Wtime();

        /* flush once for all time records */
        WAIT_ALL_REQS

        cfg.flush_time += MPI_Wtime() - timing;
    }
    MPI_Barrier(cfg.sub_comm); /*---------------------------------------*/
    timing = MPI_Wtime();

    if (fix_int_buf != NULL) free(fix_int_buf);
    if (fix_dbl_buf != NULL) free(fix_dbl_buf);
    if (rec_int_buf != NULL) free(rec_int_buf);
    if (rec_txt_buf != NULL) free(rec_txt_buf);
    if (rec_dbl_buf != NULL) free(rec_dbl_buf);
    if (rec_buf     != NULL) free(rec_buf);
    if (varids      != NULL) free(varids);

    /* obtain the write amount made by far */
    if (cfg.api == pnetcdf) INQ_PUT_SIZE(total_size)

    FILE_CLOSE

    cfg.close_time = MPI_Wtime() - timing;

    /* for hdf5 blob I/O, write amount is calculated after file is closed */
    if (cfg.api == hdf5) total_size = cfg.amount_WR;

    total_size -= previous_size;

    cfg.num_flushes     = nflushes;
    cfg.num_decomp      = decom.num_decomp;
    cfg.num_decomp_vars = num_decomp_vars;
    cfg.my_nreqs        = my_nreqs;
    cfg.metadata_WR     = metadata_size;
    cfg.amount_WR       = total_size;
    cfg.end2end_time    = MPI_Wtime() - cfg.end2end_time;

    /* check if there is any PnetCDF internal malloc residue */
    check_malloc(&cfg, &driver);

    /* report timing breakdowns */
    report_timing_WR(&cfg, &driver, base_name);

    /* print MPI-IO hints actually used */
    if (cfg.verbose && global_rank == 0) print_info(&info_used);

err_out:
    if (err < 0 && ncid >= 0)
#ifdef USE_PNETCDF_DIRECTLY
        ncmpi_close(ncid);
#else
        driver.close(ncid);
#endif

    if (info_used != MPI_INFO_NULL) MPI_Info_free(&info_used);
    if (!cfg.keep_outfile && sub_rank == 0) unlink(outfile);
    fflush(stdout);

    return err;
}

