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

#include <e3sm_io_case_F.hpp>
#include <e3sm_io.h>
#include <e3sm_io_err.h>
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

/*----< write_small_vars_F_case() >------------------------------------------*/
static
int write_small_vars_F_case(e3sm_io_driver &driver,
                            int             ncid,
                            int             vid, /* starting variable ID */
                            int            *varids,
                            int             rec_no,
                            int             gap,
                            MPI_Offset      lev,
                            MPI_Offset      ilev,
                            MPI_Offset      nbnd,
                            MPI_Offset      nchars,
                            int           **int_buf,
                            char          **txt_buf,
                            double        **dbl_buf,
                            int            *nreqs)
{
    int i, err, my_nreqs=0;
    MPI_Offset start[2], count[2];

    /* There are 27 variables in total written in this function.
     * double type:
     *        3 [lev], 3 [ilev], 1 [nbnd], 8 [1]
     * int type:
     *        10 [1]
     * char type:
     *        2 [nchars]
     */

    /* scalar and small variables are written by sub_rank 0 only */
    i = vid;

    if (rec_no == 0) {
        /* write 7 fixed-size variables */
        IPUT_VAR_DBL(*dbl_buf,  lev + gap) /* lev */
        IPUT_VAR_DBL(*dbl_buf,  lev + gap) /* hyam */
        IPUT_VAR_DBL(*dbl_buf,  lev + gap) /* hybm */
        IPUT_VAR_DBL(*dbl_buf,    1 + gap) /* P0 */
        IPUT_VAR_DBL(*dbl_buf, ilev + gap) /* ilev */
        IPUT_VAR_DBL(*dbl_buf, ilev + gap) /* hyai */
        IPUT_VAR_DBL(*dbl_buf, ilev + gap) /* hybi */
    } else
        i += 7;

    /* below are all record variables */
    start[0] = rec_no;
    start[1] = 0;
    count[0] = 1;

    IPUT_VAR1_DBL(*dbl_buf,       1 + gap) /* time */
    IPUT_VAR1_INT(*int_buf,       1 + gap) /* date */
    IPUT_VAR1_INT(*int_buf,       1 + gap) /* datesec */
    count[1] = nbnd;
    IPUT_VARA_DBL(*dbl_buf,    nbnd + gap) /* time_bnds */
    count[1] = nchars;
    IPUT_VARA_TXT(*txt_buf, nchars + gap) /* date_written */
    count[1] = nchars;
    IPUT_VARA_TXT(*txt_buf, nchars + gap) /* time_written */

    if (rec_no == 0) {
        IPUT_VAR_INT(*int_buf, 1) /* ndbase */
        IPUT_VAR_INT(*int_buf, 1) /* nsbase */
        IPUT_VAR_INT(*int_buf, 1) /* nbdate */
        IPUT_VAR_INT(*int_buf, 1) /* nbsec */
        IPUT_VAR_INT(*int_buf, 1) /* mdt */
    } else
        i += 5;

    IPUT_VAR1_INT(*int_buf, 1 + gap) /* ndcur */
    IPUT_VAR1_INT(*int_buf, 1 + gap) /* nscur */
    IPUT_VAR1_DBL(*dbl_buf, 1 + gap) /* co2vmr */
    IPUT_VAR1_DBL(*dbl_buf, 1 + gap) /* ch4vmr */
    IPUT_VAR1_DBL(*dbl_buf, 1 + gap) /* n2ovmr */
    IPUT_VAR1_DBL(*dbl_buf, 1 + gap) /* f11vmr */
    IPUT_VAR1_DBL(*dbl_buf, 1 + gap) /* f12vmr */
    IPUT_VAR1_DBL(*dbl_buf, 1 + gap) /* sol_tsi */
    IPUT_VAR1_INT(*int_buf, 1 + gap) /* nsteph */

    *nreqs = my_nreqs;

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
    char outfile[1040], *ext;;
    const char *hist;
    int i, j, err, sub_rank, global_rank, ncid=-1, nflushes=0, *varids;
    int rec_no, gap = 0, my_nreqs, num_decomp_vars;
    int contig_nreqs[MAX_NUM_DECOMP], nvars_D[MAX_NUM_DECOMP];
    double pre_timing, open_timing, post_timing, wait_timing, close_timing;
    double def_timing, timing, total_timing;
    MPI_Offset metadata_size, put_size, total_size, max_nreqs, total_nreqs;
    MPI_Offset start[2], count[2];
    MPI_Offset start_D2[2], count_D2[2], start_D3[2], count_D3[2];
    MPI_Offset blob_start[MAX_NUM_DECOMP], blob_count[MAX_NUM_DECOMP];
    MPI_Offset previous_size, sum_size, m_alloc=0, max_alloc;
    MPI_Info info_used = MPI_INFO_NULL;
    size_t ii, txt_buflen, int_buflen, dbl_buflen, rec_buflen;
    itype  *rec_buf=NULL, *rec_buf_ptr; /* buffer for rec float var */
    double *dbl_buf=NULL, *dbl_buf_ptr; /* buffer for fixed double var */
    char   *txt_buf=NULL, *txt_buf_ptr; /* buffer for char var */
    int    *int_buf=NULL, *int_buf_ptr; /* buffer for int var */

    MPI_Barrier(cfg.io_comm); /*-----------------------------------------*/
    total_timing = pre_timing = MPI_Wtime();

    post_timing  = 0.0;
    open_timing  = 0.0;
    def_timing   = 0.0;
    wait_timing  = 0.0;
    close_timing = 0.0;

    if (cfg.api == hdf5) /* I/O amount from previous I/O */
        INQ_PUT_SIZE(previous_size)
    else
        previous_size = 0;

    MPI_Comm_rank(cfg.io_comm,  &global_rank);
    MPI_Comm_rank(cfg.sub_comm, &sub_rank);

    if (cfg.non_contig_buf) gap = 10;

    /* allocate write buffer for small climate variables */
    dbl_buflen = decom.count[1] * 2          /* lat[ncol] and lon[ncol] */
               + decom.count[0]              /* area[ncol] */
               + 3 * gap
               + 3 * decom.dims[2][0]        /* 3 [lev] */
               + 3 * (decom.dims[2][0] + 1)  /* 3 [ilev] */
               + 2                           /* 1 [nbnd] */
               + 8                           /* 8 single-element variables */
               + 15 * gap;

    txt_buflen = 2 * 8     /* 2 [nchars] */
               + 10 * gap;
    int_buflen = 10        /* 10 [1] */
               + 10 * gap;

    /* allocate and initialize write buffer for large variables */
    if (cfg.nvars == 414)
        rec_buflen = decom.count[1] * 323
                   + decom.count[2] * 63
                   + (323 + 63) * gap;
    else
        rec_buflen = decom.count[1] * 22
                   + decom.count[2]
                   + (22 + 1) * gap;

#define FLUSH_ALL_RECORDS_AT_ONCE
#ifdef FLUSH_ALL_RECORDS_AT_ONCE
    dbl_buflen *= cfg.nrec;
    txt_buflen *= cfg.nrec;
    int_buflen *= cfg.nrec;
    rec_buflen *= cfg.nrec;
#else
    if (cfg.api == hdf5) {
        printf("Error in %s:%d: %s() FLUSH_ALL_RECORDS_AT_ONCE must be enabled when using HDF5 blob I/O",
               __FILE__, __LINE__, __func__);
        err = -1;
        goto err_out;
    }
#endif

    dbl_buf = (double*) malloc(dbl_buflen * sizeof(double));
    txt_buf = (char*)   malloc(txt_buflen * sizeof(char));
    int_buf = (int*)    malloc(int_buflen * sizeof(int));
    rec_buf = (itype*)  malloc(rec_buflen * sizeof(itype));

    for (ii=0; ii<dbl_buflen; ii++) dbl_buf[ii] = global_rank;
    for (ii=0; ii<txt_buflen; ii++) txt_buf[ii] = 'a' + global_rank;
    for (ii=0; ii<rec_buflen; ii++) rec_buf[ii] = global_rank;
    for (ii=0; ii<int_buflen; ii++) int_buf[ii] = global_rank;

    /* there are num_decomp_vars number of decomposition variables */
    num_decomp_vars = decom.num_decomp * NVARS_DECOMP;

    /* allocate space for all variable IDs */
    varids = (int*) malloc((cfg.nvars + num_decomp_vars) * sizeof(int));

    pre_timing = MPI_Wtime() - pre_timing;

    MPI_Barrier(cfg.sub_comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    /* construct h0/h1 file name */
    hist = (cfg.nvars == 414) ? "_h0" :  "_h1";
    ext = strrchr(cfg.out_path, '.');
    if (ext == NULL || (strcmp(ext, ".nc") && strcmp(ext, ".h5")))
        sprintf(outfile, "%s%s.%04d", cfg.out_path, hist, cfg.subfile_ID);
    else {
        sprintf(outfile, "%s", cfg.out_path);
        sprintf(outfile + (ext - cfg.out_path), "%s%s.%04d", hist, ext, cfg.subfile_ID);
    }

    /* set output subfile name */
    if (cfg.verbose && sub_rank == 0)
        printf("global_rank=%d sub_rank=%d outfile=%s\n",global_rank,sub_rank,outfile);

    /* create the output file */
    FILE_CREATE(outfile)

    open_timing += MPI_Wtime() - timing;

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

    def_timing += MPI_Wtime() - timing;

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
    dbl_buf_ptr = dbl_buf;

    /* lat */
    start[0] = decom.start[1];
    count[0] = decom.count[1];
    IPUT_VARA_DBL(dbl_buf_ptr, decom.count[1])
    nvars_D[1]++;

    /* lon */
    IPUT_VARA_DBL(dbl_buf_ptr, decom.count[1])
    nvars_D[1]++;

    /* area */
    start[0] = decom.start[0];
    count[0] = decom.count[0];
    IPUT_VARA_DBL(dbl_buf_ptr, decom.count[0])
    nvars_D[0]++;

#ifdef FLUSH_ALL_RECORDS_AT_ONCE
    rec_buf_ptr = rec_buf;
    int_buf_ptr = int_buf;
    txt_buf_ptr = txt_buf;
#endif

    for (rec_no=0; rec_no<cfg.nrec; rec_no++) {
#ifndef FLUSH_ALL_RECORDS_AT_ONCE
        dbl_buf_ptr = dbl_buf + decom.count[1] * 2 + decom.count[0] + gap * 3;
        int_buf_ptr = int_buf;
        txt_buf_ptr = txt_buf;
        rec_buf_ptr = rec_buf;
#endif
        i = 3 + num_decomp_vars;   /* 3 is lat, lon, and area */

        /* next 27 small variables are written by sub_rank 0 only */
        if (sub_rank == 0) {
            /* post nonblocking requests for small datasets */
	    err = write_small_vars_F_case(driver, ncid, i, varids, rec_no, gap,
					  decom.dims[2][0], decom.dims[2][0] +
					  1, 2, 8, &int_buf_ptr, &txt_buf_ptr,
                                          &dbl_buf_ptr, &my_nreqs);
            CHECK_ERR
        }
        i += 27;

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
#ifndef FLUSH_ALL_RECORDS_AT_ONCE
        post_timing += MPI_Wtime() - timing;

        MPI_Barrier(cfg.sub_comm); /*---------------------------------------*/
        timing = MPI_Wtime();

        /* flush once per time record */
        WAIT_ALL_REQS

        wait_timing += MPI_Wtime() - timing;

        timing = MPI_Wtime();
#endif
        if (cfg.nvars == 414)
            break; /* h0 file stores only one time stamp */
    }
#ifdef FLUSH_ALL_RECORDS_AT_ONCE
    post_timing += MPI_Wtime() - timing;

    MPI_Barrier(cfg.sub_comm); /*---------------------------------------*/
    timing = MPI_Wtime();

    /* flush once for all time records */
    WAIT_ALL_REQS

    wait_timing += MPI_Wtime() - timing;
#endif
    MPI_Barrier(cfg.sub_comm); /*---------------------------------------*/
    timing = MPI_Wtime();

    if (rec_buf != NULL) free(rec_buf);
    if (dbl_buf != NULL) free(dbl_buf);
    if (txt_buf != NULL) free(txt_buf);
    if (int_buf != NULL) free(int_buf);
    free(varids);

    if (cfg.api == pnetcdf) INQ_PUT_SIZE(total_size)

    FILE_CLOSE

    close_timing += MPI_Wtime() - timing;

    if (cfg.api == hdf5) INQ_PUT_SIZE(total_size)

    total_size -= previous_size;
    put_size = total_size - metadata_size;

    total_timing = MPI_Wtime() - total_timing;

    /* check if there is any PnetCDF internal malloc residue */
    if (cfg.api == pnetcdf) {
        MPI_Offset malloc_size;
        err = driver.inq_malloc_size(&malloc_size);
        if (err == NC_NOERR) {
            MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, cfg.io_comm);
            if (global_rank == 0 && sum_size > 0) {
                printf("-----------------------------------------------------------\n");
                printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                        sum_size);
            }
        }
        driver.inq_malloc_max_size(&m_alloc);
    }

    MPI_Offset off_tmp[3], sum_off[3];
    off_tmp[0] = my_nreqs;
    off_tmp[1] = put_size;
    off_tmp[2] = total_size;
    MPI_Reduce(off_tmp, sum_off, 3, MPI_OFFSET, MPI_SUM, 0, cfg.io_comm);
    total_nreqs = sum_off[0];
    put_size    = sum_off[1];
    total_size  = sum_off[2];

    MPI_Offset max_off[2];
    off_tmp[0] = my_nreqs;
    off_tmp[1] = m_alloc;
    MPI_Reduce(off_tmp, max_off, 2, MPI_OFFSET, MPI_MAX, 0, cfg.io_comm);
    max_nreqs = max_off[0];
    max_alloc = max_off[1];

    double dbl_tmp[7], max_dbl[7];
    dbl_tmp[0] =   pre_timing;
    dbl_tmp[1] =  open_timing;
    dbl_tmp[2] =   def_timing;
    dbl_tmp[3] =  post_timing;
    dbl_tmp[4] =  wait_timing;
    dbl_tmp[5] = close_timing;
    dbl_tmp[6] = total_timing;
    MPI_Reduce(dbl_tmp, max_dbl, 7, MPI_DOUBLE, MPI_MAX, 0, cfg.io_comm);
      pre_timing = max_dbl[0];
     open_timing = max_dbl[1];
      def_timing = max_dbl[2];
     post_timing = max_dbl[3];
     wait_timing = max_dbl[4];
    close_timing = max_dbl[5];
    total_timing = max_dbl[6];

    if (global_rank == 0) {
        int nvars_noD = cfg.nvars;
        for (i = 0; i < 3; i++) nvars_noD -= nvars_D[i];
        printf("History output file                = %s\n", outfile);
        printf("No. decomposition variables        = %3d\n", num_decomp_vars);
        printf("No. variables use no decomposition = %3d\n", nvars_noD);
        printf("No. variables use decomposition D1 = %3d\n", nvars_D[0]);
        printf("No. variables use decomposition D2 = %3d\n", nvars_D[1]);
        printf("No. variables use decomposition D3 = %3d\n", nvars_D[2]);
        printf("Total number of variables          = %3d\n", cfg.nvars + num_decomp_vars);
        if (cfg.api == pnetcdf)
            printf("MAX heap memory allocated by PnetCDF internally is %.2f MiB\n",
                   (float)max_alloc / 1048576);
        printf("Total write amount                 = %.2f MiB = %.2f GiB\n",
               (double)total_size / 1048576, (double)total_size / 1073741824);
        printf("Total no. noncontiguous requests   = %lld\n", total_nreqs);
        printf("Max   no. noncontiguous requests   = %lld\n", max_nreqs);
        if (cfg.api == pnetcdf)
            printf("No. I/O flush calls                = %d\n", nflushes);
        printf("Max Time of I/O preparing          = %.4f sec\n",   pre_timing);
        printf("Max Time of file open/create       = %.4f sec\n",  open_timing);
        printf("Max Time of define variables       = %.4f sec\n",   def_timing);
        printf("Max Time of posting iput requests  = %.4f sec\n",  post_timing);
        if (cfg.api == pnetcdf)
            printf("Max Time of write flushing         = %.4f sec\n",  wait_timing);
        printf("Max Time of close                  = %.4f sec\n", close_timing);
        printf("Max Time of TOTAL                  = %.4f sec\n", total_timing);
        printf("I/O bandwidth (open-to-close)      = %.4f MiB/sec\n",
               (double)total_size / 1048576.0 / total_timing);
        if (cfg.api == pnetcdf)
            printf("I/O bandwidth (write-only)         = %.4f MiB/sec\n",
                   (double)put_size / 1048576.0 / wait_timing);
        if (cfg.verbose) print_info(&info_used);
        printf("-----------------------------------------------------------\n");
    }

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

