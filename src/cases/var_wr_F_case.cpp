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
    err = driver.create(filename, comm, cfg.info, &ncid); \
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
#define WAIT_ALL_REQS { \
    err = driver.wait(ncid); \
    CHECK_ERR \
    nflushes++; \
}

#define FIX_IPUT_DBL(dp, varid) {                                         \
    if (cfg.strategy == canonical) {                                      \
        err = driver.put_varn(ncid, varid, MPI_DOUBLE,                    \
                              decom.contig_nreqs[dp],                     \
                              decom.w_startx[dp], decom.w_countx[dp],     \
                              fix_dbl_buf_ptr, nb);                       \
        my_nreqs += decom.contig_nreqs[dp];                               \
    }                                                                     \
    else {                                                                \
        err = driver.put_vara(ncid, varid, MPI_DOUBLE, &decom.start[dp],  \
                              &decom.count[dp], fix_dbl_buf_ptr, nb);     \
        my_nreqs++;                                                       \
    }                                                                     \
    CHECK_VAR_ERR(varid)                                                  \
    fix_dbl_buf_ptr += decom.count[dp] + gap;                             \
    nvars_D[dp]++;                                                        \
}

#define REC_IPUT(dp, numVars, vid) {                                      \
    int varid = vid + num_decomp_vars;                                    \
    if (cfg.strategy == blob) {                                           \
        for (j=0; j<numVars; j++) {                                       \
            err = driver.put_vara(ncid, varid, REC_ITYPE, starts[dp],     \
                                  counts[dp], rec_buf_ptr, nb);           \
            CHECK_VAR_ERR(varid)                                          \
            rec_buf_ptr += decom.count[dp] + gap;                         \
            my_nreqs++;                                                   \
            varid++;                                                      \
        }                                                                 \
    }                                                                     \
    else {                                                                \
        for (j=0; j<numVars; j++) {                                       \
            err = driver.put_varn(ncid, varid, REC_ITYPE,                 \
                                  decom.contig_nreqs[dp],                 \
                                  decom.w_starts[dp], decom.w_counts[dp], \
                                  rec_buf_ptr, nb);                       \
            CHECK_VAR_ERR(varid)                                          \
            rec_buf_ptr += decom.count[dp] + gap;                         \
            my_nreqs += decom.contig_nreqs[dp];                           \
            varid++;                                                      \
        }                                                                 \
    }                                                                     \
    if (rec_no == 0) nvars_D[dp] += numVars;                              \
}

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
    /* small fixed-size variable IDs (relative) */
    int varids[12] = {3, 4, 5, 6, 7, 8, 9,
                      16, 17, 18, 19, 20};

    int i, err, my_nreqs=0;

    for (i=0; i<12; i++) varids[i] += start_varid;

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

    if (nreqs != NULL) *nreqs += my_nreqs;

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
    /* small record variable IDs (relative) */
    int varids[15] = {10, 11, 12, 13, 14, 15,
                      21, 22, 23, 24, 25, 26, 27, 28, 29};

    int i, err, my_nreqs=0;
    MPI_Offset start[2], count[2];

    for (i=0; i<15; i++) varids[i] += start_varid;

    start[0] = rec_no;
    start[1] = 0;
    count[0] = 1;
    i = 0;

    /* 8 of type double, 5 of type int, and 2 of type char */
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

    if (nreqs != NULL) *nreqs += my_nreqs;

err_out:
    return err;
}

/*----< var_wr_F_case() >----------------------------------------------------*/
int var_wr_F_case(e3sm_io_config &cfg,
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
    MPI_Offset blob_start[MAX_NUM_DECOMP], blob_count[MAX_NUM_DECOMP];
    MPI_Offset starts[MAX_NUM_DECOMP][2], counts[MAX_NUM_DECOMP][2];
    MPI_Info info_used = MPI_INFO_NULL;
    MPI_Comm comm;

    size_t ii, rec_buflen;
    size_t                 fix_int_buflen, fix_dbl_buflen;
    size_t rec_txt_buflen, rec_int_buflen, rec_dbl_buflen;
    int    *fix_int_buf=NULL, *fix_int_buf_ptr;
    double *fix_dbl_buf=NULL, *fix_dbl_buf_ptr;
    int    *rec_int_buf=NULL, *rec_int_buf_ptr;
    char   *rec_txt_buf=NULL, *rec_txt_buf_ptr;
    double *rec_dbl_buf=NULL, *rec_dbl_buf_ptr;
    vtype  *rec_buf=NULL, *rec_buf_ptr;

    MPI_Barrier(cfg.io_comm); /*-----------------------------------------*/
    cfg.end2end_time = cfg.pre_time = MPI_Wtime();

    cfg.post_time  = 0.0;
    cfg.flush_time = 0.0;

    comm = (cfg.strategy == blob) ? cfg.sub_comm : cfg.io_comm;

    MPI_Comm_rank(cfg.io_comm,  &global_rank);
    MPI_Comm_rank(comm,         &sub_rank);

    /* I/O amount from previous I/O */
    previous_size = cfg.amount_WR;

#define FLUSH_ALL_RECORDS_AT_ONCE
#ifdef FLUSH_ALL_RECORDS_AT_ONCE
    one_flush = 1;
#else
    one_flush = 0;
    if (cfg.strategy == blob && cfg.api != pnetcdf) {
        printf("Error in %s:%d: %s() FLUSH_ALL_RECORDS_AT_ONCE must be enabled for all blob I/O except PnetCDF",
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
                   + 10 * gap;

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
        rec_buflen = 321 * decom.count[1]
                   +  63 * decom.count[2]
                   + (321 + 63) * gap;
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

    /* allocate and initialize write buffers */
    fix_dbl_buf = (double*) malloc(fix_dbl_buflen * sizeof(double));
    fix_int_buf = (int*)    malloc(fix_int_buflen * sizeof(int));
    rec_dbl_buf = (double*) malloc(rec_dbl_buflen * sizeof(double));
    rec_txt_buf = (char*)   malloc(rec_txt_buflen * sizeof(char));
    rec_int_buf = (int*)    malloc(rec_int_buflen * sizeof(int));
    rec_buf     = (vtype*)  malloc(rec_buflen     * sizeof(vtype));

    for (ii=0; ii<fix_dbl_buflen; ii++) fix_dbl_buf[ii] = global_rank;
    for (ii=0; ii<fix_int_buflen; ii++) fix_int_buf[ii] = global_rank;
    for (ii=0; ii<rec_dbl_buflen; ii++) rec_dbl_buf[ii] = global_rank;
    for (ii=0; ii<rec_txt_buflen; ii++) rec_txt_buf[ii] = 'a' + global_rank;
    for (ii=0; ii<rec_int_buflen; ii++) rec_int_buf[ii] = global_rank;
    for (ii=0; ii<rec_buflen;     ii++) rec_buf[ii]     = global_rank;

    /* there are num_decomp_vars number of decomposition variables */
    if (cfg.strategy == blob)
        num_decomp_vars = decom.num_decomp * NVARS_DECOMP;
    else
        num_decomp_vars = 0;

    /* allocate space for all variable IDs */
    varids = (int*) malloc((cfg.nvars + num_decomp_vars) * sizeof(int));

    cfg.pre_time = MPI_Wtime() - cfg.pre_time;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    /* construct h0/h1 file name */
    hist = (cfg.nvars == 414) ? "_h0" :  "_h1";
    ext = strrchr(cfg.out_path, '.');
    if (ext == NULL || (strcmp(ext, ".nc") && strcmp(ext, ".h5"))) {
        sprintf(base_name, "%s%s", cfg.out_path, hist);
        if (cfg.strategy == blob) /* set output subfile name */
            sprintf(outfile, "%s%s.%04d", cfg.out_path, hist, cfg.subfile_ID);
        else
            sprintf(outfile, "%s%s", cfg.out_path, hist);
    }
    else {
        sprintf(outfile, "%s", cfg.out_path);
        if (cfg.strategy == blob) /* set output subfile name */
            sprintf(outfile + (ext - cfg.out_path), "%s%s.%04d", hist, ext,
                    cfg.subfile_ID);
        else
            sprintf(outfile + (ext - cfg.out_path), "%s%s", hist, ext);
        sprintf(base_name, "%s", cfg.out_path);
        sprintf(base_name + (ext - cfg.out_path), "%s%s", hist, ext);
    }

    if (cfg.verbose) {
        if (cfg.strategy == blob && sub_rank == 0)
            printf("global_rank=%d sub_rank=%d outfile=%s\n",
                   global_rank,sub_rank,outfile);
        else if (global_rank == 0)
            printf("global_rank=%d outfile=%s\n",global_rank,outfile);
    }

    /* create the output file */
    FILE_CREATE(outfile)

    cfg.open_time = MPI_Wtime() - timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
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
    if (cfg.api != pnetcdf) metadata_size -= previous_size;

    INQ_FILE_INFO(info_used)

    cfg.def_time = MPI_Wtime() - timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    for (j=0; j<decom.num_decomp; j++)
        nvars_D[j] = 0; /* number of variables using decomposition j */

    i = 0;        /* index used in varids[] */
    my_nreqs = 0; /* number of noncontiguous requests written by this proc */

    if (cfg.strategy == blob) {
        /* First, write decomposition variables */
        for (j=0; j<decom.num_decomp; j++) {
            MPI_Offset start[2], count[2];

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

            starts[j][0] = 0;
            starts[j][1] = decom.start[j];
            counts[j][0] = 1;
            counts[j][1] = decom.count[j];
        }
    }

    /* Now, write climate variables */
    fix_dbl_buf_ptr = fix_dbl_buf;

    /* lat */
    FIX_IPUT_DBL(1, i)

    /* lon */
    FIX_IPUT_DBL(1, i+1)

    /* area */
    FIX_IPUT_DBL(0, i+2)

    /* post iput for the remaining fixed-size variables */
    fix_int_buf_ptr = fix_int_buf;
    if (sub_rank == 0) {
        err = write_small_fix_vars_F_case(driver, ncid, num_decomp_vars, gap,
                                          decom.dims[2][0],
                                          decom.dims[2][0] + 1,
                                          &fix_int_buf_ptr, &fix_dbl_buf_ptr,
                                          &my_nreqs);
        CHECK_ERR
    }

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

        /* set the start index for the next record */
        for (i=0; i<decom.num_decomp; i++) {
            if (cfg.strategy == blob)
                starts[i][0] = rec_no;
            else {
                for (j=0; j<decom.contig_nreqs[i]; j++)
                    decom.w_starts[i][j][0] = rec_no;
            }
        }

        if (cfg.nvars == 414) { /* h0 file */
            if (cfg.two_buf) {
                /* write 2D variables */
                REC_IPUT(1,   1,  30)  /* AEROD_v */
                REC_IPUT(1,  19,  33)  /* AODABS ... AODVIS */
                REC_IPUT(1,   6,  54)  /* AQ_DMS ... AQ_SOAG */
                REC_IPUT(1,   4,  64)  /* BURDEN1 ... BURDEN4 */
                REC_IPUT(1,   2,  69)  /* CDNUMC and CLDHGH */
                REC_IPUT(1,   3,  73)  /* CLDLOW ... CLDTOT */
                REC_IPUT(1,  11,  80)  /* DF_DMS ... DSTSFMBL */
                REC_IPUT(1,   2,  92)  /* DTENDTH and DTENDTQ */
                REC_IPUT(1,   7,  96)  /* FLDS ... FLUTC */
                REC_IPUT(1,  15, 107)  /* FSDS ... ICEFRAC */
                REC_IPUT(1,   2, 125)  /* LANDFRAC and LHFLX */
                REC_IPUT(1,   1, 131)  /* LINOZ_SFCSINK */
                REC_IPUT(1,   3, 133)  /* LINOZ_SZA ... LWCF */
                REC_IPUT(1,   2, 148)  /* O3_SRF and OCNFRAC */
                REC_IPUT(1,   1, 151)  /* OMEGA500 */
                REC_IPUT(1,   8, 153)  /* PBLH ... PSL */
                REC_IPUT(1,   2, 162)  /* QFLX and QREFHT */
                REC_IPUT(1,   1, 167)  /* RAM1 */
                REC_IPUT(1,  37, 169)  /* SFDMS ... SNOWHLND */
                REC_IPUT(1,  10, 208)  /* SO2_CLXF ... SWCF */
                REC_IPUT(1,  19, 219)  /* TAUGWX ... TVQ */
                REC_IPUT(1,   1, 239)  /* U10 */
                REC_IPUT(1,   3, 246)  /* WD_H2O2 ... WD_SO2 */
                REC_IPUT(1,  32, 252)  /* airFV ... dst_c3SFWET */
                REC_IPUT(1, 129, 285)  /* mlip ... soa_c3SFWET */
                /* write 3D variables */
                REC_IPUT(2,   2,  31)  /* ANRAIN and ANSNOW */
                REC_IPUT(2,   2,  52)  /* AQRAIN and AQSNOW */
                REC_IPUT(2,   4,  60)  /* AREI ... AWNI */
                REC_IPUT(2,   1,  68)  /* CCN3 */
                REC_IPUT(2,   2,  71)  /* CLDICE and CLDLIQ */
                REC_IPUT(2,   4,  76)  /* CLOUD ... DCQ */
                REC_IPUT(2,   1,  91)  /* DTCOND */
                REC_IPUT(2,   2,  94)  /* EXTINCT and FICE */
                REC_IPUT(2,   4, 103)  /* FREQI ... FREQS */
                REC_IPUT(2,   3, 122)  /* ICIMR ... IWC */
                REC_IPUT(2,   4, 127)  /* LINOZ_DO3 ... LINOZ_O3COL */
                REC_IPUT(2,   1, 132)  /* LINOZ_SSO3 */
                REC_IPUT(2,  12, 136)  /* Mass_bc ... O3 */
                REC_IPUT(2,   1, 150)  /* OMEGA */
                REC_IPUT(2,   1, 152)  /* OMEGAT */
                REC_IPUT(2,   1, 161)  /* Q */
                REC_IPUT(2,   3, 164)  /* QRL ... RAINQM */
                REC_IPUT(2,   1, 168)  /* RELHUM */
                REC_IPUT(2,   2, 206)  /* SNOWQM and SO2 */
                REC_IPUT(2,   1, 218)  /* T */
                REC_IPUT(2,   1, 238)  /* U */
                REC_IPUT(2,   6, 240)  /* UU ... VV */
                REC_IPUT(2,   3, 249)  /* WSUB ... aero_water */
                REC_IPUT(2,   1, 284)  /* hstobie_linoz */
            } else {
                /* write variables in the same order as they defined */
                REC_IPUT(1,   1,  30)  /* AEROD_v */
                REC_IPUT(2,   2,  31)  /* ANRAIN and ANSNOW */
                REC_IPUT(1,  19,  33)  /* AODABS ... AODVIS */
                REC_IPUT(2,   2,  52)  /* AQRAIN and AQSNOW */
                REC_IPUT(1,   6,  54)  /* AQ_DMS ... AQ_SOAG */
                REC_IPUT(2,   4,  60)  /* AREI ... AWNI */
                REC_IPUT(1,   4,  64)  /* BURDEN1 ... BURDEN4 */
                REC_IPUT(2,   1,  68)  /* CCN3 */
                REC_IPUT(1,   2,  69)  /* CDNUMC and CLDHGH */
                REC_IPUT(2,   2,  71)  /* CLDICE and CLDLIQ */
                REC_IPUT(1,   3,  73)  /* CLDLOW ... CLDTOT */
                REC_IPUT(2,   4,  76)  /* CLOUD ... DCQ */
                REC_IPUT(1,  11,  80)  /* DF_DMS ... DSTSFMBL */
                REC_IPUT(2,   1,  91)  /* DTCOND */
                REC_IPUT(1,   2,  92)  /* DTENDTH and DTENDTQ */
                REC_IPUT(2,   2,  94)  /* EXTINCT and FICE */
                REC_IPUT(1,   7,  96)  /* FLDS ... FLUTC */
                REC_IPUT(2,   4, 103)  /* FREQI ... FREQS */
                REC_IPUT(1,  15, 107)  /* FSDS ... ICEFRAC */
                REC_IPUT(2,   3, 122)  /* ICIMR ... IWC */
                REC_IPUT(1,   2, 125)  /* LANDFRAC and LHFLX */
                REC_IPUT(2,   4, 127)  /* LINOZ_DO3 ... LINOZ_O3COL */
                REC_IPUT(1,   1, 131)  /* LINOZ_SFCSINK */
                REC_IPUT(2,   1, 132)  /* LINOZ_SSO3 */
                REC_IPUT(1,   3, 133)  /* LINOZ_SZA ... LWCF */
                REC_IPUT(2,  12, 136)  /* Mass_bc ... O3 */
                REC_IPUT(1,   2, 148)  /* O3_SRF and OCNFRAC */
                REC_IPUT(2,   1, 150)  /* OMEGA */
                REC_IPUT(1,   1, 151)  /* OMEGA500 */
                REC_IPUT(2,   1, 152)  /* OMEGAT */
                REC_IPUT(1,   8, 153)  /* PBLH ... PSL */
                REC_IPUT(2,   1, 161)  /* Q */
                REC_IPUT(1,   2, 162)  /* QFLX and QREFHT */
                REC_IPUT(2,   3, 164)  /* QRL ... RAINQM */
                REC_IPUT(1,   1, 167)  /* RAM1 */
                REC_IPUT(2,   1, 168)  /* RELHUM */
                REC_IPUT(1,  37, 169)  /* SFDMS ... SNOWHLND */
                REC_IPUT(2,   2, 206)  /* SNOWQM and SO2 */
                REC_IPUT(1,  10, 208)  /* SO2_CLXF ... SWCF */
                REC_IPUT(2,   1, 218)  /* T */
                REC_IPUT(1,  19, 219)  /* TAUGWX ... TVQ */
                REC_IPUT(2,   1, 238)  /* U */
                REC_IPUT(1,   1, 239)  /* U10 */
                REC_IPUT(2,   6, 240)  /* UU ... VV */
                REC_IPUT(1,   3, 246)  /* WD_H2O2 ... WD_SO2 */
                REC_IPUT(2,   3, 249)  /* WSUB ... aero_water */
                REC_IPUT(1,  32, 252)  /* airFV ... dst_c3SFWET */
                REC_IPUT(2,   1, 284)  /* hstobie_linoz */
                REC_IPUT(1, 129, 285)  /* mlip ... soa_c3SFWET */
            }
        } else { /* h1 file */
            if (cfg.two_buf) {
                /* write 2D variables followed by 3D variables */
                REC_IPUT(1, 13, 30)  /* CLDHGH ... T5 */
                REC_IPUT(1,  7, 44)  /* U250 ... Z500 */
                REC_IPUT(2,  1, 43)  /* U */
            } else {
                /* write variables in the same order as they defined */
                REC_IPUT(1, 13, 30)  /* CLDHGH ... T5 */
                REC_IPUT(2,  1, 43)  /* U */
                REC_IPUT(1,  7, 44)  /* U250 ... Z500 */
            }
        }

        if (!one_flush) {
            cfg.post_time += MPI_Wtime() - timing;

            MPI_Barrier(comm); /*-----------------------------------*/
            timing = MPI_Wtime();

            /* flush once per time record */
            WAIT_ALL_REQS

            cfg.flush_time += MPI_Wtime() - timing;

            timing = MPI_Wtime();
        }
    }

    if (one_flush) {
        cfg.post_time += MPI_Wtime() - timing;

        MPI_Barrier(comm); /*---------------------------------------*/
        timing = MPI_Wtime();

        /* flush once for all time records */
        WAIT_ALL_REQS

        cfg.flush_time += MPI_Wtime() - timing;
    }
    MPI_Barrier(comm); /*---------------------------------------*/
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
    if (cfg.api != pnetcdf) total_size = cfg.amount_WR - previous_size;

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
    report_timing_WR(&cfg, &driver, &decom, base_name);

    /* print MPI-IO hints actually used */
    if (cfg.verbose && global_rank == 0) print_info(&info_used);

err_out:
    if (err < 0 && ncid >= 0) driver.close(ncid);
    if (info_used != MPI_INFO_NULL) MPI_Info_free(&info_used);
    if (!cfg.keep_outfile && sub_rank == 0) unlink(outfile);
    fflush(stdout);

    return err;
}

