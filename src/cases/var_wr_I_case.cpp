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
#define INQ_DIM_LEN(name, size) {                 \
    int dimid;                                    \
    err = driver.inq_dim(ncid, name, &dimid);     \
    CHECK_ERR                                     \
    err = driver.inq_dimlen(ncid, dimid, &size);  \
    CHECK_ERR                                     \
}
#define INQ_PUT_SIZE(size) {                 \
    err = driver.inq_put_size(ncid, &size);  \
    CHECK_ERR                                \
}
#define INQ_FILE_INFO(info) {                \
    err = driver.inq_file_info(ncid, &info); \
    CHECK_ERR                                \
}
#define IPUT_VAR_VTYPE(varid, buf, adv) { \
    err = driver.put_vara(ncid, varid, REC_ITYPE, NULL, NULL, buf, nb); \
    CHECK_VAR_ERR(varid) \
    buf += (adv); \
    my_nreqs++; \
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
#define IPUT_VAR1_VTYPE(varid, buf, adv) { \
    err = driver.put_vara(ncid, varid, REC_ITYPE, start, NULL, buf, nb); \
    CHECK_VAR_ERR(varid) \
    buf += (adv); \
    my_nreqs++; \
}
#define IPUT_VAR1_INT(varid, buf, adv) { \
    err = driver.put_vara(ncid, varid, MPI_INT, start, NULL, buf, nb); \
    CHECK_VAR_ERR(varid) \
    buf += (adv); \
    my_nreqs++; \
}
#define IPUT_VARA_INT_NOADV(buf) { \
    err = driver.put_vara(ncid, varid, MPI_INT, start, count, buf, nb); \
    CHECK_VAR_ERR(varid) \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VARA_INT64_NOADV(buf) { \
    err = driver.put_vara(ncid, varid, MPI_LONG_LONG, start, count, buf, nb); \
    CHECK_VAR_ERR(varid) \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VARA_DBL(varid, buf, adv) { \
    err = driver.put_vara(ncid, varid, MPI_DOUBLE, start, count, buf, nb); \
    CHECK_VAR_ERR(varid) \
    buf += (adv); \
    my_nreqs++; \
}
#define IPUT_VARA_TXT(varid, buf, adv) { \
    err = driver.put_vara(ncid, varid, MPI_CHAR, start, count, buf, nb); \
    CHECK_VAR_ERR(varid) \
    buf += (adv); \
    my_nreqs++; \
}
#define WAIT_ALL_REQS { \
    err = driver.wait(ncid); \
    CHECK_ERR \
    nflushes++; \
}

#define FIX_VAR_IPUT(varid, itype, buf) {                               \
    int vid = varid + num_decomp_vars;                                  \
    int dp = v_decomp_ids[vid];                                         \
    if (cfg.strategy == canonical) {                                    \
        err = driver.put_varn(ncid, vid, itype, decom.contig_nreqs[dp], \
                              decom.w_startx[dp], decom.w_countx[dp],   \
                              buf, nb);                                 \
        my_nreqs += decom.contig_nreqs[dp];                             \
    }                                                                   \
    else {                                                              \
        err = driver.put_vara(ncid, vid, itype, &decom.start[dp],       \
                              &decom.count[dp], buf, nb);               \
        my_nreqs++;                                                     \
    }                                                                   \
    CHECK_VAR_ERR(vid)                                                  \
    buf += decom.count[dp] + gap;                                       \
    nvars_D[dp]++;                                                      \
}

#define REC_VAR_IPUT(varid, itype, buf) {                               \
    int vid = varid + num_decomp_vars;                                  \
    int dp = v_decomp_ids[vid];                                         \
    assert(dp >= 0 && dp < decom.num_decomp);                           \
    if (cfg.strategy == canonical) {                                    \
        err = driver.put_varn(ncid, vid, itype, decom.contig_nreqs[dp], \
                              decom.w_starts[dp], decom.w_counts[dp],   \
                              buf, nb);                                 \
        my_nreqs += decom.contig_nreqs[dp];                             \
    }                                                                   \
    else {                                                              \
        err = driver.put_vara(ncid, vid, itype, starts[dp],             \
                              counts[dp], buf, nb);                     \
        my_nreqs++;                                                     \
    }                                                                   \
    CHECK_VAR_ERR(vid)                                                  \
    buf += decom.count[dp] + gap;                                       \
    if (rec_no == 0) nvars_D[dp]++;                                     \
}

/*----< write_small_fix_vars_F_case() >--------------------------------------*/
static
int write_small_fix_vars_I_case(e3sm_io_driver  &driver,
                                int              ncid,
                                int              start_varid,
                                int              gap,
                                vtype         **fix_buf,
                                int             *nreqs)
{
    /* small fixed-size variable IDs (relative) */
    int varids[5] = {0, 1, 2, 12, 13};

    int i, err, my_nreqs=0;
    MPI_Offset lat, lon, levgrnd, levdcmp, levlak;

    for (i=0; i<18; i++) varids[i] += start_varid;

    INQ_DIM_LEN("lat",     lat)
    INQ_DIM_LEN("lon",     lon)
    INQ_DIM_LEN("levgrnd", levgrnd)
    INQ_DIM_LEN("levdcmp", levdcmp)
    INQ_DIM_LEN("levlak",  levlak)

    /* 5 of type vtype */
    IPUT_VAR_VTYPE(varids[0], *fix_buf,  levgrnd + gap) /* levgrnd */
    IPUT_VAR_VTYPE(varids[1], *fix_buf,  levlak  + gap) /* levlak */
    IPUT_VAR_VTYPE(varids[2], *fix_buf,  levdcmp + gap) /* levdcmp */
    IPUT_VAR_VTYPE(varids[3], *fix_buf,  lon     + gap) /* lon */
    IPUT_VAR_VTYPE(varids[4], *fix_buf,  lat     + gap) /* lat */

    if (nreqs != NULL) *nreqs += my_nreqs;

err_out:
    return err;
}

/*----< write_small_rec_vars_I_case() >--------------------------------------*/
static
int write_small_rec_vars_I_case(e3sm_io_driver  &driver,
                                int              ncid,
                                int              start_varid,
                                int              rec_no,
                                int              gap,
                                char           **txt_buf,
                                int            **int_buf,
                                double         **dbl_buf,
                                vtype          **buf,
                                int             *nreqs)
{
    /* small record variable IDs (relative) */
    int varids[9] = {3, 4, 5, 6, 7, 8, 9, 10, 11};

    int i, err, my_nreqs=0;
    MPI_Offset hist_interval, string_length, start[2], count[2];

    for (i=0; i<9; i++) varids[i] += start_varid;

    INQ_DIM_LEN("hist_interval", hist_interval)
    INQ_DIM_LEN("string_length", string_length)

    start[0] = rec_no;
    start[1] = 0;
    count[0] = 1;

    /* 1 of type vtype */
    IPUT_VAR1_VTYPE(varids[0], *buf,  1 + gap) /* time */

    /* 5 of type int */
    IPUT_VAR1_INT(varids[1], *int_buf, 1 + gap) /* mcdate */
    IPUT_VAR1_INT(varids[2], *int_buf, 1 + gap) /* mcsec */
    IPUT_VAR1_INT(varids[3], *int_buf, 1 + gap) /* mdcur */
    IPUT_VAR1_INT(varids[4], *int_buf, 1 + gap) /* mscur */
    IPUT_VAR1_INT(varids[5], *int_buf, 1 + gap) /* nstep */

    /* 1 of type double */
    count[1] = hist_interval;
    IPUT_VARA_DBL(varids[6], *dbl_buf, hist_interval + gap) /* time_bounds */

    /* 2 of type char */
    count[1] = string_length;
    IPUT_VARA_TXT(varids[7], *txt_buf, string_length + gap) /* date_written */
    IPUT_VARA_TXT(varids[8], *txt_buf, string_length + gap) /* time_written */

    if (nreqs != NULL) *nreqs += my_nreqs;

err_out:
    return err;
}

void wr_buf_init(io_buffers &buf)
{
    buf.fix_int_buflen  = 0;
    buf.fix_buflen      = 0;
    buf.rec_txt_buflen  = 0;
    buf.rec_int_buflen  = 0;
    buf.rec_dbl_buflen  = 0;
    buf.rec_buflen      = 0;

    buf.fix_int_buf  = NULL;
    buf.fix_buf      = NULL;
    buf.rec_txt_buf  = NULL;
    buf.rec_int_buf  = NULL;
    buf.rec_dbl_buf  = NULL;
    buf.rec_buf      = NULL;
}

int wr_buf_malloc(e3sm_io_config &cfg,
                  int             one_flush,
                  io_buffers     &buf)
{
    int rank;
    size_t j;

    MPI_Comm_rank(cfg.io_comm, &rank);

    if (one_flush && cfg.api == pnetcdf) {
        /* write buffers should not be touched when using PnetCDF iput before
         * ncmpi_wait_all is called. For HDF5 and ADIOS blob I/O, write data
         * will be copied and cached into internally allocated buffers and user
         * buffers can be reused after put call returned.
         */
        buf.rec_dbl_buflen *= cfg.nrecs;
        buf.rec_int_buflen *= cfg.nrecs;
        buf.rec_txt_buflen *= cfg.nrecs;
        buf.rec_buflen     *= cfg.nrecs;
    }

    /* allocate and initialize write buffers */
    buf.fix_buf     = (vtype*)  malloc(buf.fix_buflen     * sizeof(vtype));
    buf.fix_int_buf = (int*)    malloc(buf.fix_int_buflen * sizeof(int));
    buf.rec_dbl_buf = (double*) malloc(buf.rec_dbl_buflen * sizeof(double));
    buf.rec_int_buf = (int*)    malloc(buf.rec_int_buflen * sizeof(int));
    buf.rec_txt_buf = (char*)   malloc(buf.rec_txt_buflen * sizeof(char));
    buf.rec_buf     = (vtype*)  malloc(buf.rec_buflen     * sizeof(vtype));

    for (j=0; j<buf.fix_buflen;     j++) buf.fix_buf[j]     = rank;
    for (j=0; j<buf.fix_int_buflen; j++) buf.fix_int_buf[j] = rank;
    for (j=0; j<buf.rec_dbl_buflen; j++) buf.rec_dbl_buf[j] = rank;
    for (j=0; j<buf.rec_int_buflen; j++) buf.rec_int_buf[j] = rank;
    for (j=0; j<buf.rec_txt_buflen; j++) buf.rec_txt_buf[j] = 'a' + rank;
    for (j=0; j<buf.rec_buflen;     j++) buf.rec_buf[j]     = rank;

    return 0;
}

void wr_buf_free(io_buffers &buf)
{
    if (buf.fix_int_buf  != NULL) free(buf.fix_int_buf);
    if (buf.fix_buf      != NULL) free(buf.fix_buf);
    if (buf.rec_txt_buf  != NULL) free(buf.rec_txt_buf);
    if (buf.rec_int_buf  != NULL) free(buf.rec_int_buf);
    if (buf.rec_dbl_buf  != NULL) free(buf.rec_dbl_buf);
    if (buf.rec_buf      != NULL) free(buf.rec_buf);

    buf.fix_int_buf  = NULL;
    buf.fix_buf      = NULL;
    buf.rec_txt_buf  = NULL;
    buf.rec_int_buf  = NULL;
    buf.rec_dbl_buf  = NULL;
    buf.rec_buf      = NULL;

    buf.fix_int_buflen  = 0;
    buf.fix_buflen      = 0;
    buf.rec_txt_buflen  = 0;
    buf.rec_int_buflen  = 0;
    buf.rec_dbl_buflen  = 0;
    buf.rec_buflen      = 0;
}

/*----< var_wr_I_case() >----------------------------------------------------*/
int var_wr_I_case(e3sm_io_config &cfg,
                  e3sm_io_decom  &decom,
                  e3sm_io_driver &driver)
{
    char outfile[1040], base_name[1024], *ext, *rec_txt_buf_ptr;
    const char *hist;
    int i, j, err, sub_rank, global_rank, ncid=-1, nflushes=0, one_flush;
    int rec_no, gap=0, my_nreqs, num_decomp_vars, *v_decomp_ids=NULL;
    int varid, contig_nreqs[MAX_NUM_DECOMP], *nvars_D=cfg.nvars_D;
    int *fix_int_buf_ptr, *rec_int_buf_ptr;
    double *rec_dbl_buf_ptr, timing;
    MPI_Offset previous_size, metadata_size, total_size;
    MPI_Offset blob_start[MAX_NUM_DECOMP], blob_count[MAX_NUM_DECOMP];
    MPI_Offset starts[MAX_NUM_DECOMP][2], counts[MAX_NUM_DECOMP][2];
    MPI_Info info_used = MPI_INFO_NULL;
    MPI_Comm comm;
    vtype  *fix_buf_ptr, *rec_buf_ptr;
    io_buffers wr_buf;

    MPI_Barrier(cfg.io_comm); /*-----------------------------------------*/
    cfg.end2end_time = timing = MPI_Wtime();

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

    /* construct h0/h1 file name */
    hist = (cfg.nvars == 560) ? "_h0" :  "_h1";
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

    /* there are num_decomp_vars number of decomposition variables */
    if (cfg.strategy == blob)
        num_decomp_vars = decom.num_decomp * NVARS_DECOMP;
    else
        num_decomp_vars = 0;

    /* allocate space for all variable IDs */
    v_decomp_ids = (int*) malloc((cfg.nvars + num_decomp_vars) * sizeof(int));
    for (i=0; i<cfg.nvars+num_decomp_vars; i++) v_decomp_ids[i] = -1;

    wr_buf_init(wr_buf);
 
    /* define dimensions, variables, and attributes.
     * 560 climate variables
     *     18 are fixed-size,      542 are record variables
     *     14 are not partitioned, 546 are partitioned.
     * Thus, 14 not partitioned variable (also small) are written by root only:
     *     5 fixed-size, 9 record variables
     */

    /* h0 file contains 560 variablesi, while h1 has 8 less */
    err = def_I_case(cfg, decom, driver, ncid, v_decomp_ids, &wr_buf);
    CHECK_ERR

    /* exit define mode and enter data mode */
    ENDDEF

    /* I/O amount so far */
    INQ_PUT_SIZE(metadata_size)
    if (cfg.api != pnetcdf) metadata_size -= previous_size;

    INQ_FILE_INFO(info_used)

    cfg.def_time = MPI_Wtime() - timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    if (cfg.non_contig_buf) gap = 10;

    /* allocate write buffers */
    wr_buf_malloc(cfg, one_flush, wr_buf);

    for (j=0; j<decom.num_decomp; j++)
        nvars_D[j] = 0; /* number of variables using decomposition j */

    cfg.pre_time = MPI_Wtime() - timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    my_nreqs = 0; /* number of noncontiguous requests written by this proc */

    if (cfg.strategy == blob) {
        varid=0;
        /* write decomposition variables, they are defined first in file */
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

    /* Now, write climate variables, fixed-size first */
    fix_buf_ptr     = wr_buf.fix_buf;
    fix_int_buf_ptr = wr_buf.fix_int_buf;

    if (sub_rank == 0) {
        /* write 5 small fixed-size variables */
        err = write_small_fix_vars_I_case(driver, ncid, num_decomp_vars, gap,
                                          &fix_buf_ptr, &my_nreqs);
        CHECK_ERR
    }

    /* area, topo, landfrac */
    for (varid=14; varid<17; varid++)
        FIX_VAR_IPUT(varid, REC_ITYPE, fix_buf_ptr)

    /* landmask, pftmask */
    for (varid=17; varid<19; varid++)
        FIX_VAR_IPUT(varid, MPI_INT, fix_int_buf_ptr)

    if (cfg.nvars == 560) {
        /* ZSOI, DZSOI, WATSAT, SUCSAT, BSW, HKSAT, ZLAKE, DZLAKE */
        for (varid=19; varid<27; varid++)
            FIX_VAR_IPUT(varid, REC_ITYPE, fix_buf_ptr)
    }

    /* now write record variables */
    rec_txt_buf_ptr = wr_buf.rec_txt_buf;
    rec_int_buf_ptr = wr_buf.rec_int_buf;
    rec_dbl_buf_ptr = wr_buf.rec_dbl_buf;
    rec_buf_ptr     = wr_buf.rec_buf;

    for (rec_no=0; rec_no<cfg.nrecs; rec_no++) {
        if (!one_flush || (cfg.strategy == blob && cfg.api == hdf5)) {
            /* reset the pointers to the beginning of the buffers */
            rec_txt_buf_ptr = wr_buf.rec_txt_buf;
            rec_int_buf_ptr = wr_buf.rec_int_buf;
            rec_dbl_buf_ptr = wr_buf.rec_dbl_buf;
            rec_buf_ptr     = wr_buf.rec_buf;
        }

        if (sub_rank == 0) {
            /* write 15 small record variables by rank 0 only */
            err = write_small_rec_vars_I_case(driver, ncid, num_decomp_vars,
                                              rec_no, gap,
                                              &rec_txt_buf_ptr,
                                              &rec_int_buf_ptr,
                                              &rec_dbl_buf_ptr,
                                              &rec_buf_ptr,
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

        if (cfg.nvars == 560) { /* h0 file */
            /* first partitioned record variable is ACTUAL_IMMOB, id = 27 */
            for (varid=27; varid<cfg.nvars; varid++)
                REC_VAR_IPUT(varid, REC_ITYPE, rec_buf_ptr)

        } else { /* h1 file */
            /* first partitioned record variable is ACTUAL_IMMOB, id = 19 */
            for (varid=19; varid<cfg.nvars; varid++)
                REC_VAR_IPUT(varid, REC_ITYPE, rec_buf_ptr)
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

    wr_buf_free(wr_buf);
    if (v_decomp_ids != NULL) free(v_decomp_ids);

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

