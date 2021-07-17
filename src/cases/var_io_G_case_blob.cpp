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
#include <stdlib.h> /* strtoll() */
#include <string.h> /* strcpy(), strncpy() */
#include <string>

#include <unistd.h> /* getopt() unlink() */

#include <mpi.h>

#include <e3sm_io_case_G.hpp>
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
#define IPUT_VAR_DBL(adv) { \
    err = ncmpi_iput_var_double(ncid, *varid, dbl_buf_ptr, NULL); \
    CHECK_VAR_ERR(*varid) \
    dbl_buf_ptr += (adv); \
    varid++; \
}
#define IPUT_VAR_INT(adv) { \
    err = ncmpi_iput_var_int(ncid, *varid, int_buf_ptr, NULL); \
    CHECK_VAR_ERR(*varid) \
    int_buf_ptr += (adv); \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VAR1_DBL(adv) { \
    err = ncmpi_iput_var1_double(ncid, *varid, start, dbl_buf_ptr, NULL); \
    CHECK_VAR_ERR(*varid) \
    dbl_buf_ptr += (adv); \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VAR1_INT(adv) { \
    err = ncmpi_iput_var1_int(ncid, *varid, start, int_buf_ptr, NULL); \
    CHECK_VAR_ERR(*varid) \
    int_buf_ptr += (adv); \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VARA_INT(adv) { \
    err = ncmpi_iput_vara_int(ncid, *varid, start, count, int_buf_ptr, NULL); \
    CHECK_VAR_ERR(*varid) \
    int_buf_ptr += (adv); \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VARA_INT_NOADV(buf) { \
    err = ncmpi_iput_vara_int(ncid, *varid, start, count, buf, NULL); \
    CHECK_VAR_ERR(*varid) \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VARA_INT64_NOADV(buf) { \
    err = ncmpi_iput_vara_longlong(ncid, *varid, start, count, buf, NULL); \
    CHECK_VAR_ERR(*varid) \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VARA_DBL(adv) { \
    err = ncmpi_iput_vara_double(ncid, *varid, start, count, dbl_buf_ptr, NULL); \
    CHECK_VAR_ERR(*varid) \
    dbl_buf_ptr += (adv); \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VARA_CHR(adv) { \
    err = ncmpi_iput_vara_text(ncid, *varid, start, count, chr_buf_ptr, NULL); \
    CHECK_VAR_ERR(*varid) \
    chr_buf_ptr += (adv); \
    my_nreqs++; \
    varid++; \
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
#define POST_VARA_DBL(decomp) {                                       \
    err = ncmpi_iput_vara_double(ncid, *varid, start_D[decomp],       \
                                 count_D[decomp], dbl_buf_ptr, NULL); \
    CHECK_VAR_ERR(*varid)                                             \
    dbl_buf_ptr += count_D[decomp][1] + gap;                          \
    my_nreqs++;                                                       \
    varid++;                                                          \
    if (rec_no == 0) nvars_D[decomp]++;                               \
}
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
#define IPUT_VAR_DBL(adv) { \
    err = driver.put_vara(ncid, *varid, MPI_DOUBLE, NULL, NULL, dbl_buf_ptr, nb); \
    CHECK_VAR_ERR(*varid) \
    dbl_buf_ptr += (adv); \
    varid++; \
}
#define IPUT_VAR_INT(adv) { \
    err = driver.put_vara(ncid, *varid, MPI_INT, NULL, NULL, int_buf_ptr, nb); \
    CHECK_VAR_ERR(*varid) \
    int_buf_ptr += (adv); \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VAR1_DBL(adv) { \
    err = driver.put_vara(ncid, *varid, MPI_DOUBLE, start, NULL, dbl_buf_ptr, nb); \
    CHECK_VAR_ERR(*varid) \
    dbl_buf_ptr += (adv); \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VAR1_INT(adv) { \
    err = driver.put_vara(ncid, *varid, MPI_INT, start, NULL, int_buf_ptr, nb); \
    CHECK_VAR_ERR(*varid) \
    int_buf_ptr += (adv); \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VARA_INT(adv) { \
    err = driver.put_vara(ncid, *varid, MPI_INT, start, count, int_buf_ptr, nb); \
    CHECK_VAR_ERR(*varid) \
    int_buf_ptr += (adv); \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VARA_INT_NOADV(buf) { \
    err = driver.put_vara(ncid, *varid, MPI_INT, start, count, buf, nb); \
    CHECK_VAR_ERR(*varid) \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VARA_INT64_NOADV(buf) { \
    err = driver.put_vara(ncid, *varid, MPI_LONG_LONG, start, count, buf, nb); \
    CHECK_VAR_ERR(*varid) \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VARA_DBL(adv) { \
    err = driver.put_vara(ncid, *varid, MPI_DOUBLE, start, count, dbl_buf_ptr, nb); \
    CHECK_VAR_ERR(*varid) \
    dbl_buf_ptr += (adv); \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VARA_CHR(adv) { \
    err = driver.put_vara(ncid, *varid, MPI_CHAR, start, count, chr_buf_ptr, nb); \
    CHECK_VAR_ERR(*varid) \
    chr_buf_ptr += (adv); \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VARA(varid, start, count, buf) { \
    err = driver.put_vara(ncid, varid, REC_ITYPE, start, count, buf, nb); \
    CHECK_VAR_ERR(varid) \
    my_nreqs++; \
    varid++; \
}
#define POST_VARA_DBL(decomp) {                                       \
    err = driver.put_vara(ncid, *varid, MPI_DOUBLE, start_D[decomp],  \
                          count_D[decomp], dbl_buf_ptr, nb);          \
    CHECK_VAR_ERR(*varid)                                             \
    dbl_buf_ptr += count_D[decomp][1] + gap;                          \
    my_nreqs++;                                                       \
    varid++;                                                          \
    if (rec_no == 0) nvars_D[decomp]++;                               \
}
#define WAIT_ALL_REQS { \
    err = driver.wait(ncid); \
    CHECK_ERR \
    nflushes++; \
}
#endif

/*----< blob_G_case() >------------------------------------------------------*/
int blob_G_case(e3sm_io_config &cfg,
                e3sm_io_decom  &decom,
                e3sm_io_driver &driver)
{
    char outfile[1040];
    int i, err, sub_rank, global_rank, ncid=-1, nflushes=0, *varids;
    int rec_no, gap=0, my_nreqs, num_decomp_vars;
    int contig_nreqs[MAX_NUM_DECOMP], nvars_D[MAX_NUM_DECOMP];
    double pre_timing, open_timing, post_timing, wait_timing, close_timing;
    double def_timing, timing, total_timing;
    MPI_Offset metadata_size, put_size, total_size, max_nreqs, total_nreqs;
    MPI_Offset start[2], count[2];
    MPI_Offset start_D[MAX_NUM_DECOMP][2], count_D[MAX_NUM_DECOMP][2];
    MPI_Offset blob_start[MAX_NUM_DECOMP], blob_count[MAX_NUM_DECOMP];
    MPI_Offset previous_size, sum_size, m_alloc=0, max_alloc;
    MPI_Info info_used = MPI_INFO_NULL;
    size_t j, chr_buflen, int_buflen, dbl_buflen;
    double *dbl_buf=NULL, *dbl_buf_ptr; /* buffer for fixed double var */
    int    *int_buf=NULL, *int_buf_ptr; /* buffer for int var */
    char   *chr_buf=NULL, *chr_buf_ptr; /* buffer for char var */

    int dec_varids[NVARS_DECOMP*MAX_NUM_DECOMP], *varid;
    int fix_varids[11] = {8, 9, 10, 11, 12, 13, 14, 33, 35, 36, 37};
    int rec_varids[41] = { 0,  1,  2,  3,  4,  5,  6,  7,
                                              15, 16, 17, 18, 19,
                          20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                          30, 31, 32,     34,             38, 39,
                          40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
                          50, 51};

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    total_timing = pre_timing = MPI_Wtime ();

    open_timing  = 0.0;
    def_timing   = 0.0;
    post_timing  = 0.0;
    wait_timing  = 0.0;
    close_timing = 0.0;

    if (cfg.api == hdf5) /* I/O amount from previous I/O */
        INQ_PUT_SIZE(previous_size)
    else
        previous_size = 0;

    MPI_Comm_rank(cfg.io_comm,  &global_rank);
    MPI_Comm_rank(cfg.sub_comm, &sub_rank);

    /* there are num_decomp_vars number of decomposition variables */
    num_decomp_vars = decom.num_decomp * NVARS_DECOMP;

    for (i=0; i<num_decomp_vars; i++) dec_varids[i] = i;
    for (i=0; i<11; i++) fix_varids[i] += num_decomp_vars;
    for (i=0; i<41; i++) rec_varids[i] += num_decomp_vars;

    if (cfg.non_contig_buf) gap = 10;

    /* allocate write buffer for small climate variables */
    dbl_buflen = decom.dims[2][1] * 4
               + decom.count[0] * 5
               + decom.count[2] * 24
               + decom.count[3]
               + decom.count[4]
               + decom.count[5] * 4
               + 6                     /* 6 single-element variables */
               + 45 * gap;

    int_buflen = decom.count[0]
               + decom.count[1] * 2
               + decom.count[2]
               + decom.count[3]
               + decom.count[4]
               + 6 * gap;

    chr_buflen = 64 + gap;

    /* allocate and initialize write buffer for large variables */
#define FLUSH_ALL_RECORDS_AT_ONCE
#ifdef FLUSH_ALL_RECORDS_AT_ONCE
    dbl_buflen *= cfg.nrec;
    int_buflen *= cfg.nrec;
    chr_buflen *= cfg.nrec;
#else
    if (cfg.api == hdf5) {
        printf("Error in %s:%d: %s() FLUSH_ALL_RECORDS_AT_ONCE must be enabled when using HDF5 blob I/O",
               __FILE__, __LINE__, __func__);
        err = -1;
        goto err_out;
    }
#endif

    dbl_buf = (double*) malloc(dbl_buflen * sizeof(double));
    int_buf = (int*)    malloc(int_buflen * sizeof(int));
    chr_buf = (char*)   malloc(chr_buflen * sizeof(char));

    for (j=0; j<dbl_buflen; j++) dbl_buf[j] = global_rank;
    for (j=0; j<int_buflen; j++) int_buf[j] = global_rank;
    for (j=0; j<chr_buflen; j++) chr_buf[j] = 'a' + global_rank;

    pre_timing = MPI_Wtime() - pre_timing;

    MPI_Barrier(cfg.sub_comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    sprintf(outfile, "%s.%04d", cfg.out_path, cfg.subfile_ID);

    /* set output subfile name */
    if (cfg.verbose && sub_rank == 0)
        printf("global_rank=%d sub_rank=%d outfile=%s\n",global_rank,sub_rank,outfile);

    /* create the output file */
    FILE_CREATE(outfile)

    open_timing += MPI_Wtime() - timing;

    MPI_Barrier(cfg.sub_comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    /* allocate space for all variable IDs */
    varids = (int*) malloc((cfg.nvars + num_decomp_vars) * sizeof(int));

    /* define dimensions, variables, and attributes */
    err = def_G_case(cfg, decom, driver, ncid, varids);
    CHECK_ERR

    free(varids);

    /* exit define mode and enter data mode */
    ENDDEF

    /* I/O amount so far */
    INQ_PUT_SIZE(metadata_size)
    metadata_size -= previous_size;

    INQ_FILE_INFO(info_used)

    def_timing += MPI_Wtime() - timing;

    MPI_Barrier(cfg.sub_comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    my_nreqs = 0;

    /* First, write decomposition variables */
    varid = dec_varids;
    for (i=0; i<decom.num_decomp; i++) {
        nvars_D[i] = 0; /* number of variables using decomposition i */
        start[0] = sub_rank;
        count[0] = 1;
        start[1] = 0;
        count[1] = decom.contig_nreqs[i];

        /* write to D*.nreqs, 1D array */
        contig_nreqs[i] = decom.contig_nreqs[i];
        IPUT_VARA_INT_NOADV(contig_nreqs+i)

        /* write to D*.blob_start, 1D array */
        blob_start[i] = decom.start[i];
        IPUT_VARA_INT64_NOADV(blob_start+i)

        /* write to D*.blob_count, 1D array */
        blob_count[i] = decom.count[i];
        IPUT_VARA_INT64_NOADV(blob_count+i);

        /* write to D*.offsets, 2D array */
        IPUT_VARA_INT_NOADV(decom.disps[i])

        /* write to D*.lengths, 2D array */
        IPUT_VARA_INT_NOADV(decom.blocklens[i])
    }

    /* Now, write climate variables */
    dbl_buf_ptr = dbl_buf;
    int_buf_ptr = int_buf;
    chr_buf_ptr = chr_buf;

    /* write all 11 fixed-size variables first */
    varid = fix_varids;

    /* int maxLevelEdgeTop(nEdges) */
    start[0] = decom.start[1];
    count[0] = decom.count[1];
    IPUT_VARA_INT(count[0] + gap)
    nvars_D[1]++;

    /* double vertCoordMovementWeights(nVertLevels) */
    if (sub_rank == 0)
        IPUT_VAR_DBL(decom.dims[2][1] + gap)
    else
        varid++;

    /* int edgeMask(nEdges, nVertLevels) */
    start[0] = decom.start[3];
    count[0] = decom.count[3];
    IPUT_VARA_INT(count[0] + gap)
    nvars_D[3]++;

    /* int cellMask(nCells, nVertLevels) */
    start[0] = decom.start[2];
    count[0] = decom.count[2];
    IPUT_VARA_INT(count[0] + gap)
    nvars_D[2]++;

    /* int vertexMask(nVertices, nVertLevels) */
    start[0] = decom.start[4];
    count[0] = decom.count[4];
    IPUT_VARA_INT(count[0] + gap)
    nvars_D[4]++;

    /* double refZMid(nVertLevels) */
    if (sub_rank == 0)
        IPUT_VAR_DBL(decom.dims[2][1] + gap)
    else
        varid++;

    /* double refLayerThickness(nVertLevels) */
    if (sub_rank == 0)
        IPUT_VAR_DBL(decom.dims[2][1] + gap)
    else
        varid++;

    /* double refBottomDepth(nVertLevels) */
    if (sub_rank == 0)
        IPUT_VAR_DBL(decom.dims[2][1] + gap)
    else
        varid++;

    /* double bottomDepth(nCells) */
    start[0] = decom.start[0];
    count[0] = decom.count[0];
    IPUT_VARA_DBL(count[0] + gap)
    nvars_D[0]++;

    /* int maxLevelCell(nCells) */
    start[0] = decom.start[0];
    count[0] = decom.count[0];
    IPUT_VARA_INT(count[0] + gap)
    nvars_D[0]++;

    /* int maxLevelEdgeBot(nEdges) */
    start[0] = decom.start[1];
    count[0] = decom.count[1];
    IPUT_VARA_INT(count[0] + gap)
    nvars_D[1]++;

    /* Now write record variables */

    for (i=0; i<decom.num_decomp; i++) {
        start_D[i][0] = 0;
        start_D[i][1] = decom.start[i];
        count_D[i][0] = 1;
        count_D[i][1] = decom.count[i];
    }

    for (rec_no=0; rec_no<cfg.nrec; rec_no++) {
#ifndef FLUSH_ALL_RECORDS_AT_ONCE
        timing = MPI_Wtime();
        dbl_buf_ptr = dbl_buf
                    + decom.dims[2][1] * 4
                    + decom.count[0] + gap * 5;
        int_buf_ptr = int_buf
                    + decom.count[0]
                    + decom.count[1] * 2
                    + decom.count[3]
                    + decom.count[2]
                    + decom.count[4] + gap * 6;
        chr_buf_ptr = chr_buf;
#endif

        for (i=0; i<decom.num_decomp; i++)
            start_D[i][0] = rec_no;

        varid = rec_varids;

        /* write 41 record variables */

        /* double salinitySurfaceRestoringTendency(Time, nCells) */
        POST_VARA_DBL(0)
        /* double vertTransportVelocityTop(Time, nCells, nVertLevelsP1) */
        POST_VARA_DBL(5)
        /* double vertGMBolusVelocityTop(Time, nCells, nVertLevelsP1) */
        POST_VARA_DBL(5)
        /* double vertAleTransportTop(Time, nCells, nVertLevelsP1) */
        POST_VARA_DBL(5)
        /* double tendSSH(Time, nCells) */
        POST_VARA_DBL(0)
        /* double layerThickness(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double normalVelocity(Time, nEdges, nVertLevels) */
        POST_VARA_DBL(3)
        /* double ssh(Time, nCells) */
        POST_VARA_DBL(0)
        /* char xtime(Time, StrLen) */
        if (sub_rank == 0) {
            start[0] = rec_no;
            start[1] = 0;
            count[0] = 1;
            count[1] = 64;
            IPUT_VARA_CHR(64 + gap)
        }
        else
            varid++;
        /* double kineticEnergyCell(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double relativeVorticityCell(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double relativeVorticity(Time, nVertices, nVertLevels) */
        POST_VARA_DBL(4)
        /* double divergence(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        if (sub_rank == 0) {
            start[0] = rec_no;
            count[0] = 1;
            /* double areaCellGlobal(Time) */
            IPUT_VAR1_DBL(1 + gap)
            /* double areaEdgeGlobal(Time) */
            IPUT_VAR1_DBL(1 + gap)
            /* double areaTriangleGlobal(Time) */
            IPUT_VAR1_DBL(1 + gap)
            /* double volumeCellGlobal(Time) */
            IPUT_VAR1_DBL(1 + gap)
            /* double volumeEdgeGlobal(Time) */
            IPUT_VAR1_DBL(1 + gap)
            /* double CFLNumberGlobal(Time) */
            IPUT_VAR1_DBL(1 + gap)
        }
        else
            varid += 6;
        /* double BruntVaisalaFreqTop(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double vertVelocityTop(Time, nCells, nVertLevelsP1) */
        POST_VARA_DBL(5)
        /* double velocityZonal(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double velocityMeridional(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double displacedDensity(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double potentialDensity(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double pressure(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double zMid(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double columnIntegratedSpeed(Time, nCells) */
        POST_VARA_DBL(0)
        /* double temperatureHorizontalAdvectionTendency(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double salinityHorizontalAdvectionTendency(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double temperatureVerticalAdvectionTendency(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double salinityVerticalAdvectionTendency(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double temperatureVertMixTendency(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double salinityVertMixTendency(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double temperatureSurfaceFluxTendency(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double salinitySurfaceFluxTendency(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double temperatureShortWaveTendency(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double temperatureNonLocalTendency(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double salinityNonLocalTendency(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double temperature(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)
        /* double salinity(Time, nCells, nVertLevels) */
        POST_VARA_DBL(2)

#ifndef FLUSH_ALL_RECORDS_AT_ONCE
        post_timing += MPI_Wtime() - timing;

        MPI_Barrier(cfg.sub_comm); /*---------------------------------------*/
        timing = MPI_Wtime();

        /* flush once per time record */
        WAIT_ALL_REQS

        wait_timing += MPI_Wtime() - timing;

        timing = MPI_Wtime();
#endif
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

    if (dbl_buf != NULL) free(dbl_buf);
    if (chr_buf != NULL) free(chr_buf);
    if (int_buf != NULL) free(int_buf);

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
        for (i=0; i<decom.num_decomp; i++) nvars_noD -= nvars_D[i];
        printf("History output file                = %s\n", outfile);
        printf("No. decomposition variables        = %3d\n", num_decomp_vars);
        printf("No. variables use no decomposition = %3d\n", nvars_noD);
        printf("No. variables use decomposition D1 = %3d\n", nvars_D[0]);
        printf("No. variables use decomposition D2 = %3d\n", nvars_D[1]);
        printf("No. variables use decomposition D3 = %3d\n", nvars_D[2]);
        printf("No. variables use decomposition D4 = %3d\n", nvars_D[3]);
        printf("No. variables use decomposition D5 = %3d\n", nvars_D[4]);
        printf("No. variables use decomposition D6 = %3d\n", nvars_D[5]);
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

