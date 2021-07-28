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
#include <assert.h>
#include <unistd.h> /* getopt() unlink() */

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
#define IPUT_VAR_DBL(varid, adv) { \
    err = driver.put_vara(ncid, varid, MPI_DOUBLE, NULL, NULL, \
                          dbl_buf_ptr, nb); \
    CHECK_VAR_ERR(varid) \
    dbl_buf_ptr += (adv); \
    my_nreqs++; \
}
#define IPUT_VAR1_DBL(varid, adv) { \
    err = driver.put_vara(ncid, varid, MPI_DOUBLE, start, NULL, \
                          dbl_buf_ptr, nb); \
    CHECK_VAR_ERR(varid) \
    dbl_buf_ptr += (adv); \
    my_nreqs++; \
}
#define FIX_IPUT_INT(dp, varid) {                                         \
    if (cfg.strategy == canonical) {                                      \
        err = driver.put_varn(ncid, varid, MPI_INT,                       \
                              decom.contig_nreqs[dp],                     \
                              decom.w_startx[dp], decom.w_countx[dp],     \
                              int_buf_ptr, nb);                           \
        my_nreqs += decom.contig_nreqs[dp];                               \
    }                                                                     \
    else {                                                                \
        err = driver.put_vara(ncid, varid, MPI_INT, &decom.start[dp],     \
                              &decom.count[dp], int_buf_ptr, nb);         \
        my_nreqs++;                                                       \
    }                                                                     \
    CHECK_VAR_ERR(varid)                                                  \
    int_buf_ptr += decom.count[dp] + gap;                                 \
    nvars_D[dp]++;                                                        \
}

#define FIX_IPUT_DBL(dp, varid) {                                         \
    if (cfg.strategy == canonical) {                                      \
        err = driver.put_varn(ncid, varid, MPI_DOUBLE,                    \
                              decom.contig_nreqs[dp],                     \
                              decom.w_startx[dp], decom.w_countx[dp],     \
                              dbl_buf_ptr, nb);                           \
        my_nreqs += decom.contig_nreqs[dp];                               \
    }                                                                     \
    else {                                                                \
        err = driver.put_vara(ncid, varid, MPI_DOUBLE, &decom.start[dp],  \
                              &decom.count[dp], dbl_buf_ptr, nb);         \
        my_nreqs++;                                                       \
    }                                                                     \
    CHECK_VAR_ERR(varid)                                                  \
    dbl_buf_ptr += decom.count[dp] + gap;                                 \
    nvars_D[dp]++;                                                        \
}

#define IPUT_VARA_INT_NOADV(varid, buf) { \
    err = driver.put_vara(ncid, varid, MPI_INT, start, count, buf, nb); \
    CHECK_VAR_ERR(varid) \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VARA_INT64_NOADV(varid, buf) { \
    err = driver.put_vara(ncid, varid, MPI_LONG_LONG, start, count, buf, nb); \
    CHECK_VAR_ERR(varid) \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VARA_CHR(varid, adv) { \
    err = driver.put_vara(ncid, varid, MPI_CHAR, start, count, \
                          chr_buf_ptr, nb); \
    CHECK_VAR_ERR(varid) \
    chr_buf_ptr += (adv); \
    my_nreqs++; \
}
#define WAIT_ALL_REQS { \
    err = driver.wait(ncid); \
    CHECK_ERR \
    nflushes++; \
}

#define REC_IPUT(dp, varid) {                                         \
    if (cfg.strategy == blob) {                                       \
        err = driver.put_vara(ncid, varid, MPI_DOUBLE, starts[dp],    \
                              counts[dp], dbl_buf_ptr, nb);           \
        CHECK_VAR_ERR(varid)                                          \
        dbl_buf_ptr += decom.count[dp] + gap;                         \
        my_nreqs++;                                                   \
    }                                                                 \
    else {                                                            \
        err = driver.put_varn(ncid, varid, MPI_DOUBLE,                \
                              decom.contig_nreqs[dp],                 \
                              decom.w_starts[dp], decom.w_counts[dp], \
                              dbl_buf_ptr, nb);                       \
        CHECK_VAR_ERR(varid)                                          \
        dbl_buf_ptr += decom.count[dp] + gap;                         \
        my_nreqs += decom.contig_nreqs[dp];                           \
    }                                                                 \
    if (rec_no == 0) nvars_D[dp]++;                                   \
}

/*----< var_wr_G_case() >----------------------------------------------------*/
int var_wr_G_case(e3sm_io_config &cfg,
                  e3sm_io_decom  &decom,
                  e3sm_io_driver &driver)
{
    char outfile[1040];
    int i, j, err, sub_rank, global_rank, ncid=-1, nflushes=0, *varids;
    int rec_no, gap=0, my_nreqs, num_decomp_vars, one_flush;
    int contig_nreqs[MAX_NUM_DECOMP], *nvars_D=cfg.nvars_D;
    double timing;
    MPI_Offset metadata_size, total_size;
    MPI_Offset start[2], count[2];
    MPI_Offset starts[MAX_NUM_DECOMP][2], counts[MAX_NUM_DECOMP][2];
    MPI_Offset blob_start[MAX_NUM_DECOMP], blob_count[MAX_NUM_DECOMP];
    MPI_Offset previous_size;
    MPI_Info info_used = MPI_INFO_NULL;
    MPI_Comm comm;

    size_t ii, chr_buflen, int_buflen, dbl_buflen;
    double *dbl_buf=NULL, *dbl_buf_ptr; /* buffer for fixed double var */
    int    *int_buf=NULL, *int_buf_ptr; /* buffer for int var */
    char   *chr_buf=NULL, *chr_buf_ptr; /* buffer for char var */

    int fix_varids[11] = {8, 9, 10, 11, 12, 13, 14, 33, 35, 36, 37};
    int rec_varids[41] = { 0,  1,  2,  3,  4,  5,  6,  7,
                                              15, 16, 17, 18, 19,
                          20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                          30, 31, 32,     34,             38, 39,
                          40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
                          50, 51};

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
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

    /* allocate write buffer for climate variables */

    /* 45 variables of type double */
    dbl_buflen = decom.dims[2][1] * 4
               + decom.count[0] * 5
               + decom.count[2] * 24
               + decom.count[3]
               + decom.count[4]
               + decom.count[5] * 4
               + 6                     /* 6 single-element variables */
               + 45 * gap;

    /* 6 variables of type int */
    int_buflen = decom.count[0]
               + decom.count[1] * 2
               + decom.count[2]
               + decom.count[3]
               + decom.count[4]
               + 6 * gap;

    /* 1 variable of type char */
    chr_buflen = 64 + gap;

    /* allocate and initialize write buffer for large variables */
    if (one_flush && cfg.api == pnetcdf) {
        /* write buffers should not be touched when using PnetCDF iput before
         * ncmpi_wait_all is called. For HDF5 and ADIOS blob I/O, write data
         * will be copied and cached into internally allocated buffers and user
         * buffers can be reused after put call returned.
         */
        dbl_buflen *= cfg.nrecs;
        int_buflen *= cfg.nrecs;
        chr_buflen *= cfg.nrecs;
    }

    dbl_buf = (double*) malloc(dbl_buflen * sizeof(double));
    int_buf = (int*)    malloc(int_buflen * sizeof(int));
    chr_buf = (char*)   malloc(chr_buflen * sizeof(char));

    for (ii=0; ii<dbl_buflen; ii++) dbl_buf[ii] = global_rank;
    for (ii=0; ii<int_buflen; ii++) int_buf[ii] = global_rank;
    for (ii=0; ii<chr_buflen; ii++) chr_buf[ii] = 'a' + global_rank;

    /* there are num_decomp_vars number of decomposition variables */
    if (cfg.strategy == blob) {
        num_decomp_vars = decom.num_decomp * NVARS_DECOMP;
        for (i=0; i<11; i++) fix_varids[i] += num_decomp_vars;
        for (i=0; i<41; i++) rec_varids[i] += num_decomp_vars;
    }
    else
        num_decomp_vars = 0;

    /* allocate space for all variable IDs */
    varids = (int*) malloc((cfg.nvars + num_decomp_vars) * sizeof(int));

    cfg.pre_time = MPI_Wtime() - cfg.pre_time;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    if (cfg.strategy == blob) /* set output subfile name */
        sprintf(outfile, "%s.%04d", cfg.out_path, cfg.subfile_ID);
    else
        strcpy(outfile, cfg.out_path);

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
    err = def_G_case(cfg, decom, driver, ncid, varids);
    CHECK_ERR

    free(varids); /* no longer need varids[] */

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

    my_nreqs = 0;

    if (cfg.strategy == blob) {
        /* First, write decomposition variables */
        int varid=0;
        for (j=0; j<decom.num_decomp; j++) {
            start[0] = sub_rank;
            count[0] = 1;
            start[1] = 0;
            count[1] = decom.contig_nreqs[j];

            /* write to D*.nreqs, 1D array */
            contig_nreqs[j] = decom.contig_nreqs[j];
            IPUT_VARA_INT_NOADV(varid, contig_nreqs+j)

            /* write to D*.blob_start, 1D array */
            blob_start[j] = decom.start[j];
            IPUT_VARA_INT64_NOADV(varid, blob_start+j)

            /* write to D*.blob_count, 1D array */
            blob_count[j] = decom.count[j];
            IPUT_VARA_INT64_NOADV(varid, blob_count+j);

            /* write to D*.offsets, 2D array */
            IPUT_VARA_INT_NOADV(varid, decom.disps[j])

            /* write to D*.lengths, 2D array */
            IPUT_VARA_INT_NOADV(varid, decom.blocklens[j])

            starts[j][0] = 0;
            starts[j][1] = decom.start[j];
            counts[j][0] = 1;
            counts[j][1] = decom.count[j];
        }
        assert(varid == num_decomp_vars);
    }

    /* Now, write climate variables */
    dbl_buf_ptr = dbl_buf;
    int_buf_ptr = int_buf;
    chr_buf_ptr = chr_buf;

    /* write all 11 fixed-size variables first */

    /* int maxLevelEdgeTop(nEdges) */
    FIX_IPUT_INT(1, fix_varids[0])

    /* double vertCoordMovementWeights(nVertLevels) */
    if (sub_rank == 0)
        IPUT_VAR_DBL(fix_varids[1], decom.dims[2][1] + gap)

    /* int edgeMask(nEdges, nVertLevels) */
    FIX_IPUT_INT(3, fix_varids[2])

    /* int cellMask(nCells, nVertLevels) */
    FIX_IPUT_INT(2, fix_varids[3])

    /* int vertexMask(nVertices, nVertLevels) */
    FIX_IPUT_INT(4, fix_varids[4])

    if (sub_rank == 0) {
        /* double refZMid(nVertLevels) */
        IPUT_VAR_DBL(fix_varids[5], decom.dims[2][1] + gap)

        /* double refLayerThickness(nVertLevels) */
        IPUT_VAR_DBL(fix_varids[6], decom.dims[2][1] + gap)

        /* double refBottomDepth(nVertLevels) */
        IPUT_VAR_DBL(fix_varids[7], decom.dims[2][1] + gap)
    }

    /* double bottomDepth(nCells) */
    FIX_IPUT_DBL(0, fix_varids[8])

    /* int maxLevelCell(nCells) */
    FIX_IPUT_INT(0, fix_varids[9])

    /* int maxLevelEdgeBot(nEdges) */
    FIX_IPUT_INT(1, fix_varids[10])

    /* Now write record variables */

    for (rec_no=0; rec_no<cfg.nrecs; rec_no++) {
        if (!one_flush || (cfg.strategy == blob && cfg.api == hdf5)) {
            /* reset the pointers to the beginning of the buffers */
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

        /* write 41 record variables */

        /* double salinitySurfaceRestoringTendency(Time, nCells) */
        REC_IPUT(0, rec_varids[0])
        /* double vertTransportVelocityTop(Time, nCells, nVertLevelsP1) */
        REC_IPUT(5, rec_varids[1])
        /* double vertGMBolusVelocityTop(Time, nCells, nVertLevelsP1) */
        REC_IPUT(5, rec_varids[2])
        /* double vertAleTransportTop(Time, nCells, nVertLevelsP1) */
        REC_IPUT(5, rec_varids[3])
        /* double tendSSH(Time, nCells) */
        REC_IPUT(0, rec_varids[4])
        /* double layerThickness(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[5])
        /* double normalVelocity(Time, nEdges, nVertLevels) */
        REC_IPUT(3, rec_varids[6])
        /* double ssh(Time, nCells) */
        REC_IPUT(0, rec_varids[7])
        /* char xtime(Time, StrLen) */
        if (sub_rank == 0) {
            start[0] = rec_no;
            start[1] = 0;
            count[0] = 1;
            count[1] = 64;
            IPUT_VARA_CHR(rec_varids[8], 64 + gap)
        }
        /* double kineticEnergyCell(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[9])
        /* double relativeVorticityCell(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[10])
        /* double relativeVorticity(Time, nVertices, nVertLevels) */
        REC_IPUT(4, rec_varids[11])
        /* double divergence(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[12])
        if (sub_rank == 0) {
            start[0] = rec_no;
            /* double areaCellGlobal(Time) */
            IPUT_VAR1_DBL(rec_varids[13], 1 + gap)
            /* double areaEdgeGlobal(Time) */
            IPUT_VAR1_DBL(rec_varids[14], 1 + gap)
            /* double areaTriangleGlobal(Time) */
            IPUT_VAR1_DBL(rec_varids[15], 1 + gap)
            /* double volumeCellGlobal(Time) */
            IPUT_VAR1_DBL(rec_varids[16], 1 + gap)
            /* double volumeEdgeGlobal(Time) */
            IPUT_VAR1_DBL(rec_varids[17], 1 + gap)
            /* double CFLNumberGlobal(Time) */
            IPUT_VAR1_DBL(rec_varids[18], 1 + gap)
        }
        /* double BruntVaisalaFreqTop(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[19])
        /* double vertVelocityTop(Time, nCells, nVertLevelsP1) */
        REC_IPUT(5, rec_varids[20])
        /* double velocityZonal(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[21])
        /* double velocityMeridional(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[22])
        /* double displacedDensity(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[23])
        /* double potentialDensity(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[24])
        /* double pressure(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[25])
        /* double zMid(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[26])
        /* double columnIntegratedSpeed(Time, nCells) */
        REC_IPUT(0, rec_varids[27])
        /* double temperatureHorizontalAdvectionTendency(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[28])
        /* double salinityHorizontalAdvectionTendency(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[29])
        /* double temperatureVerticalAdvectionTendency(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[30])
        /* double salinityVerticalAdvectionTendency(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[31])
        /* double temperatureVertMixTendency(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[32])
        /* double salinityVertMixTendency(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[33])
        /* double temperatureSurfaceFluxTendency(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[34])
        /* double salinitySurfaceFluxTendency(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[35])
        /* double temperatureShortWaveTendency(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[36])
        /* double temperatureNonLocalTendency(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[37])
        /* double salinityNonLocalTendency(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[38])
        /* double temperature(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[39])
        /* double salinity(Time, nCells, nVertLevels) */
        REC_IPUT(2, rec_varids[40])

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

    if (dbl_buf != NULL) free(dbl_buf);
    if (chr_buf != NULL) free(chr_buf);
    if (int_buf != NULL) free(int_buf);

    if (cfg.api == pnetcdf) INQ_PUT_SIZE(total_size)

    FILE_CLOSE

    cfg.close_time = MPI_Wtime() - timing;

    /* for hdf5 blob I/O, write amount is calculated after file is closed */
    if (cfg.api != pnetcdf) total_size = cfg.amount_WR - previous_size;

    cfg.num_flushes = nflushes;
    cfg.num_decomp = decom.num_decomp;
    cfg.num_decomp_vars = num_decomp_vars;
    cfg.my_nreqs = my_nreqs;
    cfg.metadata_WR = metadata_size;
    cfg.amount_WR = total_size;
    cfg.end2end_time = MPI_Wtime() - cfg.end2end_time;

    /* check if there is any PnetCDF internal malloc residue */
    check_malloc(&cfg, &driver);

    /* report timing breakdowns */
    report_timing_WR(&cfg, &driver, cfg.out_path);

    /* print MPI-IO hints actually used */
    if (cfg.verbose && global_rank == 0) print_info(&info_used);

err_out:
    if (err < 0 && ncid >= 0) driver.close(ncid);
    if (info_used != MPI_INFO_NULL) MPI_Info_free(&info_used);
    if (!cfg.keep_outfile && sub_rank == 0) unlink(outfile);
    fflush(stdout);

    return err;
}

