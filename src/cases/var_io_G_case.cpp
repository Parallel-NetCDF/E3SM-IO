/*********************************************************************
 *
 * Copyright (C) 2020, Northwestern University
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

#define FIX_1D_VAR_STARTS_COUNTS(starts, counts, nreqs, disps, blocklens) { \
    starts = (MPI_Offset**) malloc(2 * nreqs * sizeof(MPI_Offset*)); \
    counts = starts + nreqs; \
    starts[0] = (MPI_Offset*) malloc(2 * nreqs * sizeof(MPI_Offset)); \
    counts[0] = starts[0] + nreqs; \
    \
    for (i=1; i<nreqs; i++) { \
        starts[i] = starts[i-1] + 1; \
        counts[i] = counts[i-1] + 1; \
    } \
    \
    for (i=0; i<nreqs; i++) { \
        starts[i][0] = disps[i]; \
        counts[i][0] = blocklens[i]; \
    } \
}

#define FIX_2D_VAR_STARTS_COUNTS(starts, counts, nreqs, disps, blocklens, last_dimlen) { \
    starts = (MPI_Offset**) malloc(2 * nreqs * sizeof(MPI_Offset*)); \
    counts = starts + nreqs; \
    starts[0] = (MPI_Offset*) malloc(2 * nreqs * 2 * sizeof(MPI_Offset)); \
    counts[0] = starts[0] + nreqs * 2; \
    \
    for (i=1; i<nreqs; i++) { \
        starts[i] = starts[i-1] + 2; \
        counts[i] = counts[i-1] + 2; \
    } \
    \
    j = 0; \
    starts[0][0] = disps[0] / last_dimlen; \
    starts[0][1] = disps[0] % last_dimlen; /* decomposition is 2D */ \
    counts[0][0] = 1; \
    counts[0][1] = blocklens[0]; /* each blocklens[i] is no bigger than last_dimlen */ \
    for (i=1; i<nreqs; i++) { \
        MPI_Offset _start[2]; \
        _start[0] = disps[i] / last_dimlen; \
        _start[1] = disps[i] % last_dimlen; \
        if (_start[0] == starts[j][0] + counts[j][0] && \
            _start[1] == starts[j][1] && blocklens[i] == counts[j][1]) \
            counts[j][0]++; \
        else { \
            j++; \
            starts[j][0] = _start[0]; \
            starts[j][1] = _start[1]; \
            counts[j][0] = 1; \
            counts[j][1] = blocklens[i]; /* each blocklens[i] is no bigger than last_dimlen */ \
        } \
    } \
}

#define REC_2D_VAR_STARTS_COUNTS(rec, starts, counts, nreqs, disps, blocklens) { \
    starts = (MPI_Offset**) malloc(2 * nreqs * sizeof(MPI_Offset*)); \
    counts = starts + nreqs; \
    starts[0] = (MPI_Offset*) malloc(2 * nreqs * 2 * sizeof(MPI_Offset)); \
    counts[0] = starts[0] + nreqs * 2; \
    \
    for (i=1; i<nreqs; i++) { \
        starts[i] = starts[i-1] + 2; \
        counts[i] = counts[i-1] + 2; \
    } \
    \
    for (i=0; i<nreqs; i++) { \
        starts[i][1] = disps[i]; /* decomposition is 1D */ \
        counts[i][1] = blocklens[i]; \
        \
        starts[i][0] = rec; /* record ID */ \
        counts[i][0] = 1;   /* one record only */ \
    } \
}

#define REC_3D_VAR_STARTS_COUNTS(rec, starts, counts, nreqs, disps, blocklens, last_dimlen) { \
    starts = (MPI_Offset**) malloc(2 * nreqs * sizeof(MPI_Offset*)); \
    counts = starts + nreqs; \
    starts[0] = (MPI_Offset*) malloc(2 * nreqs * 3 * sizeof(MPI_Offset)); \
    counts[0] = starts[0] + nreqs * 3; \
    \
    for (i=1; i<nreqs; i++) { \
        starts[i] = starts[i-1] + 3; \
        counts[i] = counts[i-1] + 3; \
    } \
    \
    j = 0; \
    starts[0][0] = rec; /* record ID */ \
    starts[0][1] = disps[0] / last_dimlen; \
    starts[0][2] = disps[0] % last_dimlen; /* decomposition is 2D */ \
    counts[0][0] = 1;   /* one record only */ \
    counts[0][1] = 1; \
    counts[0][2] = blocklens[0]; /* each blocklens[i] is no bigger than last_dimlen */ \
    for (i=1; i<nreqs; i++) { \
        MPI_Offset _start[2]; \
        _start[0] = disps[i] / last_dimlen; \
        _start[1] = disps[i] % last_dimlen; \
        if (starts[j][0] == rec && _start[0] == starts[j][1] + counts[j][1] && \
            _start[1] == starts[j][2] && blocklens[i] == counts[j][2]) \
            /* this request can be combined into the previous one */ \
            counts[j][1]++; \
        else { \
            /* this request cannot be combined into the previous one */ \
            j++; \
            starts[j][0] = rec; \
            starts[j][1] = _start[0]; \
            starts[j][2] = _start[1]; \
            counts[j][0] = 1; \
            counts[j][1] = 1; \
            counts[j][2] = blocklens[i]; /* each blocklens[i] is no bigger than last_dimlen */ \
        } \
    } \
    nreqs = j+1; \
}

#define FREE_N_NULL(A) if (A##p == NULL) free(A)
#define ASSIGN_BUF(A) *(A##p) = A

/*----< run_varn_G_case() >--------------------------------------------------*/
int
run_varn_G_case(MPI_Comm io_comm,         /* MPI communicator that includes all the tasks involved in IO */
                const char *out_dir,      /* output folder name */
                const char *outfile,      /* output file name */
                int         nvars,        /* number of variables 51 */
                int         num_recs,     /* number of records */
                MPI_Info    info,
                MPI_Offset  dims[6][2],   /* dimension lengths decomposition 1-6 */
                const int   nreqs[6],     /* no. request in decomposition 1-6 */
                int* const  disps[6],     /* request's displacements */
                int* const  blocklens[6], /* request's block lengths */
                int *D1_fix_int_bufp,     /* D1 fix int buffer */
                int *D2_fix_int_bufp,     /* D2 fix int buffer */
                int *D3_fix_int_bufp,     /* D3 fix int buffer */
                int *D4_fix_int_bufp,     /* D4 fix int buffer */
                int *D5_fix_int_bufp,     /* D5 fix int buffer */
                double *D1_rec_dbl_bufp,  /* D1 rec double buffer */
                double *D3_rec_dbl_bufp,  /* D3 rec double buffer */
                double *D4_rec_dbl_bufp,  /* D4 rec double buffer */
                double *D5_rec_dbl_bufp,  /* D5 rec double buffer */
                double *D6_rec_dbl_bufp,  /* D6 rec double buffer */
                double *D1_fix_dbl_bufp)  /* D1 fix double buffer */
{
    char outfname[512];
    int i, j, k, err, nerrs=0, rank, ncid, cmode, *varids;
    int rec_no, my_nreqs, nvars_D[6];
    size_t rec_buflen, nelems[6];
    double *D1_rec_dbl_buf, *D3_rec_dbl_buf, *D4_rec_dbl_buf, *D5_rec_dbl_buf, *D6_rec_dbl_buf, *rec_buf_ptr;
    int *D1_fix_int_buf, *D2_fix_int_buf, *D3_fix_int_buf, *D4_fix_int_buf, *D5_fix_int_buf;
    double *D1_fix_dbl_buf;
    double pre_timing, open_timing, post_timing, wait_timing, close_timing;
    double timing, total_timing, max_timing;
    MPI_Offset tmp, metadata_size, put_size, total_size, total_nreqs, max_nreqs;
    MPI_Offset **fix_starts_D1, **fix_counts_D1;
    MPI_Offset **fix_starts_D2, **fix_counts_D2;
    MPI_Offset **fix_starts_D3, **fix_counts_D3;
    MPI_Offset **fix_starts_D4, **fix_counts_D4;
    MPI_Offset **fix_starts_D5, **fix_counts_D5;
    MPI_Offset **starts_D1, **counts_D1;
    MPI_Offset **starts_D3, **counts_D3;
    MPI_Offset **starts_D4, **counts_D4;
    MPI_Offset **starts_D5, **counts_D5;
    MPI_Offset **starts_D6, **counts_D6;
    MPI_Info info_used=MPI_INFO_NULL;
    int nD1_rec_2d_vars = 4;
    int D1_rec_2d_varids[4] = {0, 4, 7, 38};
    int nD3_rec_3d_vars = 24;
    int D3_rec_3d_varids[24] = {5, 16, 17, 19, 26, 28, 29, 30, 31, 32,
                               34, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                               48, 49, 50, 51};
    int nD4_rec_3d_vars = 1;
    int D4_rec_3d_varids[1] = {6};
    int nD5_rec_3d_vars = 1;
    int D5_rec_3d_varids[1] = {18};
    int nD6_rec_3d_vars = 4;
    int D6_rec_3d_varids[4] = {1, 2, 3, 27};
    MPI_Offset stride[2] = {1, 1};
    MPI_Offset start[2] = {0, 0};
    MPI_Offset count[2] = {1, 1};
    double *dummy_double_buf=NULL;
    char dummy_char_buf[64];
    int xnreqs[6]; /* number of requests after combination */

    for (i=0; i<6; i++) {
        xnreqs[i] = nreqs[i];
        nvars_D[i] = 0;  /* number of variables using decomposition i */
    }

    MPI_Barrier(io_comm); /*-----------------------------------------*/
    total_timing = pre_timing = MPI_Wtime();

    open_timing = 0.0;
    post_timing = 0.0;
    wait_timing = 0.0;
    close_timing = 0.0;

    MPI_Comm_rank(io_comm, &rank);

    /* number of variable elements from 6 decompositions */
    total_nreqs = max_nreqs = 0;
    my_nreqs = 0;
    for (i=0; i<6; i++) {
        for (nelems[i]=0, k=0; k<xnreqs[i]; k++)
            nelems[i] += blocklens[i][k];
    }

    if (verbose && rank == 0)
        printf("nelems=%zd %zd %zd %zd %zd %zd\n",
               nelems[0],nelems[1],nelems[2],nelems[3],nelems[4],nelems[5]);

    /* construct varn API arguments starts[][] and counts[][] */
    if (xnreqs[0] > 0) {
        FIX_1D_VAR_STARTS_COUNTS(fix_starts_D1, fix_counts_D1, xnreqs[0], disps[0], blocklens[0])
        REC_2D_VAR_STARTS_COUNTS(0, starts_D1, counts_D1, xnreqs[0], disps[0], blocklens[0])
    }
    else {
        fix_starts_D1 = NULL;
        fix_counts_D1 = NULL;
        starts_D1 = NULL;
        counts_D1 = NULL;
    }

    if (xnreqs[1] > 0) {
        FIX_1D_VAR_STARTS_COUNTS(fix_starts_D2, fix_counts_D2, xnreqs[1], disps[1], blocklens[1])
    }
    else {
        fix_starts_D2 = NULL;
        fix_counts_D2 = NULL;
    }

    if (xnreqs[2] > 0) {
        FIX_2D_VAR_STARTS_COUNTS(fix_starts_D3, fix_counts_D3, xnreqs[2], disps[2], blocklens[2], dims[2][1])
        REC_3D_VAR_STARTS_COUNTS(0, starts_D3, counts_D3, xnreqs[2], disps[2], blocklens[2], dims[2][1])
    }
    else {
        fix_starts_D3 = NULL;
        fix_counts_D3 = NULL;
        starts_D3 = NULL;
        counts_D3 = NULL;
    }

    if (xnreqs[3] > 0) {
        FIX_2D_VAR_STARTS_COUNTS(fix_starts_D4, fix_counts_D4, xnreqs[3], disps[3], blocklens[3], dims[3][1])
        REC_3D_VAR_STARTS_COUNTS(0, starts_D4, counts_D4, xnreqs[3], disps[3], blocklens[3], dims[3][1])
    }
    else {
        fix_starts_D4 = NULL;
        fix_counts_D4 = NULL;
        starts_D4 = NULL;
        counts_D4 = NULL;
    }

    if (xnreqs[4] > 0) {
        FIX_2D_VAR_STARTS_COUNTS(fix_starts_D5, fix_counts_D5, xnreqs[4], disps[4], blocklens[4], dims[4][1])
        REC_3D_VAR_STARTS_COUNTS(0, starts_D5, counts_D5, xnreqs[4], disps[4], blocklens[4], dims[4][1])
    }
    else {
        fix_starts_D5 = NULL;
        fix_counts_D5 = NULL;
        starts_D5 = NULL;
        counts_D5 = NULL;
    }

    if (xnreqs[5] > 0) {
        REC_3D_VAR_STARTS_COUNTS(0, starts_D6, counts_D6, xnreqs[5], disps[5], blocklens[5], dims[5][1])
    }
    else {
        starts_D6 = NULL;
        counts_D6 = NULL;
    }

    /* allocate and initialize write buffer for 7 fixed-size variables */
    /* int (nCells): maxLevelCell */
    if (nelems[0] > 0) {        
        if (D1_fix_int_bufp != NULL){
            D1_fix_int_buf = D1_fix_int_bufp;
        }
        else{
            D1_fix_int_buf = (int*) malloc(nelems[0] * sizeof(int));
            for (i = 0; i < nelems[0]; i++) D1_fix_int_buf[i] = rank + i;
        }
    }
    else
        D1_fix_int_buf = NULL;

    /* int (nEdges): maxLevelEdgeTop and maxLevelEdgeBot */
    if (nelems[1] > 0) {
        if (D2_fix_int_bufp != NULL){
            D2_fix_int_buf = D2_fix_int_bufp;
        }
        else{
            D2_fix_int_buf = (int*) malloc(2 * nelems[1] * sizeof(int));
            for (i = 0; i < 2 * nelems[1]; i++) D2_fix_int_buf[i] = rank + i;
        }
    }
    else
        D2_fix_int_buf = NULL;

    /* int (nCells, nVertLevels): cellMask */
    if (nelems[2] > 0) {
        if (D3_fix_int_bufp != NULL){
            D3_fix_int_buf = D3_fix_int_bufp;
        }
        else{
            D3_fix_int_buf = (int*) malloc(nelems[2] * sizeof(int));
            for (i = 0; i < nelems[2]; i++) D3_fix_int_buf[i] = rank + i;
        }
    }
    else
        D3_fix_int_buf = NULL;

    /* int (nEdges, nVertLevels): edgeMask */
    if (nelems[3] > 0) {
        if (D4_fix_int_bufp != NULL){
            D4_fix_int_buf = D4_fix_int_bufp;
        }
        else{
            D4_fix_int_buf = (int*) malloc(nelems[3] * sizeof(int));
            for (i = 0; i < nelems[3]; i++) D4_fix_int_buf[i] = rank + i;
        }
    }
    else
        D4_fix_int_buf = NULL;

    /* int (nVertices, nVertLevels): vertexMask */
    if (nelems[4] > 0) {
        if (D5_fix_int_bufp != NULL){
            D5_fix_int_buf = D5_fix_int_bufp;
        }
        else{
            D5_fix_int_buf = (int*) malloc(nelems[4] * sizeof(int));
            for (i = 0; i < nelems[4]; i++) D5_fix_int_buf[i] = rank + i;
        }
    }
    else
        D5_fix_int_buf = NULL;

    /* double (nCells): bottomDepth */
    if (nelems[0] > 0) {
        if (D1_fix_dbl_bufp != NULL){
            D1_fix_dbl_buf = D1_fix_dbl_bufp;
        }
        else{
            D1_fix_dbl_buf = (double*) malloc(nelems[0] * sizeof(double));
            for (i = 0; i < nelems[0]; i++) D1_fix_dbl_buf[i] = rank + i;
        }
    }
    else
        D1_fix_dbl_buf = NULL;

    /* allocate and initialize write buffer for 34 record variables */
    if (nelems[0] > 0) {
        rec_buflen = nelems[0] * nD1_rec_2d_vars;
        if (D1_rec_dbl_bufp != NULL){
            D1_rec_dbl_buf = D1_rec_dbl_bufp;
        }
        else{
            D1_rec_dbl_buf = (double*) malloc(rec_buflen * sizeof(double));
            for (i = 0; i < rec_buflen; i++) D1_rec_dbl_buf[i] = rank + i;
        }
    }
    else
        D1_rec_dbl_buf = NULL;

    if (nelems[2] > 0) {
        rec_buflen = nelems[2] * nD3_rec_3d_vars;
        if (D3_rec_dbl_bufp != NULL){
            D3_rec_dbl_buf = D3_rec_dbl_bufp;
        }
        else{
            D3_rec_dbl_buf = (double*) malloc(rec_buflen * sizeof(double));
            for (i = 0; i < rec_buflen; i++) D3_rec_dbl_buf[i] = rank + i;
        }
    }
    else
        D3_rec_dbl_buf = NULL;

    if (nelems[3] > 0) {
        rec_buflen = nelems[3] * nD4_rec_3d_vars;
        if (D4_rec_dbl_bufp != NULL){
            D4_rec_dbl_buf = D4_rec_dbl_bufp;
        }
        else{
            D4_rec_dbl_buf = (double*) malloc(rec_buflen * sizeof(double));
            for (i = 0; i < rec_buflen; i++) D4_rec_dbl_buf[i] = rank + i;
        }
    }
    else
        D4_rec_dbl_buf = NULL;

    if (nelems[4] > 0) {
        rec_buflen = nelems[4] * nD5_rec_3d_vars;
        if (D5_rec_dbl_bufp != NULL){
            D5_rec_dbl_buf = D5_rec_dbl_bufp;
        }
        else{
            D5_rec_dbl_buf = (double*) malloc(rec_buflen * sizeof(double));
            for (i = 0; i < rec_buflen; i++) D5_rec_dbl_buf[i] = rank + i;
        }
    }
    else
        D5_rec_dbl_buf = NULL;

    if (nelems[5] > 0) {
        rec_buflen = nelems[5] * nD6_rec_3d_vars;
        if (D6_rec_dbl_bufp != NULL){
            D6_rec_dbl_buf = D6_rec_dbl_bufp;
        }
        else{
            D6_rec_dbl_buf = (double*) malloc(rec_buflen * sizeof(double));
            for (i = 0; i < rec_buflen; i++) D6_rec_dbl_buf[i] = rank + i;
        }
    }
    else
        D6_rec_dbl_buf = NULL;

    /* initialize write buffer for 11 small variables */
    dummy_double_buf = (double*) malloc(dims[2][1] * sizeof(double));
    for (i = 0; i < dims[2][1]; i++) dummy_double_buf[i] = rank + i;
    for (i = 0; i < 64; i++) dummy_char_buf[i] = 'a' + rank + i;

    varids = (int*) malloc(nvars * sizeof(int));

    pre_timing = MPI_Wtime() - pre_timing;

    MPI_Barrier(io_comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    /* set output file name */
    sprintf(outfname, "%s/%s",out_dir, outfile);

    /* create a new CDF-5 file for writing */
    cmode = NC_CLOBBER | NC_64BIT_DATA;
    err = ncmpi_create(io_comm, outfname, cmode, info, &ncid); ERR

    MPI_Offset put_buffer_size_limit = 10485760;
    err = ncmpi_buffer_attach(ncid, put_buffer_size_limit); ERR

    /* define dimensions, variables, and attributes */
    err = def_G_case_h0(ncid, dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], nvars, varids); ERR

    /* exit define mode and enter data mode */
    err = ncmpi_enddef(ncid); ERR

    /* I/O amount so far */
    err = ncmpi_inq_put_size(ncid, &metadata_size); ERR
    err = ncmpi_inq_file_info(ncid, &info_used); ERR
    open_timing += MPI_Wtime() - timing;

    MPI_Barrier(io_comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    /* write 7 fixed-size variables */

    /* int maxLevelEdgeTop(nEdges) */
    err = ncmpi_iput_varn(ncid, 8, xnreqs[1], fix_starts_D2, fix_counts_D2,
                          D2_fix_int_buf, nelems[1], MPI_INT, NULL); ERR
    my_nreqs += xnreqs[1];
    nvars_D[1]++;

    /* int maxLevelEdgeBot(nEdges) */
    err = ncmpi_iput_varn(ncid, 37, xnreqs[1], fix_starts_D2, fix_counts_D2,
                          D2_fix_int_buf + nelems[1], nelems[1], MPI_INT, NULL); ERR
    my_nreqs += xnreqs[1];
    nvars_D[1]++;

    /* int edgeMask(nEdges, nVertLevels) */
    err = ncmpi_iput_varn(ncid, 10, xnreqs[3], fix_starts_D4, fix_counts_D4,
                          D4_fix_int_buf, nelems[3], MPI_INT, NULL); ERR
    my_nreqs += xnreqs[3];
    nvars_D[3]++;

    /* int cellMask(nCells, nVertLevels) */
    err = ncmpi_iput_varn(ncid, 11, xnreqs[2], fix_starts_D3, fix_counts_D3,
                          D3_fix_int_buf, nelems[2], MPI_INT, NULL); ERR
    my_nreqs += xnreqs[2];
    nvars_D[2]++;

    /* int vertexMask(nVertices, nVertLevels) */
    err = ncmpi_iput_varn(ncid, 12, xnreqs[4], fix_starts_D5, fix_counts_D5,
                          D5_fix_int_buf, nelems[4], MPI_INT, NULL); ERR
    my_nreqs += xnreqs[4];
    nvars_D[4]++;

    /* double bottomDepth(nCells)  */
    err = ncmpi_iput_varn(ncid, 35, xnreqs[0], fix_starts_D1, fix_counts_D1,
                          D1_fix_dbl_buf, nelems[0], MPI_DOUBLE, NULL); ERR
    my_nreqs += xnreqs[0];
    nvars_D[0]++;

    /* int maxLevelCell(nCells) */
    err = ncmpi_iput_varn(ncid, 36, xnreqs[0], fix_starts_D1, fix_counts_D1,
                          D1_fix_int_buf, nelems[0], MPI_INT, NULL); ERR
    my_nreqs += xnreqs[0];
    nvars_D[0]++;

    /* next 11 small variables are written by rank 0 only */
    if (rank == 0) {
        count[0] = dims[2][1]; /* dimension nVertLevels */

        /* double vertCoordMovementWeights(nVertLevels) */
        err = ncmpi_bput_vars_double(ncid, 9, start, count, stride, dummy_double_buf, NULL); ERR

        /* double refZMid(nVertLevels) */
        err = ncmpi_bput_vars_double(ncid, 13, start, count, stride, dummy_double_buf, NULL); ERR

        /* double refLayerThickness(nVertLevels) */
        err = ncmpi_bput_vars_double(ncid, 14, start, count, stride, dummy_double_buf, NULL); ERR

        /* double refBottomDepth(nVertLevels) */
        err = ncmpi_bput_vars_double(ncid, 33, start, count, stride, dummy_double_buf, NULL); ERR

        my_nreqs += 4; /* 4 non-record variables */

        for (rec_no = 0; rec_no < num_recs; rec_no++) {
            start[0] = rec_no; count[0] = 1;

            count[1] = 64; /* dimension StrLen */

            /* char xtime(Time, StrLen) */
            err = ncmpi_bput_vars_text(ncid, 15, start, count, stride, dummy_char_buf, NULL); ERR

            /* double areaCellGlobal(Time) */
            err = ncmpi_bput_vars_double(ncid, 20, start, count, stride, dummy_double_buf, NULL); ERR

            /* double areaEdgeGlobal(Time) */
            err = ncmpi_bput_vars_double(ncid, 21, start, count, stride, dummy_double_buf, NULL); ERR

            /* double areaTriangleGlobal(Time) */
            err = ncmpi_bput_vars_double(ncid, 22, start, count, stride, dummy_double_buf, NULL); ERR

            /* double volumeCellGlobal(Time) */
            err = ncmpi_bput_vars_double(ncid, 23, start, count, stride, dummy_double_buf, NULL); ERR

            /* double volumeEdgeGlobal(Time) */
            err = ncmpi_bput_vars_double(ncid, 24, start, count, stride, dummy_double_buf, NULL); ERR

            /* double CFLNumberGlobal(Time) */
            err = ncmpi_bput_vars_double(ncid, 25, start, count, stride, dummy_double_buf, NULL); ERR

            my_nreqs += 7; /* 7 record variables */
        }
    }

    /* write 34 record variables */

    /* 4 D1 record variables: double (Time, nCells) */
    for (rec_no = 0; rec_no < num_recs; rec_no++) {
        for (j = 0; j < xnreqs[0]; j++) starts_D1[j][0] = rec_no;

        rec_buf_ptr = D1_rec_dbl_buf;
        for (j = 0; j < nD1_rec_2d_vars; j++) {
            err = ncmpi_iput_varn(ncid, D1_rec_2d_varids[j], xnreqs[0], starts_D1,
                                  counts_D1, rec_buf_ptr, nelems[0], MPI_DOUBLE, NULL); ERR
            rec_buf_ptr += nelems[0];
            my_nreqs += xnreqs[0];
            if (rec_no == 0) nvars_D[0]++;
        }
    }

    /* 4 D6 record variables: double (Time, nCells, nVertLevelsP1) */
    for (rec_no = 0; rec_no < num_recs; rec_no++) {
        for (j = 0; j < xnreqs[5]; j++) starts_D6[j][0] = rec_no;

        rec_buf_ptr = D6_rec_dbl_buf;
        for (j = 0; j < nD6_rec_3d_vars; j++) {
            err = ncmpi_iput_varn(ncid, D6_rec_3d_varids[j], xnreqs[5], starts_D6,
                                  counts_D6, rec_buf_ptr, nelems[5], MPI_DOUBLE, NULL); ERR
            rec_buf_ptr += nelems[5];
            my_nreqs += xnreqs[5];
            if (rec_no == 0) nvars_D[5]++;
        }
    }

    /* 24 D3 record variables: double (Time, nCells, nVertLevels) */
    for (rec_no = 0; rec_no < num_recs; rec_no++) {
        for (j = 0; j < xnreqs[2]; j++) starts_D3[j][0] = rec_no;

        rec_buf_ptr = D3_rec_dbl_buf;
        for (j = 0; j < nD3_rec_3d_vars; j++) {
            err = ncmpi_iput_varn(ncid, D3_rec_3d_varids[j], xnreqs[2], starts_D3,
                                  counts_D3, rec_buf_ptr, nelems[2], MPI_DOUBLE, NULL); ERR
            rec_buf_ptr += nelems[2];
            my_nreqs += xnreqs[2];
            if (rec_no == 0) nvars_D[2]++;
        }
    }

    /* 1 D4 record variable: double (Time, nEdges, nVertLevels) */
    for (rec_no = 0; rec_no < num_recs; rec_no++) {
        for (j = 0; j < xnreqs[3]; j++) starts_D4[j][0] = rec_no;

        rec_buf_ptr = D4_rec_dbl_buf;
        for (j = 0; j < nD4_rec_3d_vars; j++) {
            err = ncmpi_iput_varn(ncid, D4_rec_3d_varids[j], xnreqs[3], starts_D4,
                                  counts_D4, rec_buf_ptr, nelems[3], MPI_DOUBLE, NULL); ERR
            rec_buf_ptr += nelems[3];
            my_nreqs += xnreqs[3];
            if (rec_no == 0) nvars_D[3]++;
        }
    }

    /* 1 D5 record variable: double (Time, nVertices, nVertLevels) */
    for (rec_no = 0; rec_no < num_recs; rec_no++) {
        for (j = 0; j < xnreqs[4]; j++) starts_D5[j][0] = rec_no;

        rec_buf_ptr = D5_rec_dbl_buf;
        for (j = 0; j < nD5_rec_3d_vars; j++) {
            err = ncmpi_iput_varn(ncid, D5_rec_3d_varids[j], xnreqs[4], starts_D5,
                                  counts_D5, rec_buf_ptr, nelems[4], MPI_DOUBLE, NULL); ERR
            rec_buf_ptr += nelems[4];
            my_nreqs += xnreqs[4];
            if (rec_no == 0) nvars_D[4]++;
        }
    }
    total_nreqs += my_nreqs;

    post_timing += MPI_Wtime() - timing;

    MPI_Barrier(io_comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); ERR

    wait_timing += MPI_Wtime() - timing;

    MPI_Barrier(io_comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    err = ncmpi_inq_put_size(ncid, &total_size); ERR
    put_size = total_size - metadata_size;
    err = ncmpi_buffer_detach(ncid); ERR
    err = ncmpi_close(ncid); ERR
    close_timing += MPI_Wtime() - timing;

    if (dummy_double_buf != NULL) free(dummy_double_buf);

    if (xnreqs[0] > 0) {
        free(fix_starts_D1[0]); free(fix_starts_D1);
        FREE_N_NULL(D1_fix_int_buf); FREE_N_NULL(D1_fix_dbl_buf);
        free(starts_D1[0]); free(starts_D1);
        FREE_N_NULL(D1_rec_dbl_buf);
    }

    if (xnreqs[1] > 0) {
        free(fix_starts_D2[0]); free(fix_starts_D2);
        FREE_N_NULL(D2_fix_int_buf);
    }

    if (xnreqs[2] > 0) {
        free(fix_starts_D3[0]); free(fix_starts_D3);
        FREE_N_NULL(D3_fix_int_buf);
        free(starts_D3[0]); free(starts_D3);
        FREE_N_NULL(D3_rec_dbl_buf);
    }

    if (xnreqs[3] > 0) {
        free(fix_starts_D4[0]); free(fix_starts_D4);
        FREE_N_NULL(D4_fix_int_buf);
        free(starts_D4[0]); free(starts_D4);
        FREE_N_NULL(D4_rec_dbl_buf);
    }

    if (xnreqs[4] > 0) {
        free(fix_starts_D5[0]); free(fix_starts_D5);
        FREE_N_NULL(D5_fix_int_buf);
        free(starts_D5[0]); free(starts_D5);
        FREE_N_NULL(D5_rec_dbl_buf);
    }

    if (xnreqs[5] > 0) {
        free(starts_D6[0]); free(starts_D6);
        FREE_N_NULL(D6_rec_dbl_buf);
    }

    free(varids);

    total_timing = MPI_Wtime() - total_timing;

    MPI_Reduce(&total_nreqs,   &max_nreqs,  1, MPI_OFFSET, MPI_MAX, 0, io_comm);
    MPI_Reduce(&total_nreqs,   &tmp,        1, MPI_OFFSET, MPI_SUM, 0, io_comm);
    total_nreqs = tmp;
    MPI_Reduce(&put_size,      &tmp,        1, MPI_OFFSET, MPI_SUM, 0, io_comm);
    put_size = tmp;
    MPI_Reduce(&total_size,    &tmp,        1, MPI_OFFSET, MPI_SUM, 0, io_comm);
    total_size = tmp;
    MPI_Reduce(&open_timing,   &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, io_comm);
    open_timing = max_timing;
    MPI_Reduce(&pre_timing,    &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, io_comm);
    pre_timing = max_timing;
    MPI_Reduce(&post_timing,   &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, io_comm);
    post_timing = max_timing;
    MPI_Reduce(&wait_timing,   &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, io_comm);
    wait_timing = max_timing;
    MPI_Reduce(&close_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, io_comm);
    close_timing = max_timing;
    MPI_Reduce(&total_timing,  &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, io_comm);
    total_timing = max_timing;

    /* check if there is any PnetCDF internal malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, io_comm);
        if (rank == 0 && sum_size > 0) {
            printf("-----------------------------------------------------------\n");
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
        }
    }
    MPI_Offset m_alloc=0, max_alloc;
    ncmpi_inq_malloc_max_size(&m_alloc);
    MPI_Reduce(&m_alloc, &max_alloc, 1, MPI_OFFSET, MPI_MAX, 0, io_comm);
    if (rank == 0) {
        int nvars_noD = nvars;
        for (i=0; i<6; i++) nvars_noD -= nvars_D[i];
        printf("History output file                = %s\n", outfile);
        printf("No. variables use no decomposition = %3d\n", nvars_noD);
        printf("No. variables use decomposition D1 = %3d\n", nvars_D[0]);
        printf("No. variables use decomposition D2 = %3d\n", nvars_D[1]);
        printf("No. variables use decomposition D3 = %3d\n", nvars_D[2]);
        printf("No. variables use decomposition D4 = %3d\n", nvars_D[3]);
        printf("No. variables use decomposition D5 = %3d\n", nvars_D[4]);
        printf("No. variables use decomposition D6 = %3d\n", nvars_D[5]);
        printf("Total number of variables          = %3d\n",nvars);
        printf("MAX heap memory allocated by PnetCDF internally is %.2f MiB\n",
               (float)max_alloc/1048576);
        printf("Total write amount                 = %.2f MiB = %.2f GiB\n",
               (double)total_size/1048576,(double)total_size/1073741824);
        printf("Total number of requests           = %lld\n",total_nreqs);
        printf("Max number of requests             = %lld\n",max_nreqs);
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
    MPI_Barrier(io_comm);
    return nerrs;
}

/*----< run_varn_G_case_rd() >--------------------------------------------------*/
int
run_varn_G_case_rd( MPI_Comm io_comm,         /* MPI communicator that includes all the tasks involved in IO */
                    char       *out_dir,      /* output folder name */
                    char       *outfile,      /* output file name */
                    int         nvars,        /* number of variables 51 */
                    int         num_recs,     /* number of records */
                    MPI_Info    info,
                    MPI_Offset  dims[6][2],   /* dimension lengths decomposition 1-6 */
                    int         nreqs[6],     /* no. request in decomposition 1-6 */
                    int        *disps[6],     /* request's displacements */
                    int        *blocklens[6], /* request's block lengths */
                    int **D1_fix_int_bufp,     /* D1 fix int buffer */
                    int **D2_fix_int_bufp,     /* D2 fix int buffer */
                    int **D3_fix_int_bufp,     /* D3 fix int buffer */
                    int **D4_fix_int_bufp,     /* D4 fix int buffer */
                    int **D5_fix_int_bufp,     /* D5 fix int buffer */
                    double **D1_rec_dbl_bufp,  /* D1 rec double buffer */
                    double **D3_rec_dbl_bufp,  /* D3 rec double buffer */
                    double **D4_rec_dbl_bufp,  /* D4 rec double buffer */
                    double **D5_rec_dbl_bufp,  /* D5 rec double buffer */
                    double **D6_rec_dbl_bufp,  /* D6 rec double buffer */
                    double **D1_fix_dbl_bufp)  /* D1 fix double buffer */
{
    char outfname[512];
    int i, j, k, err, nerrs=0, rank, ncid, cmode, *varids;
    int rec_no, my_nreqs;
    size_t rec_buflen, nelems[6];
    double *D1_rec_dbl_buf, *D3_rec_dbl_buf, *D4_rec_dbl_buf, *D5_rec_dbl_buf, *D6_rec_dbl_buf, *rec_buf_ptr;
    int *D1_fix_int_buf, *D2_fix_int_buf, *D3_fix_int_buf, *D4_fix_int_buf, *D5_fix_int_buf;
    double *D1_fix_dbl_buf;
    double pre_timing, open_timing, post_timing, wait_timing, close_timing;
    double timing, total_timing, max_timing;
    MPI_Offset tmp, metadata_size, get_size, total_size, total_nreqs, max_nreqs;
    MPI_Offset **fix_starts_D1, **fix_counts_D1;
    MPI_Offset **fix_starts_D2, **fix_counts_D2;
    MPI_Offset **fix_starts_D3, **fix_counts_D3;
    MPI_Offset **fix_starts_D4, **fix_counts_D4;
    MPI_Offset **fix_starts_D5, **fix_counts_D5;
    MPI_Offset **starts_D1, **counts_D1;
    MPI_Offset **starts_D3, **counts_D3;
    MPI_Offset **starts_D4, **counts_D4;
    MPI_Offset **starts_D5, **counts_D5;
    MPI_Offset **starts_D6, **counts_D6;
    MPI_Comm comm=io_comm;
    MPI_Info info_used=MPI_INFO_NULL;
    int nD1_rec_2d_vars = 4;
    int D1_rec_2d_varids[4] = {0, 4, 7, 38};
    int nD3_rec_3d_vars = 24;
    int D3_rec_3d_varids[24] = {5, 16, 17, 19, 26, 28, 29, 30, 31, 32,
                               34, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                               48, 49, 50, 51};
    int nD4_rec_3d_vars = 1;
    int D4_rec_3d_varids[1] = {6};
    int nD5_rec_3d_vars = 1;
    int D5_rec_3d_varids[1] = {18};
    int nD6_rec_3d_vars = 4;
    int D6_rec_3d_varids[4] = {1, 2, 3, 27};
    MPI_Offset stride[2] = {1, 1};
    MPI_Offset start[2] = {0, 0};
    MPI_Offset count[2] = {1, 1};
    double dummy_double_buf[80];
    char dummy_char_buf[64];

    MPI_Barrier(comm); /*-----------------------------------------*/
    total_timing = pre_timing = MPI_Wtime();

    open_timing = 0.0;
    post_timing = 0.0;
    wait_timing = 0.0;
    close_timing = 0.0;

    MPI_Comm_rank(comm, &rank);

    /* number of variable elements from 6 decompositions */
    total_nreqs = max_nreqs = 0;
    my_nreqs = 0;
    for (i=0; i<6; i++) {
        for (nelems[i]=0, k=0; k<nreqs[i]; k++)
            nelems[i] += blocklens[i][k];
    }

    if (verbose && rank == 0)
        printf("nelems=%zd %zd %zd %zd %zd %zd\n",
               nelems[0],nelems[1],nelems[2],nelems[3],nelems[4],nelems[5]);

    /* construct varn API arguments starts[][] and counts[][] */
    if (nreqs[0] > 0) {
        FIX_1D_VAR_STARTS_COUNTS(fix_starts_D1, fix_counts_D1, nreqs[0], disps[0], blocklens[0])
        REC_2D_VAR_STARTS_COUNTS(0, starts_D1, counts_D1, nreqs[0], disps[0], blocklens[0])
    }
    else {
        fix_starts_D1 = NULL;
        fix_counts_D1 = NULL;
        starts_D1 = NULL;
        counts_D1 = NULL;
    }

    if (nreqs[1] > 0) {
        FIX_1D_VAR_STARTS_COUNTS(fix_starts_D2, fix_counts_D2, nreqs[1], disps[1], blocklens[1])
    }
    else {
        fix_starts_D2 = NULL;
        fix_counts_D2 = NULL;
    }

    if (nreqs[2] > 0) {
        FIX_2D_VAR_STARTS_COUNTS(fix_starts_D3, fix_counts_D3, nreqs[2], disps[2], blocklens[2], dims[2][1])
        REC_3D_VAR_STARTS_COUNTS(0, starts_D3, counts_D3, nreqs[2], disps[2], blocklens[2], dims[2][1])
    }
    else {
        fix_starts_D3 = NULL;
        fix_counts_D3 = NULL;
        starts_D3 = NULL;
        counts_D3 = NULL;
    }

    if (nreqs[3] > 0) {
        FIX_2D_VAR_STARTS_COUNTS(fix_starts_D4, fix_counts_D4, nreqs[3], disps[3], blocklens[3], dims[3][1])
        REC_3D_VAR_STARTS_COUNTS(0, starts_D4, counts_D4, nreqs[3], disps[3], blocklens[3], dims[3][1])
    }
    else {
        fix_starts_D4 = NULL;
        fix_counts_D4 = NULL;
        starts_D4 = NULL;
        counts_D4 = NULL;
    }

    if (nreqs[4] > 0) {
        FIX_2D_VAR_STARTS_COUNTS(fix_starts_D5, fix_counts_D5, nreqs[4], disps[4], blocklens[4], dims[4][1])
        REC_3D_VAR_STARTS_COUNTS(0, starts_D5, counts_D5, nreqs[4], disps[4], blocklens[4], dims[4][1])
    }
    else {
        fix_starts_D5 = NULL;
        fix_counts_D5 = NULL;
        starts_D5 = NULL;
        counts_D5 = NULL;
    }

    if (nreqs[5] > 0) {
        REC_3D_VAR_STARTS_COUNTS(0, starts_D6, counts_D6, nreqs[5], disps[5], blocklens[5], dims[5][1])
    }
    else {
        starts_D6 = NULL;
        counts_D6 = NULL;
    }

    //printf("Rank: %d, nelems = [%lld, %lld, %lld, %lld, %lld, %lld], rec_buflen = %lld\n", rank, nelems[0], nelems[1], nelems[2], nelems[3], nelems[4], nelems[5], rec_buflen); fflush(stdout);

    /* allocate and initialize write buffer for 7 fixed-size variables */
    /* int (nCells): maxLevelCell */
    if (nelems[0] > 0) {
        D1_fix_int_buf = (int*) malloc(nelems[0] * sizeof(int));
        for (i = 0; i < nelems[0]; i++) D1_fix_int_buf[i] = rank + i;
        *D1_fix_int_bufp = D1_fix_int_buf;
    }
    else
        D1_fix_int_buf = NULL;

    /* int (nEdges): maxLevelEdgeTop and maxLevelEdgeBot */
    if (nelems[1] > 0) {
        D2_fix_int_buf = (int*) malloc(2 * nelems[1] * sizeof(int));
        for (i = 0; i < 2 * nelems[1]; i++) D2_fix_int_buf[i] = rank + i;
        *D2_fix_int_bufp = D2_fix_int_buf;
    }
    else
        D2_fix_int_buf = NULL;

    /* int (nCells, nVertLevels): cellMask */
    if (nelems[2] > 0) {
        D3_fix_int_buf = (int*) malloc(nelems[2] * sizeof(int));
        for (i = 0; i < nelems[2]; i++) D3_fix_int_buf[i] = rank + i;
        *D3_fix_int_bufp = D3_fix_int_buf;
    }
    else
        D3_fix_int_buf = NULL;

    /* int (nEdges, nVertLevels): edgeMask */
    if (nelems[3] > 0) {
        D4_fix_int_buf = (int*) malloc(nelems[3] * sizeof(int));
        for (i = 0; i < nelems[3]; i++) D4_fix_int_buf[i] = rank + i;
        *D4_fix_int_bufp = D4_fix_int_buf;
    }
    else
        D4_fix_int_buf = NULL;

    /* int (nVertices, nVertLevels): vertexMask */
    if (nelems[4] > 0) {
        D5_fix_int_buf = (int*) malloc(nelems[4] * sizeof(int));
        for (i = 0; i < nelems[4]; i++) D5_fix_int_buf[i] = rank + i;
        *D5_fix_int_bufp = D5_fix_int_buf;
    }
    else
        D5_fix_int_buf = NULL;

    /* double (nCells): bottomDepth */
    if (nelems[0] > 0) {
        D1_fix_dbl_buf = (double*) malloc(nelems[0] * sizeof(double));
        for (i = 0; i < nelems[0]; i++) D1_fix_dbl_buf[i] = rank + i;
        *D1_fix_dbl_bufp = D1_fix_dbl_buf;
    }
    else
        D1_fix_dbl_buf = NULL;

    /* allocate and initialize write buffer for 34 record variables */
    if (nelems[0] > 0) {
        rec_buflen = nelems[0] * nD1_rec_2d_vars;
        D1_rec_dbl_buf = (double*) malloc(rec_buflen * sizeof(double));
        for (i = 0; i < rec_buflen; i++) D1_rec_dbl_buf[i] = rank + i;
        *D1_rec_dbl_bufp = D1_rec_dbl_buf;
    }
    else
        D1_rec_dbl_buf = NULL;

    if (nelems[2] > 0) {
        rec_buflen = nelems[2] * nD3_rec_3d_vars;
        D3_rec_dbl_buf = (double*) malloc(rec_buflen * sizeof(double));
        for (i = 0; i < rec_buflen; i++) D3_rec_dbl_buf[i] = rank + i;
        *D3_rec_dbl_bufp = D3_rec_dbl_buf;
    }
    else
        D3_rec_dbl_buf = NULL;

    if (nelems[3] > 0) {
        rec_buflen = nelems[3] * nD4_rec_3d_vars;
        D4_rec_dbl_buf = (double*) malloc(rec_buflen * sizeof(double));
        for (i = 0; i < rec_buflen; i++) D4_rec_dbl_buf[i] = rank + i;
        *D4_rec_dbl_bufp = D4_rec_dbl_buf;
    }
    else
        D4_rec_dbl_buf = NULL;

    if (nelems[4] > 0) {
        rec_buflen = nelems[4] * nD5_rec_3d_vars;
        D5_rec_dbl_buf = (double*) malloc(rec_buflen * sizeof(double));
        for (i = 0; i < rec_buflen; i++) D5_rec_dbl_buf[i] = rank + i;
        *D5_rec_dbl_bufp = D5_rec_dbl_buf;
    }
    else
        D5_rec_dbl_buf = NULL;

    if (nelems[5] > 0) {
        rec_buflen = nelems[5] * nD6_rec_3d_vars;
        D6_rec_dbl_buf = (double*) malloc(rec_buflen * sizeof(double));
        for (i = 0; i < rec_buflen; i++) D6_rec_dbl_buf[i] = rank + i;
        *D6_rec_dbl_bufp = D6_rec_dbl_buf;
    }
    else
        D6_rec_dbl_buf = NULL;

    /* initialize write buffer for 11 small variables */
    for (i = 0; i < 80; i++) dummy_double_buf[i] = rank + i;
    for (i = 0; i < 64; i++) dummy_char_buf[i] = 'a' + rank + i;

    varids = (int*) malloc(nvars * sizeof(int));

    pre_timing = MPI_Wtime() - pre_timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    /* set output file name */
    sprintf(outfname, "%s/%s",out_dir, outfile);

    /* create a new CDF-5 file for writing */
    cmode = NC_64BIT_DATA;
    err = ncmpi_open(comm, outfname, cmode, info, &ncid); ERR

    MPI_Offset put_buffer_size_limit = 10485760;
    err = ncmpi_buffer_attach(ncid, put_buffer_size_limit); ERR

    /* define dimensions, variables, and attributes */
    err = inq_G_case_h0(ncid, dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], nvars, varids); ERR

    /* exit define mode and enter data mode */
    //err = ncmpi_enddef(ncid); ERR

    /* I/O amount so far */
    err = ncmpi_inq_get_size(ncid, &metadata_size); ERR
    err = ncmpi_inq_file_info(ncid, &info_used); ERR
    open_timing += MPI_Wtime() - timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    /* read 7 fixed-size variables */

    /* int maxLevelEdgeTop(nEdges) */
    err = ncmpi_iget_varn(ncid, 8, nreqs[1], fix_starts_D2, fix_counts_D2,
                          D2_fix_int_buf, nelems[1], MPI_INT, NULL); ERR
    my_nreqs += nreqs[1];

    /* int maxLevelEdgeBot(nEdges) */
    err = ncmpi_iget_varn(ncid, 37, nreqs[1], fix_starts_D2, fix_counts_D2,
                          D2_fix_int_buf + nelems[1], nelems[1], MPI_INT, NULL); ERR
    my_nreqs += nreqs[1];

    /* int edgeMask(nEdges, nVertLevels) */
    err = ncmpi_iget_varn(ncid, 10, nreqs[3], fix_starts_D4, fix_counts_D4,
                          D4_fix_int_buf, nelems[3], MPI_INT, NULL); ERR
    my_nreqs += nreqs[3];

    /* int cellMask(nCells, nVertLevels) */
    err = ncmpi_iget_varn(ncid, 11, nreqs[2], fix_starts_D3, fix_counts_D3,
                          D3_fix_int_buf, nelems[2], MPI_INT, NULL); ERR
    my_nreqs += nreqs[2];

    /* int vertexMask(nVertices, nVertLevels) */
    err = ncmpi_iget_varn(ncid, 12, nreqs[4], fix_starts_D5, fix_counts_D5,
                          D5_fix_int_buf, nelems[4], MPI_INT, NULL); ERR
    my_nreqs += nreqs[4];

    /* double bottomDepth(nCells)  */
    err = ncmpi_iget_varn(ncid, 35, nreqs[0], fix_starts_D1, fix_counts_D1,
                          D1_fix_dbl_buf, nelems[0], MPI_DOUBLE, NULL); ERR
    my_nreqs += nreqs[0];

    /* int maxLevelCell(nCells) */
    err = ncmpi_iget_varn(ncid, 36, nreqs[0], fix_starts_D1, fix_counts_D1,
                          D1_fix_int_buf, nelems[0], MPI_INT, NULL); ERR
    my_nreqs += nreqs[0];

    /* next 11 small variables are written by rank 0 only */
    if (rank == 0) {
        count[0] = dims[2][1]; /* dimension nVertLevels */

        /* double vertCoordMovementWeights(nVertLevels) */
        err = ncmpi_iget_vars_double(ncid, 9, start, count, stride, dummy_double_buf, NULL); ERR

        /* double refZMid(nVertLevels) */
        err = ncmpi_iget_vars_double(ncid, 13, start, count, stride, dummy_double_buf, NULL); ERR

        /* double refLayerThickness(nVertLevels) */
        err = ncmpi_iget_vars_double(ncid, 14, start, count, stride, dummy_double_buf, NULL); ERR

        /* double refBottomDepth(nVertLevels) */
        err = ncmpi_iget_vars_double(ncid, 33, start, count, stride, dummy_double_buf, NULL); ERR

        my_nreqs += 4; /* 4 non-record variables */

        for (rec_no = 0; rec_no < num_recs; rec_no++) {
            start[0] = rec_no; count[0] = 1;

            count[1] = 64; /* dimension StrLen */

            /* char xtime(Time, StrLen) */
            err = ncmpi_iget_vars_text(ncid, 15, start, count, stride, dummy_char_buf, NULL); ERR

            /* double areaCellGlobal(Time) */
            err = ncmpi_iget_vars_double(ncid, 20, start, count, stride, dummy_double_buf, NULL); ERR

            /* double areaEdgeGlobal(Time) */
            err = ncmpi_iget_vars_double(ncid, 21, start, count, stride, dummy_double_buf, NULL); ERR

            /* double areaTriangleGlobal(Time) */
            err = ncmpi_iget_vars_double(ncid, 22, start, count, stride, dummy_double_buf, NULL); ERR

            /* double volumeCellGlobal(Time) */
            err = ncmpi_iget_vars_double(ncid, 23, start, count, stride, dummy_double_buf, NULL); ERR

            /* double volumeEdgeGlobal(Time) */
            err = ncmpi_iget_vars_double(ncid, 24, start, count, stride, dummy_double_buf, NULL); ERR

            /* double CFLNumberGlobal(Time) */
            err = ncmpi_iget_vars_double(ncid, 25, start, count, stride, dummy_double_buf, NULL); ERR

            my_nreqs += 7; /* 7 record variables */
        }
    }

    /* read 34 record variables */

    /* 4 D1 record variables: double (Time, nCells) */
    for (rec_no = 0; rec_no < num_recs; rec_no++) {
        for (j = 0; j < nreqs[0]; j++) starts_D1[j][0] = rec_no;

        rec_buf_ptr = D1_rec_dbl_buf;
        for (j = 0; j < nD1_rec_2d_vars; j++) {
            err = ncmpi_iget_varn(ncid, D1_rec_2d_varids[j], nreqs[0], starts_D1,
                                  counts_D1, rec_buf_ptr, nelems[0], MPI_DOUBLE, NULL); ERR
            /*for(int k=0;k<nD1_rec_2d_vars;k++){
                printf("Rank: %d, start=[%lld, %lld], count = [%lld, %lld]\n", rank, starts_D1[k][0], starts_D1[k][1], counts_D1[k][0], counts_D1[k][1]); fflush(stdout);
            }*/
            rec_buf_ptr += nelems[0];
            my_nreqs += nreqs[0];
        }
    }

    /* 4 D6 record variables: double (Time, nCells, nVertLevelsP1) */
    for (rec_no = 0; rec_no < num_recs; rec_no++) {
        for (j = 0; j < nreqs[5]; j++) starts_D6[j][0] = rec_no;

        rec_buf_ptr = D6_rec_dbl_buf;
        for (j = 0; j < nD6_rec_3d_vars; j++) {
            err = ncmpi_iget_varn(ncid, D6_rec_3d_varids[j], nreqs[5], starts_D6,
                                  counts_D6, rec_buf_ptr, nelems[5], MPI_DOUBLE, NULL); ERR
            rec_buf_ptr += nelems[5];
            my_nreqs += nreqs[5];
        }
    }

    /* 24 D3 record variables: double (Time, nCells, nVertLevels) */
    for (rec_no = 0; rec_no < num_recs; rec_no++) {
        for (j = 0; j < nreqs[2]; j++) starts_D3[j][0] = rec_no;

        rec_buf_ptr = D3_rec_dbl_buf;
        for (j = 0; j < nD3_rec_3d_vars; j++) {
            err = ncmpi_iget_varn(ncid, D3_rec_3d_varids[j], nreqs[2], starts_D3,
                                  counts_D3, rec_buf_ptr, nelems[2], MPI_DOUBLE, NULL); ERR
            rec_buf_ptr += nelems[2];
            my_nreqs += nreqs[2];
        }
    }

    /* 1 D4 record variable: double (Time, nEdges, nVertLevels) */
    for (rec_no = 0; rec_no < num_recs; rec_no++) {
        for (j = 0; j < nreqs[3]; j++) starts_D4[j][0] = rec_no;

        rec_buf_ptr = D4_rec_dbl_buf;
        for (j = 0; j < nD4_rec_3d_vars; j++) {
            err = ncmpi_iget_varn(ncid, D4_rec_3d_varids[j], nreqs[3], starts_D4,
                                  counts_D4, rec_buf_ptr, nelems[3], MPI_DOUBLE, NULL); ERR
            rec_buf_ptr += nelems[3];
            my_nreqs += nreqs[3];
        }
    }

    /* 1 D5 record variable: double (Time, nVertices, nVertLevels) */
    for (rec_no = 0; rec_no < num_recs; rec_no++) {
        for (j = 0; j < nreqs[4]; j++) starts_D5[j][0] = rec_no;

        rec_buf_ptr = D5_rec_dbl_buf;
        for (j = 0; j < nD5_rec_3d_vars; j++) {
            err = ncmpi_iget_varn(ncid, D5_rec_3d_varids[j], nreqs[4], starts_D5,
                                  counts_D5, rec_buf_ptr, nelems[4], MPI_DOUBLE, NULL); ERR
            rec_buf_ptr += nelems[4];
            my_nreqs += nreqs[4];
        }
    }
    total_nreqs += my_nreqs;

    post_timing += MPI_Wtime() - timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    err = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL); ERR

    wait_timing += MPI_Wtime() - timing;

    MPI_Barrier(comm); /*-----------------------------------------*/
    timing = MPI_Wtime();

    err = ncmpi_inq_get_size(ncid, &total_size); ERR
    get_size = total_size - metadata_size;
    err = ncmpi_buffer_detach(ncid); ERR
    err = ncmpi_close(ncid); ERR
    close_timing += MPI_Wtime() - timing;

    if (nreqs[0] > 0) {
        free(fix_starts_D1[0]); free(fix_starts_D1);
        ASSIGN_BUF(D1_fix_int_buf); ASSIGN_BUF(D1_fix_dbl_buf);
        free(starts_D1[0]); free(starts_D1);
        ASSIGN_BUF(D1_rec_dbl_buf);
    }

    if (nreqs[1] > 0) {
        free(fix_starts_D2[0]); free(fix_starts_D2);
        ASSIGN_BUF(D2_fix_int_buf);
    }

    if (nreqs[2] > 0) {
        free(fix_starts_D3[0]); free(fix_starts_D3);
        ASSIGN_BUF(D3_fix_int_buf);
        free(starts_D3[0]); free(starts_D3);
        ASSIGN_BUF(D3_rec_dbl_buf);
    }

    if (nreqs[3] > 0) {
        free(fix_starts_D4[0]); free(fix_starts_D4);
        ASSIGN_BUF(D4_fix_int_buf);
        free(starts_D4[0]); free(starts_D4);
        ASSIGN_BUF(D4_rec_dbl_buf);
    }

    if (nreqs[4] > 0) {
        free(fix_starts_D5[0]); free(fix_starts_D5);
        ASSIGN_BUF(D5_fix_int_buf);
        free(starts_D5[0]); free(starts_D5);
        ASSIGN_BUF(D5_rec_dbl_buf);
    }

    if (nreqs[5] > 0) {
        free(starts_D6[0]); free(starts_D6);
        ASSIGN_BUF(D6_rec_dbl_buf);
    }

    free(varids);

    total_timing = MPI_Wtime() - total_timing;

    MPI_Reduce(&total_nreqs,   &max_nreqs,  1, MPI_OFFSET, MPI_MAX, 0, comm);
    MPI_Reduce(&total_nreqs,   &tmp,        1, MPI_OFFSET, MPI_SUM, 0, comm);
    total_nreqs = tmp;
    MPI_Reduce(&get_size,      &tmp,        1, MPI_OFFSET, MPI_SUM, 0, comm);
    get_size = tmp;
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
        printf("Total read amount                 = %.2f MiB = %.2f GiB\n",
               (double)total_size/1048576,(double)total_size/1073741824);
        printf("Total number of requests           = %lld\n",total_nreqs);
        printf("Max number of requests             = %lld\n",max_nreqs);
        printf("Max Time of open + metadata define = %.4f sec\n",open_timing);
        printf("Max Time of I/O preparing          = %.4f sec\n",pre_timing);
        printf("Max Time of ncmpi_iget_varn        = %.4f sec\n",post_timing);
        printf("Max Time of ncmpi_wait_all         = %.4f sec\n",wait_timing);
        printf("Max Time of close                  = %.4f sec\n",close_timing);
        printf("Max Time of TOTAL                  = %.4f sec\n",total_timing);
        printf("I/O bandwidth (open-to-close)      = %.4f MiB/sec\n",
               (double)total_size/1048576.0/total_timing);
        printf("I/O bandwidth (read-only)         = %.4f MiB/sec\n",
               (double)get_size/1048576.0/wait_timing);
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
