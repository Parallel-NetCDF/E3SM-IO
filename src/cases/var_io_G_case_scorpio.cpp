/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 *********************************************************************/

/* Mimic the I/O kernel of E3SM's scorpio module on G case
 * Compared to the canonical method used by PnetCDF, the blob strategy used by
 *
 * The file written by scorpio blob strategy differs in the canonical layout
 * in the following points:
 *
 * Scorpio saves the decomposition maps as global variables.
 * Each decomposition map contains 3 attributes - "dimlen" (size of the dimensions),
 * "ndims" (dimensionality), and "piotype" (datatype in PIO type enum).
 * The decomposition variables and their attributes are written by all processes.
 *
 * E3SM dimensions are stored as scalar variables of type uint64_t.
 * The dimensions are "chars", "ilev", "lev", "nbnd", "ncol", and "time".
 * The time dimension is not used and is always 0.
 *
 * All E3SM variables are stored as 1-dimensional ADIOS2 local variables.
 *
 * For each E3SM variable that is associated with a decomposition map, Scorpio
 * adds additional scalar variables and attributes.
 * The variables are "decomp_id" (ID of the associated decomposition map),
 * and frame_id (time steps). If the variable is of double type, there is an
 * additional variable "fillval_id" attached.
 * When a process writes a block of data to the E3SM variable it will also write
 * a block (cell) to the associated scalar variables.
 * The attributes are "decomp" (decomposition map ID in string representation), "dims"
 * (name of the dimensions as string array), "ncop" (string representation of the
 * type of NetCDF operation performed to write the variable), "nctype" (datatype
 * in NetCDF type enum), "ndims" (dimensionality). The attributes are only written
 * by rank 0.
 *
 * For each small E3SM variable that is not associated with a decomposition map,
 * Scorpio does not introduce additional variables but attaches attributes to them.
 * The attributes are "adiostype" (datatype in ADIOS2 type enum), "ncop" (string
 * representation of the type of NetCDF operation performed to write the variable),
 * "nctype" (datatype in NetCDF type enum), "ndims" (dimensionality), and, if the
 * variable is not a scalar, "dims" (name of the dimensions as string array).
 * The attributes are only written by rank 0.
 * Unlike variables associated with a decomposition map, small variables are stored
 * as byte stream (uint_8) together with their metadata (start and count). The actual
 * type of the data is stored as an attribute. An exception is scalar variables which
 * retains its original type.
 *
 * There is one scalar global variable "Nproc" (number of processes writing the file)
 * written by rank 0.
 *
 * There is one scalar global attribute "fillmode" written by all processes.
 *
 * When sub-filling is enabled, ADIOS2 records data objects in a subfile only when one of
 * the process writing to the subfile writes those data objects. As a result, the data
 * objects written only by rank 0 will only appear in the first subfile.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//
#include <stdio.h>
#include <stdlib.h> /* strtoll() */
#include <string.h> /* strcpy(), strncpy() */

#include <string>
//
#include <unistd.h> /* getopt() unlink() */
//
#include <mpi.h>
//
#include <e3sm_io.h>
#include <e3sm_io_err.h>

#include <e3sm_io_case_scorpio.hpp>
#include <e3sm_io_driver.hpp>
#include <e3sm_io_driver_pnc.hpp>
#ifdef ENABLE_ADIOS2
#include <e3sm_io_driver_adios2.hpp>
#endif

#define IPUT_VAR_DOUBLE(F, D, B, R) \
    e3sm_io_scorpio_write_var (driver, rec_no, F, D, MPI_DOUBLE, B, nb);
#define IPUT_VAR_FLOAT(F, D, B, R) \
    e3sm_io_scorpio_write_var (driver, rec_no, F, D, MPI_FLOAT, B, nb);
#define IPUT_VAR_INT(F, D, B, R)  e3sm_io_scorpio_write_var (driver, rec_no, F, D, MPI_INT, B, nb);
#define IPUT_VAR_CHAR(F, D, B, R) e3sm_io_scorpio_write_var (driver, rec_no, F, D, MPI_CHAR, B, nb);
#define IPUT_VAR1_DOUBLE(F, D, S, B, R) \
    e3sm_io_scorpio_write_var (driver, rec_no, F, D, MPI_DOUBLE, B, nb);
#define IPUT_VAR1_FLOAT(F, D, S, B, R) \
    e3sm_io_scorpio_write_var (driver, rec_no, F, D, MPI_FLOAT, B, nb);
#define IPUT_VAR1_INT(F, D, S, B, R) \
    e3sm_io_scorpio_write_var (driver, rec_no, F, D, MPI_INT, B, nb);
#define IPUT_VAR1_CHAR(F, D, S, B, R) \
    e3sm_io_scorpio_write_var (driver, rec_no, F, D, MPI_CHAR, B, nb);
#define IPUT_VARA_DOUBLE(F, D, S, C, B, R) \
    e3sm_io_scorpio_write_var (driver, rec_no, F, D, MPI_DOUBLE, B, nb);
#define IPUT_VARA_FLOAT(F, D, S, C, B, R) \
    e3sm_io_scorpio_write_var (driver, rec_no, F, D, MPI_FLOAT, B, nb);
#define IPUT_VARA_INT(F, D, S, C, B, R) \
    e3sm_io_scorpio_write_var (driver, rec_no, F, D, MPI_INT, B, nb);
#define IPUT_VARA_CHAR(F, D, S, C, B, R) \
    e3sm_io_scorpio_write_var (driver, rec_no, F, D, MPI_CHAR, B, nb);
#define IPUT_VARN(F, D, N, S, C, B, BC, BT, R) \
    e3sm_io_scorpio_write_var (driver, rec_no, F, D, BT, B, nb);
#define BPUT_VARS_DOUBLE(F, D, S, C, ST, B, R) \
    e3sm_io_scorpio_write_var (driver, rec_no, F, D, MPI_FLOAT, B, nb);
#define BPUT_VARS_CHAR(F, D, S, C, ST, B, R) \
    e3sm_io_scorpio_write_var (driver, rec_no, F, D, MPI_CHAR, B, nb);

#define WAIT_ALL_REQS(F, D, B, R) driver.wait (F);

#define FIX_1D_VAR_STARTS_COUNTS(starts, counts, contig_nreqs, disps, blocklens)      \
    {                                                                                 \
        starts    = (MPI_Offset **)malloc (2 * contig_nreqs * sizeof (MPI_Offset *)); \
        counts    = starts + contig_nreqs;                                            \
        starts[0] = (MPI_Offset *)malloc (2 * contig_nreqs * sizeof (MPI_Offset));    \
        counts[0] = starts[0] + contig_nreqs;                                         \
                                                                                      \
        for (i = 1; i < contig_nreqs; i++) {                                          \
            starts[i] = starts[i - 1] + 1;                                            \
            counts[i] = counts[i - 1] + 1;                                            \
        }                                                                             \
                                                                                      \
        for (i = 0; i < contig_nreqs; i++) {                                          \
            starts[i][0] = disps[i];                                                  \
            counts[i][0] = blocklens[i];                                              \
        }                                                                             \
    }

#define FIX_2D_VAR_STARTS_COUNTS(starts, counts, contig_nreqs, disps, blocklens, last_dimlen)      \
    {                                                                                              \
        starts    = (MPI_Offset **)malloc (2 * contig_nreqs * sizeof (MPI_Offset *));              \
        counts    = starts + contig_nreqs;                                                         \
        starts[0] = (MPI_Offset *)malloc (2 * contig_nreqs * 2 * sizeof (MPI_Offset));             \
        counts[0] = starts[0] + contig_nreqs * 2;                                                  \
                                                                                                   \
        for (i = 1; i < contig_nreqs; i++) {                                                       \
            starts[i] = starts[i - 1] + 2;                                                         \
            counts[i] = counts[i - 1] + 2;                                                         \
        }                                                                                          \
                                                                                                   \
        j            = 0;                                                                          \
        starts[0][0] = disps[0] / last_dimlen;                                                     \
        starts[0][1] = disps[0] % last_dimlen; /* decomposition is 2D */                           \
        counts[0][0] = 1;                                                                          \
        counts[0][1] = blocklens[0]; /* each blocklens[i] is no bigger than last_dimlen */         \
        for (i = 1; i < contig_nreqs; i++) {                                                       \
            MPI_Offset _start[2];                                                                  \
            _start[0] = disps[i] / last_dimlen;                                                    \
            _start[1] = disps[i] % last_dimlen;                                                    \
            if (_start[0] == starts[j][0] + counts[j][0] && _start[1] == starts[j][1] &&           \
                blocklens[i] == counts[j][1])                                                      \
                counts[j][0]++;                                                                    \
            else {                                                                                 \
                j++;                                                                               \
                starts[j][0] = _start[0];                                                          \
                starts[j][1] = _start[1];                                                          \
                counts[j][0] = 1;                                                                  \
                counts[j][1] = blocklens[i]; /* each blocklens[i] is no bigger than last_dimlen */ \
            }                                                                                      \
        }                                                                                          \
    }

#define REC_2D_VAR_STARTS_COUNTS(rec, starts, counts, contig_nreqs, disps, blocklens)  \
    {                                                                                  \
        starts    = (MPI_Offset **)malloc (2 * contig_nreqs * sizeof (MPI_Offset *));  \
        counts    = starts + contig_nreqs;                                             \
        starts[0] = (MPI_Offset *)malloc (2 * contig_nreqs * 2 * sizeof (MPI_Offset)); \
        counts[0] = starts[0] + contig_nreqs * 2;                                      \
                                                                                       \
        for (i = 1; i < contig_nreqs; i++) {                                           \
            starts[i] = starts[i - 1] + 2;                                             \
            counts[i] = counts[i - 1] + 2;                                             \
        }                                                                              \
                                                                                       \
        for (i = 0; i < contig_nreqs; i++) {                                           \
            starts[i][1] = disps[i]; /* decomposition is 1D */                         \
            counts[i][1] = blocklens[i];                                               \
                                                                                       \
            starts[i][0] = rec; /* record ID */                                        \
            counts[i][0] = 1;   /* one record only */                                  \
        }                                                                              \
    }

#define REC_3D_VAR_STARTS_COUNTS(rec, starts, counts, contig_nreqs, disps, blocklens, last_dimlen) \
    {                                                                                              \
        starts    = (MPI_Offset **)malloc (2 * contig_nreqs * sizeof (MPI_Offset *));              \
        counts    = starts + contig_nreqs;                                                         \
        starts[0] = (MPI_Offset *)malloc (2 * contig_nreqs * 3 * sizeof (MPI_Offset));             \
        counts[0] = starts[0] + contig_nreqs * 3;                                                  \
                                                                                                   \
        for (i = 1; i < contig_nreqs; i++) {                                                       \
            starts[i] = starts[i - 1] + 3;                                                         \
            counts[i] = counts[i - 1] + 3;                                                         \
        }                                                                                          \
                                                                                                   \
        j            = 0;                                                                          \
        starts[0][0] = rec; /* record ID */                                                        \
        starts[0][1] = disps[0] / last_dimlen;                                                     \
        starts[0][2] = disps[0] % last_dimlen; /* decomposition is 2D */                           \
        counts[0][0] = 1;                      /* one record only */                               \
        counts[0][1] = 1;                                                                          \
        counts[0][2] = blocklens[0]; /* each blocklens[i] is no bigger than last_dimlen */         \
        for (i = 1; i < contig_nreqs; i++) {                                                       \
            MPI_Offset _start[2];                                                                  \
            _start[0] = disps[i] / last_dimlen;                                                    \
            _start[1] = disps[i] % last_dimlen;                                                    \
            if (starts[j][0] == rec && _start[0] == starts[j][1] + counts[j][1] &&                 \
                _start[1] == starts[j][2] && blocklens[i] == counts[j][2])                         \
                /* this request can be combined into the previous one */                           \
                counts[j][1]++;                                                                    \
            else {                                                                                 \
                /* this request cannot be combined into the previous one */                        \
                j++;                                                                               \
                starts[j][0] = rec;                                                                \
                starts[j][1] = _start[0];                                                          \
                starts[j][2] = _start[1];                                                          \
                counts[j][0] = 1;                                                                  \
                counts[j][1] = 1;                                                                  \
                counts[j][2] = blocklens[i]; /* each blocklens[i] is no bigger than last_dimlen */ \
            }                                                                                      \
        }                                                                                          \
        contig_nreqs = j + 1;                                                                      \
    }

#define FREE_N_NULL(A) \
    if (A##p == NULL) free (A)
#define ASSIGN_BUF(A) *(A##p) = A

/*----< run_varn_G_case() >--------------------------------------------------*/
int run_varn_G_case_scorpio (e3sm_io_config &cfg,
                     e3sm_io_decom &decom,
                     e3sm_io_driver &driver,
                     char           *outfile,
                     int *D1_fix_int_bufp,    /* D1 fix int buffer */
                     int *D2_fix_int_bufp,    /* D2 fix int buffer */
                     int *D3_fix_int_bufp,    /* D3 fix int buffer */
                     int *D4_fix_int_bufp,    /* D4 fix int buffer */
                     int *D5_fix_int_bufp,    /* D5 fix int buffer */
                     double *D1_rec_dbl_bufp, /* D1 rec double buffer */
                     double *D3_rec_dbl_bufp, /* D3 rec double buffer */
                     double *D4_rec_dbl_bufp, /* D4 rec double buffer */
                     double *D5_rec_dbl_bufp, /* D5 rec double buffer */
                     double *D6_rec_dbl_bufp, /* D6 rec double buffer */
                     double *D1_fix_dbl_bufp) /* D1 fix double buffer */
{
    int i, j, err, rank, ncid;
    e3sm_io_scorpio_var *varids;
    int scorpiovars[7];
    int nrecs, rec_no = -1, my_nreqs, *nvars_D;
    size_t ii, rec_buflen, nelems[6];
    double *D1_fix_dbl_buf, *D1_rec_dbl_buf, *D3_rec_dbl_buf, *D4_rec_dbl_buf;
    double *D5_rec_dbl_buf, *D6_rec_dbl_buf, *rec_buf_ptr;
    int    *D1_fix_int_buf, *D2_fix_int_buf, *D3_fix_int_buf, *D4_fix_int_buf;
    int    *D5_fix_int_buf;
    MPI_Offset metadata_size=0, total_size;
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
    MPI_Info info_used       = MPI_INFO_NULL;
    int nD1_rec_2d_vars      = 4;
    int D1_rec_2d_varids[4]  = {0, 4, 7, 38};
    int nD3_rec_3d_vars      = 24;
    int D3_rec_3d_varids[24] = {5,  16, 17, 19, 26, 28, 29, 30, 31, 32, 34, 39,
                                40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51};
    int nD4_rec_3d_vars      = 1;
    int D4_rec_3d_varids[1]  = {6};
    int nD5_rec_3d_vars      = 1;
    int D5_rec_3d_varids[1]  = {18};
    int nD6_rec_3d_vars      = 4;
    int D6_rec_3d_varids[4]  = {1, 2, 3, 27};
    double *dummy_double_buf = NULL;
    char dummy_char_buf[64 + 64];
    int xnreqs[6]; /* number of requests after combination */
    MPI_Offset previous_size;
    std::vector<int> decomids;
    case_meta *pr;

    if (cfg.run_case == F) {
        if (cfg.hist == h0) pr = &cfg.F_case_h0;
        else                pr = &cfg.F_case_h1;
    }
    else if (cfg.run_case == G)
        pr = &cfg.G_case;
    else if (cfg.run_case == I) {
        if (cfg.hist == h0) pr = &cfg.I_case_h0;
        else                pr = &cfg.I_case_h1;
    }
    else return -1;

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    pr->end2end_time = pr->pre_time = MPI_Wtime();

    /* write amount from previous I/O */
    driver.inq_put_size(&previous_size);

    MPI_Comm_rank (cfg.io_comm, &rank);

    nvars_D = pr->nvars_D;

    for (i = 0; i < 6; i++) {
        xnreqs[i]  = decom.contig_nreqs[i];
        nvars_D[i] = 0; /* number of variables using decomposition i */
    }

    /* number of variable elements from 6 decompositions */
    my_nreqs = 0;
    for (i = 0; i < 6; i++) {
        nelems[i] = decom.raw_nreqs[i];
    }

    if (cfg.verbose && rank == 0)
        printf ("nelems=%zd %zd %zd %zd %zd %zd\n", nelems[0], nelems[1], nelems[2], nelems[3],
                nelems[4], nelems[5]);

    /* construct varn API arguments starts[][] and counts[][] */
    if (xnreqs[0] > 0) {
        FIX_1D_VAR_STARTS_COUNTS (fix_starts_D1, fix_counts_D1, xnreqs[0], decom.disps[0],
                                  decom.blocklens[0])
        REC_2D_VAR_STARTS_COUNTS (0, starts_D1, counts_D1, xnreqs[0], decom.disps[0],
                                  decom.blocklens[0])
    } else {
        fix_starts_D1 = NULL;
        fix_counts_D1 = NULL;
        starts_D1     = NULL;
        counts_D1     = NULL;
    }

    if (xnreqs[1] > 0) {
        FIX_1D_VAR_STARTS_COUNTS (fix_starts_D2, fix_counts_D2, xnreqs[1], decom.disps[1],
                                  decom.blocklens[1])
    } else {
        fix_starts_D2 = NULL;
        fix_counts_D2 = NULL;
    }

    if (xnreqs[2] > 0) {
        FIX_2D_VAR_STARTS_COUNTS (fix_starts_D3, fix_counts_D3, xnreqs[2], decom.disps[2],
                                  decom.blocklens[2], decom.dims[2][1])
        REC_3D_VAR_STARTS_COUNTS (0, starts_D3, counts_D3, xnreqs[2], decom.disps[2],
                                  decom.blocklens[2], decom.dims[2][1])
    } else {
        fix_starts_D3 = NULL;
        fix_counts_D3 = NULL;
        starts_D3     = NULL;
        counts_D3     = NULL;
    }

    if (xnreqs[3] > 0) {
        FIX_2D_VAR_STARTS_COUNTS (fix_starts_D4, fix_counts_D4, xnreqs[3], decom.disps[3],
                                  decom.blocklens[3], decom.dims[3][1])
        REC_3D_VAR_STARTS_COUNTS (0, starts_D4, counts_D4, xnreqs[3], decom.disps[3],
                                  decom.blocklens[3], decom.dims[3][1])
    } else {
        fix_starts_D4 = NULL;
        fix_counts_D4 = NULL;
        starts_D4     = NULL;
        counts_D4     = NULL;
    }

    if (xnreqs[4] > 0) {
        FIX_2D_VAR_STARTS_COUNTS (fix_starts_D5, fix_counts_D5, xnreqs[4], decom.disps[4],
                                  decom.blocklens[4], decom.dims[4][1])
        REC_3D_VAR_STARTS_COUNTS (0, starts_D5, counts_D5, xnreqs[4], decom.disps[4],
                                  decom.blocklens[4], decom.dims[4][1])
    } else {
        fix_starts_D5 = NULL;
        fix_counts_D5 = NULL;
        starts_D5     = NULL;
        counts_D5     = NULL;
    }

    if (xnreqs[5] > 0) {
        REC_3D_VAR_STARTS_COUNTS (0, starts_D6, counts_D6, xnreqs[5], decom.disps[5],
                                  decom.blocklens[5], decom.dims[5][1])
    } else {
        starts_D6 = NULL;
        counts_D6 = NULL;
    }

    /* allocate and initialize write buffer for 7 fixed-size variables */
    /* int (nCells): maxLevelCell */
    if (nelems[0] + 64 > 0) {
        if (D1_fix_int_bufp != NULL) {
            D1_fix_int_buf = D1_fix_int_bufp;
        } else {
            D1_fix_int_buf = (int *)malloc (nelems[0] * sizeof (int) + 64);
            for (ii = 0; ii < nelems[0]; ii++) D1_fix_int_buf[ii] = rank + ii;
        }
    } else
        D1_fix_int_buf = NULL;

    /* int (nEdges): maxLevelEdgeTop and maxLevelEdgeBot */
    if (nelems[1] + 64 > 0) {
        if (D2_fix_int_bufp != NULL) {
            D2_fix_int_buf = D2_fix_int_bufp;
        } else {
            D2_fix_int_buf = (int *)malloc (2 * nelems[1] * sizeof (int) + 64);
            for (ii = 0; ii < 2 * nelems[1]; ii++) D2_fix_int_buf[ii] = rank + ii;
        }
    } else
        D2_fix_int_buf = NULL;

    /* int (nCells, nVertLevels): cellMask */
    if (nelems[2] + 64 > 0) {
        if (D3_fix_int_bufp != NULL) {
            D3_fix_int_buf = D3_fix_int_bufp;
        } else {
            D3_fix_int_buf = (int *)malloc (nelems[2] * sizeof (int) + 64);
            for (ii = 0; ii < nelems[2]; ii++) D3_fix_int_buf[ii] = rank + ii;
        }
    } else
        D3_fix_int_buf = NULL;

    /* int (nEdges, nVertLevels): edgeMask */
    if (nelems[3] + 64 > 0) {
        if (D4_fix_int_bufp != NULL) {
            D4_fix_int_buf = D4_fix_int_bufp;
        } else {
            D4_fix_int_buf = (int *)malloc (nelems[3] * sizeof (int) + 64);
            for (ii = 0; ii < nelems[3]; ii++) D4_fix_int_buf[ii] = rank + ii;
        }
    } else
        D4_fix_int_buf = NULL;

    /* int (nVertices, nVertLevels): vertexMask */
    if (nelems[4] + 64 > 0) {
        if (D5_fix_int_bufp != NULL) {
            D5_fix_int_buf = D5_fix_int_bufp;
        } else {
            D5_fix_int_buf = (int *)malloc (nelems[4] * sizeof (int) + 64);
            for (ii = 0; ii < nelems[4]; ii++) D5_fix_int_buf[ii] = rank + ii;
        }
    } else
        D5_fix_int_buf = NULL;

    /* double (nCells): bottomDepth */
    if (nelems[0] > 0) {
        if (D1_fix_dbl_bufp != NULL) {
            D1_fix_dbl_buf = D1_fix_dbl_bufp;
        } else {
            D1_fix_dbl_buf = (double *)malloc (nelems[0] * sizeof (double));
            for (ii = 0; ii < nelems[0]; ii++) D1_fix_dbl_buf[ii] = rank + ii;
        }
    } else
        D1_fix_dbl_buf = NULL;

    /* allocate and initialize write buffer for 34 record variables */
    if (nelems[0] > 0) {
        rec_buflen = nelems[0] * nD1_rec_2d_vars;
        if (D1_rec_dbl_bufp != NULL) {
            D1_rec_dbl_buf = D1_rec_dbl_bufp;
        } else {
            D1_rec_dbl_buf = (double *)malloc (rec_buflen * sizeof (double));
            for (ii = 0; ii < rec_buflen; ii++) D1_rec_dbl_buf[ii] = rank + ii;
        }
    } else
        D1_rec_dbl_buf = NULL;

    if (nelems[2] > 0) {
        rec_buflen = nelems[2] * nD3_rec_3d_vars;
        if (D3_rec_dbl_bufp != NULL) {
            D3_rec_dbl_buf = D3_rec_dbl_bufp;
        } else {
            D3_rec_dbl_buf = (double *)malloc (rec_buflen * sizeof (double));
            for (ii = 0; ii < rec_buflen; ii++) D3_rec_dbl_buf[ii] = rank + ii;
        }
    } else
        D3_rec_dbl_buf = NULL;

    if (nelems[3] > 0) {
        rec_buflen = nelems[3] * nD4_rec_3d_vars;
        if (D4_rec_dbl_bufp != NULL) {
            D4_rec_dbl_buf = D4_rec_dbl_bufp;
        } else {
            D4_rec_dbl_buf = (double *)malloc (rec_buflen * sizeof (double));
            for (ii = 0; ii < rec_buflen; ii++) D4_rec_dbl_buf[ii] = rank + ii;
        }
    } else
        D4_rec_dbl_buf = NULL;

    if (nelems[4] > 0) {
        rec_buflen = nelems[4] * nD5_rec_3d_vars;
        if (D5_rec_dbl_bufp != NULL) {
            D5_rec_dbl_buf = D5_rec_dbl_bufp;
        } else {
            D5_rec_dbl_buf = (double *)malloc (rec_buflen * sizeof (double));
            for (ii = 0; ii < rec_buflen; ii++) D5_rec_dbl_buf[ii] = rank + ii;
        }
    } else
        D5_rec_dbl_buf = NULL;

    if (nelems[5] > 0) {
        rec_buflen = nelems[5] * nD6_rec_3d_vars;
        if (D6_rec_dbl_bufp != NULL) {
            D6_rec_dbl_buf = D6_rec_dbl_bufp;
        } else {
            D6_rec_dbl_buf = (double *)malloc (rec_buflen * sizeof (double));
            for (ii = 0; ii < rec_buflen; ii++) D6_rec_dbl_buf[ii] = rank + ii;
        }
    } else
        D6_rec_dbl_buf = NULL;

    /* initialize write buffer for 11 small variables */
    dummy_double_buf = (double *)malloc ((decom.dims[2][1] + 8) * sizeof (double));
    for (i = 0; i < decom.dims[2][1] + 8; i++) dummy_double_buf[i] = rank + i;
    for (i = 0; i < 64 + 64; i++) dummy_char_buf[i] = 'a' + rank + i;

    varids = (e3sm_io_scorpio_var *)malloc (cfg.nvars * sizeof (e3sm_io_scorpio_var));
    decomids.resize (cfg.nvars);

    // Record decomids
    // For each varialbe associated with a decomposition map,
    // scorpio record its decomposition map ID in an associated scalar variable
    /* int maxLevelEdgeTop(nEdges) */
    decomids[8] = 1;
    /* int maxLevelEdgeBot(nEdges) */
    decomids[37] = 1;
    /* int edgeMask(nEdges, nVertLevels) */
    decomids[10] = 3;
    /* int cellMask(nCells, nVertLevels) */
    decomids[11] = 2;
    /* int vertexMask(nVertices, nVertLevels) */
    decomids[12] = 4;
    /* double bottomDepth(nCells)  */
    decomids[35] = 0;
    /* int maxLevelCell(nCells) */
    decomids[36] = 0;

    /* 11 small variables has no decom id*/
    /* double vertCoordMovementWeights(nVertLevels) */
    decomids[9] = -1;
    /* double refZMid(nVertLevels) */
    decomids[13] = -1;
    /* double refLayerThickness(nVertLevels) */
    decomids[14] = -1;
    /* double refBottomDepth(nVertLevels) */
    decomids[33] = -1;
    /* char xtime(Time, StrLen) */
    decomids[15] = -1;
    /* double areaCellGlobal(Time) */
    decomids[20] = -1;
    /* double areaEdgeGlobal(Time) */
    decomids[21] = -1;
    /* double areaTriangleGlobal(Time) */
    decomids[22] = -1;
    /* double volumeCellGlobal(Time) */
    decomids[23] = -1;
    /* double volumeEdgeGlobal(Time) */
    decomids[24] = -1;
    /* double CFLNumberGlobal(Time) */
    decomids[25] = -1;

    /* 34 record variables */
    /* 4 D1 record variables: double (Time, nCells) */
    for (j = 0; j < nD1_rec_2d_vars; j++) { decomids[D1_rec_2d_varids[j]] = 0; }
    /* 4 D6 record variables: double (Time, nCells, nVertLevelsP1) */
    for (j = 0; j < nD6_rec_3d_vars; j++) { decomids[D6_rec_3d_varids[j]] = 5; }
    /* 24 D3 record variables: double (Time, nCells, nVertLevels) */
    for (j = 0; j < nD3_rec_3d_vars; j++) { decomids[D3_rec_3d_varids[j]] = 2; }
    /* 1 D4 record variable: double (Time, nEdges, nVertLevels) */
    for (j = 0; j < nD4_rec_3d_vars; j++) { decomids[D4_rec_3d_varids[j]] = 3; }
    /* 1 D5 record variable: double (Time, nVertices, nVertLevels) */
    for (j = 0; j < nD5_rec_3d_vars; j++) { decomids[D5_rec_3d_varids[j]] = 4; }

    pr->pre_time = MPI_Wtime() - pr->pre_time;

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    pr->open_time = MPI_Wtime ();

    /* create a new CDF-5 file for writing */
    err = driver.create (outfile, cfg.io_comm, cfg.info, &ncid);
    CHECK_ERR

    pr->open_time = MPI_Wtime() - pr->open_time;

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    pr->def_time = MPI_Wtime ();

    /* define dimensions, variables, and attributes */
    err = def_G_case_scorpio (cfg, decom, driver, ncid,decomids, varids, scorpiovars);
    CHECK_ERR

    /* exit define mode and enter data mode */
    err = driver.enddef (ncid);
    CHECK_ERR

    /* I/O amount so far */
    err = driver.inq_put_size (&metadata_size);
    CHECK_ERR
    err = driver.inq_file_info (ncid, &info_used);
    CHECK_ERR

    pr->def_time = MPI_Wtime() - pr->def_time;

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    pr->post_time = MPI_Wtime ();

    // Write pio scalar vars (one time)

    // Nproc only written by rank 0
    if (rank == 0) {
        err = driver.put_varl (ncid, scorpiovars[6], MPI_INT, &(cfg.np), nb);
        CHECK_ERR
    }

    /* write 7 fixed-size variables */

    /* int maxLevelEdgeTop(nEdges) */
    err = IPUT_VARN (ncid, varids[8], xnreqs[1], fix_starts_D2, fix_counts_D2, D2_fix_int_buf,
                     nelems[1], MPI_INT, NULL);
    CHECK_ERR
    my_nreqs += xnreqs[1];
    nvars_D[1]++;

    /* int maxLevelEdgeBot(nEdges) */
    err = IPUT_VARN (ncid, varids[37], xnreqs[1], fix_starts_D2, fix_counts_D2,
                     D2_fix_int_buf + nelems[1], nelems[1], MPI_INT, NULL);
    CHECK_ERR
    my_nreqs += xnreqs[1];
    nvars_D[1]++;

    /* int edgeMask(nEdges, nVertLevels) */
    err = IPUT_VARN (ncid, varids[10], xnreqs[3], fix_starts_D4, fix_counts_D4, D4_fix_int_buf,
                     nelems[3], MPI_INT, NULL);
    CHECK_ERR
    my_nreqs += xnreqs[3];
    nvars_D[3]++;

    /* int cellMask(nCells, nVertLevels) */
    err = IPUT_VARN (ncid, varids[11], xnreqs[2], fix_starts_D3, fix_counts_D3, D3_fix_int_buf,
                     nelems[2], MPI_INT, NULL);
    CHECK_ERR
    my_nreqs += xnreqs[2];
    nvars_D[2]++;

    /* int vertexMask(nVertices, nVertLevels) */
    err = IPUT_VARN (ncid, varids[12], xnreqs[4], fix_starts_D5, fix_counts_D5, D5_fix_int_buf,
                     nelems[4], MPI_INT, NULL);
    CHECK_ERR
    my_nreqs += xnreqs[4];
    nvars_D[4]++;

    /* double bottomDepth(nCells)  */
    err = IPUT_VARN (ncid, varids[35], xnreqs[0], fix_starts_D1, fix_counts_D1, D1_fix_dbl_buf,
                     nelems[0], MPI_DOUBLE, NULL);
    CHECK_ERR
    my_nreqs += xnreqs[0];
    nvars_D[0]++;

    /* int maxLevelCell(nCells) */
    err = IPUT_VARN (ncid, varids[36], xnreqs[0], fix_starts_D1, fix_counts_D1, D1_fix_int_buf,
                     nelems[0], MPI_INT, NULL);
    CHECK_ERR
    my_nreqs += xnreqs[0];
    nvars_D[0]++;

    /* next 11 small variables are written by rank 0 only */
    if (rank == 0) {
        // count[0] = decom.dims[2][1]; /* dimension nVertLevels */

        /* double vertCoordMovementWeights(nVertLevels) */
        err = BPUT_VARS_DOUBLE (ncid, varids[9], start, count, stride, dummy_double_buf, NULL);
        CHECK_ERR

        /* double refZMid(nVertLevels) */
        err = BPUT_VARS_DOUBLE (ncid, varids[13], start, count, stride, dummy_double_buf, NULL);
        CHECK_ERR

        /* double refLayerThickness(nVertLevels) */
        err = BPUT_VARS_DOUBLE (ncid, varids[14], start, count, stride, dummy_double_buf, NULL);
        CHECK_ERR

        /* double refBottomDepth(nVertLevels) */
        err = BPUT_VARS_DOUBLE (ncid, varids[33], start, count, stride, dummy_double_buf, NULL);
        CHECK_ERR

        my_nreqs += 4; /* 4 non-record variables */
    }

    // Write PIO decom vars, assume there are 6
    for (j = 0; j < 6; j++) {
        int piodecomid[] = {0, 1, 2, 3, 4, 5};

        err = driver.put_varl (ncid, scorpiovars[j], MPI_LONG_LONG,
                                decom.raw_offsets[piodecomid[j]], nb);
        CHECK_ERR
    }

    if (cfg.run_case == F) {
        if (cfg.hist == h0) nrecs = cfg.F_case_h0.nrecs;
        else                nrecs = cfg.F_case_h1.nrecs;
    }
    else if (cfg.run_case == G) nrecs = cfg.G_case.nrecs;
    else if (cfg.run_case == I) {
        if (cfg.hist == h0) nrecs = cfg.I_case_h0.nrecs;
        else                nrecs = cfg.I_case_h1.nrecs;
    }

    for (rec_no = 0; rec_no < nrecs; rec_no++) {
        // start[0] = rec_no;
        // count[0] = 1;
        // count[1] = 64; /* dimension StrLen */

        if (rank == 0) {
            /* char xtime(Time, StrLen) */
            err = BPUT_VARS_CHAR (ncid, varids[15], start, count, stride, dummy_char_buf, NULL);
            CHECK_ERR

            /* double areaCellGlobal(Time) */
            err = BPUT_VARS_DOUBLE (ncid, varids[20], start, count, stride, dummy_double_buf, NULL);
            CHECK_ERR

            /* double areaEdgeGlobal(Time) */
            err = BPUT_VARS_DOUBLE (ncid, varids[21], start, count, stride, dummy_double_buf, NULL);
            CHECK_ERR

            /* double areaTriangleGlobal(Time) */
            err = BPUT_VARS_DOUBLE (ncid, varids[22], start, count, stride, dummy_double_buf, NULL);
            CHECK_ERR

            /* double volumeCellGlobal(Time) */
            err = BPUT_VARS_DOUBLE (ncid, varids[23], start, count, stride, dummy_double_buf, NULL);
            CHECK_ERR

            /* double volumeEdgeGlobal(Time) */
            err = BPUT_VARS_DOUBLE (ncid, varids[24], start, count, stride, dummy_double_buf, NULL);
            CHECK_ERR

            /* double CFLNumberGlobal(Time) */
            err = BPUT_VARS_DOUBLE (ncid, varids[25], start, count, stride, dummy_double_buf, NULL);
            CHECK_ERR

            my_nreqs += 7; /* 7 record variables */
        }

        /* write 34 record variables */

        /* 4 D1 record variables: double (Time, nCells) */
        for (j = 0; j < xnreqs[0]; j++) starts_D1[j][0] = rec_no;
        rec_buf_ptr = D1_rec_dbl_buf;
        for (j = 0; j < nD1_rec_2d_vars; j++) {
            err = IPUT_VARN (ncid, varids[D1_rec_2d_varids[j]], xnreqs[0], starts_D1, counts_D1,
                             rec_buf_ptr, nelems[0], MPI_DOUBLE, NULL);
            CHECK_ERR
            rec_buf_ptr += nelems[0];
            my_nreqs += xnreqs[0];
            if (rec_no == 0) nvars_D[0]++;
        }

        /* 4 D6 record variables: double (Time, nCells, nVertLevelsP1) */
        for (j = 0; j < xnreqs[5]; j++) starts_D6[j][0] = rec_no;
        rec_buf_ptr = D6_rec_dbl_buf;
        for (j = 0; j < nD6_rec_3d_vars; j++) {
            err = IPUT_VARN (ncid, varids[D6_rec_3d_varids[j]], xnreqs[5], starts_D6, counts_D6,
                             rec_buf_ptr, nelems[5], MPI_DOUBLE, NULL);
            CHECK_ERR
            rec_buf_ptr += nelems[5];
            my_nreqs += xnreqs[5];
            if (rec_no == 0) nvars_D[5]++;
        }

        /* 24 D3 record variables: double (Time, nCells, nVertLevels) */
        for (j = 0; j < xnreqs[2]; j++) starts_D3[j][0] = rec_no;
        rec_buf_ptr = D3_rec_dbl_buf;
        for (j = 0; j < nD3_rec_3d_vars; j++) {
            err = IPUT_VARN (ncid, varids[D3_rec_3d_varids[j]], xnreqs[2], starts_D3, counts_D3,
                             rec_buf_ptr, nelems[2], MPI_DOUBLE, NULL);
            CHECK_ERR
            rec_buf_ptr += nelems[2];
            my_nreqs += xnreqs[2];
            if (rec_no == 0) nvars_D[2]++;
        }

        /* 1 D4 record variable: double (Time, nEdges, nVertLevels) */
        for (j = 0; j < xnreqs[3]; j++) starts_D4[j][0] = rec_no;
        rec_buf_ptr = D4_rec_dbl_buf;
        for (j = 0; j < nD4_rec_3d_vars; j++) {
            err = IPUT_VARN (ncid, varids[D4_rec_3d_varids[j]], xnreqs[3], starts_D4, counts_D4,
                             rec_buf_ptr, nelems[3], MPI_DOUBLE, NULL);
            CHECK_ERR
            rec_buf_ptr += nelems[3];
            my_nreqs += xnreqs[3];
            if (rec_no == 0) nvars_D[3]++;
        }

        /* 1 D5 record variable: double (Time, nVertices, nVertLevels) */
        for (j = 0; j < xnreqs[4]; j++) starts_D5[j][0] = rec_no;
        rec_buf_ptr = D5_rec_dbl_buf;
        for (j = 0; j < nD5_rec_3d_vars; j++) {
            err = IPUT_VARN (ncid, varids[D5_rec_3d_varids[j]], xnreqs[4], starts_D5, counts_D5,
                             rec_buf_ptr, nelems[4], MPI_DOUBLE, NULL);
            CHECK_ERR
            rec_buf_ptr += nelems[4];
            my_nreqs += xnreqs[4];
            if (rec_no == 0) nvars_D[4]++;
        }

        pr->post_time = MPI_Wtime() - pr->post_time;

        MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
        pr->flush_time = MPI_Wtime ();

        /* in ADIOS,  data is actually flushed at file close */
        err = driver.wait (ncid);
        CHECK_ERR

        pr->flush_time = MPI_Wtime() - pr->flush_time;
    }

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    pr->close_time = MPI_Wtime ();

    err = driver.close (ncid);
    CHECK_ERR

    if (dummy_double_buf != NULL) free (dummy_double_buf);

    if (xnreqs[0] > 0) {
        free (fix_starts_D1[0]);
        free (fix_starts_D1);
        FREE_N_NULL (D1_fix_int_buf);
        FREE_N_NULL (D1_fix_dbl_buf);
        free (starts_D1[0]);
        free (starts_D1);
        FREE_N_NULL (D1_rec_dbl_buf);
    }

    if (xnreqs[1] > 0) {
        free (fix_starts_D2[0]);
        free (fix_starts_D2);
        FREE_N_NULL (D2_fix_int_buf);
    }

    if (xnreqs[2] > 0) {
        free (fix_starts_D3[0]);
        free (fix_starts_D3);
        FREE_N_NULL (D3_fix_int_buf);
        free (starts_D3[0]);
        free (starts_D3);
        FREE_N_NULL (D3_rec_dbl_buf);
    }

    if (xnreqs[3] > 0) {
        free (fix_starts_D4[0]);
        free (fix_starts_D4);
        FREE_N_NULL (D4_fix_int_buf);
        free (starts_D4[0]);
        free (starts_D4);
        FREE_N_NULL (D4_rec_dbl_buf);
    }

    if (xnreqs[4] > 0) {
        free (fix_starts_D5[0]);
        free (fix_starts_D5);
        FREE_N_NULL (D5_fix_int_buf);
        free (starts_D5[0]);
        free (starts_D5);
        FREE_N_NULL (D5_rec_dbl_buf);
    }

    if (xnreqs[5] > 0) {
        free (starts_D6[0]);
        free (starts_D6);
        FREE_N_NULL (D6_rec_dbl_buf);
    }

    free (varids);

    pr->close_time = MPI_Wtime() - pr->close_time;

    /* obtain the write amount tracked by the driver */
    driver.inq_put_size(&total_size);
    total_size -= previous_size;

    driver.inq_file_size(outfile, &pr->file_size);

    pr->nvars           = cfg.nvars;
    pr->num_flushes     = 1;
    pr->num_decomp_vars = decom.num_decomp;
    pr->my_nreqs        = my_nreqs;
    pr->metadata_WR     = metadata_size;
    pr->amount_WR       = total_size;
    pr->end2end_time    = MPI_Wtime() - pr->end2end_time;

    /* check if there is any PnetCDF internal malloc residue */
    check_malloc(&cfg, &driver);

    /* print MPI-IO hints actually used */
    if (cfg.verbose && rank == 0) print_info(&info_used);

err_out:
    if (info_used != MPI_INFO_NULL) MPI_Info_free (&info_used);
    if (!cfg.keep_outfile) unlink (outfile);
    return err;
}
