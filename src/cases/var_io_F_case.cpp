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
//
#include <stdio.h>
#include <stdlib.h>
#include <string>
//
#include <unistd.h> /* unlink() */
//
#include <mpi.h>
//
#include <e3sm_io_case_F.hpp>
#include <e3sm_io.h>
#include <e3sm_io_err.h>
#include <e3sm_io_driver.hpp>

#define IPUT_VAR_DOUBLE(F, D, B, R)            driver.put_vara (F, D, MPI_DOUBLE, NULL, NULL, B, nb);
#define IPUT_VAR_FLOAT(F, D, B, R)             driver.put_vara (F, D, MPI_FLOAT, NULL, NULL, B, nb);
#define IPUT_VAR_INT(F, D, B, R)               driver.put_vara (F, D, MPI_INT, NULL, NULL, B, nb);
#define IPUT_VAR_CHAR(F, D, B, R)              driver.put_vara (F, D, MPI_CHAR, NULL, NULL, B, nb);
#define IPUT_VAR1_DOUBLE(F, D, S, B, R)        driver.put_vara (F, D, MPI_DOUBLE, S, NULL, B, nb);
#define IPUT_VAR1_FLOAT(F, D, S, B, R)         driver.put_vara (F, D, MPI_FLOAT, S, NULL, B, nb);
#define IPUT_VAR1_INT(F, D, S, B, R)           driver.put_vara (F, D, MPI_INT, S, NULL, B, nb);
#define IPUT_VAR1_CHAR(F, D, S, B, R)          driver.put_vara (F, D, MPI_CHAR, S, NULL, B, nb);
#define IPUT_VARA_DOUBLE(F, D, S, C, B, R)     driver.put_vara (F, D, MPI_DOUBLE, S, C, B, nb);
#define IPUT_VARA_FLOAT(F, D, S, C, B, R)      driver.put_vara (F, D, MPI_FLOAT, S, C, B, nb);
#define IPUT_VARA_INT(F, D, S, C, B, R)        driver.put_vara (F, D, MPI_INT, S, C, B, nb);
#define IPUT_VARA_CHAR(F, D, S, C, B, R)       driver.put_vara (F, D, MPI_CHAR, S, C, B, nb);
#define IPUT_VARN(F, D, N, S, C, B, BC, BT, R) driver.put_varn (F, D, BT, N, S, C, B, nb);

#define PUT_VARD_ALL(F, D, FT, B, N, T) driver.put_vard (F, D, T, N, FT, B, coll);
#define PUT_VARD(F, D, FT, B, N, T)     driver.put_vard (F, D, T, N, FT, B, indep);

#define IGET_VAR_DOUBLE(F, D, B, R)            driver.get_vara (F, D, MPI_DOUBLE, NULL, NULL, B, nb);
#define IGET_VAR_FLOAT(F, D, B, R)             driver.get_vara (F, D, MPI_FLOAT, NULL, NULL, B, nb);
#define IGET_VAR_INT(F, D, B, R)               driver.get_vara (F, D, MPI_INT, NULL, NULL, B, nb);
#define IGET_VAR_CHAR(F, D, B, R)              driver.get_vara (F, D, MPI_CHAR, NULL, NULL, B, nb);
#define IGET_VAR1_DOUBLE(F, D, S, B, R)        driver.get_vara (F, D, MPI_DOUBLE, S, NULL, B, nb);
#define IGET_VAR1_FLOAT(F, D, S, B, R)         driver.get_vara (F, D, MPI_FLOAT, S, NULL, B, nb);
#define IGET_VAR1_INT(F, D, S, B, R)           driver.get_vara (F, D, MPI_INT, S, NULL, B, nb);
#define IGET_VAR1_CHAR(F, D, S, B, R)          driver.get_vara (F, D, MPI_CHAR, S, NULL, B, nb);
#define IGET_VARA_DOUBLE(F, D, S, C, B, R)     driver.get_vara (F, D, MPI_DOUBLE, S, C, B, nb);
#define IGET_VARA_FLOAT(F, D, S, C, B, R)      driver.get_vara (F, D, MPI_FLOAT, S, C, B, nb);
#define IGET_VARA_INT(F, D, S, C, B, R)        driver.get_vara (F, D, MPI_INT, S, C, B, nb);
#define IGET_VARA_CHAR(F, D, S, C, B, R)       driver.get_vara (F, D, MPI_CHAR, S, C, B, nb);
#define IGET_VARN(F, D, N, S, C, B, BC, BT, R) driver.get_varn (F, D, BT, N, S, C, B, nb);

#define GET_VARD_ALL(F, D, FT, B, N, T) driver.put_vard (F, D, T, N, FT, B, coll);
#define GET_VARD(F, D, FT, B, N, T)     driver.put_vard (F, D, T, N, FT, B, indep);

#define WAIT_ALL_REQS(F, D, B, R) driver.wait (F);

/*----< write_small_vars_F_case() >------------------------------------------*/
static int write_small_vars_F_case (e3sm_io_driver &driver,
                                    int ncid,
                                    int vid, /* starting variable ID */
                                    int *varids,
                                    int rec_no,
                                    int gap,
                                    MPI_Offset lev,
                                    MPI_Offset ilev,
                                    MPI_Offset nbnd,
                                    MPI_Offset nchars,
                                    int **int_buf,
                                    char **txt_buf,
                                    double **dbl_buf) {
    int i, err, nerrs = 0;
    MPI_Offset start[2], count[2];

    /* scalar and small variables are written by rank 0 only */
    i = vid;

    if (rec_no == 0) {
        /* lev */
        err = IPUT_VAR_DOUBLE (ncid, varids[i++], *dbl_buf, NULL);
        CHECK_ERR
        *dbl_buf += lev + gap;
        /* hyam */
        err = IPUT_VAR_DOUBLE (ncid, varids[i++], *dbl_buf, NULL);
        CHECK_ERR
        *dbl_buf += lev + gap;
        /* hybm */
        err = IPUT_VAR_DOUBLE (ncid, varids[i++], *dbl_buf, NULL);
        CHECK_ERR
        *dbl_buf += lev + gap;
        /* P0 */
        err = IPUT_VAR_DOUBLE (ncid, varids[i++], *dbl_buf, NULL);
        CHECK_ERR
        *dbl_buf += 1 + gap;
        /* ilev */
        err = IPUT_VAR_DOUBLE (ncid, varids[i++], *dbl_buf, NULL);
        CHECK_ERR
        *dbl_buf += ilev + gap;
        /* hyai */
        err = IPUT_VAR_DOUBLE (ncid, varids[i++], *dbl_buf, NULL);
        CHECK_ERR
        *dbl_buf += ilev + gap;
        /* hybi */
        err = IPUT_VAR_DOUBLE (ncid, varids[i++], *dbl_buf, NULL);
        CHECK_ERR
        *dbl_buf += ilev + gap;
    } else
        i += 7;

    /* time */
    start[0] = rec_no;
    err      = IPUT_VAR1_DOUBLE (ncid, varids[i++], start, *dbl_buf, NULL);
    CHECK_ERR
    *dbl_buf += 1 + gap;
    /* date */
    start[0] = rec_no;
    err      = IPUT_VAR1_INT (ncid, varids[i++], start, *int_buf, NULL);
    CHECK_ERR
    *int_buf += 1;
    /* datesec */
    start[0] = rec_no;
    err      = IPUT_VAR1_INT (ncid, varids[i++], start, *int_buf, NULL);
    CHECK_ERR
    *int_buf += 1;
    /* time_bnds */
    start[0] = rec_no;
    start[1] = 0;
    count[0] = 1;
    count[1] = nbnd;
    err      = IPUT_VARA_DOUBLE (ncid, varids[i++], start, count, *dbl_buf, NULL);
    CHECK_ERR
    *dbl_buf += nbnd + gap;
    /* date_written */
    start[0] = rec_no;
    start[1] = 0;
    count[0] = 1;
    count[1] = nchars;
    err      = IPUT_VARA_CHAR (ncid, varids[i++], start, count, *txt_buf, NULL);
    CHECK_ERR
    *txt_buf += nchars;
    /* time_written */
    start[0] = rec_no;
    start[1] = 0;
    count[0] = 1;
    count[1] = nchars;
    err      = IPUT_VARA_CHAR (ncid, varids[i++], start, count, *txt_buf, NULL);
    CHECK_ERR
    *txt_buf += nchars;

    if (rec_no == 0) {
        /* ndbase */
        err = IPUT_VAR_INT (ncid, varids[i++], *int_buf, NULL);
        CHECK_ERR
        *int_buf += 1;
        /* nsbase */
        err = IPUT_VAR_INT (ncid, varids[i++], *int_buf, NULL);
        CHECK_ERR
        *int_buf += 1;
        /* nbdate */
        err = IPUT_VAR_INT (ncid, varids[i++], *int_buf, NULL);
        CHECK_ERR
        *int_buf += 1;
        /* nbsec */
        err = IPUT_VAR_INT (ncid, varids[i++], *int_buf, NULL);
        CHECK_ERR
        *int_buf += 1;
        /* mdt */
        err = IPUT_VAR_INT (ncid, varids[i++], *int_buf, NULL);
        CHECK_ERR
        *int_buf += 1;
    } else
        i += 5;

    /* ndcur */
    start[0] = rec_no;
    err      = IPUT_VAR1_INT (ncid, varids[i++], start, *int_buf, NULL);
    CHECK_ERR
    *int_buf += 1;
    /* nscur */
    start[0] = rec_no;
    err      = IPUT_VAR1_INT (ncid, varids[i++], start, *int_buf, NULL);
    CHECK_ERR
    *int_buf += 1;
    /* co2vmr */
    start[0] = rec_no;
    err      = IPUT_VAR1_DOUBLE (ncid, varids[i++], start, *dbl_buf, NULL);
    CHECK_ERR
    *dbl_buf += 1 + gap;
    /* ch4vmr */
    start[0] = rec_no;
    err      = IPUT_VAR1_DOUBLE (ncid, varids[i++], start, *dbl_buf, NULL);
    CHECK_ERR
    *dbl_buf += 1 + gap;
    /* n2ovmr */
    start[0] = rec_no;
    err      = IPUT_VAR1_DOUBLE (ncid, varids[i++], start, *dbl_buf, NULL);
    CHECK_ERR
    *dbl_buf += 1 + gap;
    /* f11vmr */
    start[0] = rec_no;
    err      = IPUT_VAR1_DOUBLE (ncid, varids[i++], start, *dbl_buf, NULL);
    CHECK_ERR
    *dbl_buf += 1 + gap;
    /* f12vmr */
    start[0] = rec_no;
    err      = IPUT_VAR1_DOUBLE (ncid, varids[i++], start, *dbl_buf, NULL);
    CHECK_ERR
    *dbl_buf += 1 + gap;
    /* sol_tsi */
    start[0] = rec_no;
    err      = IPUT_VAR1_DOUBLE (ncid, varids[i++], start, *dbl_buf, NULL);
    CHECK_ERR
    *dbl_buf += 1 + gap;
    /* nsteph */
    start[0] = rec_no;
    err      = IPUT_VAR1_INT (ncid, varids[i++], start, *int_buf, NULL);
    CHECK_ERR
    *int_buf += 1;
err_out:
    return err;
}

/*----< read_small_vars_F_case() >------------------------------------------*/
static int read_small_vars_F_case (e3sm_io_driver &driver,
                                   int ncid,
                                   int vid, /* starting variable ID */
                                   int *varids,
                                   int rec_no,
                                   int gap,
                                   MPI_Offset lev,
                                   MPI_Offset ilev,
                                   MPI_Offset nbnd,
                                   MPI_Offset nchars,
                                   int **int_buf,
                                   char **txt_buf,
                                   double **dbl_buf) {
    int i, err, nerrs = 0;
    MPI_Offset start[2], count[2];

    /* scalar and small variables are written by rank 0 only */
    i = vid;

    if (rec_no == 0) {
        /* lev */
        err = IGET_VAR_DOUBLE (ncid, varids[i++], *dbl_buf, NULL);
        CHECK_ERR
        *dbl_buf += lev + gap;
        /* hyam */
        err = IGET_VAR_DOUBLE (ncid, varids[i++], *dbl_buf, NULL);
        CHECK_ERR
        *dbl_buf += lev + gap;
        /* hybm */
        err = IGET_VAR_DOUBLE (ncid, varids[i++], *dbl_buf, NULL);
        CHECK_ERR
        *dbl_buf += lev + gap;
        /* P0 */
        err = IGET_VAR_DOUBLE (ncid, varids[i++], *dbl_buf, NULL);
        CHECK_ERR
        *dbl_buf += 1 + gap;
        /* ilev */
        err = IGET_VAR_DOUBLE (ncid, varids[i++], *dbl_buf, NULL);
        CHECK_ERR
        *dbl_buf += ilev + gap;
        /* hyai */
        err = IGET_VAR_DOUBLE (ncid, varids[i++], *dbl_buf, NULL);
        CHECK_ERR
        *dbl_buf += ilev + gap;
        /* hybi */
        err = IGET_VAR_DOUBLE (ncid, varids[i++], *dbl_buf, NULL);
        CHECK_ERR
        *dbl_buf += ilev + gap;
    } else
        i += 7;

    /* time */
    start[0] = rec_no;
    err      = IGET_VAR1_DOUBLE (ncid, varids[i++], start, *dbl_buf, NULL);
    CHECK_ERR
    *dbl_buf += 1 + gap;
    /* date */
    start[0] = rec_no;
    err      = IGET_VAR1_INT (ncid, varids[i++], start, *int_buf, NULL);
    CHECK_ERR
    *int_buf += 1;
    /* datesec */
    start[0] = rec_no;
    err      = IGET_VAR1_INT (ncid, varids[i++], start, *int_buf, NULL);
    CHECK_ERR
    *int_buf += 1;
    /* time_bnds */
    start[0] = rec_no;
    start[1] = 0;
    count[0] = 1;
    count[1] = nbnd;
    err      = IGET_VARA_DOUBLE (ncid, varids[i++], start, count, *dbl_buf, NULL);
    CHECK_ERR
    *dbl_buf += nbnd + gap;
    /* date_written */
    start[0] = rec_no;
    start[1] = 0;
    count[0] = 1;
    count[1] = nchars;
    err      = IGET_VARA_CHAR (ncid, varids[i++], start, count, *txt_buf, NULL);
    CHECK_ERR
    *txt_buf += nchars;
    /* time_written */
    start[0] = rec_no;
    start[1] = 0;
    count[0] = 1;
    count[1] = nchars;
    err      = IGET_VARA_CHAR (ncid, varids[i++], start, count, *txt_buf, NULL);
    CHECK_ERR
    *txt_buf += nchars;

    if (rec_no == 0) {
        /* ndbase */
        err = IGET_VAR_INT (ncid, varids[i++], *int_buf, NULL);
        CHECK_ERR
        *int_buf += 1;
        /* nsbase */
        err = IGET_VAR_INT (ncid, varids[i++], *int_buf, NULL);
        CHECK_ERR
        *int_buf += 1;
        /* nbdate */
        err = IGET_VAR_INT (ncid, varids[i++], *int_buf, NULL);
        CHECK_ERR
        *int_buf += 1;
        /* nbsec */
        err = IGET_VAR_INT (ncid, varids[i++], *int_buf, NULL);
        CHECK_ERR
        *int_buf += 1;
        /* mdt */
        err = IGET_VAR_INT (ncid, varids[i++], *int_buf, NULL);
        CHECK_ERR
        *int_buf += 1;
    } else
        i += 5;

    /* ndcur */
    start[0] = rec_no;
    err      = IGET_VAR1_INT (ncid, varids[i++], start, *int_buf, NULL);
    CHECK_ERR
    *int_buf += 1;
    /* nscur */
    start[0] = rec_no;
    err      = IGET_VAR1_INT (ncid, varids[i++], start, *int_buf, NULL);
    CHECK_ERR
    *int_buf += 1;
    /* co2vmr */
    start[0] = rec_no;
    err      = IGET_VAR1_DOUBLE (ncid, varids[i++], start, *dbl_buf, NULL);
    CHECK_ERR
    *dbl_buf += 1 + gap;
    /* ch4vmr */
    start[0] = rec_no;
    err      = IGET_VAR1_DOUBLE (ncid, varids[i++], start, *dbl_buf, NULL);
    CHECK_ERR
    *dbl_buf += 1 + gap;
    /* n2ovmr */
    start[0] = rec_no;
    err      = IGET_VAR1_DOUBLE (ncid, varids[i++], start, *dbl_buf, NULL);
    CHECK_ERR
    *dbl_buf += 1 + gap;
    /* f11vmr */
    start[0] = rec_no;
    err      = IGET_VAR1_DOUBLE (ncid, varids[i++], start, *dbl_buf, NULL);
    CHECK_ERR
    *dbl_buf += 1 + gap;
    /* f12vmr */
    start[0] = rec_no;
    err      = IGET_VAR1_DOUBLE (ncid, varids[i++], start, *dbl_buf, NULL);
    CHECK_ERR
    *dbl_buf += 1 + gap;
    /* sol_tsi */
    start[0] = rec_no;
    err      = IGET_VAR1_DOUBLE (ncid, varids[i++], start, *dbl_buf, NULL);
    CHECK_ERR
    *dbl_buf += 1 + gap;
    /* nsteph */
    start[0] = rec_no;
    err      = IGET_VAR1_INT (ncid, varids[i++], start, *int_buf, NULL);
    CHECK_ERR
    *int_buf += 1;
err_out:
    return err;
}

#define SET_TYPE(kind)                                                                \
    {                                                                                 \
        var_types[i] = type[kind];                                                    \
        err          = driver.inq_var_off (ncid, varids[i], &var_offset);             \
        CHECK_ERR                                                                     \
        var_disps[i] = var_offset - offset_rec;                                       \
        if (kind == 2) {                                                              \
            my_nreqs += xnreqs[1];                                                    \
            if (i < cfg.nvars - 1)                                                    \
                buf_disps[i + 1] = buf_disps[i] + (nelems[1] + gap) * sizeof (itype); \
            buf_blocklens[i] = nelems[1];                                             \
        } else { /* kind == 3 */                                                      \
            my_nreqs += xnreqs[2];                                                    \
            if (i < cfg.nvars - 1)                                                    \
                buf_disps[i + 1] = buf_disps[i] + (nelems[2] + gap) * sizeof (itype); \
            buf_blocklens[i] = nelems[2];                                             \
        }                                                                             \
        i++;                                                                          \
    }

#define SET_TYPES(kind, num)                                                          \
    for (j = 0; j < num; j++) {                                                       \
        var_types[i] = type[kind];                                                    \
        err          = driver.inq_var_off (ncid, varids[i], &var_offset);             \
        CHECK_ERR                                                                     \
        var_disps[i] = var_offset - offset_rec;                                       \
        if (kind == 2) {                                                              \
            my_nreqs += xnreqs[1];                                                    \
            if (i < cfg.nvars - 1)                                                    \
                buf_disps[i + 1] = buf_disps[i] + (nelems[1] + gap) * sizeof (itype); \
            buf_blocklens[i] = nelems[1];                                             \
        } else { /* kind == 3 */                                                      \
            my_nreqs += xnreqs[2];                                                    \
            if (i < cfg.nvars - 1)                                                    \
                buf_disps[i + 1] = buf_disps[i] + (nelems[2] + gap) * sizeof (itype); \
            buf_blocklens[i] = nelems[2];                                             \
        }                                                                             \
        i++;                                                                          \
    }

/*----< run_vard_F_case() >--------------------------------------------------*/
int run_vard_F_case (e3sm_io_config &cfg,
                     e3sm_io_decom &decom,
                     e3sm_io_driver &driver,
                     std::string outfile) /* request's block lengths */
{
    std::string targetfname;
    char txt_buf[16], *txt_buf_ptr;
    int i, j, k, err, nerrs                                   = 0, rank, ncid, cmode, *varids;
    int *var_blocklens, *buf_blocklens, my_nreqs, rec_no, gap = 0;
    int int_buf[10], *int_buf_ptr, xnreqs[3];
    size_t fix_buflen, dbl_buflen, rec_buflen;
    size_t nelems[3];
    itype *rec_buf;
    double *dbl_buf, *dbl_buf_ptr;
    double pre_timing, open_timing, io_timing, close_timing;
    double timing, total_timing, max_timing;
    MPI_Aint *var_disps, *buf_disps;
    MPI_Offset tmp, metadata_size, rec_size, put_size, total_size;
    MPI_Offset offset_fix, offset_rec, var_offset, max_nreqs, total_nreqs;
    MPI_Datatype *var_types, type[4], *filetype_rec, filetype_dbl;
    MPI_Datatype buftype_rec, buftype_dbl;
    MPI_Info info_used = MPI_INFO_NULL;

    MPI_Offset malloc_size, sum_size;
    MPI_Offset m_alloc = 0, max_alloc;

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    total_timing = pre_timing = MPI_Wtime ();

    MPI_Comm_rank (cfg.io_comm, &rank);

    if (cfg.non_contig_buf) gap = 10;

    for (i = 0; i < 3; i++) xnreqs[i] = decom.contig_nreqs[i];

    varids = (int *)malloc (cfg.nvars * sizeof (int));

    /* allocate arrays for constructing fileview and buffer type */
    var_types     = (MPI_Datatype *)malloc (cfg.nvars * sizeof (MPI_Datatype));
    var_blocklens = (int *)malloc (cfg.nvars * 2 * sizeof (int));
    buf_blocklens = var_blocklens + cfg.nvars;
    var_disps     = (MPI_Aint *)malloc (cfg.nvars * 2 * sizeof (MPI_Aint));
    buf_disps     = var_disps + cfg.nvars;

    /* define MPI datatypes for 4 kinds from 3 decompositions */
    MPI_Type_indexed (xnreqs[0], decom.blocklens[0], decom.disps[0], MPI_DOUBLE, &type[0]);
    MPI_Type_commit (&type[0]);
    MPI_Type_indexed (xnreqs[1], decom.blocklens[1], decom.disps[1], MPI_DOUBLE, &type[1]);
    MPI_Type_commit (&type[1]);
    MPI_Type_indexed (xnreqs[1], decom.blocklens[1], decom.disps[1], REC_ITYPE, &type[2]);
    MPI_Type_commit (&type[2]);
    MPI_Type_indexed (xnreqs[2], decom.blocklens[2], decom.disps[2], REC_ITYPE, &type[3]);
    MPI_Type_commit (&type[3]);

    /* number of variable elements from 3 decompositions */
    for (i = 0; i < cfg.nvars; i++) var_blocklens[i] = 1;
    for (i = 0; i < 3; i++)
        for (nelems[i] = 0, j = 0; j < xnreqs[i]; j++) nelems[i] += decom.blocklens[i][j];

    if (cfg.verbose && rank == 0) printf ("nelems=%zd %zd %zd\n", nelems[0], nelems[1], nelems[2]);

    /* allocate and initialize write buffer for small variables */
    dbl_buflen = nelems[1] * 2 + nelems[0] + gap * 3 + 3 * decom.dims[2][0] +
                 3 * (decom.dims[2][0] + 1) + 8 + 2 + 20 * gap;

    dbl_buf = (double *)malloc (dbl_buflen * sizeof (double));
    for (i = 0; i < dbl_buflen; i++) dbl_buf[i] = rank + i;

    for (i = 0; i < 10; i++) int_buf[i] = rank + i;
    for (i = 0; i < 16; i++) txt_buf[i] = 'a' + rank + i;

    /* allocate and initialize write buffer for large variables */
    if (cfg.nvars == 414)
        rec_buflen = nelems[1] * 321 + nelems[2] * 63 + (321 + 63) * gap;
    else
        rec_buflen = nelems[1] * 20 + nelems[2] + (20 + 1) * gap;

    rec_buf = (itype *)malloc (rec_buflen * sizeof (itype));
    for (i = 0; i < rec_buflen; i++) rec_buf[i] = rank + i;

    pre_timing = MPI_Wtime () - pre_timing;

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    open_timing = MPI_Wtime ();

    /* set output file name */
    targetfname = std::string (cfg.targetdir) + '/' + outfile;

    /* create a new CDF-5 file for writing */
    err = driver.create (targetfname, cfg.io_comm, cfg.info, &ncid);
    CHECK_ERR

    /* define dimensions, variables, and attributes */
    if (cfg.nvars == 414) {
        /* for h0 file */
        err = def_F_case_h0 (driver, ncid, decom.dims[2], cfg.nvars, varids);
        CHECK_ERR
    } else {
        /* for h1 file */
        err = def_F_case_h1 (driver, ncid, decom.dims[2], cfg.nvars, varids);
        CHECK_ERR
    }

    /* exit define mode and enter data mode */
    err = driver.enddef (ncid);
    CHECK_ERR

    /* I/O amount so far */
    err = driver.inq_put_size (ncid, &metadata_size);
    CHECK_ERR
    err = driver.inq_file_info (ncid, &info_used);
    CHECK_ERR
    open_timing = MPI_Wtime () - open_timing;

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    timing = MPI_Wtime ();

    /* the first 3 variables are of type MPI_DOUBLE -------------------*/
    i        = 0;
    my_nreqs = 0;

    /* lat */
    var_types[i] = type[1];
    err          = driver.inq_var_off (ncid, varids[i], &offset_fix);
    CHECK_ERR
    var_disps[i]     = 0;
    buf_disps[0]     = 0;
    buf_blocklens[0] = nelems[1];
    i++;
    my_nreqs += xnreqs[1];

    /* lon */
    var_types[i] = type[1];
    err          = driver.inq_var_off (ncid, varids[i], &var_offset);
    CHECK_ERR
    var_disps[i]     = var_offset - offset_fix;
    buf_disps[i]     = buf_disps[i - 1] + (nelems[1] + gap) * sizeof (double);
    buf_blocklens[1] = nelems[1];
    i++;
    my_nreqs += xnreqs[1];

    /* area */
    var_types[i] = type[0];
    err          = driver.inq_var_off (ncid, varids[i], &var_offset);
    CHECK_ERR
    var_disps[i]     = var_offset - offset_fix;
    buf_disps[i]     = buf_disps[i - 1] + (nelems[1] + gap) * sizeof (double);
    buf_blocklens[2] = nelems[0];
    i++;
    my_nreqs += xnreqs[0];
    fix_buflen = nelems[1] * 2 + nelems[0] + gap * 3;

    /* skip next 27 small variables */
    i += 27;

    /* concatenate 3 var_types[] into filetype_dbl */
    MPI_Type_create_struct (3, var_blocklens, var_disps, var_types, &filetype_dbl);
    MPI_Type_commit (&filetype_dbl);

    if (cfg.non_contig_buf) {
        /* construct buffer type for 3 variables */
        MPI_Type_create_hindexed (3, buf_blocklens, buf_disps, MPI_DOUBLE, &buftype_dbl);
        MPI_Type_commit (&buftype_dbl);
    } else {
        /* buffer type is contiguous */
        buftype_dbl = MPI_DOUBLE;
    }

    err = driver.inq_var_off (ncid, varids[i], &offset_rec);
    CHECK_ERR
    err = driver.inq_rec_size (ncid, &rec_size);
    CHECK_ERR
    buf_disps[i] = 0;

    if (cfg.nvars == 414) {
        SET_TYPE (2)       /* AEROD_v */
        SET_TYPES (3, 2)   /* ANRAIN and ANSNOW */
        SET_TYPES (2, 19)  /* AODABS ... AODVIS */
        SET_TYPES (3, 2)   /* AQRAIN and AQSNOW */
        SET_TYPES (2, 6)   /* AQ_DMS ... AQ_SOAG */
        SET_TYPES (3, 4)   /* AREI ... AWNI */
        SET_TYPES (2, 4)   /* BURDEN1 ... BURDEN4 */
        SET_TYPE (3)       /* CCN3 */
        SET_TYPES (2, 2)   /* CDNUMC and CLDHGH */
        SET_TYPES (3, 2)   /* CLDICE and CLDLIQ */
        SET_TYPES (2, 3)   /* CLDLOW ... CLDTOT */
        SET_TYPES (3, 4)   /* CLOUD ... DCQ */
        SET_TYPES (2, 11)  /* DF_DMS ... DSTSFMBL */
        SET_TYPE (3)       /* DTCOND */
        SET_TYPES (2, 2)   /* DTENDTH and DTENDTQ */
        SET_TYPES (3, 2)   /* EXTINCT and FICE */
        SET_TYPES (2, 7)   /* FLDS ... FLUTC */
        SET_TYPES (3, 4)   /* FREQI ... FREQS */
        SET_TYPES (2, 15)  /* FSDS ... ICEFRAC */
        SET_TYPES (3, 3)   /* ICIMR ... IWC */
        SET_TYPES (2, 2)   /* LANDFRAC and LHFLX */
        SET_TYPES (3, 4)   /* LINOZ_DO3 ... LINOZ_O3COL */
        SET_TYPE (2)       /* LINOZ_SFCSINK */
        SET_TYPE (3)       /* LINOZ_SSO3 */
        SET_TYPES (2, 3)   /* LINOZ_SZA ... LWCF */
        SET_TYPES (3, 12)  /* Mass_bc ... O3 */
        SET_TYPES (2, 2)   /* O3_SRF and OCNFRAC */
        SET_TYPE (3)       /* OMEGA */
        SET_TYPE (2)       /* OMEGA500 */
        SET_TYPE (3)       /* OMEGAT */
        SET_TYPES (2, 8)   /* PBLH ... PSL */
        SET_TYPE (3)       /* Q */
        SET_TYPES (2, 2)   /* QFLX and QREFHT */
        SET_TYPES (3, 3)   /* QRL ... RAINQM */
        SET_TYPE (2)       /* RAM1 */
        SET_TYPE (3)       /* RELHUM */
        SET_TYPES (2, 37)  /* SFDMS ... SNOWHLND */
        SET_TYPES (3, 2)   /* SNOWQM and SO2 */
        SET_TYPES (2, 10)  /* SO2_CLXF ... SWCF */
        SET_TYPE (3)       /* T */
        SET_TYPES (2, 19)  /* TAUGWX ... TVQ */
        SET_TYPE (3)       /* U */
        SET_TYPE (2)       /* U10 */
        SET_TYPES (3, 6)   /* UU ... VV */
        SET_TYPES (2, 3)   /* WD_H2O2 ... WD_SO2 */
        SET_TYPES (3, 3)   /* WSUB ... aero_water */
        SET_TYPES (2, 32)  /* airFV ... dst_c3SFWET */
        SET_TYPE (3)       /* hstobie_linoz */
        SET_TYPES (2, 129) /* mlip ... soa_c3SFWET */
    } else {
        SET_TYPES (2, 13) /* CLDHGH ... T5 */
        SET_TYPE (3)      /* U */
        SET_TYPES (2, 7)  /* U250 ... Z500 */
    }

    if (cfg.non_contig_buf) {
        /* construct buffer type for record variables */
        MPI_Type_create_hindexed (cfg.nvars - 30, buf_blocklens + 30, buf_disps + 30, REC_ITYPE,
                                  &buftype_rec);
        MPI_Type_commit (&buftype_rec);
    } else {
        /* all record variables are in a single contiguous buffer */
        buftype_rec = REC_ITYPE;
    }

    filetype_rec = (MPI_Datatype *)malloc (cfg.nrec * sizeof (MPI_Datatype));
    for (j = 0; j < cfg.nrec; j++) {
        if (j > 0) {
            for (k = 30; k < cfg.nvars; k++) var_disps[k] += rec_size;
        }
        /* concatenate cfg.nvars-30 var_types[] into filetype_rec[j] */
        MPI_Type_create_struct (cfg.nvars - 30, var_blocklens + 30, var_disps + 30, var_types + 30,
                                filetype_rec + j);
        MPI_Type_commit (filetype_rec + j);
    }

    for (j = 0; j < 4; j++) MPI_Type_free (&type[j]);
    free (var_types);
    free (var_disps);
    free (var_blocklens);

    pre_timing += MPI_Wtime () - timing;

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    io_timing = MPI_Wtime ();

    if (cfg.non_contig_buf)
        fix_buflen = rec_buflen = 1;
    else {
        if (cfg.nvars == 414)
            rec_buflen = nelems[1] * 321 + nelems[2] * 63;
        else
            rec_buflen = nelems[1] * 20 + nelems[2];
    }

    /* write first 3 MPI_DOUBLE fixed-size variables in one vard call */
    err = PUT_VARD_ALL (ncid, varids[0], filetype_dbl, dbl_buf, fix_buflen, buftype_dbl);
    CHECK_ERR

    for (rec_no = 0; rec_no < cfg.nrec; rec_no++) {
        i           = 3;
        dbl_buf_ptr = dbl_buf + nelems[1] * 2 + nelems[0] + gap * 3;
        int_buf_ptr = int_buf;
        txt_buf_ptr = txt_buf;

        /* next 27 small variables are written by rank 0 only */
        if (rank == 0) {
            /* post nonblocking requests using IPUT_VARN() */
            err = write_small_vars_F_case (driver, ncid, i, varids, rec_no, gap, decom.dims[2][0],
                                           decom.dims[2][0] + 1, 2, 8, &int_buf_ptr, &txt_buf_ptr,
                                           &dbl_buf_ptr);
            CHECK_ERR
            my_nreqs += 27;
        }
        i += 27;

        err = WAIT_ALL_REQS (ncid, NC_REQ_ALL, NULL, NULL);
        CHECK_ERR

        /* write remaining record variables in one vard call */
        err =
            PUT_VARD_ALL (ncid, varids[30], filetype_rec[rec_no], rec_buf, rec_buflen, buftype_rec);
        CHECK_ERR
    }
    io_timing = MPI_Wtime () - io_timing;

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    close_timing = MPI_Wtime ();

    err = driver.inq_put_size (ncid, &total_size);
    CHECK_ERR
    put_size = total_size - metadata_size;
    err      = driver.close (ncid);
    CHECK_ERR
    close_timing = MPI_Wtime () - close_timing;

    for (j = 0; j < cfg.nrec; j++) MPI_Type_free (filetype_rec + j);
    free (filetype_rec);
    MPI_Type_free (&filetype_dbl);

    if (cfg.non_contig_buf) {
        MPI_Type_free (&buftype_rec);
        MPI_Type_free (&buftype_dbl);
    }

    free (rec_buf);
    free (dbl_buf);
    free (varids);

    total_timing = MPI_Wtime () - total_timing;

    tmp = my_nreqs;
    MPI_Reduce (&tmp, &max_nreqs, 1, MPI_OFFSET, MPI_MAX, 0, cfg.io_comm);
    MPI_Reduce (&tmp, &total_nreqs, 1, MPI_OFFSET, MPI_MAX, 0, cfg.io_comm);
    MPI_Reduce (&put_size, &tmp, 1, MPI_OFFSET, MPI_SUM, 0, cfg.io_comm);
    put_size = tmp;
    MPI_Reduce (&total_size, &tmp, 1, MPI_OFFSET, MPI_SUM, 0, cfg.io_comm);
    total_size = tmp;
    MPI_Reduce (&open_timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, cfg.io_comm);
    open_timing = max_timing;
    MPI_Reduce (&pre_timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, cfg.io_comm);
    pre_timing = max_timing;
    MPI_Reduce (&io_timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, cfg.io_comm);
    io_timing = max_timing;
    MPI_Reduce (&close_timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, cfg.io_comm);
    close_timing = max_timing;
    MPI_Reduce (&total_timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, cfg.io_comm);
    total_timing = max_timing;

    /* check if there is any PnetCDF internal malloc residue */

    err = driver.inq_malloc_size (&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce (&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, cfg.io_comm);
        if (rank == 0 && sum_size > 0) {
            printf ("-----------------------------------------------------------\n");
            printf ("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                    sum_size);
        }
    }

    driver.inq_malloc_max_size (&m_alloc);
    MPI_Reduce (&m_alloc, &max_alloc, 1, MPI_OFFSET, MPI_MAX, 0, cfg.io_comm);

    if (rank == 0) {
        printf ("History output file                = %s\n", outfile.c_str ());
        printf ("MAX heap memory allocated by PnetCDF internally is %.2f MiB\n",
                (float)max_alloc / 1048576);
        printf ("Total number of variables          = %d\n", cfg.nvars);
        printf ("Total write amount                 = %.2f MiB = %.2f GiB\n",
                (double)total_size / 1048576, (double)total_size / 1073741824);
        printf ("Total number of requests           = %lld\n", total_nreqs);
        printf ("Max number of requests             = %lld\n", max_nreqs);
        printf ("Max Time of open + metadata define = %.4f sec\n", open_timing);
        printf ("Max Time of I/O preparing          = %.4f sec\n", pre_timing);
        printf ("Max Time of driver.put_vard        = %.4f sec\n", io_timing);
        printf ("Max Time of close                  = %.4f sec\n", close_timing);
        printf ("Max Time of TOTAL                  = %.4f sec\n", total_timing);
        printf ("I/O bandwidth (open-to-close)      = %.4f MiB/sec\n",
                (double)total_size / 1048576.0 / total_timing);
        printf ("I/O bandwidth (write-only)         = %.4f MiB/sec\n",
                (double)put_size / 1048576.0 / io_timing);
        if (cfg.verbose) print_info (&info_used);
        printf ("-----------------------------------------------------------\n");
    }
err_out:
    if (info_used != MPI_INFO_NULL) MPI_Info_free (&info_used);
    if (!cfg.keep_outfile) unlink (targetfname.c_str ());
    fflush (stdout);
    MPI_Barrier (cfg.io_comm);
    return nerrs;
}

static inline void FIX_1D_VAR_STARTS_COUNTS (
    MPI_Offset **&starts, MPI_Offset **&counts, int nreq, int *disps, int *blocklens) {
    int j;

    starts    = (MPI_Offset **)malloc (2 * nreq * sizeof (MPI_Offset *));
    counts    = starts + nreq;
    starts[0] = (MPI_Offset *)malloc (2 * nreq * sizeof (MPI_Offset));
    counts[0] = starts[0] + nreq;

    for (j = 1; j < nreq; j++) {
        starts[j] = starts[j - 1] + 1;
        counts[j] = counts[j - 1] + 1;
    }

    for (j = 0; j < nreq; j++) {
        starts[j][0] = disps[j];
        counts[j][0] = blocklens[j];
    }
}

static inline void FIX_2D_VAR_STARTS_COUNTS (MPI_Offset **&starts,
                                             MPI_Offset **&counts,
                                             int &nreq,
                                             int *disps,
                                             int *blocklens,
                                             MPI_Offset last_dimlen) {
    int j, k;

    starts    = (MPI_Offset **)malloc (2 * nreq * sizeof (MPI_Offset *));
    counts    = starts + nreq;
    starts[0] = (MPI_Offset *)malloc (2 * nreq * 2 * sizeof (MPI_Offset));
    counts[0] = starts[0] + nreq * 2;

    for (j = 1; j < nreq; j++) {
        starts[j] = starts[j - 1] + 2;
        counts[j] = counts[j - 1] + 2;
    }

    k            = 0;
    starts[0][0] = disps[0] / last_dimlen;
    starts[0][1] = disps[0] % last_dimlen; /* decomposition is 2D */
    counts[0][0] = 1;
    counts[0][1] = blocklens[0]; /* each blocklens[j] is no bigger than last_dimlen */
    for (j = 1; j < nreq; j++) {
        MPI_Offset _start[2];
        _start[0] = disps[j] / last_dimlen;
        _start[1] = disps[j] % last_dimlen;
        if (_start[0] == starts[k][0] + counts[k][0] && _start[1] == starts[k][1] &&
            blocklens[j] == counts[k][1])
            counts[k][0]++;
        else {
            k++;
            starts[k][0] = _start[0];
            starts[k][1] = _start[1];
            counts[k][0] = 1;
            counts[k][1] = blocklens[j]; /* each blocklens[j] is no bigger than last_dimlen */
        }
    }
    nreq = k + 1;
}

static inline void REC_2D_VAR_STARTS_COUNTS (MPI_Offset rec,
                                             MPI_Offset **&starts,
                                             MPI_Offset **&counts,
                                             int nreq,
                                             int *disps,
                                             int *blocklens) {
    int j;

    starts    = (MPI_Offset **)malloc (2 * nreq * sizeof (MPI_Offset *));
    counts    = starts + nreq;
    starts[0] = (MPI_Offset *)malloc (2 * nreq * 2 * sizeof (MPI_Offset));
    counts[0] = starts[0] + nreq * 2;

    for (j = 1; j < nreq; j++) {
        starts[j] = starts[j - 1] + 2;
        counts[j] = counts[j - 1] + 2;
    }

    for (j = 0; j < nreq; j++) {
        starts[j][1] = disps[j]; /* decomposition is 1D */
        counts[j][1] = blocklens[j];

        starts[j][0] = rec; /* record ID */
        counts[j][0] = 1;   /* one record only */
    }
}

static inline void REC_3D_VAR_STARTS_COUNTS (MPI_Offset rec,
                                             MPI_Offset **&starts,
                                             MPI_Offset **&counts,
                                             int &nreq,
                                             int *disps,
                                             int *blocklens,
                                             MPI_Offset last_dimlen) {
    int j, k;

    starts    = (MPI_Offset **)malloc (2 * nreq * sizeof (MPI_Offset *));
    counts    = starts + nreq;
    starts[0] = (MPI_Offset *)malloc (2 * nreq * 3 * sizeof (MPI_Offset));
    counts[0] = starts[0] + nreq * 3;

    for (j = 1; j < nreq; j++) {
        starts[j] = starts[j - 1] + 3;
        counts[j] = counts[j - 1] + 3;
    }

    k            = 0;
    starts[0][0] = rec; /* record ID */
    starts[0][1] = disps[0] / last_dimlen;
    starts[0][2] = disps[0] % last_dimlen; /* decomposition is 2D */
    counts[0][0] = 1;                      /* one record only */
    counts[0][1] = 1;
    counts[0][2] = blocklens[0]; /* each blocklens[j] is no bigger than last_dimlen */
    for (j = 1; j < nreq; j++) {
        MPI_Offset _start[2];
        _start[0] = disps[j] / last_dimlen;
        _start[1] = disps[j] % last_dimlen;
        if (starts[k][0] == rec && _start[0] == starts[k][1] + counts[k][1] &&
            _start[1] == starts[k][2] && blocklens[j] == counts[k][2])
            counts[k][1]++;
        else {
            k++;
            starts[k][0] = rec;
            starts[k][1] = _start[0];
            starts[k][2] = _start[1];
            counts[k][0] = 1;
            counts[k][1] = 1;
            counts[k][2] = blocklens[j]; /* each blocklens[j] is no bigger than last_dimlen */
        }
    }
    nreq = k + 1;
}

#define POST_VARN(k, num, vid)                                                                    \
    for (j = 0; j < num; j++) {                                                                   \
        err = IPUT_VARN (ncid, vid + j, xnreqs[k - 1], starts_D##k, counts_D##k, rec_buf_ptr, -1, \
                         REC_ITYPE, NULL);                                                        \
        CHECK_ERR                                                                                 \
        rec_buf_ptr += nelems[k - 1] + gap;                                                       \
        my_nreqs += xnreqs[k - 1];                                                                \
        if (rec_no == 0) nvars_D[k - 1]++;                                                        \
    }

/*----< run_varn_F_case() >--------------------------------------------------*/
int run_varn_F_case (e3sm_io_config &cfg,
                     e3sm_io_decom &decom,
                     e3sm_io_driver &driver,
                     std::string outfile, /* output file name */
                     double *dbl_bufp,    /* buffer for fixed size double var */
                     itype *rec_bufp,     /* buffer for rec floating point var */
                     char *txt_buf,       /* buffer for char var */
                     int *int_buf)        /* buffer for int var */
{
    std::string targetfname;
    char *txt_buf_ptr;
    int i, j, k, err, nerrs = 0, rank, ncid, cmode, *varids, nvars_D[3];
    int rec_no, gap = 0, my_nreqs, *int_buf_ptr, xnreqs[3];
    size_t dbl_buflen, rec_buflen, nelems[3];
    itype *rec_buf  = NULL, *rec_buf_ptr;
    double *dbl_buf = NULL, *dbl_buf_ptr;
    double pre_timing, open_timing, post_timing, wait_timing, close_timing;
    double timing, total_timing, max_timing;
    MPI_Offset tmp, metadata_size, put_size, total_size, max_nreqs, total_nreqs;
    MPI_Offset **starts_D2 = NULL, **counts_D2 = NULL;
    MPI_Offset **starts_D3 = NULL, **counts_D3 = NULL;
    MPI_Info info_used = MPI_INFO_NULL;
    MPI_Offset malloc_size, sum_size;
    MPI_Offset m_alloc = 0, max_alloc;

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    total_timing = pre_timing = MPI_Wtime ();

    open_timing  = 0.0;
    post_timing  = 0.0;
    wait_timing  = 0.0;
    close_timing = 0.0;

    MPI_Comm_rank (cfg.io_comm, &rank);

    if (cfg.non_contig_buf) gap = 10;

    varids = (int *)malloc (cfg.nvars * sizeof (int));

    for (i = 0; i < 3; i++) {
        xnreqs[i]  = decom.contig_nreqs[i];
        nvars_D[i] = 0; /* number of variables using decomposition i */
    }

    /* calculate number of variable elements from 3 decompositions */
    my_nreqs = 0;
    for (i = 0; i < 3; i++) {
        for (nelems[i] = 0, k = 0; k < xnreqs[i]; k++) nelems[i] += decom.blocklens[i][k];
    }
    if (cfg.verbose && rank == 0) printf ("nelems=%zd %zd %zd\n", nelems[0], nelems[1], nelems[2]);

    /* allocate and initialize write buffer for small variables */
    dbl_buflen = nelems[1] * 2 + nelems[0] + 3 * decom.dims[2][0] + 3 * (decom.dims[2][0] + 1) + 8 +
                 2 + 20 * gap;
    if (dbl_bufp != NULL) {
        dbl_buf = dbl_bufp;
    } else {
        dbl_buf = (double *)malloc (dbl_buflen * sizeof (double));
        for (i = 0; i < dbl_buflen; i++) dbl_buf[i] = rank;
    }

    /* allocate and initialize write buffer for large variables */
    if (cfg.nvars == 414)
        rec_buflen = nelems[1] * 323 + nelems[2] * 63 + (323 + 63) * gap;
    else
        rec_buflen = nelems[1] * 22 + nelems[2] + (20 + 1) * gap;

    if (rec_bufp != NULL) {
        rec_buf = rec_bufp;
    } else {
        rec_buf = (itype *)malloc (rec_buflen * sizeof (itype));

        for (i = 0; i < rec_buflen; i++) rec_buf[i] = rank;
        for (i = 0; i < 10; i++) int_buf[i] = rank;
        for (i = 0; i < 16; i++) txt_buf[i] = 'a' + rank;
    }

    pre_timing = MPI_Wtime () - pre_timing;

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    timing = MPI_Wtime ();

    /* set output file name */
    targetfname = std::string (cfg.targetdir) + '/' + outfile;

    /* create a new CDF-5 file for writing */
    err = driver.create (targetfname, cfg.io_comm, cfg.info, &ncid);
    CHECK_ERR

    /* define dimensions, variables, and attributes */
    if (cfg.nvars == 414) {
        /* for h0 file */
        err = def_F_case_h0 (driver, ncid, decom.dims[2], cfg.nvars, varids);
        CHECK_ERR
    } else {
        /* for h1 file */
        err = def_F_case_h1 (driver, ncid, decom.dims[2], cfg.nvars, varids);
        CHECK_ERR
    }

    /* exit define mode and enter data mode */
    err = driver.enddef (ncid);
    CHECK_ERR

    /* I/O amount so far */
    err = driver.inq_put_size (ncid, &metadata_size);
    CHECK_ERR
    err = driver.inq_file_info (ncid, &info_used);
    CHECK_ERR
    open_timing += MPI_Wtime () - timing;

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    timing = MPI_Wtime ();

    i           = 0;
    dbl_buf_ptr = dbl_buf;

    if (xnreqs[1] > 0) {
        /* lat */
        MPI_Offset **fix_starts_D2, **fix_counts_D2;

        /* construct varn API arguments starts[][] and counts[][] */
        int num = xnreqs[1];
        FIX_1D_VAR_STARTS_COUNTS (fix_starts_D2, fix_counts_D2, num, decom.disps[1],
                                  decom.blocklens[1]);

        REC_2D_VAR_STARTS_COUNTS (0, starts_D2, counts_D2, xnreqs[1], decom.disps[1],
                                  decom.blocklens[1]);

        err = IPUT_VARN (ncid, varids[i++], xnreqs[1], fix_starts_D2, fix_counts_D2, dbl_buf_ptr,
                         nelems[1], MPI_DOUBLE, NULL);
        CHECK_ERR
        dbl_buf_ptr += nelems[1] + gap;
        my_nreqs += xnreqs[1];
        nvars_D[1]++;

        /* lon */
        err = IPUT_VARN (ncid, varids[i++], xnreqs[1], fix_starts_D2, fix_counts_D2, dbl_buf_ptr,
                         nelems[1], MPI_DOUBLE, NULL);
        CHECK_ERR
        dbl_buf_ptr += nelems[1] + gap;
        my_nreqs += xnreqs[1];
        nvars_D[1]++;

        free (fix_starts_D2[0]);
        free (fix_starts_D2);
    } else
        i += 2;

    /* area */
    if (xnreqs[0] > 0) {
        MPI_Offset **fix_starts_D1, **fix_counts_D1;

        /* construct varn API arguments starts[][] and counts[][] */
        FIX_1D_VAR_STARTS_COUNTS (fix_starts_D1, fix_counts_D1, xnreqs[0], decom.disps[0],
                                  decom.blocklens[0]);

        err = IPUT_VARN (ncid, varids[i++], xnreqs[0], fix_starts_D1, fix_counts_D1, dbl_buf_ptr,
                         nelems[0], MPI_DOUBLE, NULL);
        CHECK_ERR
        dbl_buf_ptr += nelems[0] + gap;
        my_nreqs += xnreqs[0];
        nvars_D[0]++;

        free (fix_starts_D1[0]);
        free (fix_starts_D1);
    } else
        i++;

    /* construct varn API arguments starts[][] and counts[][] */
    if (xnreqs[2] > 0)
        REC_3D_VAR_STARTS_COUNTS (0, starts_D3, counts_D3, xnreqs[2], decom.disps[2],
                                  decom.blocklens[2], decom.dims[2][1]);

    post_timing += MPI_Wtime () - timing;

    for (rec_no = 0; rec_no < cfg.nrec; rec_no++) {
        MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
        timing = MPI_Wtime ();

        i           = 3;
        dbl_buf_ptr = dbl_buf + nelems[1] * 2 + nelems[0] + gap * 3;
        int_buf_ptr = int_buf;
        txt_buf_ptr = txt_buf;

        /* next 27 small variables are written by rank 0 only */
        if (rank == 0) {
            my_nreqs += 27;
            /* post nonblocking requests using IPUT_VARN() */
            err = write_small_vars_F_case (driver, ncid, i, varids, rec_no, gap, decom.dims[2][0],
                                           decom.dims[2][0] + 1, 2, 8, &int_buf_ptr, &txt_buf_ptr,
                                           &dbl_buf_ptr);
            CHECK_ERR
        }
        i += 27;

        post_timing += MPI_Wtime () - timing;

        MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
        timing = MPI_Wtime ();

        /* flush fixed-size and small variables */
        err = WAIT_ALL_REQS (ncid, NC_REQ_ALL, NULL, NULL);
        CHECK_ERR

        wait_timing += MPI_Wtime () - timing;

        MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
        timing = MPI_Wtime ();

        rec_buf_ptr = rec_buf;

        for (j = 0; j < xnreqs[1]; j++) starts_D2[j][0] = rec_no;
        for (j = 0; j < xnreqs[2]; j++) starts_D3[j][0] = rec_no;

        if (cfg.nvars == 414) {
            if (cfg.two_buf) {
                /* write 2D variables */
                POST_VARN (2, 1, 30)    /* AEROD_v */
                POST_VARN (2, 19, 33)   /* AODABS ... AODVIS */
                POST_VARN (2, 6, 54)    /* AQ_DMS ... AQ_SOAG */
                POST_VARN (2, 4, 64)    /* BURDEN1 ... BURDEN4 */
                POST_VARN (2, 2, 69)    /* CDNUMC and CLDHGH */
                POST_VARN (2, 3, 73)    /* CLDLOW ... CLDTOT */
                POST_VARN (2, 11, 80)   /* DF_DMS ... DSTSFMBL */
                POST_VARN (2, 2, 92)    /* DTENDTH and DTENDTQ */
                POST_VARN (2, 7, 96)    /* FLDS ... FLUTC */
                POST_VARN (2, 15, 107)  /* FSDS ... ICEFRAC */
                POST_VARN (2, 2, 125)   /* LANDFRAC and LHFLX */
                POST_VARN (2, 1, 131)   /* LINOZ_SFCSINK */
                POST_VARN (2, 3, 133)   /* LINOZ_SZA ... LWCF */
                POST_VARN (2, 2, 148)   /* O3_SRF and OCNFRAC */
                POST_VARN (2, 1, 151)   /* OMEGA500 */
                POST_VARN (2, 8, 153)   /* PBLH ... PSL */
                POST_VARN (2, 2, 162)   /* QFLX and QREFHT */
                POST_VARN (2, 1, 167)   /* RAM1 */
                POST_VARN (2, 37, 169)  /* SFDMS ... SNOWHLND */
                POST_VARN (2, 10, 208)  /* SO2_CLXF ... SWCF */
                POST_VARN (2, 19, 219)  /* TAUGWX ... TVQ */
                POST_VARN (2, 1, 239)   /* U10 */
                POST_VARN (2, 3, 246)   /* WD_H2O2 ... WD_SO2 */
                POST_VARN (2, 32, 252)  /* airFV ... dst_c3SFWET */
                POST_VARN (2, 129, 285) /* mlip ... soa_c3SFWET */
                /* write 3D variables */
                POST_VARN (3, 2, 31)   /* ANRAIN and ANSNOW */
                POST_VARN (3, 2, 52)   /* AQRAIN and AQSNOW */
                POST_VARN (3, 4, 60)   /* AREI ... AWNI */
                POST_VARN (3, 1, 68)   /* CCN3 */
                POST_VARN (3, 2, 71)   /* CLDICE and CLDLIQ */
                POST_VARN (3, 4, 76)   /* CLOUD ... DCQ */
                POST_VARN (3, 1, 91)   /* DTCOND */
                POST_VARN (3, 2, 94)   /* EXTINCT and FICE */
                POST_VARN (3, 4, 103)  /* FREQI ... FREQS */
                POST_VARN (3, 3, 122)  /* ICIMR ... IWC */
                POST_VARN (3, 4, 127)  /* LINOZ_DO3 ... LINOZ_O3COL */
                POST_VARN (3, 1, 132)  /* LINOZ_SSO3 */
                POST_VARN (3, 12, 136) /* Mass_bc ... O3 */
                POST_VARN (3, 1, 150)  /* OMEGA */
                POST_VARN (3, 1, 152)  /* OMEGAT */
                POST_VARN (3, 1, 161)  /* Q */
                POST_VARN (3, 3, 164)  /* QRL ... RAINQM */
                POST_VARN (3, 1, 168)  /* RELHUM */
                POST_VARN (3, 2, 206)  /* SNOWQM and SO2 */
                POST_VARN (3, 1, 218)  /* T */
                POST_VARN (3, 1, 238)  /* U */
                POST_VARN (3, 6, 240)  /* UU ... VV */
                POST_VARN (3, 3, 249)  /* WSUB ... aero_water */
                POST_VARN (3, 1, 284)  /* hstobie_linoz */
            } else {
                /* write variables in the same order as they defined */
                POST_VARN (2, 1, 30)    /* AEROD_v */
                POST_VARN (3, 2, 31)    /* ANRAIN and ANSNOW */
                POST_VARN (2, 19, 33)   /* AODABS ... AODVIS */
                POST_VARN (3, 2, 52)    /* AQRAIN and AQSNOW */
                POST_VARN (2, 6, 54)    /* AQ_DMS ... AQ_SOAG */
                POST_VARN (3, 4, 60)    /* AREI ... AWNI */
                POST_VARN (2, 4, 64)    /* BURDEN1 ... BURDEN4 */
                POST_VARN (3, 1, 68)    /* CCN3 */
                POST_VARN (2, 2, 69)    /* CDNUMC and CLDHGH */
                POST_VARN (3, 2, 71)    /* CLDICE and CLDLIQ */
                POST_VARN (2, 3, 73)    /* CLDLOW ... CLDTOT */
                POST_VARN (3, 4, 76)    /* CLOUD ... DCQ */
                POST_VARN (2, 11, 80)   /* DF_DMS ... DSTSFMBL */
                POST_VARN (3, 1, 91)    /* DTCOND */
                POST_VARN (2, 2, 92)    /* DTENDTH and DTENDTQ */
                POST_VARN (3, 2, 94)    /* EXTINCT and FICE */
                POST_VARN (2, 7, 96)    /* FLDS ... FLUTC */
                POST_VARN (3, 4, 103)   /* FREQI ... FREQS */
                POST_VARN (2, 15, 107)  /* FSDS ... ICEFRAC */
                POST_VARN (3, 3, 122)   /* ICIMR ... IWC */
                POST_VARN (2, 2, 125)   /* LANDFRAC and LHFLX */
                POST_VARN (3, 4, 127)   /* LINOZ_DO3 ... LINOZ_O3COL */
                POST_VARN (2, 1, 131)   /* LINOZ_SFCSINK */
                POST_VARN (3, 1, 132)   /* LINOZ_SSO3 */
                POST_VARN (2, 3, 133)   /* LINOZ_SZA ... LWCF */
                POST_VARN (3, 12, 136)  /* Mass_bc ... O3 */
                POST_VARN (2, 2, 148)   /* O3_SRF and OCNFRAC */
                POST_VARN (3, 1, 150)   /* OMEGA */
                POST_VARN (2, 1, 151)   /* OMEGA500 */
                POST_VARN (3, 1, 152)   /* OMEGAT */
                POST_VARN (2, 8, 153)   /* PBLH ... PSL */
                POST_VARN (3, 1, 161)   /* Q */
                POST_VARN (2, 2, 162)   /* QFLX and QREFHT */
                POST_VARN (3, 3, 164)   /* QRL ... RAINQM */
                POST_VARN (2, 1, 167)   /* RAM1 */
                POST_VARN (3, 1, 168)   /* RELHUM */
                POST_VARN (2, 37, 169)  /* SFDMS ... SNOWHLND */
                POST_VARN (3, 2, 206)   /* SNOWQM and SO2 */
                POST_VARN (2, 10, 208)  /* SO2_CLXF ... SWCF */
                POST_VARN (3, 1, 218)   /* T */
                POST_VARN (2, 19, 219)  /* TAUGWX ... TVQ */
                POST_VARN (3, 1, 238)   /* U */
                POST_VARN (2, 1, 239)   /* U10 */
                POST_VARN (3, 6, 240)   /* UU ... VV */
                POST_VARN (2, 3, 246)   /* WD_H2O2 ... WD_SO2 */
                POST_VARN (3, 3, 249)   /* WSUB ... aero_water */
                POST_VARN (2, 32, 252)  /* airFV ... dst_c3SFWET */
                POST_VARN (3, 1, 284)   /* hstobie_linoz */
                POST_VARN (2, 129, 285) /* mlip ... soa_c3SFWET */
            }
        } else {
            if (cfg.two_buf) {
                /* write 2D variables followed by 3D variables */
                POST_VARN (2, 13, 30) /* CLDHGH ... T5 */
                POST_VARN (2, 7, 44)  /* U250 ... Z500 */
                POST_VARN (3, 1, 43)  /* U */
            } else {
                /* write variables in the same order as they defined */
                POST_VARN (2, 13, 30) /* CLDHGH ... T5 */
                POST_VARN (3, 1, 43)  /* U */
                POST_VARN (2, 7, 44)  /* U250 ... Z500 */
            }
        }

        post_timing += MPI_Wtime () - timing;

        MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
        timing = MPI_Wtime ();

        err = WAIT_ALL_REQS (ncid, NC_REQ_ALL, NULL, NULL);
        CHECK_ERR

        wait_timing += MPI_Wtime () - timing;
    }

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    timing = MPI_Wtime ();

    err = driver.inq_put_size (ncid, &total_size);
    CHECK_ERR
    put_size = total_size - metadata_size;
    err      = driver.close (ncid);
    CHECK_ERR
    close_timing += MPI_Wtime () - timing;

    if (starts_D3 != NULL) {
        free (starts_D3[0]);
        free (starts_D3);
    }
    if (starts_D2 != NULL) {
        free (starts_D2[0]);
        free (starts_D2);
    }
    if (rec_buf != NULL) free (rec_buf);
    if (dbl_buf != NULL) free (dbl_buf);
    free (varids);

    total_timing = MPI_Wtime () - total_timing;

    tmp = my_nreqs;
    MPI_Reduce (&tmp, &max_nreqs, 1, MPI_OFFSET, MPI_MAX, 0, cfg.io_comm);
    MPI_Reduce (&tmp, &total_nreqs, 1, MPI_OFFSET, MPI_SUM, 0, cfg.io_comm);
    MPI_Reduce (&put_size, &tmp, 1, MPI_OFFSET, MPI_SUM, 0, cfg.io_comm);
    put_size = tmp;
    MPI_Reduce (&total_size, &tmp, 1, MPI_OFFSET, MPI_SUM, 0, cfg.io_comm);
    total_size = tmp;
    MPI_Reduce (&open_timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, cfg.io_comm);
    open_timing = max_timing;
    MPI_Reduce (&pre_timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, cfg.io_comm);
    pre_timing = max_timing;
    MPI_Reduce (&post_timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, cfg.io_comm);
    post_timing = max_timing;
    MPI_Reduce (&wait_timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, cfg.io_comm);
    wait_timing = max_timing;
    MPI_Reduce (&close_timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, cfg.io_comm);
    close_timing = max_timing;
    MPI_Reduce (&total_timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, cfg.io_comm);
    total_timing = max_timing;

    /* check if there is any PnetCDF internal malloc residue */
    err = driver.inq_malloc_size (&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce (&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, cfg.io_comm);
        if (rank == 0 && sum_size > 0) {
            printf ("-----------------------------------------------------------\n");
            printf ("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                    sum_size);
        }
    }
    driver.inq_malloc_max_size (&m_alloc);
    MPI_Reduce (&m_alloc, &max_alloc, 1, MPI_OFFSET, MPI_MAX, 0, cfg.io_comm);

    if (rank == 0) {
        int nvars_noD = cfg.nvars;
        for (i = 0; i < 3; i++) nvars_noD -= nvars_D[i];
        printf ("History output file                = %s\n", outfile.c_str ());
        printf ("No. variables use no decomposition = %3d\n", nvars_noD);
        printf ("No. variables use decomposition D1 = %3d\n", nvars_D[0]);
        printf ("No. variables use decomposition D2 = %3d\n", nvars_D[1]);
        printf ("No. variables use decomposition D3 = %3d\n", nvars_D[2]);
        printf ("Total number of variables          = %3d\n", cfg.nvars);
        printf ("MAX heap memory allocated by PnetCDF internally is %.2f MiB\n",
                (float)max_alloc / 1048576);
        printf ("Total write amount                 = %.2f MiB = %.2f GiB\n",
                (double)total_size / 1048576, (double)total_size / 1073741824);
        printf ("Total number of requests           = %lld\n", total_nreqs);
        printf ("Max number of requests             = %lld\n", max_nreqs);
        printf ("Max Time of open + metadata define = %.4f sec\n", open_timing);
        printf ("Max Time of I/O preparing          = %.4f sec\n", pre_timing);
        printf ("Max Time of IPUT_VARN              = %.4f sec\n", post_timing);
        printf ("Max Time of WAIT_ALL_REQS          = %.4f sec\n", wait_timing);
        printf ("Max Time of close                  = %.4f sec\n", close_timing);
        printf ("Max Time of TOTAL                  = %.4f sec\n", total_timing);
        printf ("I/O bandwidth (open-to-close)      = %.4f MiB/sec\n",
                (double)total_size / 1048576.0 / total_timing);
        printf ("I/O bandwidth (write-only)         = %.4f MiB/sec\n",
                (double)put_size / 1048576.0 / wait_timing);
        if (cfg.verbose) print_info (&info_used);
        printf ("-----------------------------------------------------------\n");
    }

err_out:
    if (info_used != MPI_INFO_NULL) MPI_Info_free (&info_used);
    if (!cfg.keep_outfile) unlink (targetfname.c_str ());
    fflush (stdout);
    MPI_Barrier (cfg.io_comm);
    return nerrs;
}

#define POST_VARN_RD(k, num, vid)                                                            \
    for (j = 0; j < num; j++) {                                                              \
        err = IGET_VARN (ncid, vid + j, decom.contig_nreqs[k - 1], starts_D##k, counts_D##k, \
                         rec_buf_ptr, -1, REC_ITYPE, NULL);                                  \
        CHECK_ERR                                                                            \
        rec_buf_ptr += nelems[k - 1] + gap;                                                  \
        my_nreqs += decom.contig_nreqs[k - 1];                                               \
    }

/*----< run_varn_F_case_rd() >--------------------------------------------------*/
int run_varn_F_case_rd (e3sm_io_config &cfg,
                        e3sm_io_decom &decom,
                        e3sm_io_driver &driver,
                        std::string outfile, /* request's block lengths */
                        double **dbl_bufp,   /* buffer for fixed size double var */
                        itype **rec_bufp,    /* buffer for rec floating point var */
                        char *txt_buf,       /* buffer for char var */
                        int *int_buf)        /* buffer for int var */
{
    std::string targetfname;
    char *txt_buf_ptr;
    int i, j, k, err, nerrs = 0, rank, ncid, cmode, *varids, nreqs_D3_merged;
    int rec_no, gap = 0, my_nreqs, *int_buf_ptr;
    size_t dbl_buflen, rec_buflen, nelems[3];
    itype *rec_buf, *rec_buf_ptr;
    double *dbl_buf, *dbl_buf_ptr;
    double pre_timing, open_timing, post_timing, wait_timing, close_timing;
    double timing, total_timing, max_timing;
    MPI_Offset tmp, metadata_size, put_size, total_size, max_nreqs, total_nreqs;
    MPI_Offset **starts_D2 = NULL, **counts_D2 = NULL;
    MPI_Offset **starts_D3 = NULL, **counts_D3 = NULL;
    MPI_Comm comm      = cfg.io_comm;
    MPI_Info info_used = MPI_INFO_NULL;
    MPI_Offset malloc_size, sum_size;
    MPI_Offset m_alloc = 0, max_alloc;

    MPI_Barrier (comm); /*-----------------------------------------*/
    total_timing = pre_timing = MPI_Wtime ();

    open_timing  = 0.0;
    post_timing  = 0.0;
    wait_timing  = 0.0;
    close_timing = 0.0;

    MPI_Comm_rank (comm, &rank);

    if (cfg.non_contig_buf) gap = 10;

    /* calculate number of variable elements from 3 decompositions */
    my_nreqs = 0;
    for (i = 0; i < 3; i++) {
        for (nelems[i] = 0, k = 0; k < decom.contig_nreqs[i]; k++)
            nelems[i] += decom.blocklens[i][k];
    }
    if (cfg.verbose && rank == 0) printf ("nelems=%zd %zd %zd\n", nelems[0], nelems[1], nelems[2]);

    /* allocate and initialize write buffer for small variables */
    dbl_buflen = nelems[1] * 2 + nelems[0] + 3 * decom.dims[2][0] + 3 * (decom.dims[2][0] + 1) + 8 +
                 2 + 20 * gap;

    dbl_buf = (double *)malloc (dbl_buflen * sizeof (double));
    if (dbl_bufp != NULL) { *dbl_bufp = dbl_buf; }

    for (i = 0; i < dbl_buflen; i++) dbl_buf[i] = rank;

    /* allocate and initialize write buffer for large variables */
    if (cfg.nvars == 414)
        rec_buflen = nelems[1] * 321 + nelems[2] * 63 + (321 + 63) * gap;
    else
        rec_buflen = nelems[1] * 20 + nelems[2] + (20 + 1) * gap;

    rec_buf = (itype *)malloc (rec_buflen * sizeof (itype));
    if (rec_bufp != NULL) { *rec_bufp = rec_buf; }

    for (i = 0; i < rec_buflen; i++) rec_buf[i] = rank;

    for (i = 0; i < 10; i++) int_buf[i] = rank;

    for (i = 0; i < 16; i++) txt_buf[i] = 'a' + rank;

    varids = (int *)malloc (cfg.nvars * sizeof (int));

    pre_timing = MPI_Wtime () - pre_timing;

    MPI_Barrier (comm); /*-----------------------------------------*/
    timing = MPI_Wtime ();

    /* set output file name */
    targetfname = std::string (cfg.targetdir) + '/' + outfile;

    /* create a new CDF-5 file for writing */
    err = driver.open (targetfname.c_str (), comm, cfg.info, &ncid);
    CHECK_ERR

    /* define dimensions, variables, and attributes */
    if (cfg.nvars == 414) {
        /* for h0 file */
        err = inq_F_case_h0 (driver, ncid, decom.dims[2], cfg.nvars, varids);
        CHECK_ERR
    } else {
        /* for h1 file */
        err = inq_F_case_h1 (driver, ncid, decom.dims[2], cfg.nvars, varids);
        CHECK_ERR
    }

    /* exit define mode and enter data mode */
    // err = driver.enddef(ncid); CHECK_ERR

    /* I/O amount so far */
    err = driver.inq_get_size (ncid, &metadata_size);
    CHECK_ERR
    err = driver.inq_file_info (ncid, &info_used);
    CHECK_ERR
    open_timing += MPI_Wtime () - timing;

    MPI_Barrier (comm); /*-----------------------------------------*/
    timing = MPI_Wtime ();

    i           = 0;
    dbl_buf_ptr = dbl_buf;

    if (decom.contig_nreqs[1] > 0) {
        /* lat */
        MPI_Offset **fix_starts_D2, **fix_counts_D2;

        /* construct varn API arguments starts[][] and counts[][] */
        int num = decom.contig_nreqs[1];
        FIX_1D_VAR_STARTS_COUNTS (fix_starts_D2, fix_counts_D2, num, decom.disps[1],
                                  decom.blocklens[1]);

        REC_2D_VAR_STARTS_COUNTS (0, starts_D2, counts_D2, decom.contig_nreqs[1], decom.disps[1],
                                  decom.blocklens[1]);

        err = IGET_VARN (ncid, varids[i++], decom.contig_nreqs[1], fix_starts_D2, fix_counts_D2,
                         dbl_buf_ptr, nelems[1], MPI_DOUBLE, NULL);
        CHECK_ERR
        dbl_buf_ptr += nelems[1] + gap;
        my_nreqs += decom.contig_nreqs[1];

        /* lon */
        err = IGET_VARN (ncid, varids[i++], decom.contig_nreqs[1], fix_starts_D2, fix_counts_D2,
                         dbl_buf_ptr, nelems[1], MPI_DOUBLE, NULL);
        CHECK_ERR
        dbl_buf_ptr += nelems[1] + gap;
        my_nreqs += decom.contig_nreqs[1];

        free (fix_starts_D2[0]);
        free (fix_starts_D2);
    } else
        i += 2;

    /* area */
    if (decom.contig_nreqs[0] > 0) {
        MPI_Offset **fix_starts_D1, **fix_counts_D1;

        /* construct varn API arguments starts[][] and counts[][] */
        FIX_1D_VAR_STARTS_COUNTS (fix_starts_D1, fix_counts_D1, decom.contig_nreqs[0],
                                  decom.disps[0], decom.blocklens[0]);

        err = IGET_VARN (ncid, varids[i++], decom.contig_nreqs[0], fix_starts_D1, fix_counts_D1,
                         dbl_buf_ptr, nelems[0], MPI_DOUBLE, NULL);
        CHECK_ERR
        dbl_buf_ptr += nelems[0] + gap;
        my_nreqs += decom.contig_nreqs[0];

        free (fix_starts_D1[0]);
        free (fix_starts_D1);
    } else
        i++;

    /* construct varn API arguments starts[][] and counts[][] */
    if (decom.contig_nreqs[2] > 0)
        REC_3D_VAR_STARTS_COUNTS (0, starts_D3, counts_D3, decom.contig_nreqs[2], decom.disps[2],
                                  decom.blocklens[2], decom.dims[2][1]);

    post_timing += MPI_Wtime () - timing;

    for (rec_no = 0; rec_no < cfg.nrec; rec_no++) {
        MPI_Barrier (comm); /*-----------------------------------------*/
        timing = MPI_Wtime ();

        i           = 3;
        dbl_buf_ptr = dbl_buf + nelems[1] * 2 + nelems[0] + gap * 3;
        int_buf_ptr = int_buf;
        txt_buf_ptr = txt_buf;

        /* next 27 small variables are written by rank 0 only */
        if (rank == 0) {
            my_nreqs += 27;
            /* post nonblocking requests using IPUT_VARN() */
            err = read_small_vars_F_case (driver, ncid, i, varids, rec_no, gap, decom.dims[2][0],
                                          decom.dims[2][0] + 1, 2, 8, &int_buf_ptr, &txt_buf_ptr,
                                          &dbl_buf_ptr);
            CHECK_ERR
        }
        i += 27;

        post_timing += MPI_Wtime () - timing;

        MPI_Barrier (comm); /*-----------------------------------------*/
        timing = MPI_Wtime ();

        /* flush fixed-size and small variables */
        err = WAIT_ALL_REQS (ncid, NC_REQ_ALL, NULL, NULL);
        CHECK_ERR

        wait_timing += MPI_Wtime () - timing;

        MPI_Barrier (comm); /*-----------------------------------------*/
        timing = MPI_Wtime ();

        rec_buf_ptr = rec_buf;

        for (j = 0; j < decom.contig_nreqs[1]; j++) starts_D2[j][0] = rec_no;
        for (j = 0; j < decom.contig_nreqs[2]; j++) starts_D3[j][0] = rec_no;

        if (cfg.nvars == 414) {
            if (cfg.two_buf) {
                /* read 2D variables */
                POST_VARN_RD (2, 1, 30)    /* AEROD_v */
                POST_VARN_RD (2, 19, 33)   /* AODABS ... AODVIS */
                POST_VARN_RD (2, 6, 54)    /* AQ_DMS ... AQ_SOAG */
                POST_VARN_RD (2, 4, 64)    /* BURDEN1 ... BURDEN4 */
                POST_VARN_RD (2, 2, 69)    /* CDNUMC and CLDHGH */
                POST_VARN_RD (2, 3, 73)    /* CLDLOW ... CLDTOT */
                POST_VARN_RD (2, 11, 80)   /* DF_DMS ... DSTSFMBL */
                POST_VARN_RD (2, 2, 92)    /* DTENDTH and DTENDTQ */
                POST_VARN_RD (2, 7, 96)    /* FLDS ... FLUTC */
                POST_VARN_RD (2, 15, 107)  /* FSDS ... ICEFRAC */
                POST_VARN_RD (2, 2, 125)   /* LANDFRAC and LHFLX */
                POST_VARN_RD (2, 1, 131)   /* LINOZ_SFCSINK */
                POST_VARN_RD (2, 3, 133)   /* LINOZ_SZA ... LWCF */
                POST_VARN_RD (2, 2, 148)   /* O3_SRF and OCNFRAC */
                POST_VARN_RD (2, 1, 151)   /* OMEGA500 */
                POST_VARN_RD (2, 8, 153)   /* PBLH ... PSL */
                POST_VARN_RD (2, 2, 162)   /* QFLX and QREFHT */
                POST_VARN_RD (2, 1, 167)   /* RAM1 */
                POST_VARN_RD (2, 37, 169)  /* SFDMS ... SNOWHLND */
                POST_VARN_RD (2, 10, 208)  /* SO2_CLXF ... SWCF */
                POST_VARN_RD (2, 19, 219)  /* TAUGWX ... TVQ */
                POST_VARN_RD (2, 1, 239)   /* U10 */
                POST_VARN_RD (2, 3, 246)   /* WD_H2O2 ... WD_SO2 */
                POST_VARN_RD (2, 32, 252)  /* airFV ... dst_c3SFWET */
                POST_VARN_RD (2, 129, 285) /* mlip ... soa_c3SFWET */
                /* write 3D variables */
                POST_VARN_RD (3, 2, 31)   /* ANRAIN and ANSNOW */
                POST_VARN_RD (3, 2, 52)   /* AQRAIN and AQSNOW */
                POST_VARN_RD (3, 4, 60)   /* AREI ... AWNI */
                POST_VARN_RD (3, 1, 68)   /* CCN3 */
                POST_VARN_RD (3, 2, 71)   /* CLDICE and CLDLIQ */
                POST_VARN_RD (3, 4, 76)   /* CLOUD ... DCQ */
                POST_VARN_RD (3, 1, 91)   /* DTCOND */
                POST_VARN_RD (3, 2, 94)   /* EXTINCT and FICE */
                POST_VARN_RD (3, 4, 103)  /* FREQI ... FREQS */
                POST_VARN_RD (3, 3, 122)  /* ICIMR ... IWC */
                POST_VARN_RD (3, 4, 127)  /* LINOZ_DO3 ... LINOZ_O3COL */
                POST_VARN_RD (3, 1, 132)  /* LINOZ_SSO3 */
                POST_VARN_RD (3, 12, 136) /* Mass_bc ... O3 */
                POST_VARN_RD (3, 1, 150)  /* OMEGA */
                POST_VARN_RD (3, 1, 152)  /* OMEGAT */
                POST_VARN_RD (3, 1, 161)  /* Q */
                POST_VARN_RD (3, 3, 164)  /* QRL ... RAINQM */
                POST_VARN_RD (3, 1, 168)  /* RELHUM */
                POST_VARN_RD (3, 2, 206)  /* SNOWQM and SO2 */
                POST_VARN_RD (3, 1, 218)  /* T */
                POST_VARN_RD (3, 1, 238)  /* U */
                POST_VARN_RD (3, 6, 240)  /* UU ... VV */
                POST_VARN_RD (3, 3, 249)  /* WSUB ... aero_water */
                POST_VARN_RD (3, 1, 284)  /* hstobie_linoz */
            } else {
                /* read variables in the same order as they defined */
                POST_VARN_RD (2, 1, 30)    /* AEROD_v */
                POST_VARN_RD (3, 2, 31)    /* ANRAIN and ANSNOW */
                POST_VARN_RD (2, 19, 33)   /* AODABS ... AODVIS */
                POST_VARN_RD (3, 2, 52)    /* AQRAIN and AQSNOW */
                POST_VARN_RD (2, 6, 54)    /* AQ_DMS ... AQ_SOAG */
                POST_VARN_RD (3, 4, 60)    /* AREI ... AWNI */
                POST_VARN_RD (2, 4, 64)    /* BURDEN1 ... BURDEN4 */
                POST_VARN_RD (3, 1, 68)    /* CCN3 */
                POST_VARN_RD (2, 2, 69)    /* CDNUMC and CLDHGH */
                POST_VARN_RD (3, 2, 71)    /* CLDICE and CLDLIQ */
                POST_VARN_RD (2, 3, 73)    /* CLDLOW ... CLDTOT */
                POST_VARN_RD (3, 4, 76)    /* CLOUD ... DCQ */
                POST_VARN_RD (2, 11, 80)   /* DF_DMS ... DSTSFMBL */
                POST_VARN_RD (3, 1, 91)    /* DTCOND */
                POST_VARN_RD (2, 2, 92)    /* DTENDTH and DTENDTQ */
                POST_VARN_RD (3, 2, 94)    /* EXTINCT and FICE */
                POST_VARN_RD (2, 7, 96)    /* FLDS ... FLUTC */
                POST_VARN_RD (3, 4, 103)   /* FREQI ... FREQS */
                POST_VARN_RD (2, 15, 107)  /* FSDS ... ICEFRAC */
                POST_VARN_RD (3, 3, 122)   /* ICIMR ... IWC */
                POST_VARN_RD (2, 2, 125)   /* LANDFRAC and LHFLX */
                POST_VARN_RD (3, 4, 127)   /* LINOZ_DO3 ... LINOZ_O3COL */
                POST_VARN_RD (2, 1, 131)   /* LINOZ_SFCSINK */
                POST_VARN_RD (3, 1, 132)   /* LINOZ_SSO3 */
                POST_VARN_RD (2, 3, 133)   /* LINOZ_SZA ... LWCF */
                POST_VARN_RD (3, 12, 136)  /* Mass_bc ... O3 */
                POST_VARN_RD (2, 2, 148)   /* O3_SRF and OCNFRAC */
                POST_VARN_RD (3, 1, 150)   /* OMEGA */
                POST_VARN_RD (2, 1, 151)   /* OMEGA500 */
                POST_VARN_RD (3, 1, 152)   /* OMEGAT */
                POST_VARN_RD (2, 8, 153)   /* PBLH ... PSL */
                POST_VARN_RD (3, 1, 161)   /* Q */
                POST_VARN_RD (2, 2, 162)   /* QFLX and QREFHT */
                POST_VARN_RD (3, 3, 164)   /* QRL ... RAINQM */
                POST_VARN_RD (2, 1, 167)   /* RAM1 */
                POST_VARN_RD (3, 1, 168)   /* RELHUM */
                POST_VARN_RD (2, 37, 169)  /* SFDMS ... SNOWHLND */
                POST_VARN_RD (3, 2, 206)   /* SNOWQM and SO2 */
                POST_VARN_RD (2, 10, 208)  /* SO2_CLXF ... SWCF */
                POST_VARN_RD (3, 1, 218)   /* T */
                POST_VARN_RD (2, 19, 219)  /* TAUGWX ... TVQ */
                POST_VARN_RD (3, 1, 238)   /* U */
                POST_VARN_RD (2, 1, 239)   /* U10 */
                POST_VARN_RD (3, 6, 240)   /* UU ... VV */
                POST_VARN_RD (2, 3, 246)   /* WD_H2O2 ... WD_SO2 */
                POST_VARN_RD (3, 3, 249)   /* WSUB ... aero_water */
                POST_VARN_RD (2, 32, 252)  /* airFV ... dst_c3SFWET */
                POST_VARN_RD (3, 1, 284)   /* hstobie_linoz */
                POST_VARN_RD (2, 129, 285) /* mlip ... soa_c3SFWET */
            }
        } else {
            if (cfg.two_buf) {
                /* read 2D variables followed by 3D variables */
                POST_VARN_RD (2, 13, 30) /* CLDHGH ... T5 */
                POST_VARN_RD (2, 7, 44)  /* U250 ... Z500 */
                POST_VARN_RD (3, 1, 43)  /* U */
            } else {
                /* read variables in the same order as they defined */
                POST_VARN_RD (2, 13, 30) /* CLDHGH ... T5 */
                POST_VARN_RD (3, 1, 43)  /* U */
                POST_VARN_RD (2, 7, 44)  /* U250 ... Z500 */
            }
        }

        post_timing += MPI_Wtime () - timing;

        MPI_Barrier (comm); /*-----------------------------------------*/
        timing = MPI_Wtime ();

        err = WAIT_ALL_REQS (ncid, NC_REQ_ALL, NULL, NULL);
        CHECK_ERR

        wait_timing += MPI_Wtime () - timing;
    }

    MPI_Barrier (comm); /*-----------------------------------------*/
    timing = MPI_Wtime ();

    err = driver.inq_get_size (ncid, &total_size);
    CHECK_ERR
    put_size = total_size - metadata_size;
    err      = driver.close (ncid);
    CHECK_ERR
    close_timing += MPI_Wtime () - timing;

    if (starts_D3 != NULL) {
        free (starts_D3[0]);
        free (starts_D3);
    }
    if (starts_D2 != NULL) {
        free (starts_D2[0]);
        free (starts_D2);
    }
    if (rec_bufp == NULL) { free (rec_buf); }
    if (dbl_bufp == NULL) { free (dbl_buf); }
    free (varids);

    total_timing = MPI_Wtime () - total_timing;

    tmp = my_nreqs;
    MPI_Reduce (&tmp, &max_nreqs, 1, MPI_OFFSET, MPI_MAX, 0, comm);
    MPI_Reduce (&tmp, &total_nreqs, 1, MPI_OFFSET, MPI_SUM, 0, comm);
    MPI_Reduce (&put_size, &tmp, 1, MPI_OFFSET, MPI_SUM, 0, comm);
    put_size = tmp;
    MPI_Reduce (&total_size, &tmp, 1, MPI_OFFSET, MPI_SUM, 0, comm);
    total_size = tmp;
    MPI_Reduce (&open_timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    open_timing = max_timing;
    MPI_Reduce (&pre_timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    pre_timing = max_timing;
    MPI_Reduce (&post_timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    post_timing = max_timing;
    MPI_Reduce (&wait_timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    wait_timing = max_timing;
    MPI_Reduce (&close_timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    close_timing = max_timing;
    MPI_Reduce (&total_timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    total_timing = max_timing;

    /* check if there is any PnetCDF internal malloc residue */
    err = driver.inq_malloc_size (&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce (&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, comm);
        if (rank == 0 && sum_size > 0) {
            printf ("-----------------------------------------------------------\n");
            printf ("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                    sum_size);
        }
    }
    driver.inq_malloc_max_size (&m_alloc);
    MPI_Reduce (&m_alloc, &max_alloc, 1, MPI_OFFSET, MPI_MAX, 0, comm);
    if (rank == 0) {
        printf ("History output file                = %s\n", outfile.c_str ());
        printf ("MAX heap memory allocated by PnetCDF internally is %.2f MiB\n",
                (float)max_alloc / 1048576);
        printf ("Total number of variables          = %d\n", cfg.nvars);
        printf ("Total read amount                  = %.2f MiB = %.2f GiB\n",
                (double)total_size / 1048576, (double)total_size / 1073741824);
        printf ("Total number of requests           = %lld\n", total_nreqs);
        printf ("Max number of requests             = %lld\n", max_nreqs);
        printf ("Max Time of open + metadata inquery = %.4f sec\n", open_timing);
        printf ("Max Time of I/O preparing          = %.4f sec\n", pre_timing);
        printf ("Max Time of IGET_VARN              = %.4f sec\n", post_timing);
        printf ("Max Time of WAIT_ALL_REQS          = %.4f sec\n", wait_timing);
        printf ("Max Time of close                  = %.4f sec\n", close_timing);
        printf ("Max Time of TOTAL                  = %.4f sec\n", total_timing);
        printf ("I/O bandwidth (open-to-close)      = %.4f MiB/sec\n",
                (double)total_size / 1048576.0 / total_timing);
        printf ("I/O bandwidth (read-only)          = %.4f MiB/sec\n",
                (double)put_size / 1048576.0 / wait_timing);
        if (cfg.verbose) print_info (&info_used);
        printf ("-----------------------------------------------------------\n");
    }
err_out:
    if (info_used != MPI_INFO_NULL) MPI_Info_free (&info_used);
    if (!cfg.keep_outfile) unlink (targetfname.c_str ());
    fflush (stdout);
    MPI_Barrier (comm);
    return nerrs;
}
