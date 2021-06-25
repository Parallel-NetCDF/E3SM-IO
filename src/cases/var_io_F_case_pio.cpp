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
#include <e3sm_io.h>
#include <e3sm_io_err.h>

#include <e3sm_io_case_F.hpp>
#include <e3sm_io_case_F_pio.hpp>
#include <e3sm_io_case_pio.hpp>
#include <e3sm_io_driver.hpp>

#define IPUT_VAR_DOUBLE(F, D, B, R) e3sm_io_pio_write_var (driver, rec_no, F, D, MPI_DOUBLE, B, nb);
#define IPUT_VAR_FLOAT(F, D, B, R)  e3sm_io_pio_write_var (driver, rec_no, F, D, MPI_FLOAT, B, nb);
#define IPUT_VAR_INT(F, D, B, R)    e3sm_io_pio_write_var (driver, rec_no, F, D, MPI_INT, B, nb);
#define IPUT_VAR_CHAR(F, D, B, R)   e3sm_io_pio_write_var (driver, rec_no, F, D, MPI_CHAR, B, nb);
#define IPUT_VAR1_DOUBLE(F, D, S, B, R) \
    e3sm_io_pio_write_var (driver, rec_no, F, D, MPI_DOUBLE, B, nb);
#define IPUT_VAR1_FLOAT(F, D, S, B, R) \
    e3sm_io_pio_write_var (driver, rec_no, F, D, MPI_FLOAT, B, nb);
#define IPUT_VAR1_INT(F, D, S, B, R)  e3sm_io_pio_write_var (driver, rec_no, F, D, MPI_INT, B, nb);
#define IPUT_VAR1_CHAR(F, D, S, B, R) e3sm_io_pio_write_var (driver, rec_no, F, D, MPI_CHAR, B, nb);
#define IPUT_VARA_DOUBLE(F, D, S, C, B, R) \
    e3sm_io_pio_write_var (driver, rec_no, F, D, MPI_DOUBLE, B, nb);
#define IPUT_VARA_FLOAT(F, D, S, C, B, R) \
    e3sm_io_pio_write_var (driver, rec_no, F, D, MPI_FLOAT, B, nb);
#define IPUT_VARA_INT(F, D, S, C, B, R) \
    e3sm_io_pio_write_var (driver, rec_no, F, D, MPI_INT, B, nb);
#define IPUT_VARA_CHAR(F, D, S, C, B, R) \
    e3sm_io_pio_write_var (driver, rec_no, F, D, MPI_CHAR, B, nb);
#define IPUT_VARN(F, D, N, S, C, B, BC, BT, R) \
    e3sm_io_pio_write_var (driver, rec_no, F, D, BT, B, nb);

#define WAIT_ALL_REQS(F, D, B, R) driver.wait (F);

/*----< write_small_vars_F_case() >------------------------------------------*/
static int write_small_vars_F_case_pio (e3sm_io_driver &driver,
                                        int ncid,
                                        int vid, /* starting variable ID */
                                        e3sm_io_pio_var *varids,
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

#define POST_VARN(k, num, vid)                                                           \
    for (j = 0; j < num; j++) {                                                          \
        err = IPUT_VARN (ncid, varids[vid + j], xnreqs[k - 1], starts_D##k, counts_D##k, \
                         rec_buf_ptr, -1, REC_ITYPE, NULL);                              \
        CHECK_ERR                                                                        \
        rec_buf_ptr += nelems[k - 1] + gap;                                              \
        my_nreqs += xnreqs[k - 1];                                                       \
        if (rec_no == 0) nvars_D[k - 1]++;                                               \
    }

#define ASSIGN_DECOMID(k, num, vid) \
    for (j = 0; j < num; j++) { decomids[vid + j] = k - 1; }

/*----< run_varn_F_case() >--------------------------------------------------*/
int run_varn_F_case_pio (e3sm_io_config &cfg,
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
    int i, j, k, err, nerrs = 0, rank, ncid, cmode, nvars_D[3];
    e3sm_io_pio_var *varids;
    int piovars[6];
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
    long long lbuf;
    std::vector<int> decomids;

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    total_timing = pre_timing = MPI_Wtime ();

    open_timing  = 0.0;
    post_timing  = 0.0;
    wait_timing  = 0.0;
    close_timing = 0.0;

    MPI_Comm_rank (cfg.io_comm, &rank);

    if (cfg.non_contig_buf) gap = 10;

    varids = (e3sm_io_pio_var *)malloc (cfg.nvars * sizeof (e3sm_io_pio_var));
    decomids.resize (cfg.nvars);

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

    // Assign decom ID
    i = 0;
    if (xnreqs[1] > 0) {
        /* lat */
        decomids[i++] = 1;
        /* lon */
        decomids[i++] = 1;
    } else
        i += 2;

    /* area */
    if (xnreqs[0] > 0) {
        decomids[i++] = 0;
    } else
        i++;

    /* 27 small variables has no decomposition map */
    for (; i < 30; i++) { decomids[i] = -1; }

    if (cfg.nvars == 414) {
        ASSIGN_DECOMID (2, 1, 30)    /* AEROD_v */
        ASSIGN_DECOMID (3, 2, 31)    /* ANRAIN and ANSNOW */
        ASSIGN_DECOMID (2, 19, 33)   /* AODABS ... AODVIS */
        ASSIGN_DECOMID (3, 2, 52)    /* AQRAIN and AQSNOW */
        ASSIGN_DECOMID (2, 6, 54)    /* AQ_DMS ... AQ_SOAG */
        ASSIGN_DECOMID (3, 4, 60)    /* AREI ... AWNI */
        ASSIGN_DECOMID (2, 4, 64)    /* BURDEN1 ... BURDEN4 */
        ASSIGN_DECOMID (3, 1, 68)    /* CCN3 */
        ASSIGN_DECOMID (2, 2, 69)    /* CDNUMC and CLDHGH */
        ASSIGN_DECOMID (3, 2, 71)    /* CLDICE and CLDLIQ */
        ASSIGN_DECOMID (2, 3, 73)    /* CLDLOW ... CLDTOT */
        ASSIGN_DECOMID (3, 4, 76)    /* CLOUD ... DCQ */
        ASSIGN_DECOMID (2, 11, 80)   /* DF_DMS ... DSTSFMBL */
        ASSIGN_DECOMID (3, 1, 91)    /* DTCOND */
        ASSIGN_DECOMID (2, 2, 92)    /* DTENDTH and DTENDTQ */
        ASSIGN_DECOMID (3, 2, 94)    /* EXTINCT and FICE */
        ASSIGN_DECOMID (2, 7, 96)    /* FLDS ... FLUTC */
        ASSIGN_DECOMID (3, 4, 103)   /* FREQI ... FREQS */
        ASSIGN_DECOMID (2, 15, 107)  /* FSDS ... ICEFRAC */
        ASSIGN_DECOMID (3, 3, 122)   /* ICIMR ... IWC */
        ASSIGN_DECOMID (2, 2, 125)   /* LANDFRAC and LHFLX */
        ASSIGN_DECOMID (3, 4, 127)   /* LINOZ_DO3 ... LINOZ_O3COL */
        ASSIGN_DECOMID (2, 1, 131)   /* LINOZ_SFCSINK */
        ASSIGN_DECOMID (3, 1, 132)   /* LINOZ_SSO3 */
        ASSIGN_DECOMID (2, 3, 133)   /* LINOZ_SZA ... LWCF */
        ASSIGN_DECOMID (3, 12, 136)  /* Mass_bc ... O3 */
        ASSIGN_DECOMID (2, 2, 148)   /* O3_SRF and OCNFRAC */
        ASSIGN_DECOMID (3, 1, 150)   /* OMEGA */
        ASSIGN_DECOMID (2, 1, 151)   /* OMEGA500 */
        ASSIGN_DECOMID (3, 1, 152)   /* OMEGAT */
        ASSIGN_DECOMID (2, 8, 153)   /* PBLH ... PSL */
        ASSIGN_DECOMID (3, 1, 161)   /* Q */
        ASSIGN_DECOMID (2, 2, 162)   /* QFLX and QREFHT */
        ASSIGN_DECOMID (3, 3, 164)   /* QRL ... RAINQM */
        ASSIGN_DECOMID (2, 1, 167)   /* RAM1 */
        ASSIGN_DECOMID (3, 1, 168)   /* RELHUM */
        ASSIGN_DECOMID (2, 37, 169)  /* SFDMS ... SNOWHLND */
        ASSIGN_DECOMID (3, 2, 206)   /* SNOWQM and SO2 */
        ASSIGN_DECOMID (2, 10, 208)  /* SO2_CLXF ... SWCF */
        ASSIGN_DECOMID (3, 1, 218)   /* T */
        ASSIGN_DECOMID (2, 19, 219)  /* TAUGWX ... TVQ */
        ASSIGN_DECOMID (3, 1, 238)   /* U */
        ASSIGN_DECOMID (2, 1, 239)   /* U10 */
        ASSIGN_DECOMID (3, 6, 240)   /* UU ... VV */
        ASSIGN_DECOMID (2, 3, 246)   /* WD_H2O2 ... WD_SO2 */
        ASSIGN_DECOMID (3, 3, 249)   /* WSUB ... aero_water */
        ASSIGN_DECOMID (2, 32, 252)  /* airFV ... dst_c3SFWET */
        ASSIGN_DECOMID (3, 1, 284)   /* hstobie_linoz */
        ASSIGN_DECOMID (2, 129, 285) /* mlip ... soa_c3SFWET */
    } else {
        ASSIGN_DECOMID (2, 13, 30) /* CLDHGH ... T5 */
        ASSIGN_DECOMID (3, 1, 43)  /* U */
        ASSIGN_DECOMID (2, 7, 44)  /* U250 ... Z500 */
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
        err = def_F_case_h0_pio (driver, decom, ncid, decom.dims[2], cfg.nvars, decomids, varids,
                                 piovars);
        CHECK_ERR
    } else {
        /* for h1 file */
        err = def_F_case_h1_pio (driver, decom, ncid, decom.dims[2], cfg.nvars, decomids, varids,
                                 piovars);
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

    // Write pio scalar vars (one time)
    // TODO: only the first subfile contain nproc
    // nproc
    err = driver.put_varl (ncid, piovars[5], MPI_LONG_LONG, &(cfg.np), nb);
    CHECK_ERR

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

        // Write PIO decom vars
        for (j = 0; j < 5; j++) {
            int piodecomid[] = {0, 1, 1, 1, 2};

            err = driver.put_varl (ncid, piovars[j], MPI_LONG_LONG,
                                   decom.raw_offsets[piodecomid[j]], nb);
            CHECK_ERR
        }

        i           = 3;
        dbl_buf_ptr = dbl_buf + nelems[1] * 2 + nelems[0] + gap * 3;
        int_buf_ptr = int_buf;
        txt_buf_ptr = txt_buf;

        /* next 27 small variables are written by rank 0 only */
        if (rank == 0) {
            my_nreqs += 27;
            /* post nonblocking requests using IPUT_VARN() */
            err = write_small_vars_F_case_pio (driver, ncid, i, varids, rec_no, gap,
                                               decom.dims[2][0], decom.dims[2][0] + 1, 2, 8,
                                               &int_buf_ptr, &txt_buf_ptr, &dbl_buf_ptr);
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
                POST_VARN (2, 37, 169)  /* SCO ... SNOWHLND */
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
            printf (
                "heap memory allocated by PnetCDF internally has %lld bytes yet to be "
                "freed\n",
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