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
#include <string>
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


#define IPUT_VARN(F, D, N, S, C, B, BC, BT, R) driver.put_varn (F, D, BT, N, S, C, B, nb);


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


#define WAIT_ALL_REQS { \
    err = driver.wait(ncid); \
    CHECK_ERR \
    nflushes++; \
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

    if (nreqs != NULL) *nreqs = my_nreqs;

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

#define POST_VARN(k, num, vid)                                              \
    for (j=0; j<num; j++) {                                                 \
        err = IPUT_VARN(ncid, vid+j, xnreqs[k-1], starts_D##k, counts_D##k, \
                        rec_buf_ptr, -1, REC_ITYPE, NULL);                  \
        CHECK_ERR                                                           \
        rec_buf_ptr += nelems[k - 1] + gap;                                 \
        my_nreqs += xnreqs[k - 1];                                          \
        if (rec_no == 0) nvars_D[k - 1]++;                                  \
    }

/*----< run_varn_F_case() >--------------------------------------------------*/
int run_varn_F_case(e3sm_io_config &cfg,
                    e3sm_io_decom &decom,
                    e3sm_io_driver &driver)
{
    char outfile[1040], *ext;
    int i, j, k, err, rank, ncid, *varids=NULL, *nvars_D=cfg.nvars_D;
    int rec_no, gap=0, my_nreqs=0, xnreqs[3];
    int one_flush, nflushes=0;;
    double timing;
    MPI_Offset metadata_size, total_size;
    MPI_Offset **starts_D2 = NULL, **counts_D2 = NULL;
    MPI_Offset **starts_D3 = NULL, **counts_D3 = NULL;
    MPI_Info info_used = MPI_INFO_NULL;
    size_t ii, rec_buflen, nelems[3];
    size_t                 fix_int_buflen, fix_dbl_buflen;
    size_t rec_txt_buflen, rec_int_buflen, rec_dbl_buflen;
    int    *fix_int_buf=NULL, *fix_int_buf_ptr;
    double *fix_dbl_buf=NULL, *fix_dbl_buf_ptr;
    int    *rec_int_buf=NULL, *rec_int_buf_ptr;
    char   *rec_txt_buf=NULL, *rec_txt_buf_ptr;
    double *rec_dbl_buf=NULL, *rec_dbl_buf_ptr;
    itype  *rec_buf=NULL, *rec_buf_ptr;

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    cfg.end2end_time = cfg.pre_time = MPI_Wtime();

    cfg.post_time  = 0.0;
    cfg.flush_time = 0.0;

    MPI_Comm_rank (cfg.io_comm, &rank);

#define FLUSH_ALL_RECORDS_AT_ONCE
#ifdef FLUSH_ALL_RECORDS_AT_ONCE
    one_flush = 1;
#else
    one_flush = 0;
#endif

    if (cfg.non_contig_buf) gap = 10;

    /* allocate space for all variable IDs */
    varids = (int *)malloc (cfg.nvars * sizeof (int));

    for (i = 0; i < 3; i++) {
        xnreqs[i]  = decom.contig_nreqs[i];
        nvars_D[i] = 0; /* number of variables using decomposition i */
    }

    /* calculate number of variable elements from 3 decompositions */
    for (i=0; i<3; i++) {
        for (nelems[i]=0, k=0; k<xnreqs[i]; k++)
             nelems[i] += decom.blocklens[i][k];
    }
    if (cfg.verbose && rank == 0)
        printf ("nelems=%zd %zd %zd\n", nelems[0], nelems[1], nelems[2]);

    /* allocate write buffer for fixed-size variables */
    fix_dbl_buflen = nelems[1] * 2               /* lat[ncol] and lon[ncol] */
                   + nelems[0]                   /* area[ncol] */
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
        rec_buflen = nelems[1] * 321
                   + nelems[2] * 63
                   + (321 + 63) * gap;
    else
        rec_buflen = 13 * nelems[1]   /* CLDHGH ... TS */
                   +  1 * nelems[2]   /* U */
                   +  7 * nelems[1]   /* U250 ... Z500 */
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
    rec_buf     = (itype*)  malloc(rec_buflen     * sizeof(itype));

    for (ii=0; ii<fix_dbl_buflen; ii++) fix_dbl_buf[ii] = rank;
    for (ii=0; ii<fix_int_buflen; ii++) fix_int_buf[ii] = rank;
    for (ii=0; ii<rec_dbl_buflen; ii++) rec_dbl_buf[ii] = rank;
    for (ii=0; ii<rec_txt_buflen; ii++) rec_txt_buf[ii] = 'a' + rank;
    for (ii=0; ii<rec_int_buflen; ii++) rec_int_buf[ii] = rank;
    for (ii=0; ii<rec_buflen;     ii++) rec_buf[ii]     = rank;

    cfg.pre_time = MPI_Wtime() - cfg.pre_time;

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    timing = MPI_Wtime ();

    /* construct h0/h1 file name */
    const char *hist = (cfg.nvars == 414) ? "_h0" :  "_h1";
    ext = strrchr(cfg.out_path, '.');
    if (ext == NULL || (strcmp(ext, ".nc") && strcmp(ext, ".h5")))
        sprintf(outfile, "%s%s", cfg.out_path, hist);
    else {
        sprintf(outfile, "%s", cfg.out_path);
        sprintf(outfile + (ext - cfg.out_path), "%s", hist);
        strcat(outfile, ext);
    }

    /* create a new CDF-5 file for writing */
    err = driver.create (outfile, cfg.io_comm, cfg.info, &ncid);
    CHECK_ERR

    cfg.open_time = MPI_Wtime() - timing;

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    timing = MPI_Wtime ();

    /* define dimensions, variables, and attributes */
    if (cfg.nvars == 414) {
        /* for h0 file */
        err = def_F_case_h0(cfg, decom, driver, ncid, varids);
        CHECK_ERR
    } else {
        /* for h1 file */
        err = def_F_case_h1(cfg, decom, driver, ncid, varids);
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

    cfg.def_time = MPI_Wtime() - timing;

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    timing = MPI_Wtime ();

    /* Now, write climate variables */
    i = 0;
    fix_dbl_buf_ptr = fix_dbl_buf;

    if (xnreqs[1] > 0) {
        /* lat */
        MPI_Offset **fix_starts_D2, **fix_counts_D2;

        /* construct varn API arguments starts[][] and counts[][] */
        int num = xnreqs[1];
        FIX_1D_VAR_STARTS_COUNTS(fix_starts_D2, fix_counts_D2, num,
                                 decom.disps[1], decom.blocklens[1]);

        REC_2D_VAR_STARTS_COUNTS(0, starts_D2, counts_D2, xnreqs[1],
                                 decom.disps[1], decom.blocklens[1]);

        err = IPUT_VARN(ncid, varids[i++], xnreqs[1], fix_starts_D2,
                        fix_counts_D2, fix_dbl_buf_ptr, nelems[1], MPI_DOUBLE,
                        NULL);
        CHECK_ERR
        fix_dbl_buf_ptr += nelems[1] + gap;
        my_nreqs += xnreqs[1];
        nvars_D[1]++;

        /* lon */
        err = IPUT_VARN(ncid, varids[i++], xnreqs[1], fix_starts_D2,
                        fix_counts_D2, fix_dbl_buf_ptr, nelems[1], MPI_DOUBLE,
                        NULL);
        CHECK_ERR
        fix_dbl_buf_ptr += nelems[1] + gap;
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
        FIX_1D_VAR_STARTS_COUNTS(fix_starts_D1, fix_counts_D1, xnreqs[0],
                                 decom.disps[0], decom.blocklens[0]);

        err = IPUT_VARN(ncid, varids[i++], xnreqs[0], fix_starts_D1,
                        fix_counts_D1, fix_dbl_buf_ptr, nelems[0], MPI_DOUBLE,
                        NULL);
        CHECK_ERR
        fix_dbl_buf_ptr += nelems[0] + gap;
        my_nreqs += xnreqs[0];
        nvars_D[0]++;

        free (fix_starts_D1[0]);
        free (fix_starts_D1);
    } else
        i++;

    /* construct varn API arguments starts[][] and counts[][] */
    if (xnreqs[2] > 0)
        REC_3D_VAR_STARTS_COUNTS(0, starts_D3, counts_D3, xnreqs[2],
                                 decom.disps[2], decom.blocklens[2],
                                 decom.dims[2][1]);

    /* post iput for the remaining fixed-size variables */
    fix_int_buf_ptr = fix_int_buf;
    err = write_small_fix_vars_F_case(driver, ncid, 0, gap, decom.dims[2][0],
                                      decom.dims[2][0] + 1,
                                      &fix_int_buf_ptr, &fix_dbl_buf_ptr, NULL);
    CHECK_ERR

    rec_int_buf_ptr = rec_int_buf;
    rec_txt_buf_ptr = rec_txt_buf;
    rec_dbl_buf_ptr = rec_dbl_buf;
    rec_buf_ptr     = rec_buf;

    /* now write record variables */
    for (rec_no = 0; rec_no < cfg.nrecs; rec_no++) {
        if (!one_flush || (cfg.strategy == blob && cfg.api == hdf5)) {
            /* reset the pointers to the beginning of the buffers */
            rec_int_buf_ptr = rec_int_buf;
            rec_txt_buf_ptr = rec_txt_buf;
            rec_dbl_buf_ptr = rec_dbl_buf;
            rec_buf_ptr     = rec_buf;
        }

        if (rank == 0) {
            /* write 15 small record variables by rank 0 only */
            err = write_small_rec_vars_F_case(driver, ncid, 0, rec_no, gap,
                                              2, 8,
                                              &rec_int_buf_ptr,
                                              &rec_txt_buf_ptr,
                                              &rec_dbl_buf_ptr, NULL);
            CHECK_ERR
        }

        for (j = 0; j < xnreqs[1]; j++) starts_D2[j][0] = rec_no;
        for (j = 0; j < xnreqs[2]; j++) starts_D3[j][0] = rec_no;

        if (cfg.nvars == 414) {
            if (cfg.two_buf) {
                /* write 2D variables */
                POST_VARN(2,   1,  30)  /* AEROD_v */
                POST_VARN(2,  19,  33)  /* AODABS ... AODVIS */
                POST_VARN(2,   6,  54)  /* AQ_DMS ... AQ_SOAG */
                POST_VARN(2,   4,  64)  /* BURDEN1 ... BURDEN4 */
                POST_VARN(2,   2,  69)  /* CDNUMC and CLDHGH */
                POST_VARN(2,   3,  73)  /* CLDLOW ... CLDTOT */
                POST_VARN(2,  11,  80)  /* DF_DMS ... DSTSFMBL */
                POST_VARN(2,   2,  92)  /* DTENDTH and DTENDTQ */
                POST_VARN(2,   7,  96)  /* FLDS ... FLUTC */
                POST_VARN(2,  15, 107)  /* FSDS ... ICEFRAC */
                POST_VARN(2,   2, 125)  /* LANDFRAC and LHFLX */
                POST_VARN(2,   1, 131)  /* LINOZ_SFCSINK */
                POST_VARN(2,   3, 133)  /* LINOZ_SZA ... LWCF */
                POST_VARN(2,   2, 148)  /* O3_SRF and OCNFRAC */
                POST_VARN(2,   1, 151)  /* OMEGA500 */
                POST_VARN(2,   8, 153)  /* PBLH ... PSL */
                POST_VARN(2,   2, 162)  /* QFLX and QREFHT */
                POST_VARN(2,   1, 167)  /* RAM1 */
                POST_VARN(2,  37, 169)  /* SFDMS ... SNOWHLND */
                POST_VARN(2,  10, 208)  /* SO2_CLXF ... SWCF */
                POST_VARN(2,  19, 219)  /* TAUGWX ... TVQ */
                POST_VARN(2,   1, 239)  /* U10 */
                POST_VARN(2,   3, 246)  /* WD_H2O2 ... WD_SO2 */
                POST_VARN(2,  32, 252)  /* airFV ... dst_c3SFWET */
                POST_VARN(2, 129, 285)  /* mlip ... soa_c3SFWET */
                /* write 3D variables */
                POST_VARN(3,   2, 31)   /* ANRAIN and ANSNOW */
                POST_VARN(3,   2, 52)   /* AQRAIN and AQSNOW */
                POST_VARN(3,   4, 60)   /* AREI ... AWNI */
                POST_VARN(3,   1, 68)   /* CCN3 */
                POST_VARN(3,   2, 71)   /* CLDICE and CLDLIQ */
                POST_VARN(3,   4, 76)   /* CLOUD ... DCQ */
                POST_VARN(3,   1, 91)   /* DTCOND */
                POST_VARN(3,   2, 94)   /* EXTINCT and FICE */
                POST_VARN(3,   4, 103)  /* FREQI ... FREQS */
                POST_VARN(3,   3, 122)  /* ICIMR ... IWC */
                POST_VARN(3,   4, 127)  /* LINOZ_DO3 ... LINOZ_O3COL */
                POST_VARN(3,   1, 132)  /* LINOZ_SSO3 */
                POST_VARN(3,  12, 136)  /* Mass_bc ... O3 */
                POST_VARN(3,   1, 150)  /* OMEGA */
                POST_VARN(3,   1, 152)  /* OMEGAT */
                POST_VARN(3,   1, 161)  /* Q */
                POST_VARN(3,   3, 164)  /* QRL ... RAINQM */
                POST_VARN(3,   1, 168)  /* RELHUM */
                POST_VARN(3,   2, 206)  /* SNOWQM and SO2 */
                POST_VARN(3,   1, 218)  /* T */
                POST_VARN(3,   1, 238)  /* U */
                POST_VARN(3,   6, 240)  /* UU ... VV */
                POST_VARN(3,   3, 249)  /* WSUB ... aero_water */
                POST_VARN(3,   1, 284)  /* hstobie_linoz */
            } else {
                /* write variables in the same order as they defined */
                POST_VARN(2,   1,  30)  /* AEROD_v */
                POST_VARN(3,   2,  31)  /* ANRAIN and ANSNOW */
                POST_VARN(2,  19,  33)  /* AODABS ... AODVIS */
                POST_VARN(3,   2,  52)  /* AQRAIN and AQSNOW */
                POST_VARN(2,   6,  54)  /* AQ_DMS ... AQ_SOAG */
                POST_VARN(3,   4,  60)  /* AREI ... AWNI */
                POST_VARN(2,   4,  64)  /* BURDEN1 ... BURDEN4 */
                POST_VARN(3,   1,  68)  /* CCN3 */
                POST_VARN(2,   2,  69)  /* CDNUMC and CLDHGH */
                POST_VARN(3,   2,  71)  /* CLDICE and CLDLIQ */
                POST_VARN(2,   3,  73)  /* CLDLOW ... CLDTOT */
                POST_VARN(3,   4,  76)  /* CLOUD ... DCQ */
                POST_VARN(2,  11,  80)  /* DF_DMS ... DSTSFMBL */
                POST_VARN(3,   1,  91)  /* DTCOND */
                POST_VARN(2,   2,  92)  /* DTENDTH and DTENDTQ */
                POST_VARN(3,   2,  94)  /* EXTINCT and FICE */
                POST_VARN(2,   7,  96)  /* FLDS ... FLUTC */
                POST_VARN(3,   4, 103)  /* FREQI ... FREQS */
                POST_VARN(2,  15, 107)  /* FSDS ... ICEFRAC */
                POST_VARN(3,   3, 122)  /* ICIMR ... IWC */
                POST_VARN(2,   2, 125)  /* LANDFRAC and LHFLX */
                POST_VARN(3,   4, 127)  /* LINOZ_DO3 ... LINOZ_O3COL */
                POST_VARN(2,   1, 131)  /* LINOZ_SFCSINK */
                POST_VARN(3,   1, 132)  /* LINOZ_SSO3 */
                POST_VARN(2,   3, 133)  /* LINOZ_SZA ... LWCF */
                POST_VARN(3,  12, 136)  /* Mass_bc ... O3 */
                POST_VARN(2,   2, 148)  /* O3_SRF and OCNFRAC */
                POST_VARN(3,   1, 150)  /* OMEGA */
                POST_VARN(2,   1, 151)  /* OMEGA500 */
                POST_VARN(3,   1, 152)  /* OMEGAT */
                POST_VARN(2,   8, 153)  /* PBLH ... PSL */
                POST_VARN(3,   1, 161)  /* Q */
                POST_VARN(2,   2, 162)  /* QFLX and QREFHT */
                POST_VARN(3,   3, 164)  /* QRL ... RAINQM */
                POST_VARN(2,   1, 167)  /* RAM1 */
                POST_VARN(3,   1, 168)  /* RELHUM */
                POST_VARN(2,  37, 169)  /* SFDMS ... SNOWHLND */
                POST_VARN(3,   2, 206)  /* SNOWQM and SO2 */
                POST_VARN(2,  10, 208)  /* SO2_CLXF ... SWCF */
                POST_VARN(3,   1, 218)  /* T */
                POST_VARN(2,  19, 219)  /* TAUGWX ... TVQ */
                POST_VARN(3,   1, 238)  /* U */
                POST_VARN(2,   1, 239)  /* U10 */
                POST_VARN(3,   6, 240)  /* UU ... VV */
                POST_VARN(2,   3, 246)  /* WD_H2O2 ... WD_SO2 */
                POST_VARN(3,   3, 249)  /* WSUB ... aero_water */
                POST_VARN(2,  32, 252)  /* airFV ... dst_c3SFWET */
                POST_VARN(3,   1, 284)  /* hstobie_linoz */
                POST_VARN(2, 129, 285)  /* mlip ... soa_c3SFWET */
            }
        } else { /* h1 file */
            if (cfg.two_buf) {
                /* write 2D variables followed by 3D variables */
                POST_VARN(2, 13, 30)  /* CLDHGH ... T5 */
                POST_VARN(2,  7, 44)  /* U250 ... Z500 */
                POST_VARN(3,  1, 43)  /* U */
            } else {
                /* write variables in the same order as they defined */
                POST_VARN(2, 13, 30)  /* CLDHGH ... T5 */
                POST_VARN(3,  1, 43)  /* U */
                POST_VARN(2,  7, 44)  /* U250 ... Z500 */
            }
        }
        if (!one_flush) {
            cfg.post_time += MPI_Wtime() - timing;

            MPI_Barrier(cfg.io_comm); /*-----------------------------------*/
            timing = MPI_Wtime();

            /* flush once per time record */
            WAIT_ALL_REQS

            cfg.flush_time += MPI_Wtime() - timing;

            timing = MPI_Wtime();
        }
    }

    if (one_flush) {
        cfg.post_time += MPI_Wtime() - timing;

        MPI_Barrier(cfg.io_comm); /*---------------------------------------*/
        timing = MPI_Wtime();

        /* flush once for all time records */
        WAIT_ALL_REQS

        cfg.flush_time += MPI_Wtime() - timing;
    }

    MPI_Barrier (cfg.io_comm); /*-----------------------------------------*/
    timing = MPI_Wtime ();

    if (starts_D3 != NULL) {
        free (starts_D3[0]);
        free (starts_D3);
    }
    if (starts_D2 != NULL) {
        free (starts_D2[0]);
        free (starts_D2);
    }
    if (fix_int_buf != NULL) free(fix_int_buf);
    if (fix_dbl_buf != NULL) free(fix_dbl_buf);
    if (rec_int_buf != NULL) free(rec_int_buf);
    if (rec_txt_buf != NULL) free(rec_txt_buf);
    if (rec_dbl_buf != NULL) free(rec_dbl_buf);
    if (rec_buf     != NULL) free(rec_buf);
    if (varids      != NULL) free(varids);

    err = driver.inq_put_size (ncid, &total_size);
    CHECK_ERR

    err = driver.close (ncid);
    CHECK_ERR

    cfg.close_time = MPI_Wtime () - timing;

    cfg.num_flushes     = nflushes;
    cfg.num_decomp      = decom.num_decomp;
    cfg.num_decomp_vars = 0;
    cfg.my_nreqs        = my_nreqs;
    cfg.metadata_WR     = metadata_size;
    cfg.amount_WR       = total_size;
    cfg.end2end_time    = MPI_Wtime() - cfg.end2end_time;

    /* check if there is any PnetCDF internal malloc residue */
    check_malloc(&cfg, &driver);

    /* report timing breakdowns */
    report_timing_WR(&cfg, &driver, outfile);

    /* print MPI-IO hints actually used */
    if (cfg.verbose && rank == 0) print_info(&info_used);

err_out:
    if (info_used != MPI_INFO_NULL) MPI_Info_free (&info_used);
    if (!cfg.keep_outfile) unlink (outfile);
    return err;
}

#define POST_VARN_RD(k, num, vid)                                                            \
    for (j = 0; j < num; j++) {                                                              \
        err = IGET_VARN (ncid, vid + j, decom.contig_nreqs[k - 1], starts_D##k, counts_D##k, \
                         rec_buf_ptr, -1, REC_ITYPE, NULL);                                  \
        CHECK_ERR                                                                            \
        rec_buf_ptr += nelems[k - 1] + gap;                                                  \
        my_nreqs += decom.contig_nreqs[k - 1];                                               \
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
    int i, err;
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

/*----< run_varn_F_case_rd() >-----------------------------------------------*/
int run_varn_F_case_rd (e3sm_io_config &cfg,
                        e3sm_io_decom &decom,
                        e3sm_io_driver &driver,
                        double **dbl_bufp,   /* buffer for fixed size double var */
                        itype **rec_bufp,    /* buffer for rec floating point var */
                        char *txt_buf,       /* buffer for char var */
                        int *int_buf)        /* buffer for int var */
{
    char *txt_buf_ptr;
    int i, j, k, err, rank, ncid, *varids, nflushes=0;
    int rec_no, gap = 0, my_nreqs, *int_buf_ptr;
    size_t dbl_buflen, rec_buflen, nelems[3];
    itype *rec_buf, *rec_buf_ptr;
    double *dbl_buf, *dbl_buf_ptr;
    double pre_timing, open_timing, post_timing, wait_timing, close_timing;
    double timing, total_timing, max_timing;
    MPI_Offset tmp, metadata_size, put_size, total_size, fsize, max_nreqs, total_nreqs;
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
    dbl_buflen = nelems[1] * 2
               + nelems[0]
               + 3 * decom.dims[2][0]
               + 3 * (decom.dims[2][0] + 1)
               + 8 + 2 + 20 * gap;

    dbl_buf = (double *)malloc (dbl_buflen * sizeof (double));
    if (dbl_bufp != NULL) { *dbl_bufp = dbl_buf; }

    /* allocate read buffer for large variables */
    if (cfg.nvars == 414)
        rec_buflen = nelems[1] * 321 + nelems[2] * 63 + (321 + 63) * gap;
    else
        rec_buflen = nelems[1] * 20 + nelems[2] + (20 + 1) * gap;

    rec_buf = (itype *)malloc (rec_buflen * sizeof (itype));
    if (rec_bufp != NULL) { *rec_bufp = rec_buf; }

    varids = (int *)malloc (cfg.nvars * sizeof (int));

    pre_timing = MPI_Wtime () - pre_timing;

    MPI_Barrier (comm); /*-----------------------------------------*/
    timing = MPI_Wtime ();

    /* open input file for reading */
    err = driver.open (cfg.in_path, comm, cfg.info, &ncid);
    CHECK_ERR

    /* inquire dimensions, variables, and attributes */
    if (cfg.nvars == 414) {
        /* for h0 file */
        err = inq_F_case_h0 (driver, ncid, decom.dims[2], cfg.nvars, varids);
        CHECK_ERR
    } else {
        /* for h1 file */
        err = inq_F_case_h1 (driver, ncid, decom.dims[2], cfg.nvars, varids);
        CHECK_ERR
    }

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

    for (rec_no = 0; rec_no < cfg.nrecs; rec_no++) {
        MPI_Barrier (comm); /*-----------------------------------------*/
        timing = MPI_Wtime ();

        i           = 3;
        dbl_buf_ptr = dbl_buf + nelems[1] * 2 + nelems[0] + gap * 3;
        int_buf_ptr = int_buf;
        txt_buf_ptr = txt_buf;

        /* next 27 small variables are read by rank 0 only */
        if (rank == 0) {
            my_nreqs += 27;
            /* post nonblocking requests using IGET_VARN() */
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
        WAIT_ALL_REQS

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

        WAIT_ALL_REQS

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

    if (cfg.rank == 0){
        err = driver.inq_file_size(cfg.in_path, &fsize);
        CHECK_ERR
    }

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
        printf ("History input file                 = %s\n", cfg.in_path);
        printf ("Input file size                    = %.2f MiB = %.2f GiB\n",
                (double)fsize / 1048576, (double)fsize / 1073741824);
        if (cfg.api == pnetcdf)
            printf ("MAX heap memory allocated by PnetCDF internally is %.2f MiB\n",
                (float)max_alloc / 1048576);
        printf ("Total number of variables          = %d\n", cfg.nvars);
        printf ("Total read amount                  = %.2f MiB = %.2f GiB\n",
                (double)total_size / 1048576, (double)total_size / 1073741824);
        printf ("Total number of requests           = %lld\n", total_nreqs);
        printf ("Max number of requests             = %lld\n", max_nreqs);
        printf ("Max Time of open + metadata inquery = %.4f sec\n", open_timing);
        printf ("Max Time of I/O preparing          = %.4f sec\n", pre_timing);
        printf ("Max Time of posting IGET_VARN      = %.4f sec\n", post_timing);
        if (cfg.api == pnetcdf)
            printf ("Max Time of read flushing          = %.4f sec\n", wait_timing);
        printf ("Max Time of close                  = %.4f sec\n", close_timing);
        printf ("Max Time of TOTAL                  = %.4f sec\n", total_timing);
        printf ("I/O bandwidth (open-to-close)      = %.4f MiB/sec\n",
                (double)total_size / 1048576.0 / total_timing);
        printf ("I/O bandwidth (read-only)          = %.4f MiB/sec\n",
                (double)put_size / 1048576.0 / wait_timing);
        if (cfg.verbose) print_info (&info_used);
        printf ("-----------------------------------------------------------\n");
    }
    fflush (stdout);

err_out:
    if (info_used != MPI_INFO_NULL) MPI_Info_free (&info_used);
    return err;
}
