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

#define VAR_ITYPE REC_ITYPE

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
    err = driver.inq_put_size(&size);        \
    CHECK_ERR                                \
}
#define INQ_FILE_INFO(info) {                \
    err = driver.inq_file_info(ncid, &info); \
    CHECK_ERR                                \
}
#define IPUT_VARA_NOADV(itype, buf) { \
    err = driver.put_vara(ncid, varid, itype, start, count, buf, nb); \
    CHECK_VAR_ERR(varid) \
    my_nreqs++; \
    varid++; \
}
#define IPUT_VARA(varid, itype, adv, buf) { \
    err = driver.put_vara(ncid, varid, itype, start, count, buf, nb); \
    CHECK_VAR_ERR(varid) \
    buf += (adv); \
    my_nreqs++; \
}
#define IPUT_VAR(varid, itype, adv, buf) { \
    err = driver.put_vara(ncid, varid, itype, NULL, NULL, buf, nb); \
    CHECK_VAR_ERR(varid) \
    buf += (adv); \
    my_nreqs++; \
}
#define WAIT_ALL_REQS { \
    err = driver.wait(ncid); \
    CHECK_ERR \
    nflushes++; \
}

#define FIX_VAR_IPUT(varid, dp, itype, buf) {                             \
    if (cfg.strategy == canonical) {                                      \
        err = driver.put_varn(ncid, varid, itype, decom.contig_nreqs[dp], \
                              decom.w_startx[dp], decom.w_countx[dp],     \
                              buf, nb);                                   \
        my_nreqs += decom.contig_nreqs[dp];                               \
    }                                                                     \
    else {                                                                \
        err = driver.put_vara(ncid, varid, itype, &decom.start[dp],       \
                              &decom.count[dp], buf, nb);                 \
        my_nreqs++;                                                       \
    }                                                                     \
    CHECK_VAR_ERR(varid)                                                  \
    buf += decom.count[dp] + gap;                                         \
    nvars_D[dp]++;                                                        \
}

#define REC_VAR_IPUT(varid, dp, itype, buf) {                             \
    if (cfg.strategy == canonical) {                                      \
        err = driver.put_varn(ncid, varid, itype, decom.contig_nreqs[dp], \
                              decom.w_starts[dp], decom.w_counts[dp],     \
                              buf, nb);                                   \
        my_nreqs += decom.contig_nreqs[dp];                               \
    }                                                                     \
    else {                                                                \
        err = driver.put_vara(ncid, varid, itype, starts[dp],             \
                              counts[dp], buf, nb);                       \
        my_nreqs++;                                                       \
    }                                                                     \
    CHECK_VAR_ERR(varid)                                                  \
    buf += decom.count[dp] + gap;                                         \
    if (rec_no == 0) nvars_D[dp]++;                                       \
}

/*----< var_wr_I_case() >----------------------------------------------------*/
int var_wr_I_case(e3sm_io_config &cfg,
                  e3sm_io_decom  &decom,
                  e3sm_io_driver &driver,
                  char           *outfile)
{
    char *fix_txt_buf_ptr, *rec_txt_buf_ptr;
    int i, j, err, sub_rank, global_rank, ncid=-1, nflushes=0, one_flush;
    int nrecs, rec_no, gap=0, my_nreqs, num_decomp_vars;
    int contig_nreqs[MAX_NUM_DECOMP], *nvars_D;
    int *fix_int_buf_ptr, *rec_int_buf_ptr;
    double *fix_dbl_buf_ptr, *rec_dbl_buf_ptr, timing;
    MPI_Offset previous_size, metadata_size, total_size;
    MPI_Offset blob_start[MAX_NUM_DECOMP], blob_count[MAX_NUM_DECOMP];
    MPI_Offset starts[MAX_NUM_DECOMP][2], counts[MAX_NUM_DECOMP][2];
    MPI_Offset start[2], count[2];
    MPI_Info info_used = MPI_INFO_NULL;
    MPI_Comm comm;
    vtype  *fix_buf_ptr, *rec_buf_ptr;
    var_meta *vars;
    io_buffers wr_buf;
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

    MPI_Barrier(cfg.io_comm); /*---------------------------------------------*/
    pr->end2end_time = timing = MPI_Wtime();

    pr->post_time  = 0.0;
    pr->flush_time = 0.0;

    comm = (cfg.strategy == blob) ? cfg.sub_comm : cfg.io_comm;

    MPI_Comm_rank(cfg.io_comm,  &global_rank);
    MPI_Comm_rank(comm,         &sub_rank);

    nvars_D = pr->nvars_D;

    /* I/O amount from previous I/O */
    INQ_PUT_SIZE(previous_size)

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

    if (cfg.verbose) {
        if (cfg.strategy == blob && sub_rank == 0)
            printf("global_rank=%d sub_rank=%d outfile=%s\n",
                   global_rank,sub_rank,outfile);
        else if (global_rank == 0)
            printf("global_rank=%d outfile=%s\n",global_rank,outfile);
    }

    /* create the output file */
    FILE_CREATE(outfile)

    pr->open_time = MPI_Wtime() - timing;

    MPI_Barrier(comm); /*----------------------------------------------------*/
    timing = MPI_Wtime();

    /* there are num_decomp_vars number of decomposition variables */
    if (cfg.strategy == blob)
        num_decomp_vars = decom.num_decomp * NVARS_DECOMP;
    else
        num_decomp_vars = 0;

    /* allocate space to store variable metadata */
    vars = (var_meta*) malloc((cfg.nvars + num_decomp_vars) * sizeof(var_meta));

    if (cfg.non_contig_buf) gap = 10;

    wr_buf_init(wr_buf, gap);

    /* define dimensions, variables, and attributes.
     * 560 climate variables
     *     18 are fixed-size,      542 are record variables
     *     14 are not partitioned, 546 are partitioned.
     * Thus, 14 not partitioned variable (also small) are written by root only:
     *     5 fixed-size, 9 record variables
     */

    /* h0 file contains 560 variables, while h1 has 8 less */
    err = def_I_case(cfg, decom, driver, ncid, vars, &wr_buf);
    CHECK_ERR

    /* exit define mode and enter data mode */
    ENDDEF

    /* I/O amount so far */
    INQ_PUT_SIZE(metadata_size)
    metadata_size -= previous_size;

    INQ_FILE_INFO(info_used)

    pr->def_time = MPI_Wtime() - timing;

    MPI_Barrier(comm); /*----------------------------------------------------*/
    timing = MPI_Wtime();

    /* allocate write buffers */
    wr_buf_malloc(cfg, one_flush, wr_buf);

    for (j=0; j<decom.num_decomp; j++)
        nvars_D[j] = 0; /* number of variables using decomposition j */

    pr->pre_time = MPI_Wtime() - timing;

    MPI_Barrier(comm); /*----------------------------------------------------*/
    timing = MPI_Wtime();

    my_nreqs = 0; /* number of noncontiguous requests written by this rank */

    if (cfg.strategy == blob) {
        int varid=0;

        /* write decomposition variables, they are defined first in file */
        for (j=0; j<decom.num_decomp; j++) {
            start[0] = sub_rank;
            count[0] = 1;
            start[1] = 0;
            count[1] = decom.contig_nreqs[j];

            /* write to D*.nreqs, 1D array */
            contig_nreqs[j] = decom.contig_nreqs[j];
            IPUT_VARA_NOADV(MPI_INT, contig_nreqs+j)

            /* write to D*.blob_start, 1D array */
            blob_start[j] = decom.start[j];
            IPUT_VARA_NOADV(MPI_LONG_LONG, blob_start+j)

            /* write to D*.blob_count, 1D array */
            blob_count[j] = decom.count[j];
            IPUT_VARA_NOADV(MPI_LONG_LONG, blob_count+j);

            /* write to D*.offsets, 2D array */
            IPUT_VARA_NOADV(MPI_INT, decom.disps[j])

            /* write to D*.lengths, 2D array */
            IPUT_VARA_NOADV(MPI_INT, decom.blocklens[j])

            /* these 4 are used for record variables in blob I/O */
            starts[j][0] = 0;
            starts[j][1] = decom.start[j];
            counts[j][0] = 1;
            counts[j][1] = decom.count[j];
        }
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

    fix_txt_buf_ptr = wr_buf.fix_txt_buf;
    fix_int_buf_ptr = wr_buf.fix_int_buf;
    fix_dbl_buf_ptr = wr_buf.fix_dbl_buf;
    fix_buf_ptr     = wr_buf.fix_buf;

    for (rec_no=0; rec_no<nrecs; rec_no++) {

        if (rec_no == 0 || !one_flush ||
            (cfg.strategy == blob && cfg.api == hdf5)) {
            /* reset the pointers to the beginning of the buffers */
            rec_txt_buf_ptr = wr_buf.rec_txt_buf;
            rec_int_buf_ptr = wr_buf.rec_int_buf;
            rec_dbl_buf_ptr = wr_buf.rec_dbl_buf;
            rec_buf_ptr     = wr_buf.rec_buf;
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

        /* write all climate variables */
        for (j=num_decomp_vars; j<cfg.nvars+num_decomp_vars; j++) {
            int          varid = vars[j].vid;
            int          dp    = vars[j].decomp_id;
            MPI_Datatype itype = vars[j].itype;
            size_t       adv   = vars[j].vlen + gap;

            if (vars[j].decomp_id >= 0) { /* this variable is partitioned */
                if (vars[j].isRecVar) { /* this is a record variable */
                    REC_VAR_IPUT(varid, dp, VAR_ITYPE, rec_buf_ptr)
                }
                else if (rec_no == 0) { /* this is a fixed-size variable */
                    if (itype == VAR_ITYPE)
                        FIX_VAR_IPUT(varid, dp, VAR_ITYPE, fix_buf_ptr)
                    else if (itype == MPI_INT)
                        FIX_VAR_IPUT(varid, dp, MPI_INT, fix_int_buf_ptr)
                }
            }
            else if (sub_rank == 0) {
                /* not-partitioned variables are written by root only */

                if (vars[j].isRecVar) { /* this is a record variable */
                    start[0] = rec_no;
                    start[1] = 0;
                    count[0] = 1;
                    count[1] = vars[j].vlen;
                    if (itype == VAR_ITYPE)
                        IPUT_VARA(varid, itype, adv,     rec_buf_ptr)
                    else if (itype == MPI_INT)
                        IPUT_VARA(varid, itype, adv, rec_int_buf_ptr)
                    else if (itype == MPI_CHAR)
                        IPUT_VARA(varid, itype, adv, rec_txt_buf_ptr)
                    else if (itype == MPI_DOUBLE)
                        IPUT_VARA(varid, itype, adv, rec_dbl_buf_ptr)
                }
                else if (rec_no == 0) { /* this is a fixed-size variable */
                    if (itype == VAR_ITYPE)
                        IPUT_VAR(varid, itype, adv,     fix_buf_ptr)
                    else if (itype == MPI_INT)
                        IPUT_VAR(varid, itype, adv, fix_int_buf_ptr)
                    else if (itype == MPI_CHAR)
                        IPUT_VAR(varid, itype, adv, fix_txt_buf_ptr)
                    else if (itype == MPI_DOUBLE)
                        IPUT_VAR(varid, itype, adv, fix_dbl_buf_ptr)
                }
            }
        }

        if (!one_flush) { /* flush out for each record */
            pr->post_time += MPI_Wtime() - timing;

            MPI_Barrier(comm); /*--------------------------------------------*/
            timing = MPI_Wtime();

            /* flush once per time record */
            WAIT_ALL_REQS

            pr->flush_time += MPI_Wtime() - timing;

            timing = MPI_Wtime();
        }
    }

    if (one_flush) { /* flush out for all records */
        pr->post_time += MPI_Wtime() - timing;

        MPI_Barrier(comm); /*------------------------------------------------*/
        timing = MPI_Wtime();

        /* flush once for all time records */
        WAIT_ALL_REQS

        pr->flush_time += MPI_Wtime() - timing;
    }
    MPI_Barrier(comm); /*----------------------------------------------------*/
    timing = MPI_Wtime();

    /* free up previously allocated heap memory space */
    wr_buf_free(wr_buf);
    if (vars != NULL) free(vars);

    FILE_CLOSE

    pr->close_time = MPI_Wtime() - timing;

    /* obtain the write amount tracked by the driver */
    INQ_PUT_SIZE(total_size)
    total_size -= previous_size;

    pr->nvars           = cfg.nvars;
    pr->num_flushes     = nflushes;
    pr->num_decomp_vars = num_decomp_vars;
    pr->my_nreqs        = my_nreqs;
    pr->metadata_WR     = metadata_size;
    pr->amount_WR       = total_size;
    pr->end2end_time    = MPI_Wtime() - pr->end2end_time;

    /* check if there is any PnetCDF internal malloc residue */
    check_malloc(&cfg, &driver);

    /* print MPI-IO hints actually used */
    if (cfg.verbose && global_rank == 0) print_info(&info_used);

err_out:
    if (err < 0 && ncid >= 0) driver.close(ncid);
    if (info_used != MPI_INFO_NULL) MPI_Info_free(&info_used);
    if (!cfg.keep_outfile && sub_rank == 0) unlink(outfile);

    return err;
}

