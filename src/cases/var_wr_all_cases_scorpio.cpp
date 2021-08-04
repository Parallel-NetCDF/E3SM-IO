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
#include <e3sm_io_case_scorpio.hpp>
#include <e3sm_io_driver.hpp>

#define VAR_ITYPE REC_ITYPE

#define CHECK_VAR_ERR(varid)                                                                  \
    {                                                                                         \
        if (err != 0) {                                                                       \
            char var_name[64];                                                                \
            driver.inq_var_name (ncid, (varid.data), var_name);                               \
            printf ("Error in %s:%d: %s() var %s\n", __FILE__, __LINE__, __func__, var_name); \
            goto err_out;                                                                     \
        }                                                                                     \
    }
#define FILE_CREATE(filename)                                  \
    {                                                          \
        err = driver.create (filename, comm, cfg.info, &ncid); \
        CHECK_ERR                                              \
    }
#define FILE_CLOSE                 \
    {                              \
        err = driver.close (ncid); \
        CHECK_ERR                  \
    }
#define ENDDEF                      \
    {                               \
        err = driver.enddef (ncid); \
        CHECK_ERR                   \
    }
#define INQ_DIM_LEN(name, size)                       \
    {                                                 \
        int dimid;                                    \
        err = driver.inq_dim (ncid, name, &dimid);    \
        CHECK_ERR                                     \
        err = driver.inq_dimlen (ncid, dimid, &size); \
        CHECK_ERR                                     \
    }
#define INQ_PUT_SIZE(size)                 \
    {                                      \
        err = driver.inq_put_size (&size); \
        CHECK_ERR                          \
    }
#define INQ_FILE_INFO(info)                       \
    {                                             \
        err = driver.inq_file_info (ncid, &info); \
        CHECK_ERR                                 \
    }
#define IPUT_VARA(varid, itype, adv, buf)                                              \
    {                                                                                  \
        err = e3sm_io_scorpio_write_var (driver, rec_no, ncid, varid, itype, buf, nb); \
        CHECK_VAR_ERR (varid)                                                          \
        buf += (adv);                                                                  \
        my_nreqs++;                                                                    \
    }
#define IPUT_VAR(varid, itype, adv, buf)                                               \
    {                                                                                  \
        err = e3sm_io_scorpio_write_var (driver, rec_no, ncid, varid, itype, buf, nb); \
        CHECK_VAR_ERR (varid)                                                          \
        buf += (adv);                                                                  \
        my_nreqs++;                                                                    \
    }
#define WAIT_ALL_REQS             \
    {                             \
        err = driver.wait (ncid); \
        CHECK_ERR                 \
        nflushes++;               \
    }

#define FIX_VAR_IPUT(varid, dp, itype, buf)                                            \
    {                                                                                  \
        err = e3sm_io_scorpio_write_var (driver, -1, ncid, varid, itype, buf, nb); \
        my_nreqs += decom.contig_nreqs[dp];                                            \
        CHECK_VAR_ERR (varid)                                                          \
        buf += decom.raw_nreqs[dp] + gap;                                                  \
        nvars_D[dp]++;                                                                 \
    }

#define REC_VAR_IPUT(varid, dp, itype, buf)                                            \
    {                                                                                  \
        err = e3sm_io_scorpio_write_var (driver, rec_no, ncid, varid, itype, buf, nb); \
        my_nreqs += decom.contig_nreqs[dp];                                            \
        CHECK_VAR_ERR (varid)                                                          \
        buf += decom.raw_nreqs[dp] + gap;                                                  \
        if (rec_no == 0) nvars_D[dp]++;                                                \
    }

/*----< var_wr_all_cases_scorpio() >-------------------------------------------------*/
int var_wr_all_cases_scorpio(e3sm_io_config &cfg,
                     e3sm_io_decom  &decom,
                     e3sm_io_driver &driver,
                     case_meta      *cmeta)
{
    char *fix_txt_buf_ptr, *rec_txt_buf_ptr;
    int i, j, err = 0, sub_rank, global_rank, ncid=-1, nflushes=0, one_flush;
    int rec_no, gap=0, my_nreqs, nvars, num_decomp_vars = 0;
    int *nvars_D;
    int *fix_int_buf_ptr, *rec_int_buf_ptr;
    double *fix_dbl_buf_ptr, *rec_dbl_buf_ptr, timing;
    MPI_Offset previous_size, metadata_size, total_size;
    MPI_Info info_used = MPI_INFO_NULL;
    MPI_Comm comm;
    vtype *fix_buf_ptr, *rec_buf_ptr;
    var_meta_scorpio *vars;
    io_buffers wr_buf;
    int scorpiovars[7];

    MPI_Barrier(cfg.io_comm); /*---------------------------------------------*/
    cmeta->end2end_time = timing = MPI_Wtime();

    cmeta->post_time  = 0.0;
    cmeta->flush_time = 0.0;

    comm = cfg.io_comm;

    MPI_Comm_rank(cfg.io_comm,  &global_rank);
    MPI_Comm_rank(comm,         &sub_rank);

    nvars_D = cmeta->nvars_D;

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
                   global_rank,sub_rank,cmeta->outfile);
        else if (global_rank == 0)
            printf("global_rank=%d outfile=%s\n",global_rank,cmeta->outfile);
    }

    /* create the output file */
    FILE_CREATE(cmeta->outfile)

    cmeta->open_time = MPI_Wtime() - timing;

    MPI_Barrier(comm); /*----------------------------------------------------*/
    timing = MPI_Wtime();

    /* there are num_decomp_vars number of decomposition variables */
    nvars = cmeta->nvars + num_decomp_vars;

    /* allocate space to store variable metadata */
    vars = (var_meta_scorpio*) malloc(nvars * sizeof(var_meta_scorpio));

    if (cfg.non_contig_buf) gap = 10;

    wr_buf_init(wr_buf, gap);

    if (cfg.run_case == F)
        ERR_OUT("Wrong function for F case scorpio")
    else if (cfg.run_case == G)
        ERR_OUT("Wrong function for G case scorpio")
    else if (cfg.run_case == I) {
        // First time only set decomp_id requried to define var in scorpio wrapper
        err = def_I_case_scorpio(cfg, decom, driver, ncid, vars, &wr_buf, scorpiovars, false);
        CHECK_ERR
        // Actual define
        err = def_I_case_scorpio(cfg, decom, driver, ncid, vars, &wr_buf, scorpiovars, true);
        CHECK_ERR
    }

    /* exit define mode and enter data mode */
    ENDDEF

    /* I/O amount so far */
    INQ_PUT_SIZE(metadata_size)
    metadata_size -= previous_size;

    INQ_FILE_INFO(info_used)

    cmeta->def_time = MPI_Wtime() - timing;

    MPI_Barrier(comm); /*----------------------------------------------------*/
    timing = MPI_Wtime();

    /* allocate write buffers */
    wr_buf.fix_int_buflen += 64;
    wr_buf.fix_dbl_buflen += 64;
    wr_buf.fix_txt_buflen += 64;
    wr_buf.fix_buflen += 64;
    wr_buf.rec_dbl_buflen += 64;
    wr_buf.rec_int_buflen += 64;
    wr_buf.rec_txt_buflen += 64;
    wr_buf.rec_buflen += 64;
    wr_buf_malloc(cfg, one_flush, wr_buf);

    for (j=0; j<decom.num_decomp; j++)
        nvars_D[j] = 0; /* number of variables using decomposition j */

    cmeta->pre_time = MPI_Wtime() - timing;

    MPI_Barrier(comm); /*----------------------------------------------------*/
    timing = MPI_Wtime();

    // Write PIO decom vars
    for (j = 0; j < decom.num_decomp; j++) {
        err = driver.put_varl (ncid, scorpiovars[j], MPI_LONG_LONG,
                                decom.raw_offsets[j], nb);
        CHECK_ERR
    }

    // Nproc only written by rank 0
    if (cfg.rank == 0) {
        err = driver.put_varl (ncid, scorpiovars[j], MPI_INT, &(cfg.np), nb);
        CHECK_ERR
    }

    my_nreqs = 0; /* number of noncontiguous requests written by this rank */

    fix_txt_buf_ptr = wr_buf.fix_txt_buf;
    fix_int_buf_ptr = wr_buf.fix_int_buf;
    fix_dbl_buf_ptr = wr_buf.fix_dbl_buf;
    fix_buf_ptr     = wr_buf.fix_buf;

    for (rec_no=0; rec_no<cmeta->nrecs; rec_no++) {
        /* reset the pointers to the beginning of the buffers */
        rec_txt_buf_ptr = wr_buf.rec_txt_buf;
        rec_int_buf_ptr = wr_buf.rec_int_buf;
        rec_dbl_buf_ptr = wr_buf.rec_dbl_buf;
        rec_buf_ptr     = wr_buf.rec_buf;

        /* set the start index for the next record */
        for (i=0; i<decom.num_decomp; i++) {
            for (j=0; j<decom.contig_nreqs[i]; j++)
                decom.w_starts[i][j][0] = rec_no;
        }

        /* write all climate variables */
        for (j=num_decomp_vars; j<nvars; j++) {
            e3sm_io_scorpio_var          varid = vars[j].vid;
            int          dp    = vars[j].decomp_id;
            MPI_Datatype itype = vars[j].itype;
            size_t       adv   = vars[j].vlen + gap;

            if (vars[j].decomp_id >= 0) { /* this variable is partitioned */
                if (vars[j].isRecVar) { /* this is a record variable */
                    REC_VAR_IPUT(varid, dp, VAR_ITYPE, rec_buf_ptr)
                    assert(rec_buf_ptr - wr_buf.rec_buf <= wr_buf.rec_buflen);
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
            cmeta->post_time += MPI_Wtime() - timing;

            MPI_Barrier(comm); /*--------------------------------------------*/
            timing = MPI_Wtime();

            /* flush once per time record */
            WAIT_ALL_REQS

            cmeta->flush_time += MPI_Wtime() - timing;

            timing = MPI_Wtime();
        }
    }

    if (one_flush) { /* flush out for all records */
        cmeta->post_time += MPI_Wtime() - timing;

        MPI_Barrier(comm); /*------------------------------------------------*/
        timing = MPI_Wtime();

        /* flush once for all time records */
        WAIT_ALL_REQS

        cmeta->flush_time += MPI_Wtime() - timing;
    }
    MPI_Barrier(comm); /*----------------------------------------------------*/
    timing = MPI_Wtime();

    /* free up previously allocated heap memory space */
    wr_buf_free(wr_buf);
    if (vars != NULL) free(vars);

    FILE_CLOSE

    cmeta->close_time = MPI_Wtime() - timing;

    /* obtain the write amount tracked by the driver */
    INQ_PUT_SIZE(total_size)
    total_size -= previous_size;

    cmeta->num_flushes     = nflushes;
    cmeta->num_decomp_vars = num_decomp_vars;
    cmeta->my_nreqs        = my_nreqs;
    cmeta->metadata_WR     = metadata_size;
    cmeta->amount_WR       = total_size;
    cmeta->end2end_time    = MPI_Wtime() - cmeta->end2end_time;

    /* check if there is any PnetCDF internal malloc residue */
    check_malloc(&cfg, &driver);

    /* print MPI-IO hints actually used */
    if (cfg.verbose && global_rank == 0) print_info(&info_used);

err_out:
    if (err < 0 && ncid >= 0) driver.close(ncid);
    if (info_used != MPI_INFO_NULL) MPI_Info_free(&info_used);
    if (!cfg.keep_outfile && sub_rank == 0) unlink(cmeta->outfile);

    return err;
}

