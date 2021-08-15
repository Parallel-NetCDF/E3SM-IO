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

#ifdef ENABLE_ADIOS2
#include <e3sm_io_driver_adios2.hpp>
#endif

#define VAR_ITYPE REC_ITYPE

#define FILE_CREATE(filename) {                                           \
    err = driver.create(filename, comm, cfg.info, &ncid);                 \
    CHECK_ERR                                                             \
}
#define FILE_CLOSE {                                                      \
    err = driver.close(ncid);                                             \
    CHECK_ERR                                                             \
}
#define ENDDEF {                                                          \
    err = driver.enddef(ncid);                                            \
    CHECK_ERR                                                             \
}
#define INQ_DIM_LEN(name, size) {                                         \
    int dimid;                                                            \
    err = driver.inq_dim(ncid, name, &dimid);                             \
    CHECK_ERR                                                             \
    err = driver.inq_dimlen(ncid, dimid, &size);                          \
    CHECK_ERR                                                             \
}
#define INQ_PUT_SIZE(size) {                                              \
    err = driver.inq_put_size(&size);                                     \
    CHECK_ERR                                                             \
}
#define INQ_FILE_INFO(info) {                                             \
    err = driver.inq_file_info(ncid, &info);                              \
    CHECK_ERR                                                             \
}
#define IPUT_VARA_NOADV(itype, buf) {                                     \
    err = driver.put_vara(ncid, varid, itype, start, count, buf, nb);     \
    CHECK_VAR_ERR(varid)                                                  \
    my_nreqs++;                                                           \
    varid++;                                                              \
}
#define IPUT_VARA(varp, itype, adv, buf) {                                \
    if (cfg.api == adios)                                                 \
        err = e3sm_io_scorpio_write_var(driver, rec_no, ncid, varp,       \
                                        itype, buf, nb);                  \
    else                                                                  \
        err = driver.put_vara(ncid, varp.vid, itype, start, count,        \
                              buf, nb);                                   \
    CHECK_VAR_ERR(varp.vid)                                               \
    buf += (adv);                                                         \
    my_nreqs++;                                                           \
}
#define IPUT_VAR(varp, itype, adv, buf) {                                 \
    if (cfg.api == adios)                                                 \
        err = e3sm_io_scorpio_write_var(driver, -1, ncid, varp,           \
                                        itype, buf, nb);                  \
    else                                                                  \
        err = driver.put_vara(ncid, varp.vid, itype, NULL, NULL,          \
                              buf, nb);                                   \
    CHECK_VAR_ERR(varp.vid)                                               \
    buf += (adv);                                                         \
    my_nreqs++;                                                           \
}
#define WAIT_ALL_REQS {                                                   \
    err = driver.wait(ncid);                                              \
    CHECK_ERR                                                             \
    nflushes++;                                                           \
}

#define FIX_VAR_IPUT(varp, dp, itype, buf) {                              \
    if (cfg.strategy == canonical) {                                      \
        err = driver.put_varn(ncid, varp.vid, itype,                      \
                              decom.contig_nreqs[dp], decom.w_startx[dp], \
                              decom.w_countx[dp], buf, nb);               \
        my_nreqs += decom.contig_nreqs[dp];                               \
        buf += decom.count[dp] + gap;                                     \
    }                                                                     \
    else if (cfg.api == adios) {                                          \
        err = e3sm_io_scorpio_write_var(driver, -1, ncid, varp, itype,    \
                                        buf, nb);                         \
        my_nreqs++;                                                       \
        buf += decom.raw_nreqs[dp] + gap;                                 \
    }                                                                     \
    else {                                                                \
        err = driver.put_vara(ncid, varp.vid, itype, &decom.start[dp],    \
                              &decom.count[dp], buf, nb);                 \
        my_nreqs++;                                                       \
        buf += decom.count[dp] + gap;                                     \
    }                                                                     \
    CHECK_VAR_ERR(varp.vid)                                               \
    nvars_D[dp]++;                                                        \
}

#define REC_VAR_IPUT(varp, dp, itype, buf) {                              \
    if (cfg.strategy == canonical) {                                      \
        err = driver.put_varn(ncid, varp.vid, itype,                      \
                              decom.contig_nreqs[dp], decom.w_starts[dp], \
                              decom.w_counts[dp], buf, nb);               \
        my_nreqs += decom.contig_nreqs[dp];                               \
        buf += decom.count[dp] + gap;                                     \
    }                                                                     \
    else if (cfg.api == adios) {                                          \
        err = e3sm_io_scorpio_write_var(driver, rec_no, ncid, varp,       \
                                        itype, buf, nb);                  \
        my_nreqs++;                                                       \
        buf += decom.raw_nreqs[dp] + gap;                                 \
    }                                                                     \
    else {                                                                \
        err = driver.put_vara(ncid, varp.vid, itype, starts[dp],          \
                              counts[dp], buf, nb);                       \
        my_nreqs++;                                                       \
        buf += decom.count[dp] + gap;                                     \
    }                                                                     \
    CHECK_VAR_ERR(varp.vid)                                               \
    if (rec_no == 0) nvars_D[dp]++;                                       \
}

/*----< var_wr_case() >------------------------------------------------------*/
int e3sm_io_case::var_wr_case(e3sm_io_config &cfg,
                              e3sm_io_decom  &decom,
                              e3sm_io_driver &driver,
                              case_meta      *cmeta)
{
    char *fix_txt_buf_ptr, *rec_txt_buf_ptr;
    int i, j, err=0, sub_rank, global_rank, ncid=-1, nflushes=0;
    int rec_no, ffreq, gap=0, my_nreqs, nvars, num_decomp_vars, *nvars_D;
    int *fix_int_buf_ptr, *rec_int_buf_ptr;
    double *fix_dbl_buf_ptr, *rec_dbl_buf_ptr, timing;
    MPI_Offset previous_size, metadata_size, total_size;
    MPI_Comm comm;
    vtype *fix_buf_ptr, *rec_buf_ptr;

    int contig_nreqs[MAX_NUM_DECOMP];
    MPI_Offset blob_start[MAX_NUM_DECOMP], blob_count[MAX_NUM_DECOMP];
    MPI_Offset starts[MAX_NUM_DECOMP][2], counts[MAX_NUM_DECOMP][2];
    MPI_Offset start[2], count[2];

    MPI_Barrier(cfg.io_comm); /*---------------------------------------------*/
    cmeta->end2end_time = timing = MPI_Wtime();

    cmeta->post_time  = 0.0;
    cmeta->flush_time = 0.0;

    if (cfg.strategy == blob && cfg.api != adios)
        comm = cfg.sub_comm;
    else
        comm = cfg.io_comm;

    MPI_Comm_rank(cfg.io_comm,  &global_rank);
    MPI_Comm_rank(comm,         &sub_rank);

    /* I/O amount from previous I/O */
    INQ_PUT_SIZE(previous_size)

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
    if (cfg.strategy == blob) {
        if (cfg.api == adios)
            num_decomp_vars = decom.num_decomp + 1;
        else
            num_decomp_vars = decom.num_decomp * NVARS_DECOMP;
    }
    else
        num_decomp_vars = 0;

    nvars = cmeta->nvars + num_decomp_vars;

    /* allocate space to store variable metadata */
    vars = (var_meta*) malloc(nvars * sizeof(var_meta));

    if (cfg.non_contig_buf) gap = BUF_GAP;

    wr_buf_init(gap);

    if (cfg.run_case == F)
        err = def_F_case(cfg, decom, driver, ncid);
    else if (cfg.run_case == G)
        err = def_G_case(cfg, decom, driver, ncid);
    else if (cfg.run_case == I)
        err = def_I_case(cfg, decom, driver, ncid);
    CHECK_ERR

    /* exit define mode and enter data mode */
    ENDDEF

    /* I/O amount so far */
    INQ_PUT_SIZE(metadata_size)
    metadata_size -= previous_size;

    INQ_FILE_INFO(cmeta->info_used)

    cmeta->def_time = MPI_Wtime() - timing;

    MPI_Barrier(comm); /*----------------------------------------------------*/
    timing = MPI_Wtime();

    /* flush drequency only affects pnetcdf API.  Note for HDF5 and ADIOS blob
     * I/O, write data is copied into their internal buffers and only flushed
     * at file close. Calling driver.wait() takes no effect.
     */
    ffreq = (cfg.api == pnetcdf) ? cmeta->ffreq : 1;

    /* allocate write buffers */
    wr_buf_malloc(cfg, ffreq);

    nvars_D = cmeta->nvars_D;
    for (j=0; j<decom.num_decomp; j++)
        nvars_D[j] = 0; /* number of variables using decomposition j */

    cmeta->pre_time = MPI_Wtime() - timing;

    MPI_Barrier(comm); /*----------------------------------------------------*/
    timing = MPI_Wtime();

    my_nreqs = 0; /* number of noncontiguous requests written by this rank */

#ifdef ENABLE_ADIOS2
    /* write decomposition maps */
    if (cfg.api == adios) {
        for (j=0; j<decom.num_decomp; j++) {
            err = driver.put_varl(ncid, vars[j].vid, MPI_LONG_LONG,
                                  decom.raw_offsets[j], nb);
            CHECK_ERR
        }
        // Nproc only written by rank 0
        if (cfg.rank == 0) {
            err = driver.put_varl(ncid, vars[j].vid, MPI_INT, &cfg.np, nb);
            CHECK_ERR
        }
    } else
#endif
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

    fix_txt_buf_ptr = wr_buf.fix_txt_buf;
    fix_int_buf_ptr = wr_buf.fix_int_buf;
    fix_dbl_buf_ptr = wr_buf.fix_dbl_buf;
    fix_buf_ptr     = wr_buf.fix_buf;

    for (rec_no=0; rec_no<cmeta->nrecs; rec_no++) {

        if (rec_no % ffreq == 0) {
            /* reset the pointers to the beginning of the buffers */
            rec_txt_buf_ptr = wr_buf.rec_txt_buf;
            rec_int_buf_ptr = wr_buf.rec_int_buf;
            rec_dbl_buf_ptr = wr_buf.rec_dbl_buf;
            rec_buf_ptr     = wr_buf.rec_buf;
        }

        /* set the start index for the next record */
        if (cfg.api != adios) {
            for (i=0; i<decom.num_decomp; i++) {
                if (cfg.strategy == blob)
                    starts[i][0] = rec_no;
                else {
                    for (j=0; j<decom.contig_nreqs[i]; j++)
                        decom.w_starts[i][j][0] = rec_no;
                }
            }
        }

        /* write all climate variables */
        for (j=num_decomp_vars; j<nvars; j++) {
            int          dp    = vars[j].decomp_id;
            MPI_Datatype itype = vars[j].itype;
            size_t       adv   = vars[j].vlen + gap;

            if (vars[j].decomp_id >= 0) { /* this variable is partitioned */
                if (vars[j].isRecVar) { /* this is a record variable */
                    REC_VAR_IPUT(vars[j], dp, VAR_ITYPE, rec_buf_ptr)
                }
                else if (rec_no == 0) { /* this is a fixed-size variable */
                    if (itype == VAR_ITYPE)
                        FIX_VAR_IPUT(vars[j], dp, VAR_ITYPE, fix_buf_ptr)
                    else if (itype == MPI_INT)
                        FIX_VAR_IPUT(vars[j], dp, MPI_INT, fix_int_buf_ptr)
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
                        IPUT_VARA(vars[j], itype, adv,     rec_buf_ptr)
                    else if (itype == MPI_INT)
                        IPUT_VARA(vars[j], itype, adv, rec_int_buf_ptr)
                    else if (itype == MPI_CHAR)
                        IPUT_VARA(vars[j], itype, adv, rec_txt_buf_ptr)
                    else if (itype == MPI_DOUBLE)
                        IPUT_VARA(vars[j], itype, adv, rec_dbl_buf_ptr)
                }
                else if (rec_no == 0) { /* this is a fixed-size variable */
                    if (itype == VAR_ITYPE)
                        IPUT_VAR(vars[j], itype, adv,     fix_buf_ptr)
                    else if (itype == MPI_INT)
                        IPUT_VAR(vars[j], itype, adv, fix_int_buf_ptr)
                    else if (itype == MPI_CHAR)
                        IPUT_VAR(vars[j], itype, adv, fix_txt_buf_ptr)
                    else if (itype == MPI_DOUBLE)
                        IPUT_VAR(vars[j], itype, adv, fix_dbl_buf_ptr)
                }
            }
        }

        /* flush out the pending iput requests */
        if ((rec_no + 1) % ffreq == 0 || (rec_no + 1) == cmeta->nrecs) {
            cmeta->post_time += MPI_Wtime() - timing;

            MPI_Barrier(comm); /*--------------------------------------------*/
            timing = MPI_Wtime();

            /* flush once per time record */
            WAIT_ALL_REQS
            cmeta->flush_time += MPI_Wtime() - timing;

            if ((rec_no + 1) < cmeta->nrecs) timing = MPI_Wtime();
        }
    }

    MPI_Barrier(comm); /*----------------------------------------------------*/
    timing = MPI_Wtime();

    /* free up previously allocated heap memory space */
    wr_buf_free();
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

    /* note inquiring file size may be expensive on some machines */
    if (cfg.verbose && global_rank == 0)
        driver.inq_file_size(cmeta->outfile, &cmeta->file_size);

err_out:
    if (err < 0 && ncid >= 0) driver.close(ncid);
    if (!cfg.keep_outfile && sub_rank == 0) unlink(cmeta->outfile);

    return err;
}

