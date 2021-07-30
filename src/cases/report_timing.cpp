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

#include <mpi.h>
#include <e3sm_io.h>
#include <e3sm_io_driver.hpp>

/*----< check_malloc() >-----------------------------------------------------*/
int check_malloc(e3sm_io_config *cfg,
                 e3sm_io_driver *driver)
{
    int err=0, global_rank;
    MPI_Offset m_alloc, s_alloc, x_alloc;

    if (!cfg->verbose || cfg->api != pnetcdf) return 0;

    MPI_Comm_rank(cfg->io_comm, &global_rank);

    /* check if there is any PnetCDF internal malloc residue */
    err = driver->inq_malloc_size(&m_alloc);
    if (err == NC_NOERR) {
        MPI_Reduce(&m_alloc, &s_alloc, 1, MPI_OFFSET, MPI_SUM, 0, cfg->io_comm);
        if (global_rank == 0 && s_alloc > 0) {
            printf("-------------------------------------------------------\n");
            printf("Residue heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   s_alloc);
        }
    }

    /* find the high water mark among all processes */
    driver->inq_malloc_max_size(&m_alloc);
    MPI_Reduce(&m_alloc, &x_alloc, 1, MPI_OFFSET, MPI_MAX, 0, cfg->io_comm);
    if (global_rank == 0)
        printf("High water mark of heap memory allocated by PnetCDF internally is %.2f MiB\n",
               (float)x_alloc / 1048576);

    return err;
}

/*---< report_timing_WR() >--------------------------------------------------*/
int report_timing_WR(e3sm_io_config *cfg,
                     e3sm_io_driver *driver,
                     e3sm_io_decom  *decom,
                     const char     *outfile)
{
    int i, err=0, global_rank;
    MPI_Offset off_msg[2], sum_off[2], max_off[2];
    MPI_Offset sum_nreqs, sum_amount_WR, max_nreqs, file_size=0;
    double pre_time, open_time, def_time, post_time, flush_time, close_time;
    double end2end_time, wTime;

    MPI_Comm_rank(cfg->io_comm, &global_rank);

    off_msg[0] = cfg->my_nreqs;
    off_msg[1] = cfg->amount_WR;
    MPI_Reduce(off_msg, sum_off, 2, MPI_OFFSET, MPI_SUM, 0, cfg->io_comm);
    sum_nreqs     = sum_off[0];
    sum_amount_WR = sum_off[1];

    off_msg[0] = cfg->my_nreqs;
    MPI_Reduce(off_msg, max_off, 1, MPI_OFFSET, MPI_MAX, 0, cfg->io_comm);
    max_nreqs = max_off[0];

    double dbl_tmp[7], max_dbl[7];
    dbl_tmp[0] = cfg->pre_time;
    dbl_tmp[1] = cfg->open_time;
    dbl_tmp[2] = cfg->def_time;
    dbl_tmp[3] = cfg->post_time;
    dbl_tmp[4] = cfg->flush_time;
    dbl_tmp[5] = cfg->close_time;
    dbl_tmp[6] = cfg->end2end_time;
    MPI_Reduce(dbl_tmp, max_dbl, 7, MPI_DOUBLE, MPI_MAX, 0, cfg->io_comm);
        pre_time = max_dbl[0];
       open_time = max_dbl[1];
        def_time = max_dbl[2];
       post_time = max_dbl[3];
      flush_time = max_dbl[4];
      close_time = max_dbl[5];
    end2end_time = max_dbl[6];

    if (global_rank == 0) {
        int nblobs = cfg->num_subfiles;
        MPI_Offset decomp_amount;
        int nvars_noD = cfg->nvars;
        for (i=0; i<cfg->num_decomp; i++) nvars_noD -= cfg->nvars_D[i];
        int nvars = cfg->nvars + nvars_noD;

        if (cfg->api == pnetcdf || cfg->api == hdf5_log || cfg->api == hdf5_md)
            wTime = flush_time;
        else /* write happens at file close for hdf5 blob and adios blob */
            wTime = close_time;

        if (cfg->strategy == blob) {
            printf("History output name base           = %s\n", cfg->out_path);
            if (cfg->api == adios) {
                err = driver->inq_file_size(outfile, &file_size);
                if (err != 0) /* ignore non-critical error */
                    printf("Error: failed inq_file_size %s\n",outfile);
                printf("History output folder name         = %s.bp.dir\n", outfile);
                printf("History output subfile names       = %s.bp.dir/%s.bp.xxxx\n",
                       outfile, outfile);
                printf("Number of subfiles                 = %d\n", cfg->num_group);
                printf("Output file size                   = %.2f MiB = %.2f GiB\n",
                    (double)file_size / 1048576, (double)file_size / 1073741824);
            }
            else {
                printf("History output subfile names       = %s.xxxx\n", outfile);
                printf("Number of subfiles                 = %3d\n", nblobs);
            }
            decomp_amount = 0;
            for (i=0; i<cfg->num_decomp; i++) {
                decomp_amount += nblobs * sizeof(int); /* D*.nreqs */
                decomp_amount += nblobs * sizeof(MPI_Offset); /* D*.blob_start */
                decomp_amount += nblobs * sizeof(MPI_Offset); /* D*.blob_count */
                decomp_amount += nblobs * decom->nelems[i] * sizeof(int); /* D*.offsets */
                decomp_amount += nblobs * decom->nelems[i] * sizeof(int); /* D*.lengths */
            }
            printf("No. decomposition variables        = %3d\n", cfg->num_decomp_vars);
            printf("Size of decomposition variables    = %.2f MiB\n", (float)decomp_amount/1048576.0);
        }
        else
            printf("History output file                = %s\n", outfile);
        printf("No. variables use no decomposition = %3d\n", nvars_noD);
        for (i=0; i<cfg->num_decomp; i++)
            printf("No. variables use decomposition D%d = %3d\n",
                   i, cfg->nvars_D[i]);
        printf("Total number of variables          = %3d\n", nvars);
        printf("Write number of records (time dim) = %3d\n", cfg->nrecs);
        printf("Total no. noncontiguous requests   = %3lld\n", sum_nreqs);
        printf("Max   no. noncontiguous requests   = %3lld\n", max_nreqs);
        printf("No. I/O flush calls                = %3d\n", cfg->num_flushes);
        printf("-----------------------------------------------------------\n");
        printf("Total write amount                 = %.2f MiB = %.2f GiB\n",
               (double)sum_amount_WR / 1048576, (double)sum_amount_WR / 1073741824);
        printf("Max Time of I/O preparing          = %.4f sec\n",   pre_time);
        printf("Max Time of file open/create       = %.4f sec\n",  open_time);
        printf("Max Time of define variables       = %.4f sec\n",   def_time);
        printf("Max Time of posting write requests = %.4f sec\n",  post_time);
        printf("Max Time of write flushing         = %.4f sec\n",  flush_time);
        printf("Max Time of close                  = %.4f sec\n", close_time);
        printf("Max end-to-end time                = %.4f sec\n", end2end_time);
        printf("I/O bandwidth (open-to-close)      = %.4f MiB/sec\n",
               (double)sum_amount_WR / 1048576.0 / end2end_time);
        printf("I/O bandwidth (write-only)         = %.4f MiB/sec\n",
                   (double)sum_amount_WR / 1048576.0 / wTime);
        printf("-----------------------------------------------------------\n");
    }
    fflush(stdout);

    return 0;
}

