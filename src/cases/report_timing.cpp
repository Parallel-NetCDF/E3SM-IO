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

/*---< print_timing_WR() >---------------------------------------------------*/
static
int print_timing_WR(e3sm_io_config *cfg,
                    e3sm_io_decom  *decom,
                    perf_report    *pr)
{
    int i, global_rank;
    MPI_Offset off_msg[2], sum_off[2], max_off[2];
    MPI_Offset sum_nreqs, sum_amount_WR, max_nreqs;
    double pre_time, open_time, def_time, post_time, flush_time, close_time;
    double end2end_time, wTime;

    MPI_Comm_rank(cfg->io_comm, &global_rank);

    off_msg[0] = pr->my_nreqs;
    off_msg[1] = pr->amount_WR;
    MPI_Reduce(off_msg, sum_off, 2, MPI_OFFSET, MPI_SUM, 0, cfg->io_comm);
    sum_nreqs     = sum_off[0];
    sum_amount_WR = sum_off[1];

    off_msg[0] = pr->my_nreqs;
    MPI_Reduce(off_msg, max_off, 1, MPI_OFFSET, MPI_MAX, 0, cfg->io_comm);
    max_nreqs = max_off[0];

    double dbl_tmp[7], max_dbl[7];
    dbl_tmp[0] = pr->pre_time;
    dbl_tmp[1] = pr->open_time;
    dbl_tmp[2] = pr->def_time;
    dbl_tmp[3] = pr->post_time;
    dbl_tmp[4] = pr->flush_time;
    dbl_tmp[5] = pr->close_time;
    dbl_tmp[6] = pr->end2end_time;
    MPI_Reduce(dbl_tmp, max_dbl, 7, MPI_DOUBLE, MPI_MAX, 0, cfg->io_comm);
        pre_time = max_dbl[0];
       open_time = max_dbl[1];
        def_time = max_dbl[2];
       post_time = max_dbl[3];
      flush_time = max_dbl[4];
      close_time = max_dbl[5];
    end2end_time = max_dbl[6];

    if (global_rank == 0) {
        MPI_Offset decomp_amount, total_raw_nreqs;
        int nblobs = cfg->num_subfiles;
        int nvars_noD = pr->nvars;
        for (i=0; i<pr->num_decomp; i++) nvars_noD -= pr->nvars_D[i];

        if (cfg->run_case == F)
            printf("==== Benchmarking F case =============================\n");
        else if (cfg->run_case == G)
            printf("==== Benchmarking G case =============================\n");
        else if (cfg->run_case == I)
            printf("==== Benchmarking I case =============================\n");
        if (cfg->strategy == blob && cfg->api != adios)
            printf("%s\n", cfg->node_info);
        printf("Total number of MPI processes      = %d\n", cfg->np);
        printf("Number of IO processes             = %d\n", cfg->num_iotasks);
        printf("Input decomposition file           = %s\n", cfg->cfg_path);
        printf("Number of decompositions           = %d\n", decom->num_decomp);
        if (cfg->rd) {
            printf ("Input file/directory               = %s\n", cfg->in_path);
            printf ("Using noncontiguous read buffer    = %s\n", cfg->non_contig_buf ? "yes" : "no");
            printf("Variable read order: same as variables are defined\n");
        }
        if (cfg->wr) {
            printf ("Output file/directory              = %s\n",cfg->out_path);
            printf ("Using noncontiguous write buffer   = %s\n",
                    cfg->non_contig_buf ? "yes" : "no");
            printf("Variable write order: same as variables are defined\n");

            if (cfg->strategy == canonical) {
                if (cfg->api == pnetcdf)
                    printf("==== PnetCDF canonical I/O using varn API ============\n");
                else if (cfg->api == hdf5_md)
                    printf("==== HDF5 canonical I/O using multi-dataset API ======\n");
            }
            else if (cfg->strategy == log) {
                if (cfg->api == hdf5_log)
                    printf("==== HDF5 using log-based VOL ========================\n");
            }
            else if (cfg->strategy == blob) {
                if (cfg->api == pnetcdf)
                    printf("==== PnetCDF blob I/O ================================\n");
                else if (cfg->api == hdf5)
                    printf("==== HDF5 blob I/O ===================================\n");
                else if (cfg->api == adios)
                    printf("==== ADIOS blob I/O ==================================\n");
            }
        }

        if (cfg->api == pnetcdf || cfg->api == hdf5_log || cfg->api == hdf5_md)
            wTime = flush_time;
        else /* write happens at file close for hdf5 blob and adios blob */
            wTime = close_time;

        if (cfg->strategy == blob) {
            printf("History output name base           = %s\n", cfg->out_path);
            if (cfg->api == adios) {
                printf("History output folder name         = %s.bp.dir\n", pr->outfile);
                printf("History output subfile names       = %s.bp.dir/%s.bp.xxxx\n",
                       pr->outfile, pr->outfile);
                printf("Number of subfiles                 = %d\n", cfg->num_group);
                printf("Output file size                   = %.2f MiB = %.2f GiB\n",
                    (double)pr->file_size / 1048576, (double)pr->file_size / 1073741824);
            }
            else {
                printf("History output subfile names       = %s.xxxx\n", pr->outfile);
                printf("Number of subfiles                 = %3d\n", nblobs);
            }
            decomp_amount = 0;
            for (i=0; i<pr->num_decomp; i++) {
                decomp_amount += nblobs * sizeof(int); /* D*.nreqs */
                decomp_amount += nblobs * sizeof(MPI_Offset); /* D*.blob_start */
                decomp_amount += nblobs * sizeof(MPI_Offset); /* D*.blob_count */
                decomp_amount += nblobs * decom->nelems[i] * sizeof(int); /* D*.offsets */
                decomp_amount += nblobs * decom->nelems[i] * sizeof(int); /* D*.lengths */
            }
            total_raw_nreqs = 0;
            for (i=0; i<pr->num_decomp; i++)
                total_raw_nreqs += decom->total_raw_nreqs[i];

            printf("No. decomposition variables        = %3d\n", pr->num_decomp_vars);
            if (cfg->api == adios)
                printf("Size of raw decomposition maps     = %.2f MiB\n", (float)total_raw_nreqs/1048576.0);
            else
                printf("Size of decomposition variables    = %.2f MiB\n", (float)decomp_amount/1048576.0);
        }
        else
            printf("History output file                = %s\n", pr->outfile);
        printf("No. variables use no decomposition = %6d\n", nvars_noD);
        for (i=0; i<pr->num_decomp; i++)
            printf("No. variables use decomposition D%d = %6d\n",
                   i, pr->nvars_D[i]);
        printf("Total no. climate variables        = %6d\n", pr->nvars);
        printf("Write no. records (time dim)       = %6d\n", pr->nrecs);
        printf("Total no. noncontiguous requests   = %6lld\n", sum_nreqs);
        printf("Max   no. noncontiguous requests   = %6lld\n", max_nreqs);
        printf("No. I/O flush calls                = %6d\n", pr->num_flushes);
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
        printf("I/O bandwidth (write-only)         = %.4f MiB/sec\n",
               (double)sum_amount_WR / 1048576.0 / wTime);
        printf("I/O bandwidth (open-to-close)      = %.4f MiB/sec\n",
               (double)sum_amount_WR / 1048576.0 / end2end_time);
        printf("-----------------------------------------------------------\n");
    }
    fflush(stdout);

    return 0;
}

/*---< report_timing_WR() >--------------------------------------------------*/
int report_timing_WR(e3sm_io_config *cfg,
                     e3sm_io_decom  *decom)
{
    if (cfg->run_case == F) {
        print_timing_WR(cfg, decom, &cfg->F_case_h0);
        print_timing_WR(cfg, decom, &cfg->F_case_h1);
    }
    else if (cfg->run_case == G)
        print_timing_WR(cfg, decom, &cfg->G_case);
    else if (cfg->run_case == I) {
        print_timing_WR(cfg, decom, &cfg->I_case_h0);
        print_timing_WR(cfg, decom, &cfg->I_case_h1);
    }
    return 0;
}
