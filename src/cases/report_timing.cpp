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
#include <libgen.h> /* basename() */
#include <assert.h>

#include <mpi.h>
#include <e3sm_io.h>

/*---< print_timing_WR() >---------------------------------------------------*/
static
int print_timing_WR(e3sm_io_config *cfg,
                    e3sm_io_decom  *decom,
                    case_meta      *cmeta)
{
    int i, ndecomp;
    MPI_Offset off_msg[3], sum_off[3];
    MPI_Offset sum_nreqs, sum_amount_WR, max_nreqs, min_nreqs, max_wr, min_wr;
    MPI_Offset vlen, sum_decomp_varlen;
    double wTime;
    MPI_Comm comm=cfg->io_comm;

    ndecomp = decom->num_decomp;

    /* calculate the space occupied by the decomposition variables */
    sum_decomp_varlen = 0;
    for (i=0; i<ndecomp; i++) {
        vlen = sizeof(int)                        /* D*.nreqs */
             + sizeof(MPI_Offset)                 /* D*.blob_start */
             + sizeof(MPI_Offset)                 /* D*.blob_count */
             + decom->max_nreqs[i] * sizeof(int)  /* D*.offsets */
             + decom->max_nreqs[i] * sizeof(int); /* D*.lengths */
        sum_decomp_varlen += vlen;
    }
    sum_decomp_varlen *= cfg->sub_nprocs;
    if (cfg->strategy == blob && cfg->api != adios && cfg->sub_rank > 0)
        sum_decomp_varlen = 0;
    /* sum_decomp_varlen is the size of all decomposition variables defined in
     * the subfile of this process. To calculate the total sizes across all
     * subfiles, only root ranks of subfiles do a parallel sum to obtain the
     * total size.
     */

    off_msg[0] = cmeta->my_nreqs;
    off_msg[1] = cmeta->amount_WR;
    off_msg[2] = sum_decomp_varlen;
    MPI_Reduce(off_msg, sum_off, 3, MPI_OFFSET, MPI_SUM, 0, comm);
    sum_nreqs         = sum_off[0];
    sum_amount_WR     = sum_off[1];
    sum_decomp_varlen = sum_off[2];

    MPI_Reduce(&cmeta->my_nreqs, &max_nreqs, 1, MPI_OFFSET, MPI_MAX, 0, comm);
    MPI_Reduce(&cmeta->my_nreqs, &min_nreqs, 1, MPI_OFFSET, MPI_MIN, 0, comm);
    MPI_Reduce(&cmeta->amount_WR, &max_wr,   1, MPI_OFFSET, MPI_MAX, 0, comm);
    MPI_Reduce(&cmeta->amount_WR, &min_wr,   1, MPI_OFFSET, MPI_MIN, 0, comm);

    double dbl_tmp[7], max_dbl[7], min_dbl[7];
    dbl_tmp[0] = cmeta->pre_time;
    dbl_tmp[1] = cmeta->open_time;
    dbl_tmp[2] = cmeta->def_time;
    dbl_tmp[3] = cmeta->post_time;
    dbl_tmp[4] = cmeta->flush_time;
    dbl_tmp[5] = cmeta->close_time;
    dbl_tmp[6] = cmeta->end2end_time;
    MPI_Reduce(dbl_tmp, max_dbl, 7, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(dbl_tmp, min_dbl, 7, MPI_DOUBLE, MPI_MIN, 0, comm);

    if (cfg->rank == 0) {
        int nvars_noD = cmeta->nvars;
        for (i=0; i<ndecomp; i++) nvars_noD -= cmeta->nvars_D[i];

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
        printf("Input decomposition file           = %s\n", cfg->decomp_path);
        printf("Number of decompositions           = %d\n", ndecomp);
        printf("Output file/directory              = %s\n", cfg->out_path);
        printf("Fill mode                          = %s\n", cfg->fill_mode ? "enabled" : "disabled");
        printf("Using noncontiguous write buffer   = %s\n", cfg->non_contig_buf ? "yes" : "no");
        printf("Variable write order: same as variables are defined\n");

        if (cfg->strategy == canonical) {
            if (cfg->api == pnetcdf)
                printf("==== PnetCDF canonical I/O using varn API ============\n");
            else if (cfg->api == hdf5)
                printf("==== HDF5 canonical I/O ==============================\n");
            else if (cfg->api == hdf5_md)
                printf("==== HDF5 canonical I/O using multi-dataset API ======\n");
            else if (cfg->api == netcdf4)
                printf("==== NetCDF-4 canonical I/O ==========================\n");
        }
        else if (cfg->strategy == log) {
            if (cfg->api == hdf5)
                printf("==== HDF5 using log-based VOL through native APIs=====\n");
            else if (cfg->api == hdf5_log)
                printf("==== HDF5 using log-based VOL APIs ===================\n");
            else if (cfg->api == netcdf4)
                printf("==== NetCDF-4 using HDF5 log-based VOL ===============\n");
        }
        else if (cfg->strategy == blob) {
            if (cfg->api == pnetcdf)
                printf("==== PnetCDF blob I/O ================================\n");
            else if (cfg->api == hdf5)
                printf("==== HDF5 blob I/O ===================================\n");
            else if (cfg->api == adios)
                printf("==== ADIOS blob I/O ==================================\n");
        }

        if (cfg->strategy == canonical && (cfg->api == netcdf4 || cfg->api == hdf5))
            wTime = max_dbl[3];
        else if (cfg->strategy == log && cfg->api == netcdf4)
            wTime = max_dbl[5];
        else if (cfg->api == pnetcdf || cfg->strategy == log || cfg->api == hdf5_md)
            wTime = max_dbl[4];
        else /* write happens at file close for hdf5 blob and adios blob */
            wTime = max_dbl[5];

        if (cfg->strategy == log && (cfg->api == hdf5_log || cfg->api == netcdf4)) {
            if (cfg->num_subfiles != 0) {
                printf("History output folder names        = %s.subfiles\n", cmeta->outfile);
                printf("History output subfile names       = %s.subfiles/%s.xxxx\n",
                       cmeta->outfile, basename(cmeta->outfile));
            }
            printf("Number of subfiles                 = %d\n", cfg->num_subfiles);
        }
        else if (cfg->strategy == blob) {
            printf("History output name base           = %s\n", cfg->out_path);
            if (cfg->api == adios) {
                printf("History output folder name         = %s\n", cmeta->outfile);
                printf("History output subfile names       = %s/dat.xxxx\n", cmeta->outfile);
                printf("Number of subfiles                 = %d\n", cfg->num_subfiles);
                if (cfg->verbose)
                    printf("Output file size                   = %.2f MiB = %.2f GiB\n",
                        (double)cmeta->file_size / 1048576, (double)cmeta->file_size / 1073741824);
            }
            else {
                printf("History output subfile names       = %s.xxxx\n", cmeta->outfile);
                printf("Number of subfiles                 = %d\n", cfg->num_subfiles);
            }
            printf("No. decomposition variables        = %d\n", cmeta->num_decomp_vars);
            if (cfg->api == adios) {
                MPI_Offset amount = 0;
                for (i=0; i<ndecomp; i++)
                    amount += decom->raw_nreqs[i];
                amount *= sizeof(int64_t);
                printf("Size of raw decomposition maps     = %.2f MiB\n",
                       (float)amount/1048576.0);
            }
            else {
                printf("Size of decomposition variables    = %.2f MiB\n",
                       (float)sum_decomp_varlen/1048576.0);
            }
        }
        else
            printf("History output file                = %s\n", cmeta->outfile);
        printf("No. variables use no decomposition = %6d\n", nvars_noD);
        for (i=0; i<ndecomp; i++)
            printf("No. variables use decomposition D%d = %6d\n",
                   i, cmeta->nvars_D[i]);
        printf("Total no. climate variables        = %6d\n", cmeta->nvars);
        printf("Total no. attributes               = %6d\n", cmeta->num_attrs);
        printf("Total no. noncontiguous requests   = %6lld\n", sum_nreqs);
        printf("Max   no. noncontiguous requests   = %6lld\n", max_nreqs);
        printf("Min   no. noncontiguous requests   = %6lld\n", min_nreqs);
        printf("Write no. records (time dim)       = %6d\n", cmeta->nrecs);
        printf("I/O flush frequency                = %6d\n", cmeta->ffreq);
        printf("No. I/O flush calls                = %6d\n", cmeta->num_flushes);
        printf("-----------------------------------------------------------\n");
        printf("Indiv write amount MAX             = %.2f MiB = %.2f GiB\n",
               (double)max_wr / 1048576, (double)max_wr / 1073741824);
        printf("Indiv write amount MIN             = %.2f MiB = %.2f GiB\n",
               (double)min_wr / 1048576, (double)min_wr / 1073741824);
        printf("Total write amount                 = %.2f MiB = %.2f GiB\n",
               (double)sum_amount_WR / 1048576, (double)sum_amount_WR / 1073741824);
        printf("-----------------------------------------------------------\n");
        printf("Time of I/O preparing              min/max = %8.4f / %8.4f\n", min_dbl[0], max_dbl[0]);
        printf("Time of file open/create           min/max = %8.4f / %8.4f\n", min_dbl[1], max_dbl[1]);
        printf("Time of define variables           min/max = %8.4f / %8.4f\n", min_dbl[2], max_dbl[2]);
        printf("Time of posting write requests     min/max = %8.4f / %8.4f\n", min_dbl[3], max_dbl[3]);
        printf("Time of write flushing             min/max = %8.4f / %8.4f\n", min_dbl[4], max_dbl[4]);
        printf("Time of close                      min/max = %8.4f / %8.4f\n", min_dbl[5], max_dbl[5]);
        printf("end-to-end time                    min/max = %8.4f / %8.4f\n", min_dbl[6], max_dbl[6]);
        printf("Emulate computation time (sleep)   min/max = %8.4f / %8.4f\n", (double)(cfg->comp_time), (double)(cfg->comp_time));
        printf("I/O bandwidth in MiB/sec (write-only)      = %.4f\n",
               (double)sum_amount_WR / 1048576.0 / wTime);
        printf("I/O bandwidth in MiB/sec (open-to-close)   = %.4f\n",
               (double)sum_amount_WR / 1048576.0 / max_dbl[6]);
        printf("-----------------------------------------------------------\n");

        /* print MPI-IO hints actually used */
        if (cfg->verbose) print_info(&cmeta->info_used);
    }
    fflush(stdout);

    if (cmeta->info_used != MPI_INFO_NULL) MPI_Info_free(&cmeta->info_used);

    return 0;
}

/*---< report_timing_WR() >--------------------------------------------------*/
int report_timing_WR(e3sm_io_config *cfg,
                     e3sm_io_decom  *decom)
{
    if (cfg->run_case == F) {
        if (cfg->hx == 0 || cfg->hx == -1) /* h0 file */
            print_timing_WR(cfg, decom, &cfg->F_case_h0);
        if (cfg->hx == 1 || cfg->hx == -1) /* h1 file */
            print_timing_WR(cfg, decom, &cfg->F_case_h1);
    }
    else if (cfg->run_case == G)
        print_timing_WR(cfg, decom, &cfg->G_case);
    else if (cfg->run_case == I) {
        if (cfg->hx == 0 || cfg->hx == -1) /* h0 file */
            print_timing_WR(cfg, decom, &cfg->I_case_h0);
        if (cfg->hx == 1 || cfg->hx == -1) /* h1 file */
            print_timing_WR(cfg, decom, &cfg->I_case_h1);
    }
    return 0;
}

/*---< print_timing_RD() >---------------------------------------------------*/
static
int print_timing_RD(e3sm_io_config *cfg,
                    e3sm_io_decom  *decom,
                    case_meta      *cmeta)
{
    int i, ndecomp;
    MPI_Offset off_msg[3], sum_off[3];
    MPI_Offset sum_nreqs, sum_amount_RD, max_nreqs, min_nreqs;
    MPI_Offset vlen, sum_decomp_varlen;
    double wTime;
    MPI_Comm comm=cfg->io_comm;

    ndecomp = decom->num_decomp;

    /* calculate the space occupied by the decomposition variables */
    sum_decomp_varlen = 0;
    for (i=0; i<ndecomp; i++) {
        vlen = sizeof(int)                        /* D*.nreqs */
             + sizeof(MPI_Offset)                 /* D*.blob_start */
             + sizeof(MPI_Offset)                 /* D*.blob_count */
             + decom->max_nreqs[i] * sizeof(int)  /* D*.offsets */
             + decom->max_nreqs[i] * sizeof(int); /* D*.lengths */
        sum_decomp_varlen += vlen;
    }
    sum_decomp_varlen *= cfg->sub_nprocs;
    if (cfg->strategy == blob && cfg->api != adios && cfg->sub_rank > 0)
        sum_decomp_varlen = 0;
    /* sum_decomp_varlen is the size of all decomposition variables defined in
     * the subfile of this process. To calculate the total sizes across all
     * subfiles, only root ranks of subfiles do a parallel sum to obtain the
     * total size.
     */

    off_msg[0] = cmeta->my_nreqs;
    off_msg[1] = cmeta->amount_RD;
    off_msg[2] = sum_decomp_varlen;
    MPI_Reduce(off_msg, sum_off, 3, MPI_OFFSET, MPI_SUM, 0, comm);
    sum_nreqs         = sum_off[0];
    sum_amount_RD     = sum_off[1];
    sum_decomp_varlen = sum_off[2];

    MPI_Reduce(&cmeta->my_nreqs, &max_nreqs, 1, MPI_OFFSET, MPI_MAX, 0, comm);
    MPI_Reduce(&cmeta->my_nreqs, &min_nreqs, 1, MPI_OFFSET, MPI_MIN, 0, comm);

    double dbl_tmp[7], max_dbl[7], min_dbl[7];
    dbl_tmp[0] = cmeta->pre_time;
    dbl_tmp[1] = cmeta->open_time;
    dbl_tmp[2] = cmeta->def_time;
    dbl_tmp[3] = cmeta->post_time;
    dbl_tmp[4] = cmeta->flush_time;
    dbl_tmp[5] = cmeta->close_time;
    dbl_tmp[6] = cmeta->end2end_time;
    MPI_Reduce(dbl_tmp, max_dbl, 7, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(dbl_tmp, min_dbl, 7, MPI_DOUBLE, MPI_MIN, 0, comm);

    if (cfg->rank == 0) {
        int nvars_noD = cmeta->nvars;
        for (i=0; i<ndecomp; i++) nvars_noD -= cmeta->nvars_D[i];

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
        printf("Input decomposition file           = %s\n", cfg->decomp_path);
        printf("Number of decompositions           = %d\n", ndecomp);
        printf("Input file/directory               = %s\n", cfg->in_path);
        printf("Using noncontiguous read buffer    = %s\n", cfg->non_contig_buf ? "yes" : "no");
        printf("Variable read order: same as variables are defined\n");

        if (cfg->strategy == canonical) {
            if (cfg->api == pnetcdf)
                printf("==== PnetCDF canonical I/O using varn API ============\n");
            else if (cfg->api == hdf5)
                printf("==== HDF5 canonical I/O ==============================\n");
            else if (cfg->api == hdf5_md)
                printf("==== HDF5 canonical I/O using multi-dataset API ======\n");
            else if (cfg->api == netcdf4)
                printf("==== NetCDF-4 canonical I/O ==========================\n");
        }
        else if (cfg->strategy == log) {
            if (cfg->api == hdf5)
                printf("==== HDF5 using log-based VOL through native APIs=====\n");
            else if (cfg->api == hdf5_log)
                printf("==== HDF5 using log-based VOL APIs ===================\n");
            else if (cfg->api == netcdf4)
                printf("==== NetCDF-4 using HDF5 log-based VOL ===============\n");
        }
        else if (cfg->strategy == blob) {
            if (cfg->api == pnetcdf)
                printf("==== PnetCDF blob I/O ================================\n");
            else if (cfg->api == hdf5)
                printf("==== HDF5 blob I/O ===================================\n");
            else if (cfg->api == adios)
                printf("==== ADIOS blob I/O ==================================\n");
        }

        if (cfg->strategy == canonical && (cfg->api == netcdf4 || cfg->api == hdf5))
            wTime = max_dbl[3];
        else if (cfg->strategy == log && cfg->api == netcdf4)
            wTime = max_dbl[5];
        else if (cfg->api == pnetcdf || cfg->strategy == log || cfg->api == hdf5_md)
            wTime = max_dbl[4];
        else /* read happens at file close for hdf5 blob and adios blob */
            wTime = max_dbl[5];

        if (cfg->strategy == log && (cfg->api == hdf5_log || cfg->api == netcdf4)) {
            if (cfg->num_subfiles != 0) {
                printf("History input folder names         = %s.subfiles\n", cmeta->outfile);
                printf("History input subfile names        = %s.subfiles/%s.xxxx\n",
                       cmeta->outfile, basename(cmeta->outfile));
            }
            printf("Number of subfiles                 = %d\n", cfg->num_subfiles);
        }
        else if (cfg->strategy == blob) {
            printf("History input name base            = %s\n", cfg->out_path);
            if (cfg->api == adios) {
                printf("History input folder name          = %s\n", cmeta->outfile);
                printf("History input subfile names        = %s/dat.xxxx\n", cmeta->outfile);
                printf("Number of subfiles                 = %d\n", cfg->num_subfiles);
                if (cfg->verbose)
                    printf("Input file size                    = %.2f MiB = %.2f GiB\n",
                        (double)cmeta->file_size / 1048576, (double)cmeta->file_size / 1073741824);
            }
            else {
                printf("History input subfile names        = %s.xxxx\n", cmeta->outfile);
                printf("Number of subfiles                 = %d\n", cfg->num_subfiles);
            }
            printf("No. decomposition variables        = %d\n", cmeta->num_decomp_vars);
            if (cfg->api == adios) {
                MPI_Offset amount = 0;
                for (i=0; i<ndecomp; i++)
                    amount += decom->raw_nreqs[i];
                amount *= sizeof(int64_t);
                printf("Size of raw decomposition maps     = %.2f MiB\n",
                       (float)amount/1048576.0);
            }
            else {
                printf("Size of decomposition variables    = %.2f MiB\n",
                       (float)sum_decomp_varlen/1048576.0);
            }
        }
        else
            printf("History input file                 = %s\n", cmeta->outfile);
        printf("No. variables use no decomposition = %6d\n", nvars_noD);
        for (i=0; i<ndecomp; i++)
            printf("No. variables use decomposition D%d = %6d\n",
                   i, cmeta->nvars_D[i]);
        printf("Total no. climate variables        = %6d\n", cmeta->nvars);
        printf("Total no. attributes               = %6d\n", cmeta->num_attrs);
        printf("Total no. noncontiguous requests   = %6lld\n", sum_nreqs);
        printf("Max   no. noncontiguous requests   = %6lld\n", max_nreqs);
        printf("Min   no. noncontiguous requests   = %6lld\n", min_nreqs);
        printf("Write no. records (time dim)       = %6d\n", cmeta->nrecs);
        printf("I/O flush frequency                = %6d\n", cmeta->ffreq);
        printf("No. I/O flush calls                = %6d\n", cmeta->num_flushes);
        printf("-----------------------------------------------------------\n");
        printf("Total read amount                          = %.2f MiB = %.2f GiB\n",
               (double)sum_amount_RD / 1048576, (double)sum_amount_RD / 1073741824);
        printf("Time of I/O preparing              min/max = %8.4f / %8.4f\n", min_dbl[0], max_dbl[0]);
        printf("Time of file open/create           min/max = %8.4f / %8.4f\n", min_dbl[1], max_dbl[1]);
        printf("Time of define variables           min/max = %8.4f / %8.4f\n", min_dbl[2], max_dbl[2]);
        printf("Time of posting read requests      min/max = %8.4f / %8.4f\n", min_dbl[3], max_dbl[3]);
        printf("Time of read flushing              min/max = %8.4f / %8.4f\n", min_dbl[4], max_dbl[4]);
        printf("Time of close                      min/max = %8.4f / %8.4f\n", min_dbl[5], max_dbl[5]);
        printf("end-to-end time                    min/max = %8.4f / %8.4f\n", min_dbl[6], max_dbl[6]);
        printf("Emulate computation time (sleep)   min/max = %8.4f / %8.4f\n", (double)(cfg->comp_time), (double)(cfg->comp_time));
        printf("I/O bandwidth in MiB/sec (read-only)       = %.4f\n",
               (double)sum_amount_RD / 1048576.0 / wTime);
        printf("I/O bandwidth in MiB/sec (open-to-close)   = %.4f\n",
               (double)sum_amount_RD / 1048576.0 / max_dbl[6]);
        printf("-----------------------------------------------------------\n");

        /* print MPI-IO hints actually used */
        if (cfg->verbose) print_info(&cmeta->info_used);
    }
    fflush(stdout);

    if (cmeta->info_used != MPI_INFO_NULL) MPI_Info_free(&cmeta->info_used);

    return 0;
}

/*---< report_timing_RD() >--------------------------------------------------*/
int report_timing_RD(e3sm_io_config *cfg,
                     e3sm_io_decom  *decom)
{
    if (cfg->run_case == F) {
        if (cfg->hx == 0 || cfg->hx == -1) /* h0 file */
            print_timing_RD(cfg, decom, &cfg->F_case_h0);
        if (cfg->hx == 1 || cfg->hx == -1) /* h1 file */
            print_timing_RD(cfg, decom, &cfg->F_case_h1);
    }
    else if (cfg->run_case == G)
        print_timing_RD(cfg, decom, &cfg->G_case);
    else if (cfg->run_case == I) {
        if (cfg->hx == 0 || cfg->hx == -1) /* h0 file */
            print_timing_RD(cfg, decom, &cfg->I_case_h0);
        if (cfg->hx == 1 || cfg->hx == -1) /* h1 file */
            print_timing_RD(cfg, decom, &cfg->I_case_h1);
    }
    return 0;
}
