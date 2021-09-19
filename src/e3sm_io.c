/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * See README.md for compile and run instructions.
 *
 *********************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */

#include <mpi.h>

#include <e3sm_io.h>
#include <e3sm_io_err.h>
#include <e3sm_io_profile.hpp>

static inline int set_info (e3sm_io_config *cfg, e3sm_io_decom *decom) {
    int err;
    MPI_Offset estimated_nc_ibuf_size;

    /* set MPI-IO hints */

    /* collective write */
    err = MPI_Info_set (cfg->info, "romio_cb_write", "enable");
    CHECK_MPIERR
    /* no independent MPI-IO */
    err = MPI_Info_set (cfg->info, "romio_no_indep_rw", "true");
    CHECK_MPIERR

    /* set PnetCDF I/O hints */

    /* no gap between variables */
    err = MPI_Info_set (cfg->info, "nc_var_align_size", "1");
    CHECK_MPIERR
    /* in-place byte swap */
    err = MPI_Info_set (cfg->info, "nc_in_place_swap", "enable");
    CHECK_MPIERR

    /* use total write amount to estimate nc_ibuf_size */
    estimated_nc_ibuf_size =
        decom->dims[2][0] * decom->dims[2][1] * sizeof (double) / cfg->num_iotasks;
    estimated_nc_ibuf_size *= cfg->nvars;
    if (estimated_nc_ibuf_size > 16777216) {
        char nc_ibuf_size_str[32];
        sprintf (nc_ibuf_size_str, "%lld", estimated_nc_ibuf_size);
        err = MPI_Info_set (cfg->info, "nc_ibuf_size", nc_ibuf_size_str);
        CHECK_MPIERR
    }

err_out:
    return err;
}

/*----< print_info() >------------------------------------------------------*/
void print_info (MPI_Info *info_used) {
    int i, nkeys;

    if (*info_used == MPI_INFO_NULL) {
        printf ("MPI File Info is NULL\n");
        return;
    }
    MPI_Info_get_nkeys (*info_used, &nkeys);
    printf ("MPI File Info: nkeys = %d\n", nkeys);
    for (i = 0; i < nkeys; i++) {
        char key[MPI_MAX_INFO_KEY], value[MPI_MAX_INFO_VAL];
        int valuelen, flag;

        MPI_Info_get_nthkey (*info_used, i, key);
        MPI_Info_get_valuelen (*info_used, key, &valuelen, &flag);
        MPI_Info_get (*info_used, key, valuelen + 1, value, &flag);
        printf ("MPI File Info: [%2d] key = %25s, value = %s\n", i, key, value);
    }
    printf("-----------------------------------------------------------\n");
}

/*----< usage() >------------------------------------------------------------*/
static void usage (char *argv0) {
    char *help =
"Usage: %s [OPTION]... FILE\n"
"       [-h] Print this help message\n"
"       [-v] Verbose mode\n"
"       [-k] Keep the output files when program exits\n"
"       [-m] Run test using noncontiguous write buffer\n"
"       [-f num] Set history output files h0 and/or h1: 0 for h0 only, 1 for\n"
"                h1 only, -1 for both (default: -1)\n"
"       [-r num] Number of records/time steps for F case h1 file (default: 1)\n"
"       [-y num] Frequency of I/O flush (affect PnetCDF API only). 1 for once\n"
"                per time step (default), -1 for once for all time steps\n"
"       [-s num] MPI rank stride for selecting processes to perform I/O tasks\n"
"                (default: 1)\n"
"       [-g num] Number of subfiles, used in blob I/O only (default: 1)\n"
"       [-i path] Enable read performance evaluation and set the input file\n"
"                 (folder) path\n"
"       [-o path] Enable write performance evaluation and set the output file\n"
"                 (folder) path\n"
"       [-a api]  I/O library name to perform write operation\n"
"           pnetcdf:   PnetCDF library (default)\n"
"           hdf5:      HDF5 library\n"
"           hdf5_log:  HDF5 library with Log-based VOL\n"
"           hdf5_md:   HDF5 library with multi-dataset APIs\n"
"           adios:     ADIOS2 library using BP3 format\n"
"       [-x strategy] I/O strategy to write\n"
"           canonical: Store E3SM variables in the canonical layout (default)\n"
"           log:       Store E3SM variables as is in log-based storage layout\n"
"           blob:      Write data is stored in a contiguous block (blob),\n"
"                      ignoring variable's canonical order\n"
"       [-c size] Data chunk size to be used when compression is enabled.\n"
"                 (default 0, i.e. no chunking)\n"
"       [-z filter] Enable data compression in write and use the supplied the\n"
"                 filter name (default: none)\n"
"       FILE: Name of input NetCDF file describing data decompositions\n\n";
    fprintf (stderr, help, argv0);
}

/*----< main() >-------------------------------------------------------------*/
int main (int argc, char **argv) {
    int i, err, nrecs=1, ffreq;
    double timing[3], max_t[3];
    e3sm_io_config cfg;
    e3sm_io_decom decom;

    MPI_Init (&argc, &argv);

    timing[0] = MPI_Wtime();

    MPI_Comm_rank (MPI_COMM_WORLD, &(cfg.rank));
    MPI_Comm_size (MPI_COMM_WORLD, &(cfg.np));

    cfg.io_comm        = MPI_COMM_WORLD;
    cfg.info           = MPI_INFO_NULL;
    cfg.num_iotasks    = cfg.np;
    cfg.num_group      = 1;
    cfg.out_path[0]    = '\0';
    cfg.in_path[0]     = '\0';
    cfg.cfg_path[0]    = '\0';
    cfg.hx             = -1;
    cfg.wr             = 0;
    cfg.rd             = 0;
    cfg.nvars          = 0;
    cfg.strategy       = undef_io;
    cfg.api            = undef_api;
    cfg.chunksize      = 0;
    cfg.filter         = none;
    cfg.verbose        = 0;
    cfg.keep_outfile   = 0;
    cfg.profiling      = 0;
    cfg.two_buf        = 0;
    cfg.non_contig_buf = 0;
    cfg.io_stride      = 1;
    cfg.sub_comm       = MPI_COMM_NULL;

    for (i = 0; i < MAX_NUM_DECOMP; i++) {
        cfg.G_case.nvars_D[i]    = 0;
        cfg.F_case_h0.nvars_D[i] = 0;
        cfg.F_case_h1.nvars_D[i] = 0;
        cfg.I_case_h0.nvars_D[i] = 0;
        cfg.I_case_h1.nvars_D[i] = 0;

        decom.blocklens[i]   = NULL;
        decom.disps[i]       = NULL;
        decom.raw_offsets[i] = NULL;
        decom.w_starts[i]    = NULL;
    }
    ffreq = 1;

    /* command-line arguments */
    while ((i = getopt (argc, argv, "vkr:s:o:i:dmf:ha:x:g:y:p")) != EOF)
        switch (i) {
            case 'v':
                cfg.verbose = 1;
                break;
            case 'k':
                cfg.keep_outfile = 1;
                break;
            case 'r':
                nrecs = atoi (optarg);
                break;
            case 'y':
                ffreq = atoi (optarg);
                break;
            case 's':
                cfg.io_stride = atoi (optarg);
                break;
            case 'a':
                if (strcmp (optarg, "pnetcdf") == 0)
                    cfg.api = pnetcdf;
                else if (strcmp (optarg, "hdf5") == 0)
                    cfg.api = hdf5;
                else if (strcmp (optarg, "hdf5_md") == 0)
                    cfg.api = hdf5_md;
                else if (strcmp (optarg, "hdf5_log") == 0)
                    cfg.api = hdf5_log;
                else if (strcmp (optarg, "adios") == 0)
                    cfg.api = adios;
#ifdef E3SM_IO_DEBUG
                /* For debug purpose only */
                else if (strcmp (optarg, "hdf5_ra") == 0)
                    cfg.api = hdf5_ra;
#endif
                else
                    ERR_OUT ("Unknown API")
                break;
                /*
            case 'l':
                if (strcmp (optarg, "contig") == 0)
                    cfg.layout = contig;
                else if (strcmp (optarg, "chunk") == 0)
                    cfg.layout = chunk;
                else
                    ERR_OUT ("Unknown layout")
                break;
                */
            case 'x':
                if (strcmp (optarg, "canonical") == 0)
                    cfg.strategy = canonical;
                else if (strcmp (optarg, "log") == 0)
                    cfg.strategy = log;
                else if (strcmp (optarg, "blob") == 0)
                    cfg.strategy = blob;
                else
                    ERR_OUT ("Unknown I/O strategy")
                break;

            case 'o':
                strncpy (cfg.out_path, optarg, E3SM_IO_MAX_PATH);
                cfg.wr = 1;
                break;
            case 'i':
                strncpy (cfg.in_path, optarg, E3SM_IO_MAX_PATH);
                cfg.rd = 1;
                break;
            case 'm':
                cfg.non_contig_buf = 1;
                break;
            case 'f':
                cfg.hx = atoi (optarg);
                break;
            case 'g':
                cfg.num_group = atoi (optarg);
                break;
            case 'c':
                cfg.chunksize = atoll (optarg);
                break;
            case 'p':
                cfg.profiling = 1;
                break;
            case 'z':
                if (strcmp (optarg, "deflate") == 0)
                    cfg.filter = deflate;
                else if (strcmp (optarg, "zlib") == 0)
                    cfg.filter = deflate;
                else
                    ERR_OUT ("Unknown filter")
                break;
            case 'h':
            default:
                if (cfg.rank == 0) usage (argv[0]);
                goto err_out;
        }

    if (optind >= argc || argv[optind] == NULL) { /* input file is mandatory */
        if (!cfg.rank) usage (argv[0]);
        ERR_OUT ("Decomposition file not provided")
    }
    strncpy (cfg.cfg_path, argv[optind], E3SM_IO_MAX_PATH);

    cfg.F_case_h0.nrecs = 1;  /* force only one record for F h0 case */
    cfg.F_case_h1.nrecs = nrecs;
    cfg.G_case.nrecs    = nrecs;
    cfg.I_case_h0.nrecs = nrecs;
    cfg.I_case_h1.nrecs = 1;  /* force only one record for I h1 case */

    cfg.F_case_h0.ffreq = ffreq;
    cfg.F_case_h1.ffreq = ffreq;
    cfg.G_case.ffreq    = ffreq;
    cfg.I_case_h0.ffreq = ffreq;
    cfg.I_case_h1.ffreq = ffreq;

    /* neither command-line option -i or -o is used */
    if (!cfg.wr && !cfg.rd)
        ERR_OUT("Error: neither command-line option -i nor -o is used")
    if (cfg.out_path[0] == '-')
        ERR_OUT ("Empty output file path")
    if (cfg.in_path[0] == '-')
        ERR_OUT ("Empty input file path")

    /* check yet to support APIs and I/O strategies */
    if (cfg.api == undef_api) cfg.api = pnetcdf;

    if (cfg.api == pnetcdf) {
        if (cfg.strategy == undef_io)
            cfg.strategy = canonical;
        else if (cfg.strategy == log)
            ERR_OUT ("PnetCDF with log I/O strategy is not supported yet")
    }

    if (cfg.api == hdf5) {
#ifndef ENABLE_HDF5
        ERR_OUT("HDF5 is not enabled at configure time")
#endif
        if (cfg.strategy == undef_io)
            cfg.strategy = canonical;
        else if (cfg.strategy == log)
            ERR_OUT ("HDF5 rearranger with log I/O strategy is not supported yet")
    }

    if (cfg.api == hdf5_md) {
#ifndef HDF5_HAVE_DWRITE_MULTI
        ERR_OUT("HDF5 does not support multi-dataset APIs")
#endif
        if (cfg.strategy == undef_io)
            cfg.strategy = canonical;
        else if (cfg.strategy == log)
            ERR_OUT ("HDF5 multi-dataset with log I/O strategy is not supported yet")
        else if (cfg.strategy == blob)
            ERR_OUT ("HDF5 multi-dataset with blob I/O strategy is not supported yet")
    }

#ifdef E3SM_IO_DEBUG
    /* For debug purpose only */
    if (cfg.api == hdf5_ra) {
        if (cfg.strategy == undef_io)
            cfg.strategy = canonical;
        else if (cfg.strategy == log)
            ERR_OUT ("HDF5 multi-dataset with log I/O strategy is not supported yet")
        else if (cfg.strategy == blob)
            ERR_OUT ("HDF5 multi-dataset with blob I/O strategy is not supported yet")
    }
#endif

    if (cfg.api == hdf5_log) {
#ifndef ENABLE_LOGVOL
        ERR_OUT("HDF5 Log VOL is not enabled at configure time")
#endif
        if (cfg.strategy == undef_io)
            cfg.strategy = log;
        else if (cfg.strategy == canonical)
            ERR_OUT ("HDF5 log-based VOL with canonical I/O strategy is not supported yet")
        else if (cfg.strategy == blob)
            ERR_OUT ("HDF5 log-based VOL with blob I/O strategy is not supported yet")
    }

    if (cfg.api == adios) {
#ifndef ENABLE_ADIOS2
        ERR_OUT("ADIOS is not enabled at configure time")
#endif
        if (cfg.strategy == undef_io)
            cfg.strategy = blob;
        else if (cfg.strategy == canonical)
            ERR_OUT ("ADIOS with canonical I/O strategy is not supported yet")
        else if (cfg.strategy == log)
            ERR_OUT ("ADIOS with log I/O strategy is not supported yet")
    }

    /* input decomposition file contains number of write requests and their
     * file access offsets (per array element) */
    PRINT_MSG (1, "Input decomposition file name = %s\n", cfg.cfg_path);

    /* print input and output file/folder names */
    PRINT_MSG (1, "Input  data file/folder name = %s\n", cfg.in_path);
    PRINT_MSG (1, "Output data file/folder name = %s\n", cfg.out_path);

    MPI_Barrier(MPI_COMM_WORLD);
    timing[1] = MPI_Wtime();

    /* read request information from decomposition file */
    err = read_decomp(&cfg, &decom);
    CHECK_ERR

    /* determine run case */
    if (decom.num_decomp == 3)
        cfg.run_case = F;
    else if (decom.num_decomp == 6)
        cfg.run_case = G;
    else if (decom.num_decomp == 5)
        cfg.run_case = I;
    else
        cfg.run_case = unknown;

    timing[1] = MPI_Wtime() - timing[1];
    MPI_Barrier(MPI_COMM_WORLD);
    timing[2] = MPI_Wtime();

    /* set MPI-IO and PnetCDF hints */
    err = MPI_Info_create (&(cfg.info));
    CHECK_MPIERR
    err += set_info (&cfg, &decom);
    if (err < 0) goto err_out;

    /* the core of this benchmark */
    err = e3sm_io_core (&cfg, &decom);
    CHECK_ERR

    timing[2] = MPI_Wtime() - timing[2];

    /* report timing breakdowns */
    report_timing_WR(&cfg, &decom);

    if (cfg.profiling) e3sm_io_print_profile(&cfg);

    timing[0] = MPI_Wtime() - timing[0];
    MPI_Reduce(timing, max_t, 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (cfg.rank == 0) {
        printf("read_decomp=%.2f e3sm_io_core=%.2f MPI init-to-finalize=%.2f\n",
               max_t[1],max_t[2],max_t[0]);
        printf("-----------------------------------------------------------\n");
        printf("\n\n");
    }

err_out:
    if (cfg.info != MPI_INFO_NULL)
        MPI_Info_free (&(cfg.info));
    if (cfg.io_comm != MPI_COMM_WORLD && cfg.io_comm != MPI_COMM_NULL)
        MPI_Comm_free (&(cfg.io_comm));

    /* Free decom */
    for (i = 0; i < MAX_NUM_DECOMP; i++) {
        if (decom.blocklens[i]) free (decom.blocklens[i]);
        if (decom.disps[i]) free (decom.disps[i]);
        if (decom.raw_offsets[i]) free (decom.raw_offsets[i]);
        if (decom.w_starts[i] != NULL) {
            free(decom.w_starts[i][0]);
            free(decom.w_starts[i]);
        }
    }

    MPI_Finalize ();

    return (err < 0) ? 1 : 0;
}
