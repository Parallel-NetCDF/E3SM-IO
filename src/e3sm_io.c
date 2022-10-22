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

    /* HDF5 may do independent I/O internally */

    if (cfg->api == pnetcdf) {
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
        estimated_nc_ibuf_size = decom->dims[2][0] * decom->dims[2][1]
                               * sizeof (double) / cfg->num_iotasks;
        estimated_nc_ibuf_size *= cfg->nvars;
        if (estimated_nc_ibuf_size > 16777216) {
            char nc_ibuf_size_str[32];
            sprintf (nc_ibuf_size_str, "%lld", estimated_nc_ibuf_size);
            err = MPI_Info_set (cfg->info, "nc_ibuf_size", nc_ibuf_size_str);
            CHECK_MPIERR
        }
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
    char *help = "Usage: %s [OPTION] FILE\n\
       [-h] Print this help message\n\
       [-v] Verbose mode\n\
       [-k] Keep the output files when program exits (default: deleted)\n\
       [-m] Run test using noncontiguous write buffer (default: contiguous)\n\
       [-f num] Output history files h0 or h1: 0 for h0 only, 1 for h1 only,\n\
                -1 for both. Affect only F and I cases. (default: -1)\n\
       [-r num] Number of time records/steps written in F case h1 file and I\n\
                case h0 file (default: 1)\n\
       [-y num] Data flush frequency. (1: flush every time step, the default,\n\
                and -1: flush once for all time steps. (No effect on ADIOS\n\
                and HDF5 blob I/O options, which always flushes at file close).\n\
       [-s num] Stride interval of ranks for selecting MPI processes to perform\n\
                I/O tasks (default: 1, i.e. all MPI processes).\n\
       [-g num] Number of subfiles, used by Log-based VOL and ADIOS I/O only,\n\
                -1 for one subfile per compute node, 0 to disable subfiling,\n\
                (default: 0).\n\
       [-t time] Add sleep time to emulate the computation in order to \n\
                 overlapping I/O when Async VOL is used.\n\
       [-i path] Input file path (folder name when subfiling is used, file\n\
                 name otherwise).\n\
       [-o path] Output file path (folder name when subfiling is used, file\n\
                 name otherwise).\n\
       [-a api]  I/O library name\n\
           pnetcdf:   PnetCDF library (default)\n\
           netcdf4:   NetCDF-4 library\n\
           hdf5:      HDF5 library\n\
           hdf5_md:   HDF5 library using multi-dataset I/O APIs\n\
           hdf5_log:  HDF5 library with Log-based VOL\n\
           adios:     ADIOS library using BP3 format\n\
       [-x strategy] I/O strategy\n\
           canonical: Store variables in the canonical layout (default).\n\
           log:       Store variables in the log-based storage layout.\n\
           blob:      Pack and store all data written locally in a contiguous\n\
                      block (blob), ignoring variable's canonical order.\n\
       FILE: Name of input file storing data decomposition maps.\n";
    fprintf (stderr, help, argv0);
}

/* command-line options -c and -z are for experimental data compression feature
 * in PnetCDF.
"       [-c size] Data chunk size used to compress data and metadata. This\n"
"                 option affects only hdf5_log. (default 0, i.e. no chunking)\n"
"       [-z filter] Filter name to compress data and metadata (default: none)\n"
*/

/*----< main() >-------------------------------------------------------------*/
int main (int argc, char **argv) {
    int i, err, nrecs=1, ffreq;
    char *env;  // HDF5_VOL_CONNECTOR environment variable
    double timing[3], max_t[3];
    e3sm_io_config cfg;
    e3sm_io_decom decom;

#ifdef E3SM_IO_THREADING
    {
        int mpi_required = 0;
        err = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_required);
        CHECK_ERR;

        if (mpi_required != MPI_THREAD_MULTIPLE) {
            ERR_OUT("MPI_Init_thread failed");
        }
    }
#else
    err = MPI_Init (&argc, &argv);
    CHECK_ERR;
#endif

    timing[0] = MPI_Wtime();

    MPI_Comm_rank (MPI_COMM_WORLD, &(cfg.rank));
    MPI_Comm_size (MPI_COMM_WORLD, &(cfg.np));

    cfg.io_comm        = MPI_COMM_WORLD;
    cfg.info           = MPI_INFO_NULL;
    cfg.num_iotasks    = cfg.np;
    cfg.num_subfiles   = 0;
    cfg.out_path[0]    = '\0';
    cfg.in_path[0]     = '\0';
    cfg.decomp_path[0] = '\0';
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
    cfg.comp_time      = 0;

    for (i = 0; i < MAX_NUM_DECOMP; i++) {
        cfg.G_case.nvars_D[i]    = 0;
        cfg.F_case_h0.nvars_D[i] = 0;
        cfg.F_case_h1.nvars_D[i] = 0;
        cfg.I_case_h0.nvars_D[i] = 0;
        cfg.I_case_h1.nvars_D[i] = 0;

        cfg.G_case.num_attrs     = 0;
        cfg.F_case_h0.num_attrs  = 0;
        cfg.F_case_h1.num_attrs  = 0;
        cfg.I_case_h0.num_attrs  = 0;
        cfg.I_case_h1.num_attrs  = 0;

        decom.blocklens[i]   = NULL;
        decom.disps[i]       = NULL;
        decom.raw_offsets[i] = NULL;
        decom.w_starts[i]    = NULL;
    }
    ffreq = 1;

    /* command-line arguments */
    while ((i = getopt (argc, argv, "vkr:s:o:i:dmf:ha:x:g:y:pt:")) != EOF)
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
                else if (strcmp (optarg, "netcdf4") == 0)
                    cfg.api = netcdf4;
                else if (strcmp (optarg, "adios") == 0)
                    cfg.api = adios;
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
                cfg.num_subfiles = atoi (optarg);
                break;
            case 'c':
                cfg.chunksize = atoll (optarg);
                break;
            case 't':
                cfg.comp_time = atoi (optarg);
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
    strncpy (cfg.decomp_path, argv[optind], E3SM_IO_MAX_PATH);

    cfg.F_case_h0.nrecs = 1;  /* force only one record for F h0 case */
    cfg.F_case_h1.nrecs = nrecs;
    cfg.G_case.nrecs    = nrecs;
    cfg.I_case_h0.nrecs = nrecs;
    cfg.I_case_h1.nrecs = 1;  /* force only one record for I h1 case */

#ifdef HDF5_HAVE_MULTI_DATASET_API
    /* HDF5 multi-dataset APIs may have the following limitations.
     * 1. Support same dataset appears twice or more times in a single call.
     * 2. Support datasets some requiring type conversion and some don't.
     *    In this case, the collective I/O mode will roll back to independent
     *    mode internally.
     * See https://github.com/HDFGroup/hdf5/issues/1859
     */
#warning TODO: HDF5 multi-dataset APIs do not support writing multiple time steps at a time. Setting flush freq to 1.
    if (cfg.api == hdf5_md) ffreq = 1;
#endif

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
    if (cfg.strategy == undef_io) {
        if (cfg.api == hdf5_log) {
            cfg.strategy = log;
        }
        else {
            cfg.strategy = canonical;
        }
    }
    env = getenv ("HDF5_VOL_CONNECTOR");
    switch (cfg.api) {
        case undef_api:;
            cfg.api = pnetcdf;
        case pnetcdf:;
            if (cfg.strategy == log)
                ERR_OUT ("PnetCDF with log I/O strategy is not supported yet")
            break;
        case netcdf4:;
#ifdef ENABLE_NETCDF4
            switch (cfg.strategy) {
                case canonical:;
                    if (env && (strncmp (env, "LOG", 3) == 0)) {
                        printf("Error: please unset HDF5_VOL_CONNECTOR to run \"-a netcdf4 -x canonical\"\n");
                        err = -1;
                        goto err_out;
                        // ERR_OUT ("The VOL set in HDF5_VOL_CONNECTOR (%s) is not compatible with NetCDF 4 canonical I/O strategy")
                    }
                    break;
                case log:;
                    if (!env || (strncmp (env, "LOG", 3) != 0)) {
                        ERR_OUT ("HDF5_VOL_CONNECTOR must be set to \"LOG\" (log-based VOL) for NetCDF 4 with log I/O strategy")
                    }
                    break;
                default:;
                    ERR_OUT ("NetCDF 4 only supports canonical and log strategies")
            }
#else
            ERR_OUT ("NetCDF 4 was not enabled in this E3SM I/O build")
#endif
            break;
        case hdf5_md:;
#ifdef HDF5_HAVE_MULTI_DATASET_API
            if (cfg.strategy != canonical){
                ERR_OUT ("HDF5 multi-dataset only support canonical strategy")
            }
#else
            ERR_OUT("HDF5 used to build E3SM I/O does not support multi-dataset APIs")
#endif
            // fall through
        case hdf5:;
#ifdef ENABLE_HDF5
            switch (cfg.strategy) {
                case blob:;
                case canonical:;
                    if (env && (strncmp (env, "LOG", 3) == 0)) {
                        ERR_OUT ("The VOL set in HDF5_VOL_CONNECTOR (%s) is not compatible with HDF5 canonical I/O strategy")
                    }
                    break;
                case log:;
                    if (env && (strncmp (env, "LOG", 3) != 0)) {
                        ERR_OUT ("HDF5_VOL_CONNECTOR must be set to \"LOG\" (log-based VOL) for HDF5 with log I/O strategy")
                    }
                    break;
                default:;
                    ERR_OUT ("NetCDF 4 only supports canonical and log strategies")
            }
#else
            ERR_OUT ("HDF5 was not enabled in this E3SM I/O build")
#endif
            break;
        case hdf5_log:;
#ifdef ENABLE_LOGVOL
            if (cfg.strategy != log){
                ERR_OUT ("HDF5-logvol only support log strategy")
            }
            if (env && (strncmp (env, "LOG", 3) != 0)) {
                ERR_OUT ("HDF5_VOL_CONNECTOR must be set to \"LOG\" (log-based VOL) for HDF5-logvol with log I/O strategy")
            }
#else
            ERR_OUT ("Log-based VOL support was not enabled in this E3SM I/O build")
#endif
            break;
        case adios:;
#ifdef ENABLE_ADIOS2
            if (cfg.strategy != blob){
                ERR_OUT ("ADIOS only support blob strategy")
            }
#else
            ERR_OUT ("ADIOS was not enabled in this E3SM I/O build")
#endif
            break;
        default:
            break;
    }

    /* input decomposition file contains number of write requests and their
     * file access offsets (per array element) */
    PRINT_MSG (1, "Input decomposition file name = %s\n", cfg.decomp_path);

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
    err = set_info (&cfg, &decom);
    if (err < 0) goto err_out;

    /* the core of this benchmark */
    err = e3sm_io_core (&cfg, &decom);
    CHECK_ERR

    timing[2] = MPI_Wtime() - timing[2];

    /* report timing breakdowns */
    if (cfg.rd) {
        report_timing_RD(&cfg, &decom);
    }
    else{
        report_timing_WR(&cfg, &decom);
    }

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
