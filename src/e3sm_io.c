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
#include <string.h> /* strcpy(), strncpy(), strstr() */
#include <unistd.h> /* getopt() */

#include <mpi.h>

#include <e3sm_io.h>
#include <e3sm_io_err.h>

#ifdef E3SM_IO_PROFILING
#include <e3sm_io_profile.hpp>
#endif

static
void check_connector_env(e3sm_io_config *cfg) {
    char *env_str;

    cfg->env_log          = 0;
    cfg->env_log_passthru = 0;
    cfg->env_log_info     = NULL;
    cfg->env_cache        = 0;
    cfg->env_async        = 0;

    env_str = getenv("H5VL_LOG_PASSTHRU");
    if (env_str != NULL && env_str[0] == '1')
        cfg->env_log_passthru = 1;

    env_str = getenv("HDF5_VOL_CONNECTOR");
    if (env_str == NULL || env_str[0] == '\0')
        /* env HDF5_VOL_CONNECTOR is not set or set with ni value */
        return;

    if (strstr(env_str, "under_vol=512") != NULL || strstr(env_str, "async ") != NULL)
        cfg->env_async = 1;
    if (strstr(env_str, "under_vol=513") != NULL || strstr(env_str, "cache_ext ") != NULL)
        cfg->env_cache = 1;

    env_str = strdup(env_str);
    char *connector = strtok(env_str, "  \t\n\r");
    if (connector == NULL) {
        free(env_str);
        return;
    }

    if (strcmp(connector, "LOG") == 0) {
        /* if LOG is set in HDF5_VOL_CONNECTOR */
        cfg->env_log = 1;
    }
    else if (strcmp(connector, "cache_ext") == 0) {
        /* if cache_ext is the 1st in HDF5_VOL_CONNECTOR */
        char *info_str = strtok(NULL, "  \t\n\r");
        cfg->env_log_info = (char*) malloc(28 + strlen(info_str));
        sprintf(cfg->env_log_info, "under_vol=513;under_info={%s}", info_str);
        cfg->env_cache = 1;
    }
    else if (strcmp(connector, "async") == 0) {
        /* if async is the 1st in HDF5_VOL_CONNECTOR */
        char *info_str = strtok(NULL, "  \t\n\r");
        cfg->env_log_info = (char*) malloc(28 + strlen(info_str));
        sprintf(cfg->env_log_info, "under_vol=512;under_info={%s}", info_str);
        cfg->env_async = 1;
    }
    free(env_str);
}

static
void parse_hint_file(e3sm_io_config  *cfg,
                     const char      *filename,
                     int             *num_hint_lines,
                     char           **hint_lines)
{
    /* read I/O hints from file e3sm_io_hints.txt, run e3sm_io_core() one per
     * line in the hint file. Each line contains hints as if it is set in the
     * environment variable PNETCDF_HINTS.
     */
    MPI_Offset fsize = 0;
    char *hint_buf=NULL;

    *num_hint_lines=0;

    if (cfg->rank == 0) {
        FILE *fptr = fopen(filename, "r");
        if (fptr != NULL) {
            fseek(fptr, 0, SEEK_END);
            fsize = ftell(fptr) + 1;
            fseek(fptr, 0, SEEK_SET);
            hint_buf = (char*) malloc(fsize);
            fread(hint_buf, 1, fsize, fptr);
            hint_buf[fsize-1] = '\0';
            fclose(fptr);
        }
    }

    MPI_Bcast(&fsize, 1, MPI_OFFSET, 0, cfg->io_comm);
    if (fsize > 0) {
        if (cfg->rank > 0) hint_buf = (char*) malloc(fsize);
        MPI_Bcast(hint_buf, fsize, MPI_BYTE, 0, cfg->io_comm);

        char *hint_str = strtok(hint_buf, "\n");
        if (hint_str != NULL && strlen(hint_str) > 0)
            hint_lines[(*num_hint_lines)++] = strdup(hint_str);

        while ((hint_str = strtok(NULL, "\n")) != NULL)
            hint_lines[(*num_hint_lines)++] = strdup(hint_str);

        free(hint_buf);
    }
}

static
int parse_hint_line(e3sm_io_config *cfg,
                    const char     *hint_str)
{
    char *warn_str="Warning: skip ill-formed hint set in PNETCDF_HINTS";
    char *hint_saved, *ptr, *key, *val, *deli;
    int err=MPI_SUCCESS;

    if (hint_str == NULL) return 1;

    /* skip blank lines */
    hint_saved = strdup(hint_str);
    if (strlen(hint_saved) == 0 || strtok(hint_saved, " \t") == NULL) {
        free(hint_saved);
        return 1;
    }
    ptr = hint_saved;

    do {
        if (*ptr == '\0') break; /* done with this line */

        key = ptr;
        deli = strchr(ptr, ';');
        if (deli != NULL) {
            *deli = '\0'; /* add terminate char */
            ptr = deli + 1;
        }
        else
            ptr += strlen(ptr); /* last hint */

        /* hint key */
        deli = strchr(key, '=');
        if (deli == NULL) {
            /* expect one token before = */
            printf("xxxx %s: '%s'\n", warn_str, key);
            break;
        }
        *deli = '\0'; /* add terminate char */

        /* hint value */
        val = deli + 1;

        /* override previouse set hint or add a new one */
        err = MPI_Info_set(cfg->info, key, val);
        CHECK_MPIERR

    } while (*ptr != '\0');

err_out:
    free(hint_saved);
    return err;
}

static
int set_info(e3sm_io_config *cfg,
             const char     *hint_str)
{
    int err;

    /* set MPI-IO hints */

    err = MPI_Info_create (&(cfg->info));
    CHECK_MPIERR

    /* collective write */
    err = MPI_Info_set (cfg->info, "romio_cb_write", "enable");
    CHECK_MPIERR

    /* HDF5 may do independent I/O internally */

    if (cfg->api == pnetcdf) {
        /* set MPI-IO hints here */
        // err = MPI_Info_set (cfg->info, "romio_no_indep_rw", "true");
        // CHECK_MPIERR

        /* set PnetCDF I/O hints */

        /* if all write buffers are in a contiguous space, then disable PnetCDF
         * internal buffering */
        if (cfg->xtype == NC_DOUBLE  && /* no type conversion is necessary */
            cfg->non_contig_buf == 0 && /* all writes are in a single buffer */
            cfg->isReqSorted) {         /* write request offsets are sorted */

            /* actually setting nc_ibuf_size to 0 is not necessary, as copying
             * write requests to an internal buffer in PnetCDF is triggered
             * only when the user buffer is non-contiguous.
             */
            err = MPI_Info_set(cfg->info, "nc_ibuf_size", "0");
            CHECK_MPIERR

            /* Actually setting nc_in_place_swap to enable is not necessary,
             * because nc_in_place_swap is automatically enabled when the write
             * request is larger than 4KB, which is the case in this E3SM-IO.
             */
            err = MPI_Info_set(cfg->info, "nc_in_place_swap", "enable");
            CHECK_MPIERR
        }

        /* no gap between variables */
        err = MPI_Info_set (cfg->info, "nc_var_align_size", "1");
        CHECK_MPIERR
        /* in-place byte swap */
        err = MPI_Info_set (cfg->info, "nc_in_place_swap", "enable");
        CHECK_MPIERR

        err = parse_hint_line(cfg, hint_str);
        if (err == 0) goto err_out;
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
       [-j] Set the external data type to NC_FLOAT. This option only affects\n\
            the F and I cases. (default: NC_DOUBLE)\n\
       [-m] Run test using noncontiguous write buffer (default: contiguous)\n\
       [-q] Do not sort write requests based on their file offsets into an\n\
            increasing order (default: yes)\n\
       [-u] Fill missing elements in decomposition maps (default: no)\n\
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
           hdf5_log:  HDF5 library with Log VOL connector\n\
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
    int i, j, err, nrecs=1, ffreq;
    double timing[5], max_t[5];
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

    timing[1] = MPI_Wtime();

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
    cfg.fill_mode      = 0;
    cfg.env_log_info   = NULL;
    cfg.xtype          = NC_DOUBLE;
    cfg.sort_reqs      = 1;
    cfg.isReqSorted    = 0;

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
        decom.max_nreqs[i]   = 0;
    }
    ffreq = 1;

    /* command-line arguments */
    while ((i = getopt (argc, argv, "vkur:s:o:i:jmqf:ha:x:g:y:pt:")) != EOF)
        switch (i) {
            case 'v':
                cfg.verbose = 1;
                break;
            case 'k':
                cfg.keep_outfile = 1;
                break;
            case 'u':
                cfg.fill_mode = 1;
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
                    ERR_OUT("Unknown API")
                break;
                /*
            case 'l':
                if (strcmp (optarg, "contig") == 0)
                    cfg.layout = contig;
                else if (strcmp (optarg, "chunk") == 0)
                    cfg.layout = chunk;
                else
                    ERR_OUT("Unknown layout")
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
                    ERR_OUT("Unknown I/O strategy")
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
            case 'q':
                cfg.sort_reqs = 0;
                break;
            case 'j':
                cfg.xtype = NC_FLOAT;
                break;
            case 'f':
                cfg.hx = atoi (optarg);
                if (cfg.hx < -1 || cfg.hx > 1) {
                    if (cfg.rank == 0) {
                        printf("Error: invalid value for option -f\n");
                        printf("       valid values are: -1, 0, 1\n");
                    }
                    goto err_out;
                }
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
                    ERR_OUT("Unknown filter")
                break;
            case 'h':
            default:
                if (cfg.rank == 0) usage (argv[0]);
                goto err_out;
        }

    if (optind >= argc || argv[optind] == NULL) { /* input file is mandatory */
        if (!cfg.rank) usage (argv[0]);
        ERR_OUT("Decomposition file not provided")
    }
    strncpy (cfg.decomp_path, argv[optind], E3SM_IO_MAX_PATH);

    cfg.F_case_h0.nrecs = 1;  /* force only one record for F h0 case */
    cfg.F_case_h1.nrecs = nrecs;
    cfg.G_case.nrecs    = nrecs;
    cfg.I_case_h0.nrecs = nrecs;
    cfg.I_case_h1.nrecs = 1;  /* force only one record for I h1 case */

#ifdef HDF5_HAVE_MULTI_DATASET_API
    /* In HDF5 1.13.3, multi-dataset APIs have the following limitations.
     * 1. Do not support same dataset appears twice or more in a single call.
     * 2. When some datasets requiring type conversion and some don't, the
     *    collective I/O mode will roll back to independent mode internally.
     *    Note F and I cases write some variables that need type conversion
     *    from double to float, while G case does not.
     * See https://github.com/HDFGroup/hdf5/issues/1859
     */
    if (cfg.api == hdf5_md && ffreq > 1) {
        if (cfg.rank == 0)
            printf("Warning: HDF5 multi-dataset APIs do not support writing multiple time steps at a time. Re-set flush freq to 1.\n");
        ffreq = 1;
    }
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
        ERR_OUT("Empty output file path")
    if (cfg.in_path[0] == '-')
        ERR_OUT("Empty input file path")

    /* check yet to support APIs and I/O strategies */
    if (cfg.strategy == undef_io) {
        if (cfg.api == hdf5_log) {
            cfg.strategy = log;
        }
        else {
            cfg.strategy = canonical;
        }
    }
    check_connector_env(&cfg);
    switch (cfg.api) {
        case undef_api:
            cfg.api = pnetcdf;
        case pnetcdf:
#ifdef ENABLE_PNC
            if (cfg.strategy == log)
                ERR_OUT("Option PnetCDF does not support log strategy")
            break;
#else
            ERR_OUT("Using PnetCDF was not enabled in this E3SM I/O build")
#endif
        case netcdf4:
#ifdef ENABLE_NETCDF4
            switch (cfg.strategy) {
                case blob:
                    ERR_OUT("Option NetCDF-4 does not support blob strategy")
                case canonical:
                    if (cfg.env_log) {
                        printf("Error: please unset HDF5_VOL_CONNECTOR to run \"-a netcdf4 -x canonical\"\n");
                        err = -1;
                        goto err_out;
                        // ERR_OUT("The VOL set in HDF5_VOL_CONNECTOR (%s) is not compatible with NetCDF 4 canonical I/O strategy")
                    }
                    break;
                case log:
                    if (cfg.env_log == 0)
                        ERR_OUT("HDF5_VOL_CONNECTOR must be set to \"LOG under_vol=0;under_info={}\" (log-based VOL) for NetCDF 4 with log I/O strategy")
                    break;
                default:
                    ERR_OUT("NetCDF 4 only supports canonical and log strategies")
            }
#else
            ERR_OUT("Using NetCDF-4 was not enabled in this E3SM I/O build")
#endif
            break;
        case hdf5_md:
#ifdef HDF5_HAVE_MULTI_DATASET_API
            if (cfg.strategy != canonical)
                ERR_OUT("HDF5 multi-dataset option only supports canonical strategy")
            if (cfg.env_log)
                unsetenv("HDF5_VOL_CONNECTOR");
            if (cfg.rank == 0 && cfg.verbose)
                printf("VERBOSE: I/O API: HDF5 multi-dataset, strategy: canonical\n");
            break;
#else
            ERR_OUT("Using HDF5 multi-dataset APIs was not enabled in this E3SM I/O build")
#endif
        case hdf5:
#ifdef ENABLE_HDF5
            switch (cfg.strategy) {
                case blob:
                    /* output file layout will be blob, but if env
                     * HDF5_VOL_CONNECTOR is set to Log VOL, then the file is
                     * also using log layout
                     */
                    if (cfg.rank == 0 && cfg.verbose)
                        printf("VERBOSE: I/O API: HDF5, strategy: blob\n");
                    break;
                case canonical:
                    /* output file layout depends on env HDF5_VOL_CONNECTOR */
                    if (cfg.rank == 0 && cfg.verbose) {
                        if (cfg.env_log)
                            /* output file will be in log layout */
                            printf("VERBOSE: I/O API: HDF5, strategy: canonical, layout: log\n");
                        else
                            /* output file will be in canonical layout */
                            printf("VERBOSE: I/O API: HDF5, strategy: canonical, layout: canonical\n");
                    }
                    break;
                case log:
#ifndef ENABLE_LOGVOL
                    ERR_OUT("Option -a hdf5 -x log required Log VOL feature enabled at the configure time")
#endif
                    /* output file layout will be log */
                    if (cfg.rank == 0 && cfg.verbose)
                        printf("VERBOSE: I/O API: HDF5, strategy: log\n");
                    /* if HDF5_VOL_CONNECTOR is not set to use Log VOL, then
                     * H5Pset_vol() will be called. Otherwise, no H5Pset_vol()
                     * is called. This requires Log VOL feature enabled at the
                     * configure time.
                     */
                    break;
                default:
                    ERR_OUT("No such I/O strategy")
                    break;
            }
#else
            ERR_OUT("Using HDF5 was not enabled in this E3SM I/O build")
#endif
            break;
        case hdf5_log:
#ifdef ENABLE_LOGVOL
            if (cfg.strategy != log)
                ERR_OUT("Option hdf5_log only support log strategy")
            if (cfg.rank == 0 && cfg.verbose)
                printf("VERBOSE: I/O API: HDF5 Log VOL connector, strategy: log\n");
#else
            ERR_OUT("Using Log VOL connector was not enabled in this E3SM I/O build")
#endif
            break;
        case adios:;
#ifdef ENABLE_ADIOS2
            if (cfg.strategy != blob)
                ERR_OUT("Option ADIOS only supports blob strategy")
#else
            ERR_OUT("Using ADIOS was not enabled in this E3SM I/O build")
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

    char *hint_lines[64];
    int num_hint_lines;
    parse_hint_file(&cfg, "e3sm_io_hints.txt", &num_hint_lines, hint_lines);

    timing[1] = MPI_Wtime() - timing[1];
    MPI_Barrier(MPI_COMM_WORLD);
    timing[2] = MPI_Wtime();

    /* read request information from decomposition file */
    err = read_decomp(&cfg, &decom);
    CHECK_ERR

    /* determine run case */
    if (decom.num_decomp == 3)
        cfg.run_case = F;
    else if (decom.num_decomp == 6) {
        cfg.run_case = G;
        /* In the G case, the external data type is NC_DOUBLE */
        cfg.xtype = NC_DOUBLE;
    }
    else if (decom.num_decomp == 5)
        cfg.run_case = I;
    else
        cfg.run_case = unknown;

    timing[2] = MPI_Wtime() - timing[2];

    j = 0;
    do {

        char *hint_str = (num_hint_lines == 0) ? NULL : hint_lines[j];

if (cfg.rank == 0) printf("\nHINTS: %s\n\n", (hint_str) ? hint_str : "");

        /* set MPI-IO and PnetCDF hints */
        err = set_info(&cfg, hint_str);
        if (err < 0) goto err_out;

        /* the core of this benchmark */
        MPI_Barrier(MPI_COMM_WORLD);
        timing[3] = MPI_Wtime();

        err = e3sm_io_core(&cfg, &decom);
        CHECK_ERR

        timing[3] = MPI_Wtime() - timing[3];

        /* report timing breakdowns */
        if (cfg.rd)
            report_timing_RD(&cfg, &decom);
        else
            report_timing_WR(&cfg, &decom);

#ifdef E3SM_IO_PROFILING
        if (cfg.profiling) e3sm_io_print_profile(&cfg);
#else
        if (cfg.profiling && cfg.rank == 0)
            printf("\nWarning: E3SM-IO internal time profiling was disabled at configure time\n\n");
#endif

        for (i = 0; i < MAX_NUM_DECOMP; i++) {
            if (decom.w_starts[i] != NULL) {
                free(decom.w_starts[i][0]);
                free(decom.w_starts[i]);
            }
        }

        if (cfg.info != MPI_INFO_NULL) MPI_Info_free (&(cfg.info));

        timing[0] = timing[1] + timing[2] + timing[3];
        MPI_Reduce(timing, max_t, 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (cfg.rank == 0) {
            printf("init=%.2f read_decomp=%.2f e3sm_io_core=%.2f end-to-end=%.2f\n",
                   max_t[1],max_t[2],max_t[3],max_t[0]);
            printf("-----------------------------------------------------------\n");
            printf("\n\n");
        }

        j++;
    } while (j < num_hint_lines);

    for (j=0; j<num_hint_lines; j++)
        free(hint_lines[j]);

err_out:
    if (cfg.io_comm != MPI_COMM_WORLD && cfg.io_comm != MPI_COMM_NULL)
        MPI_Comm_free (&(cfg.io_comm));
    if (cfg.env_log_info != NULL)
        free(cfg.env_log_info);

    /* Free decom */
    for (i = 0; i < MAX_NUM_DECOMP; i++) {
        if (decom.blocklens[i]) free (decom.blocklens[i]);
        if (decom.disps[i]) free (decom.disps[i]);
        if (decom.raw_offsets[i]) free (decom.raw_offsets[i]);
    }

    MPI_Finalize ();

    return (err < 0) ? 1 : 0;
}
