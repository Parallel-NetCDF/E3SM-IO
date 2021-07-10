/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program uses the E3SM I/O patterns recorded by the PIO library to
 * evaluate the performance of two PnetCDF APIs: ncmpi_vard_all(), and
 * ncmpi_iput_varn(). The E3SM I/O patterns consist of a large number of small,
 * noncontiguous requests on each MPI process, which presents a challenge for
 * achieving a good performance.
 *
 * See README.md for compile and run instructions.
 *
 *********************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
/**/
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
/**/
#include <unistd.h> /* getopt() */
/**/
#include <mpi.h>
/**/
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
}

/*----< usage() >------------------------------------------------------------*/
static void usage (char *argv0) {
    char *help =
        "Usage: %s [OPTION]... FILE\n"
        "       [-h] Print help\n"
        "       [-v] Verbose mode\n"
        "       [-k] Keep the output files when program exits\n"
        "       [-d] Run test that uses PnetCDF vard API\n"
        "       [-n] Run test that uses PnetCDF varn API\n"
        "       [-m] Run test using noncontiguous write buffer\n"
        "       [-t] Write 2D variables followed by 3D variables\n"
        "       [-R] Test reading performance\n"
        "       [-W] Test writing performance\n"
        "       [-f num] File number to run in F case (-1 (both) (default), 0, 1)\n"
        "       [-r num] Number of records (default 1)\n"
        "       [-s num] Stride between IO tasks (default 1)\n"
        "       [-g num] Number of IO groups (subfiles) (default 1)\n"
        "       [-o output_dir] Output directory name (default ./)\n"
        "       [-i target_dir] Path to directory containing the input files\n"
        "       [-a api] Underlying API to test (pnetcdf (default), hdf5_ra, hdf5_log, hdf5_mv, "
        "adios)\n"
        "              pnetcdf:     PnetCDF library\n"
        "              hdf5_ra:     HDF5 library's native VOL with rearranger in E3SM\n"
        "              hdf5_log:    HDF5 library with Log I/O VOL\n"
        "              hdf5_ra:     HDF5 library's experimental multi-dataset write function with rearranger in E3SM\n"
        "              adios:       ADIOS2 library using BP3 format\n"
        "       [-x strategy] I/O strategy used to write E3SM variables (canonical (default), log, blob)\n"
        "              canonical:   Store E3SM variables as is in canonical layout\n"
        "              log:         Store E3SM variables as is in log-based storage layout\n"
        "              blob:        Flatten E3SM variables into 1-dimensional data blocks. Record decomposition information in other variables and attributes\n"
        "       [-c chunk_size] Use chunked storage layout with chunk_size (0 (no chunking) "
        "(default))\n"
        "       [-z filter] Apply the filter if supported by the underlying API (none (default), "
        "deflate) \n"
        "       FILE: Name of input netCDF file describing data decompositions\n";
    "\n";
    fprintf (stderr, help, argv0);
}

/*----< main() >-------------------------------------------------------------*/
int main (int argc, char **argv) {
    int err;
    int i;
    char targetdir[E3SM_IO_MAX_PATH] = "./";
    char datadir[E3SM_IO_MAX_PATH]   = "";
    char cfgpath[E3SM_IO_MAX_PATH]   = "";
    double timing[2], max_t[2];
    e3sm_io_config cfg;
    e3sm_io_decom decom;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &(cfg.rank));
    MPI_Comm_size (MPI_COMM_WORLD, &(cfg.np));
    cfg.io_comm        = MPI_COMM_WORLD;
    cfg.info           = MPI_INFO_NULL;
    cfg.num_iotasks    = cfg.np;
    cfg.num_group      = 1;
    cfg.targetdir      = targetdir;
    cfg.datadir        = datadir;
    cfg.cfgpath        = cfgpath;
    cfg.hx             = -1;
    cfg.nrec           = 1;
    cfg.wr             = 0;
    cfg.rd             = 0;
    cfg.nvars          = 0;
    cfg.strategy       = canonical;
    cfg.api            = pnetcdf;
    cfg.chunksize      = 0;
    cfg.filter         = none;
    cfg.vard           = 0;
    cfg.verbose        = 0;
    cfg.keep_outfile   = 0;
    cfg.two_buf        = 0;
    cfg.non_contig_buf = 0;
    cfg.io_stride      = 1;
    cfg.sub_comm       = MPI_COMM_NULL;

    for (i = 0; i < MAX_NUM_DECOMP; i++) {
        decom.blocklens[i]   = NULL;
        decom.disps[i]       = NULL;
        decom.raw_offsets[i] = NULL;
    }

    /* command-line arguments */
    while ((i = getopt (argc, argv, "vkr:s:o:i:dnmtRWf:ha:x:g:")) != EOF) switch (i) {
            case 'v':
                cfg.verbose = 1;
                break;
            case 'k':
                cfg.keep_outfile = 1;
                break;
            case 'r':
                cfg.nrec = atoi (optarg);
                break;
            case 's':
                cfg.io_stride = atoi (optarg);
                break;
            case 'a':
                if (strcmp (optarg, "pnetcdf") == 0) {
                    cfg.api = pnetcdf;
                }
#ifdef ENABLE_HDF5
                else if (strcmp (optarg, "hdf5_ra") == 0) {
                    cfg.api = hdf5_ra;
                } else if (strcmp (optarg, "hdf5_mv") == 0) {
#ifdef HDF5_HAVE_DWRITE_MULTI
                    cfg.api = hdf5_mv;
#else
                    ERR_OUT ("The HDF5 used does not support multi-dataset write")
#endif
                } else if (strcmp (optarg, "hdf5_log") == 0) {
#ifdef ENABLE_LOGVOL
                    cfg.api = hdf5_log;
#else
                    ERR_OUT ("Log VOL support was not enabled in this build");
#endif
                }
#endif
#ifdef ENABLE_ADIOS2
                else if (strcmp (optarg, "adios") == 0) {
                    cfg.api = adios;
                }
#endif
                else {
                    ERR_OUT ("Unknown API")
                }
                break;
                /*
            case 'l':
                if (strcmp (optarg, "contig") == 0) {
                    cfg.layout = contig;
                } else if (strcmp (optarg, "chunk") == 0) {
                    cfg.layout = chunk;
                } else {
                    ERR_OUT ("Unknown layout")
                }
                break;
                */
            case 'x':
                if (strcmp (optarg, "canonical") == 0) {
                    cfg.strategy = canonical;
                } else if (strcmp (optarg, "log") == 0) {
                    cfg.strategy = log;
                } else if (strcmp (optarg, "blob") == 0) {
                    cfg.strategy = blob;
                } else {
                    ERR_OUT ("Unknown I/O strategy")
                }
                break;

            case 'o':
                strncpy (cfg.targetdir, optarg, E3SM_IO_MAX_PATH);
                break;
            case 'i':
                strncpy (cfg.datadir, optarg, E3SM_IO_MAX_PATH);
                break;
            case 'd':
                cfg.vard = 1;
                break;
            case 'n':
                cfg.vard = 0;
                break;
            case 'm':
                cfg.non_contig_buf = 1;
                break;
            case 't':
                cfg.two_buf = 1;
                break;
            case 'R':
                cfg.rd = 1;
                break;
            case 'W':
                cfg.wr = 1;
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
            case 'z':
                if (strcmp (optarg, "deflate") == 0) {
                    cfg.filter = deflate;
                } else if (strcmp (optarg, "zlib") == 0) {
                    cfg.filter = deflate;
                } else {
                    ERR_OUT ("Unknown filter")
                }
                break;
            case 'h':
            default:
                if (cfg.rank == 0) usage (argv[0]);
                goto err_out;
        }

    if ((optind >= argc) || (argv[optind] == NULL)) { /* input file is mandatory */
        if (!(cfg.rank)) usage (argv[0]);
        ERR_OUT ("Decomposition file not provided")
    }
    strncpy (cfg.cfgpath, argv[optind], E3SM_IO_MAX_PATH);

    if((cfg.strategy==log) && (cfg.api!=hdf5_log)){
        ERR_OUT ("Selected API does not support log-based I/O")
    }

    /* input file contains number of write requests and their file access
     * offsets (per array element) */
    PRINT_MSG (1, "input file name =%s\n", cfg.cfgpath);

    /* neither command-line option -R or -W is used, run write */
    if (!(cfg.wr || cfg.rd)) cfg.wr = 1;

    /* set the output folder name */
    PRINT_MSG (1, "Target folder name =%s\n", cfg.targetdir);
    if (cfg.datadir[0] != '\0') { PRINT_MSG (1, "Input folder name =%s\n", cfg.datadir); }

    MPI_Barrier(MPI_COMM_WORLD);
    timing[0] = MPI_Wtime();

    /* read request information from decomposition file */
    err = read_decomp(&cfg, &decom);
    CHECK_ERR

    timing[0] = MPI_Wtime() - timing[0];
    MPI_Barrier(MPI_COMM_WORLD);
    timing[1] = MPI_Wtime();

    err = MPI_Info_create (&(cfg.info));
    CHECK_MPIERR
    err += set_info (&cfg, &decom);
    if (err < 0) goto err_out;

    err = e3sm_io_core (&cfg, &decom);
    if (err < 0) goto err_out;

    timing[1] = MPI_Wtime() - timing[1];

    MPI_Reduce(timing, max_t, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (cfg.rank == 0)
        printf("read_decomp=%.2f e3sm_io_core=%.2f\n", max_t[0],max_t[1]);

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
    }

    MPI_Finalize ();

    return (err < 0) ? 1 : 0;
}
