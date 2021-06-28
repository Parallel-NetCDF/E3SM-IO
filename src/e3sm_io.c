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
    int err, nerrs = 0;
    MPI_Offset estimated_nc_ibuf_size;

    /* set MPI-IO hints */
    MPI_Info_set (cfg->info, "romio_cb_write", "enable");  /* collective write */
    MPI_Info_set (cfg->info, "romio_no_indep_rw", "true"); /* no independent MPI-IO */

    /* set PnetCDF I/O hints */
    MPI_Info_set (cfg->info, "nc_var_align_size", "1");     /* no gap between variables */
    MPI_Info_set (cfg->info, "nc_in_place_swap", "enable"); /* in-place byte swap */

    /* use total write amount to estimate nc_ibuf_size */
    estimated_nc_ibuf_size =
        decom->dims[2][0] * decom->dims[2][1] * sizeof (double) / cfg->num_iotasks;
    estimated_nc_ibuf_size *= cfg->nvars;
    if (estimated_nc_ibuf_size > 16777216) {
        char nc_ibuf_size_str[32];
        sprintf (nc_ibuf_size_str, "%lld", estimated_nc_ibuf_size);
        MPI_Info_set (cfg->info, "nc_ibuf_size", nc_ibuf_size_str);
    }

err_out:;
    return nerrs;
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
        "       [-o output_dir] Output directory name (default ./)\n"
        "       [-i target_dir] Path to directory containing the input files\n"
        "       [-a api] Underlying API to test (pnc (default), hdf5, hdf5_logvol, hdf5_multi, "
        "adios2, adios2_bp3)\n"
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
    int err, nerrs = 0;
    int i;
    char targetdir[E3SM_IO_MAX_PATH] = "./";
    char datadir[E3SM_IO_MAX_PATH]   = "";
    char cfgpath[E3SM_IO_MAX_PATH]   = "";
    e3sm_io_config cfg;
    e3sm_io_decom decom;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &(cfg.rank));
    MPI_Comm_size (MPI_COMM_WORLD, &(cfg.np));
    cfg.io_comm        = MPI_COMM_WORLD;
    cfg.info           = MPI_INFO_NULL;
    cfg.num_iotasks    = cfg.np;
    cfg.targetdir      = targetdir;
    cfg.datadir        = datadir;
    cfg.cfgpath        = cfgpath;
    cfg.hx             = -1;
    cfg.nrec           = 1;
    cfg.wr             = 0;
    cfg.rd             = 0;
    cfg.nvars          = 0;
    cfg.strate         = canonical;
    cfg.api            = pnc;
    cfg.chunksize      = 0;
    cfg.filter         = none;
    cfg.vard           = 0;
    cfg.verbose        = 0;
    cfg.keep_outfile   = 0;
    cfg.two_buf        = 0;
    cfg.non_contig_buf = 0;
    cfg.io_stride      = 1;
    cfg.filepernode    = 0;

    /* command-line arguments */
    while ((i = getopt (argc, argv, "vkr:s:o:i:dnmtRWf:ha:S:l")) != EOF) switch (i) {
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
                if (strcmp (optarg, "pnc") == 0) {
                    cfg.api = pnc;
                }
#ifdef ENABLE_HDF5
                else if (strcmp (optarg, "hdf5") == 0) {
                    cfg.api = hdf5_native;
                } else if (strcmp (optarg, "hdf5_multi") == 0) {
#ifdef HDF5_HAVE_DWRITE_MULTI
                    cfg.api = hdf5_multi;
#else
                    RET_ERR ("The HDF5 used does not support multi-dataset write")
#endif
                } else if (strcmp (optarg, "hdf5_logvol") == 0) {
#ifdef ENABLE_LOGVOL
                    cfg.api = hdf5_logvol;
#else
                    RET_ERR ("Log VOL support was not enabled in this build");
#endif
                }
#endif
#ifdef ENABLE_ADIOS2
                else if (strcmp (optarg, "adios2") == 0) {
                    cfg.api = adios2;
                } else if (strcmp (optarg, "adios2_bp3") == 0) {
                    cfg.api = adios2_bp3;
                }
#endif
                else {
                    RET_ERR ("Unknown API")
                }
                break;

            case 'l':
                cfg.filepernode = 1;
                break;
            case 'S':
                if (strcmp (optarg, "canonical") == 0) {
                    cfg.strate = canonical;
                } else if (strcmp (optarg, "log") == 0) {
                    cfg.strate = log;
                } else if (strcmp (optarg, "blob") == 0) {
                    cfg.strate = blob;
                } else {
                    RET_ERR ("Unknown I/O strategy")
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
            case 'c':
                cfg.chunksize = atoll (optarg);
                break;
            case 'z':
                if (strcmp (optarg, "deflate") == 0) {
                    cfg.filter = deflate;
                } else if (strcmp (optarg, "zlib") == 0) {
                    cfg.filter = deflate;
                } else {
                    RET_ERR ("Unknown filter")
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

    if ((cfg.strate == log) && (cfg.api != hdf5_logvol)) {
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

    /* read request information from decompositions 1, 2 and 3 */
    for (i = 0; i < MAX_NUM_DECOMP; i++) {
        decom.blocklens[i]   = NULL;
        decom.disps[i]       = NULL;
        decom.raw_offsets[i] = NULL;
    }
    if (cfg.strate == blob) {
        err = read_decomp (cfg.verbose, cfg.io_comm, cfg.cfgpath, &(decom.num_decomp), decom.dims,
                           decom.contig_nreqs, decom.ndims, decom.disps, decom.blocklens,
                           decom.raw_nreqs, decom.raw_offsets);
    } else {
        // No need to generate simulated raw offsets in canonical strategy
        err =
            read_decomp (cfg.verbose, cfg.io_comm, cfg.cfgpath, &(decom.num_decomp), decom.dims,
                         decom.contig_nreqs, decom.ndims, decom.disps, decom.blocklens, NULL, NULL);
    }
    CHECK_ERR

    // Create node-local I/O communicator
    if (cfg.filepernode) {
        int is_node_head;

        if (cfg.strate != blob) { RET_ERR ("Subfiling only available to blob strategy") }
        
        err = MPI_Comm_split_type (cfg.io_comm, MPI_COMM_TYPE_SHARED, cfg.rank, MPI_INFO_NULL,
                                   &(cfg.node_comm));
        CHECK_MPIERR

        /* Simulate 2 nodes split for debugging purpose 
            err = MPI_Comm_split (cfg.io_comm, cfg.rank & 1, cfg.rank, &(cfg.node_comm));
            CHECK_MPIERR
        */
       
        err = MPI_Comm_size (cfg.node_comm, &(cfg.node_np));
        CHECK_MPIERR
        err = MPI_Comm_rank (cfg.node_comm, &(cfg.node_rank));
        CHECK_MPIERR

        if (cfg.node_rank == 0) {
            is_node_head = 1;
        } else {
            is_node_head = 0;
        }
        cfg.node_id = 0;  // Exscan won't assign for rank 0, initialize to 0
        err         = MPI_Exscan (&is_node_head, &(cfg.node_id), 1, MPI_INT, MPI_SUM, cfg.io_comm);
        CHECK_MPIERR

        err = MPI_Bcast (&(cfg.node_id), 1, MPI_INT, 0, cfg.node_comm);
        CHECK_MPIERR

        err = MPI_Allreduce (&(cfg.node_id), &(cfg.num_node), 1, MPI_INT, MPI_MAX, cfg.node_comm);
        CHECK_MPIERR
        cfg.num_node++;
    } else {
        cfg.node_id   = 0;
        cfg.node_comm = NULL;
    }

    err = MPI_Info_create (&(cfg.info));
    CHECK_MPIERR
    nerrs += set_info (&cfg, &decom);
    CHECK_NERR

    nerrs += e3sm_io_core (&cfg, &decom);

err_out:;
    if (cfg.info != MPI_INFO_NULL) MPI_Info_free (&(cfg.info));
    if (cfg.io_comm != MPI_COMM_WORLD && cfg.io_comm != MPI_COMM_NULL) {
        MPI_Comm_free (&(cfg.io_comm));
    }
    if (cfg.node_comm != MPI_COMM_WORLD && cfg.node_comm != MPI_COMM_NULL) {
        MPI_Comm_free (&(cfg.node_comm));
    }

    // Free decom
    for (i = 0; i < MAX_NUM_DECOMP; i++) {
        if (decom.blocklens[i]) { free (decom.blocklens[i]); }
        if (decom.disps[i]) { free (decom.disps[i]); }
        if (decom.raw_offsets[i]) { free (decom.raw_offsets[i]); }
    }

    /* Non-IO tasks wait for IO tasks to complete */
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize ();

    return (nerrs > 0);
}