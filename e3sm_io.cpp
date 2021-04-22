/*********************************************************************
 *
 * Copyright (C) 2018, Northwestern University
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

#include "e3sm_io.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */

#include <string>

#include "config.h"
#include "e3sm_io.hpp"
#include "e3sm_io_case.hpp"
#include "e3sm_io_driver.hpp"
#include "e3sm_io_driver_pnc.hpp"
#include "e3sm_io_err.hpp"
#ifdef ENABLE_HDF5
#include "e3sm_io_driver_hdf5.hpp"
#endif

int verbose;      /* verbose mode to print additional messages on screen */
int keep_outfile; /* whether to keep the output files when exits */
int two_buf;

e3sm_io_config::e3sm_io_config () {
    int err, nerrs = 0;

    err = MPI_Info_create (&(this->info));
    CHECK_MPIERR

err_out:;
    if (nerrs > 0) { throw "e3sm_io_config init fail"; }
}

e3sm_io_config::~e3sm_io_config () {
    MPI_Info_free (&(this->info));
    if (this->io_comm != MPI_COMM_WORLD && this->io_comm != MPI_COMM_NULL)
        MPI_Comm_free (&(this->io_comm));
}

static inline int set_info (e3sm_io_config &cfg, e3sm_io_decom &decom) {
    int err, nerrs = 0;
    MPI_Offset estimated_nc_ibuf_size;

    /* set MPI-IO hints */
    MPI_Info_set (cfg.info, "romio_cb_write", "enable");  /* collective write */
    MPI_Info_set (cfg.info, "romio_no_indep_rw", "true"); /* no independent MPI-IO */

    /* set PnetCDF I/O hints */
    MPI_Info_set (cfg.info, "nc_var_align_size", "1");     /* no gap between variables */
    MPI_Info_set (cfg.info, "nc_in_place_swap", "enable"); /* in-place byte swap */

    /* use total write amount to estimate nc_ibuf_size */
    estimated_nc_ibuf_size =
        decom.dims[2][0] * decom.dims[2][1] * sizeof (double) / cfg.num_iotasks;
    estimated_nc_ibuf_size *= cfg.nvars;
    if (estimated_nc_ibuf_size > 16777216) {
        char nc_ibuf_size_str[32];
        sprintf (nc_ibuf_size_str, "%lld", estimated_nc_ibuf_size);
        MPI_Info_set (cfg.info, "nc_ibuf_size", nc_ibuf_size_str);
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
    std::string help =
        "Usage: %s [OPTION]... FILE\n"
        "       [-h] Print help\n"
        "       [-v] Verbose mode\n"
        "       [-k] Keep the output files when program exits\n"
        "       [-d] Run test that uses PnetCDF vard API\n"
        "       [-n] Run test that uses PnetCDF varn API\n"
        "       [-m] Run test using noncontiguous write buffer\n"
        "       [-t] Write 2D variables followed by 3D variables\n"
        "       [-r num] Number of records (default 1)\n"
        "       [-s num] Stride between IO tasks (default 1)\n"
        "       [-o output_dir] Output directory name (default ./)\n"
        "       FILE: Name of input netCDF file describing data decompositions\n";
    fprintf (stderr, help.c_str (), argv0);
}

/*----< main() >-------------------------------------------------------------*/
int main (int argc, char **argv) {
    int err, nerrs = 0;
    int i;
    int ioproc;
    int nvars;

    MPI_Init (&argc, &argv);

    {
        e3sm_io_config cfg;
        e3sm_io_decom decom;
        e3sm_io_case *tcase    = NULL;
        e3sm_io_driver *driver = NULL;

        MPI_Comm_rank (MPI_COMM_WORLD, &(cfg.rank));
        MPI_Comm_size (MPI_COMM_WORLD, &(cfg.np));

        /* command-line arguments */
        while ((i = getopt (argc, argv, "vkr:s:l:o:i:c:dntRWH:ha:")) != EOF) switch (i) {
                case 'v':
                    cfg.verbose = 1;
                    break;
                case 'k':
                    cfg.keep_outfile = true;
                    break;
                case 'r':
                    cfg.nrec = atoi (optarg);
                    break;
                case 's':
                    cfg.io_stride = atoi (optarg);
                    break;
                case 'a':
                    if (std::string (optarg) == "pnc") {
                        cfg.api = pnc;
                    } else if (std::string (optarg) == "hdf5") {
                        cfg.api = hdf5;
                    }
                    break;
                case 'l':
                    if (std::string (optarg) == "contig") {
                        cfg.layout = contig;
                    } else if (std::string (optarg) == "chunk") {
                        cfg.layout = chunk;
                    }
                    break;
                case 'o':
                    cfg.targetdir = std::string (optarg);
                    break;
                case 'i':
                    cfg.datadir = std::string (optarg);
                    break;
                case 'c':
                    cfg.cfgpath = std::string (optarg);
                    break;
                case 'd':
                    cfg.low_lvl = true;
                    break;
                case 'n':
                    cfg.non_contig_buf = true;
                    break;
                case 't':
                    cfg.two_buf = true;
                    break;
                case 'R':
                    cfg.rd = true;
                    break;
                case 'W':
                    cfg.wr = true;
                    break;
                case 'H':
                    cfg.hx = atoi (optarg);
                    break;
                case 'h':
                default:
                    if (cfg.rank == 0) usage (argv[0]);
                    MPI_Finalize ();
                    return 1;
            }

        if (cfg.cfgpath == "") { /* input file is mandatory */
            if (!cfg.rank) usage (argv[0]);
            ERR_OUT ("Decomposition file not provided")
        }
        /* input file contains number of write requests and their file access
         * offsets (per array element) */
        PRINT_MSG (1, "input file name =%s\n", cfg.cfgpath.c_str ());

        /* neither command-line option -R or -W is used, run write */
        if (!(cfg.wr || cfg.rd)) cfg.wr = true;

        /* set the output folder name */
        PRINT_MSG (1, "Target folder name =%s\n", cfg.targetdir.c_str ());

        if (cfg.io_stride < 1) cfg.io_stride = 1;
        if (cfg.io_stride == 1) {
            cfg.num_iotasks = cfg.np;
            cfg.io_comm     = MPI_COMM_WORLD;
            ioproc          = 1;
        } else {
            int *ioranks;
            MPI_Group group, iogroup;

            /* assume that the IO root (rank of the first IO task) is 0 */
            cfg.num_iotasks = (cfg.np - 1) / cfg.io_stride + 1;

            /* create an array that holds the ranks of the IO tasks */
            ioranks = (int *)malloc (cfg.num_iotasks * sizeof (int));
            ioproc  = 0;
            for (i = 0; i < cfg.num_iotasks; i++) {
                ioranks[i] = i * cfg.io_stride;
                if (ioranks[i] == cfg.rank) ioproc = 1;
            }

            /* create a group for all MPI tasks */
            MPI_Comm_group (MPI_COMM_WORLD, &group);

            /* create a sub-group for the IO tasks */
            MPI_Group_incl (group, cfg.num_iotasks, ioranks, &iogroup);

            /* create an MPI communicator for the IO tasks */
            MPI_Comm_create (MPI_COMM_WORLD, iogroup, &(cfg.io_comm));

            MPI_Group_free (&iogroup);
            MPI_Group_free (&group);
            free (ioranks);
        }

        /* only IO tasks call PnetCDF APIs */
        if (!ioproc) goto err_out;

        nerrs += set_info (cfg, decom);
        CHECK_NERR

        switch (cfg.api) {
            case pnc:
                driver = new e3sm_io_driver_pnc ();
                break;
            case hdf5:
                driver = new e3sm_io_driver_hdf5 ();
                break;
            default:
                RET_ERR ("Unknown driver")
                break;
        }

        /* read request information from decompositions 1, 2 and 3 */
        err = read_decomp (cfg.verbose, cfg.io_comm, cfg.cfgpath.c_str (), &(decom.num_decomp),
                           decom.dims, decom.contig_nreqs, decom.disps, decom.blocklens);
        CHECK_ERR

        /* F case has 3 decompositions, G case has 6 */
        if (decom.num_decomp == 3) {
            cfg.nvars = 414;
            tcase     = new e3sm_io_case_F ();
        } else if (decom.num_decomp == 6) {
            cfg.nvars = 52;
            tcase     = new e3sm_io_case_G ();
        } else {
            RET_ERR ("Unknown decom file")
        }

        if (cfg.wr) { tcase->wr_test (cfg, decom, *driver); }
        if (cfg.rd) { tcase->rd_test (cfg, decom, *driver); }

    err_out:
        delete driver;
        delete tcase;

        /* Non-IO tasks wait for IO tasks to complete */
        MPI_Barrier (MPI_COMM_WORLD);
    }

    MPI_Finalize ();

    return err;
}
