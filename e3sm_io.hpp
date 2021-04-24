/*********************************************************************
 *
 * Copyright (C) 2018, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 *********************************************************************/

#ifndef _H_E3SM_IO_
#define _H_E3SM_IO_

#include <mpi.h>
#include <pnetcdf.h>

#include <cstdio>
#include <string>

#define E3SM_IO_GLOBAL_ATTR -1

#define MAX_NUM_DECOMP 6

#define PRINT_MSG(V, ...)                                                    \
    {                                                                        \
        if ((cfg.rank == 0) && (cfg.verbose >= V)) { printf (__VA_ARGS__); } \
    }

typedef enum e3sm_io_api { pnc, hdf5, adios2 } e3sm_io_api;

typedef enum e3sm_io_layout { contig, chunk } e3sm_io_layout;

typedef struct e3sm_io_config {
    int rank;
    int np;
    MPI_Comm io_comm = MPI_COMM_WORLD;
    MPI_Info info    = MPI_INFO_NULL;
    int num_iotasks;

    std::string targetdir = "./";
    std::string datadir   = "";
    std::string cfgpath   = "";
    int hx                = 0;
    int nrec              = 0;
    bool wr               = false;
    bool rd               = false;
    int nvars             = 0;

    e3sm_io_api api;
    e3sm_io_layout layout;
    bool low_lvl = false;

    int verbose       = 0;     /* verbose mode to print additional messages on screen */
    bool keep_outfile = false; /* whether to keep the output files when exits */

    bool two_buf        = false;
    bool non_contig_buf = false;
    int io_stride       = 1;

    e3sm_io_config ();
    ~e3sm_io_config ();
} e3sm_io_config;

typedef struct e3sm_io_decom {
    int num_decomp;
    int contig_nreqs[MAX_NUM_DECOMP], *disps[MAX_NUM_DECOMP];
    int *blocklens[MAX_NUM_DECOMP];
    MPI_Offset dims[MAX_NUM_DECOMP][2];
} e3sm_io_decom;

/* In E3SM production runs, the write buffers are of type double in memory,
 * and the variables stored in NetCDF files are of type NC_FLOAT. This default
 * behavior can be changed (i.e. memory buffer of type float) by removing the
 * line "#define _DOUBLE_TYPE_" below.
 */
#define _DOUBLE_TYPE_

#ifdef _DOUBLE_TYPE_
typedef double itype; /* internal data type of buffer in memory */
#define REC_ITYPE MPI_DOUBLE
#define REC_XTYPE NC_DOUBLE
#else
typedef float itype; /* internal data type of buffer in memory */
#define REC_ITYPE MPI_FLOAT
#define REC_XTYPE NC_FLOAT
#endif

extern int read_decomp (int verbose,
                        MPI_Comm io_comm,
                        const char *infname,
                        int *num_decomp,
                        MPI_Offset dims[][2],
                        int contig_nreqs[3],
                        int *disps[3],
                        int *blocklens[3]);

extern void print_info (MPI_Info *info_used);

#endif
