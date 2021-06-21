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
#include <stdio.h>

#define E3SM_IO_MAX_PATH    1024
#define MAX_NUM_DECOMP      6
#define E3SM_IO_GLOBAL_ATTR -1

#define PRINT_MSG(V, ...)                                                    \
    {                                                                        \
        if ((cfg.rank == 0) && (cfg.verbose >= V)) { printf (__VA_ARGS__); } \
    }

#ifndef MAX
#define MAX(a, b) ((a) > (b)) ? (a) : (b)
#endif
#ifndef MIN
#define MIN(a, b) ((a) < (b)) ? (a) : (b)
#endif
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

typedef enum e3sm_io_api { pnc, hdf5_native, hdf5_multi, hdf5_logvol, adios2 } e3sm_io_api;

typedef enum e3sm_io_filter { none, deflate, bzip2 } e3sm_io_filter;

typedef struct e3sm_io_config {
    int rank;
    int np;
    MPI_Comm io_comm;
    MPI_Info info;
    int num_iotasks;

    char *targetdir;
    char *datadir;
    char *cfgpath;
    int hx;
    int nrec;
    int wr;
    int rd;
    int nvars;

    e3sm_io_api api;
    e3sm_io_filter filter;
    size_t chunksize;
    int vard;

    int verbose;      /* verbose mode to print additional messages on screen */
    int keep_outfile; /* whether to keep the output files when exits */

    int two_buf;
    int non_contig_buf;
    int io_stride;
} e3sm_io_config;

typedef struct e3sm_io_decom {
    int num_decomp;
    int contig_nreqs[MAX_NUM_DECOMP], *disps[MAX_NUM_DECOMP];
    int *blocklens[MAX_NUM_DECOMP];
    MPI_Offset dims[MAX_NUM_DECOMP][2];
} e3sm_io_decom;

#ifdef __cplusplus
extern "C" {
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
int e3sm_io_core (e3sm_io_config *cfg, e3sm_io_decom *decom);
#ifdef __cplusplus
}
#endif

#endif
