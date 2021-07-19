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

typedef enum e3sm_io_api {
    pnetcdf,
    hdf5,
    hdf5_md,
    hdf5_log,
    adios,
    undef_api
} e3sm_io_api;

typedef enum e3sm_io_strategy { canonical, blob, log, undef_io } e3sm_io_strategy;

typedef enum e3sm_io_filter { none, deflate, bzip2 } e3sm_io_filter;

typedef struct e3sm_io_config {
    int rank;
    int np;
    MPI_Comm io_comm;
    MPI_Info info;
    int num_iotasks;
    int num_group;

    char in_path[E3SM_IO_MAX_PATH];
    char out_path[E3SM_IO_MAX_PATH];
    char cfg_path[E3SM_IO_MAX_PATH];
    int hx;
    int nrec;
    int wr;
    int rd;
    int nvars;

    e3sm_io_strategy strategy;
    e3sm_io_api api;
    e3sm_io_filter filter;
    size_t chunksize;
    int vard;

    int verbose;      /* verbose mode to print additional messages on screen */
    int keep_outfile; /* whether to keep the output files when exits */

    int two_buf;
    int non_contig_buf;
    int io_stride;

    /* below 3 are used for PnetCDF blob I/O subfiling */
    int      num_subfiles; /* number of subfiles */
    int      subfile_ID;   /* unqiue file identifier for subfiles */
    MPI_Comm sub_comm;     /* communicator for a subfile */

    /* statistics */
    char *outfile;
    int num_flushes;
    int num_decomp;
    int num_decomp_vars;
    int nvars_D[MAX_NUM_DECOMP];
    MPI_Offset metadata_WR;
    MPI_Offset amount_WR;
    MPI_Offset amount_RD;
    MPI_Offset my_nreqs;

    /* timings */
    double pre_time;
    double open_time;
    double def_time;
    double post_time;
    double flush_time;
    double close_time;
    double end2end_time;

} e3sm_io_config;


/* NVARS_DECOMP is the number of new variables added to each subfile that
 * describe a decomposition. For decomposition D*, they are:
 * int   D*.nreqs(nblobs) ;
 * int64 D*.blob_start(nblobs) ;
 * int64 D*.blob_count(nblobs) ;
 * int   D*.offsets(nblobs, D*.max_nreqs) ;
 * int   D*.lengths(nblobs, D*.max_nreqs) ;
 */
#define NVARS_DECOMP 5

typedef struct e3sm_io_decom {
    int  num_decomp;                   /* number of decompositions: 3 or 6 */
    int  contig_nreqs[MAX_NUM_DECOMP]; /* number of noncontiguous requests in
                                          each decomposition */
    int  max_nreqs[MAX_NUM_DECOMP];    /* max contig_nreqs[] among processes */
    int *disps[MAX_NUM_DECOMP];        /* starting array index offset of
                                          requests to be read/written by this
                                          process */
    int *blocklens[MAX_NUM_DECOMP];    /* length (number of array elements) of
                                          individual requests */
    int  ndims[MAX_NUM_DECOMP];        /* number of dimensions in each
                                          decomposition */
    MPI_Offset dims[MAX_NUM_DECOMP][2];/* global dimension sizes of each
                                          decomposition */

    /* below 3 are used for PnetCDF blob I/O subfiling */
    MPI_Offset nelems[MAX_NUM_DECOMP]; /* total no. array elements in each
                                          decomposition */
    MPI_Offset start[MAX_NUM_DECOMP];  /* This proc's starting array index for
                                          variables using each decomposition */
    MPI_Offset count[MAX_NUM_DECOMP];  /* This proc's number of array elements
                                          for variables using each
                                          decomposition */

    /* below 3 are used for Scorpio blob I/O, which saves only offsets, but
     * no blocklens[], i.e. all blocklens are of length 1 */
    MPI_Offset raw_nreqs[MAX_NUM_DECOMP];
    MPI_Offset *raw_offsets[MAX_NUM_DECOMP];
} e3sm_io_decom;

#ifdef __cplusplus
extern "C" {
#endif
extern int read_decomp(e3sm_io_config *cfg, e3sm_io_decom *decom);
extern int blob_metadata(e3sm_io_config *cfg, e3sm_io_decom *decom);
extern void print_info (MPI_Info *info_used);
extern int e3sm_io_core (e3sm_io_config *cfg, e3sm_io_decom *decom);
#ifdef __cplusplus
}
#endif

#endif
