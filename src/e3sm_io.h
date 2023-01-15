/*********************************************************************
 *
 * Copyright (C) 2022, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 *********************************************************************/

#ifndef _H_E3SM_IO_
#define _H_E3SM_IO_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <mpi.h>

#ifdef ENABLE_NETCDF4
#include <netcdf.h>
#include <netcdf_par.h>
#endif
#ifdef ENABLE_PNC
#include <pnetcdf.h>
#endif

#define E3SM_IO_MAX_PATH    1024
#define MAX_NUM_DECOMP      6

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

#if !defined(PNETCDF_VERSION) && !defined(_NETCDF_)
#define NC_NAT          0       /**< Not A Type */
#define NC_BYTE         1       /**< signed 1 byte integer */
#define NC_CHAR         2       /**< ISO/ASCII character */
#define NC_SHORT        3       /**< signed 2 byte integer */
#define NC_INT          4       /**< signed 4 byte integer */
#define NC_LONG         NC_INT
#define NC_FLOAT        5       /**< single precision floating point number */
#define NC_DOUBLE       6       /**< double precision floating point number */
#define NC_UBYTE        7       /**< unsigned 1 byte int */
#define NC_USHORT       8       /**< unsigned 2-byte int */
#define NC_UINT         9       /**< unsigned 4-byte int */
#define NC_INT64        10      /**< signed 8-byte int */
#define NC_UINT64       11      /**< unsigned 8-byte int */
#define NC_STRING       12      /**< string */

#define NC_NOERR        0   /**< No Error */
#define NC_GLOBAL 	-1
#define NC_UNLIMITED 0L
#define _FillValue      "_FillValue"

typedef int nc_type;
#endif

/* In E3SM production runs, the write buffers are of type double in memory,
 * and the variables stored in NetCDF files are of type NC_FLOAT. This default
 * behavior can be changed (i.e. memory buffer of type float) by removing the
 * line "#define _DOUBLE_TYPE_" below.
 */
#define _DOUBLE_TYPE_

#ifdef _DOUBLE_TYPE_
typedef double vtype; /* internal data type of buffer in memory */
#define REC_ITYPE MPI_DOUBLE
#define REC_XTYPE NC_DOUBLE
#else
typedef float vtype; /* internal data type of buffer in memory */
#define REC_ITYPE MPI_FLOAT
#define REC_XTYPE NC_FLOAT
#endif

typedef enum e3sm_io_api {
    pnetcdf,
    netcdf4,
    hdf5,
    hdf5_md,
    hdf5_log,
    adios,
    undef_api
} e3sm_io_api;

typedef enum e3sm_io_strategy { canonical, blob, log, undef_io } e3sm_io_strategy;

typedef enum e3sm_io_filter { none, deflate, bzip2 } e3sm_io_filter;

typedef enum { h0, h1 } history;

#define NVARS_F_CASE_H0 414
#define NVARS_F_CASE_H1  51
#define NVARS_G_CASE     52
#define NVARS_I_CASE_H0 560
#define NVARS_I_CASE_H1 552

typedef enum { F, G, I, unknown } climate_case;

typedef struct {
    /* statistics */
    char outfile[1040];
    int nvars;        /* number of climate variables */
    int nrecs;        /* number of time records */
    int ffreq;        /* I/O flush frequency */
    int num_attrs;    /* number of attributes */
    int num_flushes;  /* number of flush called */
    int num_decomp_vars;         /* no. climate variables decomposed */
    int nvars_D[MAX_NUM_DECOMP]; /* no. climate variables per decomposition */
    MPI_Offset metadata_WR;
    MPI_Offset metadata_RD;
    MPI_Offset amount_WR;
    MPI_Offset amount_RD;
    MPI_Offset my_nreqs;
    MPI_Offset file_size;
    MPI_Info   info_used;

    /* timings */
    double pre_time;
    double open_time;
    double def_time;
    double post_time;
    double flush_time;
    double close_time;
    double end2end_time;
} case_meta;

typedef struct e3sm_io_config {
    int rank;
    int np;
    MPI_Comm io_comm;
    MPI_Info info;
    int num_iotasks;
    int num_subfiles;

    char in_path[E3SM_IO_MAX_PATH];
    char out_path[E3SM_IO_MAX_PATH];
    char decomp_path[E3SM_IO_MAX_PATH];
    int hx;
    int wr;
    int rd;
    int nvars;  /* number of climate variables */

    e3sm_io_strategy strategy;
    e3sm_io_api api;
    e3sm_io_filter filter;
    size_t chunksize;

    climate_case run_case;
    history hist;
    int verbose;      /* verbose mode to print additional messages on screen */
    int keep_outfile; /* whether to keep the output files when exits */
    int fill_mode;    /* fill missing elements in decomposition maps */

    int two_buf;
    int non_contig_buf;
    int io_stride;
    int comp_time;   /* Emulate computation time (sleep) for a time step */
    int profiling;

    /* below are used for PnetCDF blob I/O subfiling */
    int      subfile_ID;   /* unique file identifier for subfiles */
    MPI_Comm sub_comm;     /* communicator for a subfile */
    int      sub_rank;     /* rank in sub_comm */
    int      sub_nprocs;   /* nprocs in sub_comm, also nblobs in a subfile */

    char node_info[2048]; /* info about the number of compute nodes and the
                           * number of MPI processes running per node
                           */
    case_meta F_case_h0;
    case_meta F_case_h1;
    case_meta G_case;
    case_meta I_case_h0;
    case_meta I_case_h1;

    int   env_log;
    int   env_log_passthru;
    char *env_log_info;
    int   env_cache;
    int   env_async;
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

/* maximum number of dimensions of a climate variable */
#define MAX_NDIMS    4

typedef struct e3sm_io_decom {
    int  num_decomp;                   /* number of decompositions: 3 or 6 */
    int  ndims[MAX_NUM_DECOMP];        /* number of dimensions in each
                                          decomposition */
    MPI_Offset dims[MAX_NUM_DECOMP][MAX_NDIMS];
                                       /* global dimension sizes of each
                                          decomposition */

    /* below are for within a subfile, if subfile option is enabled. If not
     * then, they are for the single shared file.
     */
    int  contig_nreqs[MAX_NUM_DECOMP]; /* number of noncontiguous requests in
                                          each decomposition */
    int  max_nreqs[MAX_NUM_DECOMP];    /* max contig_nreqs[] among processes */
    int *disps[MAX_NUM_DECOMP];        /* starting array index offset of
                                          requests to be read/written by this
                                          process */
    int *blocklens[MAX_NUM_DECOMP];    /* length (number of array elements) of
                                          individual requests */

    /* the following 4 are the starts[] and counts[] used in varn APIs */
    MPI_Offset **w_starts[MAX_NUM_DECOMP]; /* [nreqs][] record variables */
    MPI_Offset **w_counts[MAX_NUM_DECOMP]; /* [nreqs][] record variables */
    MPI_Offset **w_startx[MAX_NUM_DECOMP]; /* [nreqs][] fixed-size variables */
    MPI_Offset **w_countx[MAX_NUM_DECOMP]; /* [nreqs][] fixed-size variables */

    /* below 3 are used for PnetCDF blob I/O subfiling */
    MPI_Offset nelems[MAX_NUM_DECOMP]; /* total no. array elements of each
                                          decomposition to store in a subfile */
    MPI_Offset start[MAX_NUM_DECOMP];  /* This proc's starting array index for
                                          variables using each decomposition */
    MPI_Offset count[MAX_NUM_DECOMP];  /* This proc's number of array elements
                                          for variables using each
                                          decomposition */

    /* below 2 are used for Scorpio blob I/O, which saves only offsets, but
     * no blocklens[], i.e. all blocklens are of length 1 */
    MPI_Offset  raw_nreqs[MAX_NUM_DECOMP];
    MPI_Offset *raw_offsets[MAX_NUM_DECOMP];

} e3sm_io_decom;

#ifdef __cplusplus
extern "C" {
#endif
int read_decomp(e3sm_io_config *cfg, e3sm_io_decom *decom);
#ifdef ENABLE_PNC
extern int read_decomp_pnc(e3sm_io_config *cfg, e3sm_io_decom *decom);
#endif
#ifdef ENABLE_NETCDF4
extern int read_decomp_nc4(e3sm_io_config *cfg, e3sm_io_decom *decom);
#endif
extern int calc_metadata(e3sm_io_config *cfg, e3sm_io_decom *decom);
extern void print_info (MPI_Info *info_used);
extern int e3sm_io_core (e3sm_io_config *cfg, e3sm_io_decom *decom);
extern int e3sm_io_xlen_nc_type(nc_type xtype, int *size);
extern int report_timing_WR(e3sm_io_config *cfg, e3sm_io_decom *decom);
extern int report_timing_RD(e3sm_io_config *cfg, e3sm_io_decom *decom);
#ifdef __cplusplus
}
#endif

#endif
