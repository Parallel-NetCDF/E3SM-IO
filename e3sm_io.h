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

#include <stdio.h>

#include <mpi.h>
#include <pnetcdf.h>

/* minimum PnetCDF version required is 1.10.0 */
#if (PNETCDF_VERSION_MAJOR*1000000 + PNETCDF_VERSION_MINOR*1000 + PNETCDF_VERSION_SUB < 1010000)
#error "PnetCDF 1.10.0 and later is required to build this program"
#endif

#define MAX_NUM_DECOMP 6

int verbose; /* verbose mode to print additional messages on screen */
int keep_outfile; /* whether to keep the output files when exits */
int two_buf;

#ifndef MAX
#define MAX(a,b) ((a) > (b)) ? (a) : (b)
#endif
#ifndef MIN
#define MIN(a,b) ((a) < (b)) ? (a) : (b)
#endif

#define ERR { \
    if (err != NC_NOERR) { \
        printf("Error in %s:%d: %s\n", __FILE__, __LINE__, \
               ncmpi_strerror(err)); \
        nerrs++; \
        goto fn_exit; \
    } \
}

/* In E3SM production runs, the write buffers are of type double in memory,
 * and the variables stored in NetCDF files are of type NC_FLOAT. This default
 * behavior can be changed (i.e. memory buffer of type float) by removing the
 * line "#define _DOUBLE_TYPE_" below.
 */
#define _DOUBLE_TYPE_

#ifdef _DOUBLE_TYPE_
typedef double itype;  /* internal data type of buffer in memory */
#define REC_ITYPE MPI_DOUBLE
#define REC_XTYPE NC_DOUBLE
#else
typedef float itype;   /* internal data type of buffer in memory */
#define REC_ITYPE MPI_FLOAT
#define REC_XTYPE NC_FLOAT
#endif

extern int
read_decomp(MPI_Comm io_comm, const char *infname, int *num_decomp, MPI_Offset dims[][2],
            int contig_nreqs[3], int *disps[3], int *blocklens[3]);

extern void
print_info(MPI_Info *info_used);

extern int
def_F_case_h0(int ncid, const MPI_Offset dims[2], int nvars, int *varids);
extern int 
inq_F_case_h0(  int ncid,                  /* file ID */
                MPI_Offset  dims[2],       /* dimension sizes */
                int               nvars,   /* number of variables */
                int              *varids); /* variable IDs */

extern int
def_F_case_h1(int ncid, const MPI_Offset dims[2], int nvars, int *varids);
extern int
inq_F_case_h1(  int               ncid,    /* file ID */
                MPI_Offset  dims[2],       /* dimension sizes */
                int               nvars,   /* number of variables */
                int              *varids) ;/* variable IDs */

extern int
run_vard_F_case(MPI_Comm io_comm,         /* MPI communicator that includes all the tasks involved in IO */
                const char *out_dir,      /* output folder name */
                const char *outfile,      /* output file name */
                int         nvars,        /* number of variables 408 or 51 */
                int         num_recs,     /* number of records */
                int         noncontig_buf,/* whether to us noncontiguous buffer */
                MPI_Info    info,
                MPI_Offset  dims[3][2],   /* dimension lengths */
                const int   nreqs[3],     /* no. request in decompositions 1,2,3 */
                int* const  disps[3],     /* request's displacements */
                int* const  blocklens[3]);/* request's block lengths */

extern int
run_varn_F_case(MPI_Comm io_comm,         /* MPI communicator that includes all the tasks involved in IO */
                const char *out_dir,      /* output folder name */
                const char *outfile,      /* output file name */
                int         nvars,        /* number of variables 408 or 51 */
                int         num_recs,     /* number of records */
                int         noncontig_buf,/* whether to us noncontiguous buffer */
                MPI_Info    info,
                MPI_Offset  dims[3][2],   /* dimension lengths */
                const int   nreqs[3],     /* no. request in decompositions 1,2,3 */
                int* const  disps[3],     /* request's displacements */
                int* const  blocklens[3], /* request's block lengths */
                double     *dbl_bufp,     /* buffer for fixed size double var */
                itype      *rec_bufp,     /* buffer for rec floating point var */
                char       *txt_buf,      /* buffer for char var */
                int        *int_buf);     /* buffer for int var */
extern int
run_varn_F_case_rd( MPI_Comm io_comm,         /* MPI communicator that includes all the tasks involved in IO */
                    char       *out_dir,      /* output folder name */
                    char       *outfile,      /* output file name */
                    int         nvars,        /* number of variables 408 or 51 */
                    int         num_recs,     /* number of records */
                    int         noncontig_buf,/* whether to us noncontiguous buffer */
                    MPI_Info    info,
                    MPI_Offset  dims[3][2],   /* dimension lengths */
                    int         nreqs[3],     /* no. request in decompositions 1,2,3 */
                    int        *disps[3],     /* request's displacements */
                    int        *blocklens[3], /* request's block lengths */
                    double     **dbl_bufp,    /* buffer for fixed size double var */
                    itype      **rec_bufp,    /* buffer for rec floating point var */
                    char       *txt_buf,      /* buffer for char var */
                    int        *int_buf);     /* buffer for int var */

extern int
def_G_case_h0(int               ncid,       /* file ID */
              const MPI_Offset  dims_D1[1], /* dimension sizes of decomposition 1 */
              const MPI_Offset  dims_D2[1], /* dimension sizes of decomposition 2 */
              const MPI_Offset  dims_D3[2], /* dimension sizes of decomposition 3 */
              const MPI_Offset  dims_D4[2], /* dimension sizes of decomposition 4 */
              const MPI_Offset  dims_D5[2], /* dimension sizes of decomposition 5 */
              const MPI_Offset  dims_D6[2], /* dimension sizes of decomposition 6 */
              int               nvars,      /* number of variables */
              int              *varids);    /* variable IDs */
extern int
inq_G_case_h0(int         ncid,       /* file ID */
              MPI_Offset  dims_D1[1], /* dimension sizes of decomposition 1 */
              MPI_Offset  dims_D2[1], /* dimension sizes of decomposition 2 */
              MPI_Offset  dims_D3[2], /* dimension sizes of decomposition 3 */
              MPI_Offset  dims_D4[2], /* dimension sizes of decomposition 4 */
              MPI_Offset  dims_D5[2], /* dimension sizes of decomposition 5 */
              MPI_Offset  dims_D6[2], /* dimension sizes of decomposition 6 */
              int         nvars,      /* number of variables */
              int        *varids);     /* variable IDs */

extern int
run_varn_G_case(MPI_Comm io_comm,         /* MPI communicator that includes all the tasks involved in IO */
                const char *out_dir,      /* output folder name */
                const char *outfile,      /* output file name */
                int         nvars,        /* number of variables 51 */
                int         num_recs,     /* number of records */
                MPI_Info    info,
                MPI_Offset  dims[6][2],   /* dimension lengths decomposition 1-6 */
                const int   nreqs[6],     /* no. request in decomposition 1-6 */
                int* const  disps[6],     /* request's displacements */
                int* const  blocklens[6], /* request's block lengths */
                int *D1_fix_int_bufp,     /* D1 fix int buffer */
                int *D2_fix_int_bufp,     /* D2 fix int buffer */
                int *D3_fix_int_bufp,     /* D3 fix int buffer */
                int *D4_fix_int_bufp,     /* D4 fix int buffer */
                int *D5_fix_int_bufp,     /* D5 fix int buffer */
                double *D1_rec_dbl_bufp,  /* D1 rec double buffer */
                double *D3_rec_dbl_bufp,  /* D3 rec double buffer */
                double *D4_rec_dbl_bufp,  /* D4 rec double buffer */
                double *D5_rec_dbl_bufp,  /* D5 rec double buffer */
                double *D6_rec_dbl_bufp,  /* D6 rec double buffer */
                double *D1_fix_dbl_bufp); /* D1 fix double buffer */
extern int
run_varn_G_case_rd( MPI_Comm io_comm,         /* MPI communicator that includes all the tasks involved in IO */
                    char       *out_dir,      /* output folder name */
                    char       *outfile,      /* output file name */
                    int         nvars,        /* number of variables 51 */
                    int         num_recs,     /* number of records */
                    MPI_Info    info,
                    MPI_Offset  dims[6][2],   /* dimension lengths decomposition 1-6 */
                    int         nreqs[6],     /* no. request in decomposition 1-6 */
                    int        *disps[6],     /* request's displacements */
                    int        *blocklens[6], /* request's block lengths */
                    int **D1_fix_int_bufp,     /* D1 fix int buffer */
                    int **D2_fix_int_bufp,     /* D2 fix int buffer */
                    int **D3_fix_int_bufp,     /* D3 fix int buffer */
                    int **D4_fix_int_bufp,     /* D4 fix int buffer */
                    int **D5_fix_int_bufp,     /* D5 fix int buffer */
                    double **D1_rec_dbl_bufp,  /* D1 rec double buffer */
                    double **D3_rec_dbl_bufp,  /* D3 rec double buffer */
                    double **D4_rec_dbl_bufp,  /* D4 rec double buffer */
                    double **D5_rec_dbl_bufp,  /* D5 rec double buffer */
                    double **D6_rec_dbl_bufp,  /* D6 rec double buffer */
                    double **D1_fix_dbl_bufp); /* D1 fix double buffer */

#endif

