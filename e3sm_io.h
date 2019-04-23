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
read_decomp(const char *infname, int *num_decomp, MPI_Offset dims[][2],
            int contig_nreqs[3], int *disps[3], int *blocklens[3]);

extern void
print_info(MPI_Info *info_used);

extern int
def_F_case_h0(int ncid, MPI_Offset dims[2], int nvars, int *varids);

extern int
def_F_case_h1(int ncid, MPI_Offset dims[2], int nvars, int *varids);

extern int
inq_F_case_h0(int ncid, MPI_Offset dims[2], int nvars, int *varids);

extern int
inq_F_case_h1(int ncid, MPI_Offset dims[2], int nvars, int *varids);

extern int
run_vard_F_case(char *out_dir, char *outfile, int nvars, int num_recs,
		int noncontig_buf, MPI_Info info, MPI_Offset dims[3][2],
                int nreqs[3], int *disps[3], int *blocklens[3]);
extern int
run_varn_F_case(char *out_dir, char *outfile, int nvars, int num_recs,
		int noncontig_buf, MPI_Info info, MPI_Offset dims[3][2],
                int nreqs[3], int *disps[3], int *blocklens[3], double*, itype*, char*, int*);
extern int
run_varn_F_case_rd(char *out_dir, char *outfile, int nvars, int num_recs,
		int noncontig_buf, MPI_Info info, MPI_Offset dims[3][2],
                int nreqs[3], int *disps[3], int *blocklens[3], double**, itype**, char*, int*);
extern int
def_G_case_h0(int ncid, MPI_Offset dims_D1[1], MPI_Offset dims_D2[1],
	      MPI_Offset dims_D3[2], MPI_Offset dims_D4[2],
	      MPI_Offset dims_D5[2], MPI_Offset dims_D6[2], int nvars,
              int *varids);

extern int
run_varn_G_case(char *out_dir, char *outfile, int nvars, int num_recs,
		MPI_Info info, MPI_Offset dims[6][2], int nreqs[6],
                int *disps[6], int *blocklens[6], int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*);

extern int
run_varn_G_case_rd(char *out_dir, char *outfile, int nvars, int num_recs,
		MPI_Info info, MPI_Offset dims[6][2], int nreqs[6],
                int *disps[6], int *blocklens[6], int**, int**, int**, int**, int**, double**, double**, double**, double**, double**, double**);

#endif

