/*********************************************************************
 *
 * Copyright (C) 2018, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 *********************************************************************/

#include <stdio.h>

#include <mpi.h>
#include <pnetcdf.h>

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

#ifdef _DOUBLE_TYPE_
typedef double dtype;
#define REC_DTYPE MPI_DOUBLE
#define REC_XTYPE NC_DOUBLE
#else
typedef float dtype;
#define REC_DTYPE MPI_FLOAT
#define REC_XTYPE NC_FLOAT
#endif

extern
int e3sm_io_header(int ncid, MPI_Offset dims[2], int nvars, int *varids);

extern
int e3sm_io_header1(int ncid, MPI_Offset dims[2], int nvars, int *varids);

