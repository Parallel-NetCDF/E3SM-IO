/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 * This file is extracted from PnetCDF source file ncmpio_NC.h in:
 * https://github.com/Parallel-NetCDF/PnetCDF/blob/master/src/drivers/ncmpio
 *
 */

#ifndef _BLOB_NCMPIO_H
#define _BLOB_NCMPIO_H

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>     /* size_t */

#include <mpi.h>

#if defined(__cplusplus)
extern "C" {
#endif

#define NC_NOERR 0
#define NC_UNLIMITED 0L
#define NC_GLOBAL -1

#define NCI_Malloc malloc
#define NCI_Calloc calloc
#define NCI_Realloc realloc
#define NCI_Free free

#define X_ALIGN 4
#define _RNDUP(x, unit)      ((((x) + (unit) - 1) / (unit)) * (unit))

#define X_SIZEOF_BYTE           1
#define X_SIZEOF_CHAR           1
#define X_SIZEOF_SHORT          2
#define X_SIZEOF_INT            4       /* xdr_int */
#define X_SIZEOF_FLOAT          4
#define X_SIZEOF_DOUBLE         8

/* additional data types in CDF-5 */
#define X_SIZEOF_UBYTE          1
#define X_SIZEOF_USHORT         2
#define X_SIZEOF_UINT           4
#define X_SIZEOF_LONGLONG       8
#define X_SIZEOF_ULONGLONG      8
#define X_SIZEOF_INT64          8
#define X_SIZEOF_UINT64         8

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifndef HAVE_SCHAR
typedef signed char schar;
#endif
#ifndef HAVE_UCHAR
typedef unsigned char uchar;
#endif
#ifndef HAVE_USHORT
typedef unsigned short int ushort;
#endif
#ifndef HAVE_UINT
typedef unsigned int uint;
#endif
#ifndef HAVE_LONGLONG
typedef long long longlong;
#endif
#ifndef HAVE_ULONGLONG
typedef unsigned long long ulonglong;
#endif
#ifndef HAVE_INT64
typedef long long int64;
#endif
#ifndef HAVE_UINT64
typedef unsigned long long uint64;
#endif

#define FILE_ALIGNMENT_DEFAULT 512
#define FILE_ALIGNMENT_LB      4

#ifndef NC_ARRAY_GROWBY
#define NC_ARRAY_GROWBY 64
#endif

#define MIN_NC_XSZ 32

typedef enum {
    NC_UNSPECIFIED =  0,  /* ABSENT */
    NC_DIMENSION   = 10,  /* \x00 \x00 \x00 \x0A */
    NC_VARIABLE    = 11,  /* \x00 \x00 \x00 \x0B */
    NC_ATTRIBUTE   = 12   /* \x00 \x00 \x00 \x0C */
} NC_tag;

/* netcdf file format:
     netcdf_file  = header  data
     header       = magic  numrecs  dim_list  gatt_list  var_list
     magic        = 'C'  'D'  'F'  VERSION
     VERSION      = \x01 | \x02 | \x05
     numrecs      = NON_NEG | STREAMING
     dim_list     = ABSENT | NC_DIMENSION  nelems  [dim ...]
     gatt_list    = att_list
     att_list     = ABSENT | NC_ATTRIBUTE  nelems  [attr ...]
     var_list     = ABSENT | NC_VARIABLE   nelems  [var ...]
     ABSENT       = ZERO  ZERO                  // Means list is not present
     ZERO         = \x00 \x00 \x00 \x00         // 32-bit zero

  Minimum happens when nothing is defined, i.e.
     magic              -- 4 bytes
     numrecs            -- 4 bytes for CDF-1 and CDF-2, 8 bytes for CDF-5
     dim_list = ABSENT  -- 8 bytes
     gatt_list = ABSENT -- 8 bytes
     var_list = ABSENT  -- 8 bytes
*/

typedef struct {
    MPI_Offset  size;
    size_t      name_len; /* strlen(name), for faster string compare */
    char       *name;
} NC_dim;

typedef struct NC_dimarray {
    int            ndefined;      /* number of defined dimensions */
    int            unlimited_id;  /* -1 for not defined, otherwise >= 0 */
    NC_dim       **value;
} NC_dimarray;

typedef struct {
    MPI_Offset nelems;   /* number of attribute elements */
    MPI_Offset xsz;      /* amount of space at xvalue (4-byte aligned) */
    nc_type    xtype;    /* external NC data type of the attribute */
    size_t     name_len; /* strlen(name) for faster string compare */
    char      *name;     /* name of the attributes */
    void      *xvalue;   /* the actual data, in external representation */
} NC_attr;

typedef struct {
    int            ndefined;  /* number of defined attributes */
    NC_attr      **value;
} NC_attrarray;

/*
 * NC variable: description and data
 */
typedef struct {
    int           varid;   /* variable ID */
    int           xsz;     /* byte size of 1 array element */
    nc_type       xtype;   /* variable's external NC data type */
    int           no_fill; /* whether fill mode is disabled */
    size_t        name_len;/* strlen(name) for faster string compare */
    char         *name;    /* name of the variable */
    int           ndims;   /* number of dimensions */
    int          *dimids;  /* [ndims] array of dimension IDs */
    MPI_Offset   *shape;   /* [ndims] dim->size of each dim
                              shape[0] == NC_UNLIMITED if record variable */
    MPI_Offset   *dsizes;  /* [ndims] the right to left product of shape */
    MPI_Offset    begin;   /* starting file offset of this variable */
    MPI_Offset    len;     /* this is the "vsize" defined in header format, the
                              total size in bytes of the array variable.
                              For record variable, this is the record size */
    NC_attrarray  attrs;   /* attribute array */
} NC_var;

typedef struct {
    int            ndefined;    /* number of defined variables */
    int            num_rec_vars;/* number of defined record variables */
    NC_var       **value;
} NC_vararray;

#define IS_RECVAR(vp) \
        ((vp)->shape != NULL ? (*(vp)->shape == NC_UNLIMITED) : 0 )

typedef struct {
    int           ncid;     /* file ID */
    int           format;   /* 1, 2, or 5 corresponding to CDF-1, 2, or 5 */
    MPI_Offset    xsz;      /* size of this file header, <= var[0].begin */
    MPI_Offset    recsize;  /* length of 'record': sum of single record sizes
                               of all the record variables */
    MPI_Offset    numrecs;  /* number of 'records' allocated */

    NC_dimarray   dims;     /* dimensions defined */
    NC_attrarray  attrs;    /* global attributes defined */
    NC_vararray   vars;     /* variables defined */
} NC;


typedef struct {
    MPI_Offset  offset;   /* current read/write offset in the file */
    int         size;     /* allocated size of the buffer */
    int         version;  /* 1, 2, and 5 for CDF-1, 2, and 5 respectively */
    char       *base;     /* beginning of read/write buffer */
    char       *pos;      /* current position in buffer */
} bufferinfo;

extern int blob_ncmpio_create_NC(NC *ncp);
extern int blob_ncmpio_free_NC(NC *ncp);
extern int blob_ncmpio_add_var(NC *ncp, const char *name, nc_type xtype,
                               int ndims, int *dimids, int *varidp);
extern int blob_ncmpio_add_dim(NC *ncp, const char *name, MPI_Offset size,
                               int *dimidp);
extern int blob_ncmpio_put_att(NC *ncp, int varid, const char *name,
                               nc_type xtype, MPI_Offset nelems,
                               const void *buf);
extern int blob_ncmpio_get_att(NC *ncp, int varid, const char *name, void *buf);
extern int blob_ncmpio_pack_NC(NC *ncp, size_t *buf_len, void **buf);


#if defined(__cplusplus)
}
#endif

#endif /* _BLOB_NCMPIO_H */
