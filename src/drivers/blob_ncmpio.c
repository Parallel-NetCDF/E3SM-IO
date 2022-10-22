/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 * This file is extracted from PnetCDF source files:
 * ncmpio_enddef.c, ncmpio_header_get.c, ncmpio_header_put.c, and utils.c
 * https://github.com/Parallel-NetCDF/PnetCDF/tree/master/src/drivers/ncmpio
 *
 *********************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h> /* strdup() */
#include <assert.h>
#include <mpi.h>

#include <e3sm_io.h>
#include <blob_ncmpio.h>

#define NC_MAGIC_LEN 4

#define X_SIZEOF_NC_TYPE X_SIZEOF_INT
#define X_SIZEOF_NC_TAG  X_SIZEOF_INT

/*
 * "magic number" at beginning of file: 0x43444601 (big endian)
 */
static const char ncmagic5[] = {'C', 'D', 'F', 0x05};

static const char nada[X_ALIGN] = {0, 0, 0, 0};

/*----< hdr_len_NC_dim() >---------------------------------------------------*/
static MPI_Offset
hdr_len_NC_dim(const NC_dim *dimp, int sizeof_NON_NEG)
{
    /* netCDF file format:
     *  ...
     * dim        = name dim_length
     * dim_length = NON_NEG
     * NON_NEG    = <non-negative INT> |  // CDF-1 and CDF-2
     *              <non-negative INT64>  // CDF-5
     */
    MPI_Offset sz;

    assert(dimp != NULL);

    sz  = sizeof_NON_NEG + _RNDUP(dimp->name_len, X_ALIGN); /* name */
    sz += sizeof_NON_NEG;                                   /* dim_length */

    return sz;
}

/*----< hdr_len_NC_dimarray() >----------------------------------------------*/
static MPI_Offset
hdr_len_NC_dimarray(const NC_dimarray *ncap, int sizeof_NON_NEG)
{
    /* netCDF file format:
     *  ...
     * dim_list     = ABSENT | NC_DIMENSION  nelems  [dim ...]
     * ABSENT       = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *                ZERO  ZERO64  // for CDF-5
     * ZERO         = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64       = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_DIMENSION = \x00 \x00 \x00 \x0A         // tag for list of dimensions
     * nelems       = NON_NEG       // number of elements in following sequence
     * NON_NEG      = <non-negative INT> |        // CDF-1 and CDF-2
     *                <non-negative INT64>        // CDF-5
     */
    int i;
    MPI_Offset xlen;

    xlen = X_SIZEOF_NC_TAG;           /* NC_DIMENSION */
    xlen += sizeof_NON_NEG;           /* nelems */

    if (ncap == NULL) /* ABSENT: no dimension is defined */
        return xlen;

    /* [dim ...] */
    for (i=0; i<ncap->ndefined; i++)
        xlen += hdr_len_NC_dim(ncap->value[i], sizeof_NON_NEG);

    return xlen;
}

/*----< hdr_len_NC_attr() >--------------------------------------------------*/
static MPI_Offset
hdr_len_NC_attr(const NC_attr *attrp, int sizeof_NON_NEG)
{
    /* netCDF file format:
     *  ...
     * attr    = name  nc_type  nelems  [values ...]
     * nc_type = NC_BYTE | NC_CHAR | NC_SHORT | ...
     * nelems  = NON_NEG       // number of elements in following sequence
     * values  = bytes | chars | shorts | ints | floats | doubles
     * bytes   = [BYTE ...]  padding
     * chars   = [CHAR ...]  padding
     * shorts  = [SHORT ...]  padding
     * ints    = [INT ...]
     * floats  = [FLOAT ...]
     * doubles = [DOUBLE ...]
     * padding = <0, 1, 2, or 3 bytes to next 4-byte boundary>
     * NON_NEG = <non-negative INT> |  // CDF-1 and CDF-2
     *           <non-negative INT64>  // CDF-5
     */
    MPI_Offset sz;

    assert(attrp != NULL);

    sz  = sizeof_NON_NEG + _RNDUP(attrp->name_len, X_ALIGN); /* name */
    sz += X_SIZEOF_NC_TYPE;                                  /* nc_type */
    sz += sizeof_NON_NEG;                                    /* nelems */
    sz += attrp->xsz;                                        /* [values ...] */

    return sz;
}

/*----< hdr_len_NC_attrarray() >---------------------------------------------*/
static MPI_Offset
hdr_len_NC_attrarray(const NC_attrarray *ncap, int sizeof_NON_NEG)
{
    /* netCDF file format:
     *  ...
     * att_list     = ABSENT | NC_ATTRIBUTE  nelems  [attr ...]
     * ABSENT       = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *                ZERO  ZERO64  // for CDF-5
     * ZERO         = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64       = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_ATTRIBUTE = \x00 \x00 \x00 \x0C         // tag for list of attributes
     * nelems       = NON_NEG       // number of elements in following sequence
     * NON_NEG      = <non-negative INT> |        // CDF-1 and CDF-2
     *                <non-negative INT64>        // CDF-5
     */
    int i;
    MPI_Offset xlen;

    xlen = X_SIZEOF_NC_TAG;        /* NC_ATTRIBUTE */
    xlen += sizeof_NON_NEG;        /* nelems */

    if (ncap == NULL) /* ABSENT: no attribute is defined */
        return xlen;

    for (i=0; i<ncap->ndefined; i++) /* [attr ...] */
        xlen += hdr_len_NC_attr(ncap->value[i], sizeof_NON_NEG);

    return xlen;
}

/*----< hdr_len_NC_var() >---------------------------------------------------*/
static MPI_Offset
hdr_len_NC_var(const NC_var *varp,
               int           sizeof_off_t,    /* OFFSET */
               int           sizeof_NON_NEG)  /* NON_NEG */
{
    /* netCDF file format:
     * netcdf_file = header data
     * header      = magic numrecs dim_list gatt_list var_list
     *  ...
     * var         = name nelems [dimid ...] vatt_list nc_type vsize begin
     * nelems      = NON_NEG
     * dimid       = NON_NEG
     * vatt_list   = att_list
     * nc_type     = NC_BYTE | NC_CHAR | NC_SHORT | ...
     * vsize       = NON_NEG
     * begin       = OFFSET        // Variable start location.
     * OFFSET      = <non-negative INT> |  // CDF-1
     *               <non-negative INT64>  // CDF-2 and CDF-5
     * NON_NEG     = <non-negative INT> |  // CDF-1 and CDF-2
     *               <non-negative INT64>  // CDF-5
     */
    MPI_Offset sz;

    assert(varp != NULL);

    /* for CDF-1, sizeof_off_t == 4 && sizeof_NON_NEG == 4
     * for CDF-2, sizeof_off_t == 8 && sizeof_NON_NEG == 4
     * for CDF-5, sizeof_off_t == 8 && sizeof_NON_NEG == 8
     */
    sz  = sizeof_NON_NEG + _RNDUP(varp->name_len, X_ALIGN);   /* name */
    sz += sizeof_NON_NEG;                                     /* nelems */
    sz += sizeof_NON_NEG * varp->ndims;                       /* [dimid ...] */
    sz += hdr_len_NC_attrarray(&varp->attrs, sizeof_NON_NEG); /* vatt_list */
    sz += X_SIZEOF_NC_TYPE;                                   /* nc_type */
    sz += sizeof_NON_NEG;                                     /* vsize */
    sz += sizeof_off_t;                                       /* begin */

    return sz;
}

/*----< hdr_len_NC_vararray() >----------------------------------------------*/
static MPI_Offset
hdr_len_NC_vararray(const NC_vararray *ncap,
                    int                sizeof_NON_NEG, /* NON_NEG */
                    int                sizeof_off_t)   /* OFFSET */
{
    /* netCDF file format:
     * netcdf_file = header  data
     * header      = magic  numrecs  dim_list  gatt_list  var_list
     *  ...
     * var_list    = ABSENT | NC_VARIABLE   nelems  [var ...]
     * ABSENT      = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *               ZERO  ZERO64  // for CDF-5
     * ZERO        = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64      = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_VARIABLE = \x00 \x00 \x00 \x0B         // tag for list of variables
     * nelems      = NON_NEG       // number of elements in following sequence
     * NON_NEG     = <non-negative INT> |        // CDF-1 and CDF-2
     *               <non-negative INT64>        // CDF-5
     */
    int i;
    MPI_Offset xlen;

    xlen = X_SIZEOF_NC_TAG;           /* NC_VARIABLE */
    xlen += sizeof_NON_NEG;           /* nelems */

    if (ncap == NULL) /* ABSENT: no variable is defined */
        return xlen;

    /* for CDF-1, sizeof_off_t == 4 && sizeof_NON_NEG == 4
     * for CDF-2, sizeof_off_t == 8 && sizeof_NON_NEG == 4
     * for CDF-5, sizeof_off_t == 8 && sizeof_NON_NEG == 8
     */
    for (i=0; i<ncap->ndefined; i++)  /* [var ...] */
        xlen += hdr_len_NC_var(ncap->value[i], sizeof_off_t, sizeof_NON_NEG);

    return xlen;
}

static MPI_Offset
ncmpio_hdr_len_NC(const NC *ncp)
{
    /* netCDF file format:
     * netcdf_file = header  data
     * header      = magic  numrecs  dim_list  gatt_list  var_list
     *  ...
     * numrecs     = NON_NEG | STREAMING   // length of record dimension
     * NON_NEG     = <non-negative INT> |  // CDF-1 and CDF-2
     *               <non-negative INT64>  // CDF-5
     */

    int sizeof_NON_NEG, sizeof_off_t;
    MPI_Offset xlen;

    sizeof_NON_NEG = X_SIZEOF_INT64; /* 8-byte integer for all integers */
    sizeof_off_t   = X_SIZEOF_INT64; /* 8-byte integer for var begin */

    xlen  = NC_MAGIC_LEN;                                                    /* magic */
    xlen += sizeof_NON_NEG;                                                  /* numrecs */
    xlen += hdr_len_NC_dimarray(&ncp->dims,   sizeof_NON_NEG);               /* dim_list */
    xlen += hdr_len_NC_attrarray(&ncp->attrs, sizeof_NON_NEG);               /* gatt_list */
    xlen += hdr_len_NC_vararray(&ncp->vars,   sizeof_NON_NEG, sizeof_off_t); /* var_list */

    return xlen; /* return the header size (not yet aligned) */
}

int blob_ncmpio_create_NC(NC *ncp)
{
    /* calculate the true header size (not-yet aligned) */
    ncp->xsz = ncmpio_hdr_len_NC(ncp);

    /* initialize unlimited_id as no unlimited dimension yet defined */
    ncp->dims.unlimited_id = -1;

    return 0;
}

static int
ncmpix_putn_text(void **xpp, MPI_Offset nelems, const char *tp)
{
    (void) memcpy(*xpp, tp, (size_t)nelems);
    *xpp = (void *)((char *)(*xpp) + nelems);
    return NC_NOERR;
}

static int
ncmpix_put_uint64(void **xpp, unsigned long long ip) {
    (void) memcpy(*xpp, &ip, X_SIZEOF_UINT64);
    *xpp  = (void *)((char *)(*xpp) + 8);
    return NC_NOERR;
}

static int
ncmpix_put_uint32(void **xpp, unsigned int ip) {
    (void) memcpy(*xpp, &ip, X_SIZEOF_UINT);
    *xpp  = (void *)((char *)(*xpp) + 4);
    return NC_NOERR;
}

static int
ncmpix_pad_putn_text(void **xpp, MPI_Offset nelems, const char *tp)
{
    MPI_Offset rndup = nelems % X_ALIGN;

    if (rndup) rndup = X_ALIGN - rndup;

    (void) memcpy(*xpp, tp, (size_t)nelems);
    *xpp = (void *)((char *)(*xpp) + nelems);

    if (rndup) {
        (void) memcpy(*xpp, nada, (size_t)rndup);
        *xpp = (void *)((char *)(*xpp) + rndup);
    }
    return NC_NOERR;
}

/*----< hdr_put_NC_name() >--------------------------------------------------*/
static int
hdr_put_NC_name(bufferinfo *pbp,
                const char *name)
{
    /* netCDF file format:
     *  ...
     * name       = nelems  namestring
     * nelems     = NON_NEG
     * namestring = ID1 [IDN ...] padding
     * ID1        = alphanumeric | '_'
     * IDN        = alphanumeric | special1 | special2
     * padding    = <0, 1, 2, or 3 bytes to next 4-byte boundary>
     * NON_NEG    = <non-negative INT> |  // CDF-1 and CDF-2
     *              <non-negative INT64>  // CDF-5
     */
    int err;
    size_t nchars = strlen(name);

    /* copy nelems */
    err = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)nchars);
    if (err != NC_NOERR) return err;

    /* copy namestring */
    return ncmpix_pad_putn_text((void **)(&pbp->pos), (MPI_Offset)nchars, name);
}

/*----< hdr_put_NC_dim() >---------------------------------------------------*/
static int
hdr_put_NC_dim(bufferinfo   *pbp,
               const NC_dim *dimp)
{
    /* netCDF file format:
     *  ...
     * dim        = name  dim_length
     * dim_length = NON_NEG
     * NON_NEG    = <non-negative INT> |  // CDF-1 and CDF-2
     *              <non-negative INT64>  // CDF-5
     */
    int err;

    /* copy name */
    err = hdr_put_NC_name(pbp, dimp->name);
    if (err != NC_NOERR) return err;

    err = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)dimp->size);

    return err;
}

/*----< hdr_put_NC_dimarray() >----------------------------------------------*/
static int
hdr_put_NC_dimarray(bufferinfo        *pbp,
                    const NC_dimarray *ncap)
{
    /* netCDF file format:
     *  ...
     * dim_list     = ABSENT | NC_DIMENSION  nelems  [dim ...]
     * ABSENT       = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *                ZERO  ZERO64  // for CDF-5
     * ZERO         = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64       = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_DIMENSION = \x00 \x00 \x00 \x0A         // tag for list of dimensions
     * nelems       = NON_NEG       // number of elements in following sequence
     * NON_NEG      = <non-negative INT> |        // CDF-1 and CDF-2
     *                <non-negative INT64>        // CDF-5
     */
    int i, status;

    assert(pbp != NULL);

    if (ncap == NULL || ncap->ndefined == 0) { /* ABSENT */
        status = ncmpix_put_uint32((void**)(&pbp->pos), NC_UNSPECIFIED);
        if (status != NC_NOERR) return status;

        /* put a ZERO or ZERO64 depending on which CDF format */
        status = ncmpix_put_uint64((void**)(&pbp->pos), 0);
        if (status != NC_NOERR) return status;
    }
    else {
        /* copy NC_DIMENSION */
        status = ncmpix_put_uint32((void**)(&pbp->pos), NC_DIMENSION);
        if (status != NC_NOERR) return status;

        /* copy nelems */
        status = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)ncap->ndefined);
        if (status != NC_NOERR) return status;

        /* copy [dim ...] */
        for (i=0; i<ncap->ndefined; i++) {
            status = hdr_put_NC_dim(pbp, ncap->value[i]);
            if (status != NC_NOERR) return status;
        }
    }

    return NC_NOERR;
}

/*----< hdr_put_NC_attrV() >-------------------------------------------------*/
static int
hdr_put_NC_attrV(bufferinfo    *pbp,
                 const NC_attr *attrp)
{
    /* netCDF file format:
     *  ...
     * attr    = name  nc_type  nelems  [values ...]
     *  ...
     * values  = bytes | chars | shorts | ints | floats | doubles
     * bytes   = [BYTE ...]  padding
     * chars   = [CHAR ...]  padding
     * shorts  = [SHORT ...]  padding
     * ints    = [INT ...]
     * floats  = [FLOAT ...]
     * doubles = [DOUBLE ...]
     * padding = <0, 1, 2, or 3 bytes to next 4-byte boundary>
     */
    int xsz;
    MPI_Offset padding, sz;

    /* e3sm_io_xlen_nc_type() returns the element size (unaligned) of
     * attrp->xtype attrp->xsz is the aligned total size of attribute values
     */
    e3sm_io_xlen_nc_type(attrp->xtype, &xsz);
    sz = attrp->nelems * xsz;
    padding = attrp->xsz - sz;

    memcpy(pbp->pos, attrp->xvalue, (size_t)sz);
    pbp->pos = (void *)((char *)pbp->pos + sz);

    if (padding > 0) {
        /* zero-padding is per buffer, not per element */
        memset(pbp->pos, 0, (size_t)padding);
        pbp->pos = (void *)((char *)pbp->pos + padding);
    }

    return NC_NOERR;
}

/*----< hdr_put_NC_attr() >--------------------------------------------------*/
static int
hdr_put_NC_attr(bufferinfo    *pbp,
                const NC_attr *attrp)
{
    /* netCDF file format:
     *  ...
     * attr    = name  nc_type  nelems  [values ...]
     * nc_type = NC_BYTE | NC_CHAR | NC_SHORT | ...
     * nelems  = NON_NEG       // number of elements in following sequence
     * NON_NEG = <non-negative INT> |  // CDF-1 and CDF-2
     *           <non-negative INT64>  // CDF-5
     */
    int status;

    /* copy name */
    status = hdr_put_NC_name(pbp, attrp->name);
    if (status != NC_NOERR) return status;

    /* copy nc_type */
    status = ncmpix_put_uint32((void**)(&pbp->pos), (uint)attrp->xtype);
    if (status != NC_NOERR) return status;

    /* copy nelems */
    status = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)attrp->nelems);
    if (status != NC_NOERR) return status;

    /* copy [values ...] */
    status = hdr_put_NC_attrV(pbp, attrp);
    if (status != NC_NOERR) return status;

    return NC_NOERR;
}

/*----< hdr_put_NC_attrarray() >---------------------------------------------*/
static int
hdr_put_NC_attrarray(bufferinfo         *pbp,
                     const NC_attrarray *ncap)
{
    /* netCDF file format:
     *  ...
     * att_list     = ABSENT | NC_ATTRIBUTE  nelems  [attr ...]
     * ABSENT       = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *                ZERO  ZERO64  // for CDF-5
     * ZERO         = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64       = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_ATTRIBUTE = \x00 \x00 \x00 \x0C         // tag for list of attributes
     * nelems       = NON_NEG       // number of elements in following sequence
     * NON_NEG      = <non-negative INT> |        // CDF-1 and CDF-2
     *                <non-negative INT64>        // CDF-5
     */
    int i, status;

    assert(pbp != NULL);

    if (ncap == NULL || ncap->ndefined == 0) { /* ABSENT */
        status = ncmpix_put_uint32((void**)(&pbp->pos), NC_UNSPECIFIED);
        if (status != NC_NOERR) return status;

        /* put a ZERO or ZERO64 depending on which CDF format */
        status = ncmpix_put_uint64((void**)(&pbp->pos), 0);
        if (status != NC_NOERR) return status;
    }
    else {
        /* copy NC_ATTRIBUTE */
        status = ncmpix_put_uint32((void**)(&pbp->pos), NC_ATTRIBUTE);
        if (status != NC_NOERR) return status;

        /* copy nelems */
        status = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)ncap->ndefined);
        if (status != NC_NOERR) return status;

        /* copy [attr ...] */
        for (i=0; i<ncap->ndefined; i++) {
            status = hdr_put_NC_attr(pbp, ncap->value[i]);
            if (status != NC_NOERR) return status;
        }
    }

    return NC_NOERR;
}

/*----< hdr_put_NC_var() >---------------------------------------------------*/
static int
hdr_put_NC_var(bufferinfo   *pbp,
               const NC_var *varp)
{
    /* netCDF file format:
     * netcdf_file = header data
     * header      = magic numrecs dim_list gatt_list var_list
     *  ...
     * var         = name nelems [dimid ...] vatt_list nc_type vsize begin
     * nelems      = NON_NEG
     * dimid       = NON_NEG
     * vatt_list   = att_list
     * nc_type     = NC_BYTE | NC_CHAR | NC_SHORT | ...
     * vsize       = NON_NEG
     * begin       = OFFSET        // Variable start location.
     * OFFSET      = <non-negative INT> |  // CDF-1
     *               <non-negative INT64>  // CDF-2 and CDF-5
     * NON_NEG     = <non-negative INT> |  // CDF-1 and CDF-2
     *               <non-negative INT64>  // CDF-5
     */
    int i, status;

    /* copy name */
    status = hdr_put_NC_name(pbp, varp->name);
    if (status != NC_NOERR) return status;

    /* copy nelems */
    status = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)varp->ndims);
    if (status != NC_NOERR) return status;

    /* copy [dimid ...] */
    for (i=0; i<varp->ndims; i++) {
        status = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)varp->dimids[i]);
        if (status != NC_NOERR) return status;
    }

    /* copy vatt_list */
    status = hdr_put_NC_attrarray(pbp, &varp->attrs);
    if (status != NC_NOERR) return status;

    /* copy nc_type */
    status = ncmpix_put_uint32((void**)(&pbp->pos), (uint)varp->xtype);
    if (status != NC_NOERR) return status;

    /* copy vsize */
    status = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)varp->len);
    if (status != NC_NOERR) return status;

    /* copy begin */
    status = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)varp->begin);
    if (status != NC_NOERR) return status;

    return NC_NOERR;
}

/*----< hdr_put_NC_vararray() >----------------------------------------------*/
static int
hdr_put_NC_vararray(bufferinfo        *pbp,
                    const NC_vararray *ncap)
{
    /* netCDF file format:
     * netcdf_file = header  data
     * header      = magic  numrecs  dim_list  gatt_list  var_list
     *  ...
     * var_list    = ABSENT | NC_VARIABLE   nelems  [var ...]
     * ABSENT      = ZERO  ZERO |  // list is not present for CDF-1 and 2
     *               ZERO  ZERO64  // for CDF-5
     * ZERO        = \x00 \x00 \x00 \x00                      // 32-bit zero
     * ZERO64      = \x00 \x00 \x00 \x00 \x00 \x00 \x00 \x00  // 64-bit zero
     * NC_VARIABLE = \x00 \x00 \x00 \x0B         // tag for list of variables
     * nelems      = NON_NEG       // number of elements in following sequence
     * NON_NEG     = <non-negative INT> |        // CDF-1 and CDF-2
     *               <non-negative INT64>        // CDF-5
     */
    int i, status;

    assert(pbp != NULL);

    if (ncap == NULL || ncap->ndefined == 0) { /* ABSENT */
        status = ncmpix_put_uint32((void**)(&pbp->pos), NC_UNSPECIFIED);
        if (status != NC_NOERR) return status;

        /* put a ZERO or ZERO64 depending on which CDF format */
        status = ncmpix_put_uint64((void**)(&pbp->pos), 0);
        if (status != NC_NOERR) return status;
    }
    else {
        /* copy NC_VARIABLE */
        status = ncmpix_put_uint32((void**)(&pbp->pos), NC_VARIABLE);
        if (status != NC_NOERR) return status;

        /* copy nelems */
        status = ncmpix_put_uint64((void**)(&pbp->pos), (uint64)ncap->ndefined);
        if (status != NC_NOERR) return status;

        /* copy [var ...] */
        for (i=0; i<ncap->ndefined; i++) {
            status = hdr_put_NC_var(pbp, ncap->value[i]);
            if (status != NC_NOERR) return status;
        }
    }

    return NC_NOERR;
}

/*----< ncmpio_hdr_put_NC() >-------------------------------------------*/
static int
ncmpio_hdr_put_NC(NC *ncp, void *buf)
{
    int status;
    bufferinfo putbuf;
    MPI_Offset nrecs=0;

    putbuf.offset        = 0;
    putbuf.pos           = buf;
    putbuf.base          = buf;
    putbuf.size          = ncp->xsz;

    /* netCDF file format:
     * netcdf_file  = header  data
     * header       = magic  numrecs  dim_list  gatt_list  var_list
     */

    /* copy "magic", 4 characters */
    putbuf.version = 5;
    status = ncmpix_putn_text((void **)(&putbuf.pos), sizeof(ncmagic5), ncmagic5);
    if (status != NC_NOERR) return status;

    /* copy numrecs, number of records */
    nrecs = ncp->numrecs;
    status = ncmpix_put_uint64((void**)(&putbuf.pos), (uint64)nrecs);
    if (status != NC_NOERR) return status;

    /* copy dim_list */
    status = hdr_put_NC_dimarray(&putbuf, &ncp->dims);
    if (status != NC_NOERR) return status;

    /* copy gatt_list */
    status = hdr_put_NC_attrarray(&putbuf, &ncp->attrs);
    if (status != NC_NOERR) return status;

    /* copy var_list */
    status = hdr_put_NC_vararray(&putbuf, &ncp->vars);
    if (status != NC_NOERR) return status;

    return NC_NOERR;
}

int
blob_ncmpio_pack_NC(NC      *ncp,
                    size_t  *bufLen,
                    void   **buf)
{
    int err=NC_NOERR, header_wlen;

    *bufLen = 0;
    *buf = NULL;

    /* get the true header size (not header extent) */
    ncp->xsz = ncmpio_hdr_len_NC(ncp);

    /* Do not write padding area (between ncp->xsz and ncp->begin_var) */
    header_wlen = (int) ncp->xsz;
    *bufLen = _RNDUP(header_wlen, X_ALIGN);
    *buf = NCI_Malloc(*bufLen);

    /* pack the entire local header object to buf */
    err = ncmpio_hdr_put_NC(ncp, *buf);
    if (err != NC_NOERR) { /* a fatal error */
        NCI_Free(*buf);
        *buf = NULL;
    }
    return err;
}

static void
ncmpio_free_NC_dimarray(NC_dimarray *ncap)
{
    int i;
    if (ncap->ndefined == 0) return;

    if (ncap->value != NULL) {
        for (i=0; i<ncap->ndefined; i++) {
            if (ncap->value[i] == NULL) break;
            NCI_Free(ncap->value[i]->name);
            NCI_Free(ncap->value[i]);
        }
        NCI_Free(ncap->value);
        ncap->value = NULL;
    }
    ncap->ndefined = 0;
}

static void
ncmpio_free_NC_attrarray(NC_attrarray *ncap)
{
    int i;

    if (ncap->value != NULL) {
        for (i=0; i<ncap->ndefined; i++) {
            if (ncap->value[i] == NULL) continue;
            if (ncap->value[i]->xvalue != NULL)
                NCI_Free(ncap->value[i]->xvalue);
            NCI_Free(ncap->value[i]->name);
            NCI_Free(ncap->value[i]);
        }
        NCI_Free(ncap->value);
        ncap->value = NULL;
    }
    ncap->ndefined = 0;
}

static void
ncmpio_free_NC_var(NC_var *varp)
{
    if (varp == NULL) return;

    ncmpio_free_NC_attrarray(&varp->attrs);
    NCI_Free(varp->name);
    if (varp->shape  != NULL) NCI_Free(varp->shape);
    if (varp->dsizes != NULL) NCI_Free(varp->dsizes);
    if (varp->dimids != NULL) NCI_Free(varp->dimids);

    NCI_Free(varp);
}

static void
ncmpio_free_NC_vararray(NC_vararray *ncap)
{
    int i;

    if (ncap->ndefined == 0) return;

    if (ncap->value != NULL) {
        for (i=0; i<ncap->ndefined; i++) {
            if (ncap->value[i] != NULL)
                ncmpio_free_NC_var(ncap->value[i]);
        }
        NCI_Free(ncap->value);
        ncap->value    = NULL;
    }
    ncap->ndefined = 0;
}

int blob_ncmpio_free_NC(NC *ncp) {

    /* free up space occupied by the header metadata */
    ncmpio_free_NC_dimarray(&ncp->dims);
    ncmpio_free_NC_attrarray(&ncp->attrs);
    ncmpio_free_NC_vararray(&ncp->vars);

    return 0;
}

static int
ncmpio_NC_var_shape64(NC_var            *varp,
                      const NC_dimarray *dims)
{
    int i;
    MPI_Offset product = 1;

    if (varp->ndims == 0) goto out;

    for (i=0; i<varp->ndims; i++)
        varp->shape[i] = dims->value[varp->dimids[i]]->size;

    product = 1;
    if (varp->ndims == 1) {
        if (varp->shape[0] == NC_UNLIMITED)
            varp->dsizes[0] = 1;
        else {
            varp->dsizes[0] = varp->shape[0];
            product = varp->shape[0];
        }
    }
    else { /* varp->ndims > 1 */
        varp->dsizes[varp->ndims-1] = varp->shape[varp->ndims-1];
        product = varp->shape[varp->ndims-1];
        for (i=varp->ndims-2; i>=0; i--) {
            if (varp->shape[i] != NC_UNLIMITED)
                product *= varp->shape[i];
            varp->dsizes[i] = product;
        }
    }

out:
    varp->len = product * varp->xsz;
    if (varp->len % 4 > 0)
        varp->len += 4 - varp->len % 4; /* round up */

    return NC_NOERR;
}

static NC_var*
ncmpio_new_NC_var(const char *name, int ndims)
{
    NC_var *varp = (NC_var *) NCI_Calloc(1, sizeof(NC_var));
    if (varp == NULL) return NULL;

    if (ndims > 0) {
        varp->shape  = (MPI_Offset*)NCI_Calloc(ndims, sizeof(MPI_Offset));
        varp->dsizes = (MPI_Offset*)NCI_Calloc(ndims, sizeof(MPI_Offset));
        varp->dimids = (int *)      NCI_Calloc(ndims, sizeof(int));
    }

    varp->name     = strdup(name);
    varp->name_len = strlen(name); /* name has been NULL checked */
    varp->ndims    = ndims;

    return varp;
}

int blob_ncmpio_add_var(NC         *ncp,
                        const char *name,
                        nc_type     xtype,
                        int         ndims,
                        int        *dimids,
                        int        *varidp)
{
    /* add a variable object in NC */
    int err=NC_NOERR;
    NC_var *varp=NULL;
    varp = ncmpio_new_NC_var(name, ndims);
    varp->xtype = xtype;
    e3sm_io_xlen_nc_type(xtype, &varp->xsz);

    /* copy dimids[] */
    if (ndims != 0 && dimids != NULL)
        memcpy(varp->dimids, dimids, (size_t)ndims * sizeof(int));

    /* set up array dimensional structures */
    ncmpio_NC_var_shape64(varp, &ncp->dims);

    /* allocate/expand ncp->vars.value array */
    if (ncp->vars.ndefined % NC_ARRAY_GROWBY == 0)
        ncp->vars.value = (NC_var **) NCI_Realloc(ncp->vars.value,
                          ((size_t)ncp->vars.ndefined + NC_ARRAY_GROWBY) *
                          sizeof(NC_var*));

    varp->varid = ncp->vars.ndefined; /* varid */

    /* Add a new handle to the end of an array of handles */
    ncp->vars.value[ncp->vars.ndefined] = varp;

    ncp->vars.ndefined++;

    if (varidp != NULL) *varidp = varp->varid;

    return err;
}

int blob_ncmpio_add_dim(NC         *ncp,
                        const char *name,
                        MPI_Offset  size,
                        int        *dimidp)
{
    int err=NC_NOERR, dimid;
    NC_dim *dimp=NULL;

    /* create a new dimension object (dimp->name points to nname) */
    dimp = (NC_dim*) NCI_Malloc(sizeof(NC_dim));

    dimp->size     = size;
    dimp->name     = strdup(name);
    dimp->name_len = strlen(dimp->name);

    /* allocate/expand ncp->dims.value array */
    if (ncp->dims.ndefined % NC_ARRAY_GROWBY == 0)
        ncp->dims.value = (NC_dim **) NCI_Realloc(ncp->dims.value,
                          ((size_t)ncp->dims.ndefined + NC_ARRAY_GROWBY) *
                          sizeof(NC_dim*));

    dimid = ncp->dims.ndefined;

    /* Add a new dim handle to the end of handle array */
    ncp->dims.value[dimid] = dimp;

    if (size == NC_UNLIMITED) ncp->dims.unlimited_id = dimid;

    ncp->dims.ndefined++;

    if (dimidp != NULL) *dimidp = dimid;

    return err;
}

/*----< x_len_NC_attrV() >---------------------------------------------------*/
/* How much space will 'nelems' of 'xtype' take in external representation.
 * Note the space is aligned in 4-byte boundary.
 */
static MPI_Offset
x_len_NC_attrV(nc_type    xtype,
               MPI_Offset nelems)
{
    switch(xtype) {
        case NC_BYTE:
        case NC_CHAR:
        case NC_UBYTE:  return _RNDUP(nelems, 4);
        case NC_SHORT:
        case NC_USHORT: return ((nelems + nelems%2) * 2);
        case NC_INT:
        case NC_UINT:
        case NC_FLOAT:  return (nelems * 4);
        case NC_DOUBLE:
        case NC_INT64:
        case NC_UINT64: return (nelems * 8);
        default: return -1;
    }
    return 0;
}

static int
ncmpio_new_NC_attr(const char  *name,
                   nc_type      xtype,
                   MPI_Offset   nelems,
                   NC_attr    **attrp)
{
    *attrp = (NC_attr*) NCI_Malloc(sizeof(NC_attr));

    (*attrp)->xtype    = xtype;
    (*attrp)->xsz      = 0;
    (*attrp)->nelems   = nelems;
    (*attrp)->xvalue   = NULL;
    (*attrp)->name     = strdup(name);
    (*attrp)->name_len = strlen(name);

    if (nelems > 0) {
        /* obtain 4-byte aligned size of space to store the values */
        MPI_Offset xsz = x_len_NC_attrV(xtype, nelems);
        (*attrp)->xsz    = xsz;
        (*attrp)->xvalue = NCI_Malloc((size_t)xsz);
    }
    return NC_NOERR;
}

static int
incr_NC_attrarray(NC_attrarray *ncap, NC_attr *new_attr)
{
    if (ncap->ndefined % NC_ARRAY_GROWBY == 0)
        ncap->value = (NC_attr **) NCI_Realloc(ncap->value,
                      ((size_t)ncap->ndefined + NC_ARRAY_GROWBY) *
                      sizeof(NC_attr*));

    ncap->value[ncap->ndefined++] = new_attr;

    return NC_NOERR;
}

static int
ncmpio_NC_findattr(const NC_attrarray *ncap,
                   const char         *name) /* normalized string */
{
    int i; 
    size_t nchars;

    if (ncap->ndefined == 0) return -1; /* none created yet */

    nchars = strlen(name);
    for (i=0; i<ncap->ndefined; i++) {
        if (ncap->value[i]->name_len == nchars && 
            strcmp(ncap->value[i]->name, name) == 0)
            return i;
    }
    return -1; /* the name has never been used */
}

int blob_ncmpio_put_att(NC         *ncp,
                        int         varid,
                        const char *name,
                        nc_type     xtype,
                        MPI_Offset  nelems,
                        const void *buf)
{
    /* E3SM always uses the same itype and xtype for attributes */

    int indx=0, err=NC_NOERR;
    NC_attrarray *ncap=NULL;
    NC_attr *attrp=NULL;

    /* obtain NC_attrarray object pointer, ncap */
    if (varid == NC_GLOBAL) ncap = &ncp->attrs;
    else                    ncap = &ncp->vars.value[varid]->attrs;

    /* check whether attribute already exists */
    indx = ncmpio_NC_findattr(ncap, name);
    if (indx >= 0) {
        if (varid == NC_GLOBAL)
            printf("NC_GLOBAL attribute name %s in used\n", name);
        else
            printf("var %s attribute name %s in used\n",
                   ncp->vars.value[varid]->name, name);
        return -1;
    }

    err = ncmpio_new_NC_attr(name, xtype, nelems, &attrp);
    if (err != NC_NOERR) return err;

    err = incr_NC_attrarray(ncap, attrp);
    if (err != NC_NOERR) return err;

    if (nelems == 0) return NC_NOERR; /* non-zero length attribute */

    if (buf == NULL) return NC_NOERR;

    void *xp = attrp->xvalue;
    if (xtype == NC_CHAR || xtype == NC_BYTE || xtype == NC_UBYTE)
        memcpy(xp, buf, nelems);
    else if (xtype == NC_SHORT || xtype == NC_USHORT)
        memcpy(xp, buf, nelems*sizeof(short));
    else if (xtype == NC_INT || xtype == NC_UINT)
        memcpy(xp, buf, nelems*sizeof(int));
    else if (xtype == NC_FLOAT)
        memcpy(xp, buf, nelems*sizeof(float));
    else if (xtype == NC_DOUBLE)
        memcpy(xp, buf, nelems*sizeof(double));
    else if (xtype == NC_INT64 || xtype == NC_UINT64)
        memcpy(xp, buf, nelems*sizeof(long long));
    else err = -1;

    return err;
}

int blob_ncmpio_get_att(NC         *ncp,
                        int         varid,
                        const char *name,
                        void       *buf)
{
    int indx=0, err=NC_NOERR;
    NC_attrarray *ncap=NULL;
    NC_attr *attrp=NULL;
    nc_type xtype;

    if (buf == NULL) return -1;

    if (varid == NC_GLOBAL) ncap = &ncp->attrs;
    else                    ncap = &ncp->vars.value[varid]->attrs;

    indx = ncmpio_NC_findattr(ncap, name);
    if (indx == -1) return -1;

    attrp = ncap->value[indx];

    if (attrp->nelems == 0) return NC_NOERR;

    xtype = attrp->xtype;

    if (xtype == NC_CHAR || xtype == NC_BYTE || xtype == NC_UBYTE)
        memcpy(buf, attrp->xvalue, attrp->nelems);
    else if (xtype == NC_SHORT || xtype == NC_USHORT)
        memcpy(buf, attrp->xvalue, attrp->nelems*sizeof(short));
    else if (xtype == NC_INT || xtype == NC_UINT)
        memcpy(buf, attrp->xvalue, attrp->nelems*sizeof(int));
    else if (xtype == NC_FLOAT)
        memcpy(buf, attrp->xvalue, attrp->nelems*sizeof(float));
    else if (xtype == NC_DOUBLE)
        memcpy(buf, attrp->xvalue, attrp->nelems*sizeof(double));
    else if (xtype == NC_INT64 || xtype == NC_UINT64)
        memcpy(buf, attrp->xvalue, attrp->nelems*sizeof(long long));
    else err = -1;

    return err;
}

