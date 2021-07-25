/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 *********************************************************************/
//
#pragma once
//
#include <map>
#include <vector>
//
#include <hdf5.h>
#include <mpi.h>
//
#include <e3sm_io.h>
#include <e3sm_io_err.h>
#include <e3sm_io_driver.hpp>

#define CHECK_HERR                                                    \
    {                                                                 \
        if (herr != 0) {                                              \
            printf ("Error at line %d in %s:\n", __LINE__, __FILE__); \
            H5Eprint1 (stdout);                                       \
            DEBUG_ABORT;                                              \
            err = -1;                                                 \
            goto err_out;                                             \
        }                                                             \
    }

#define CHECK_HID(A)                                                  \
    {                                                                 \
        if (A < 0) {                                                  \
            printf ("Error at line %d in %s:\n", __LINE__, __FILE__); \
            H5Eprint1 (stdout);                                       \
            DEBUG_ABORT;                                              \
            err = -1;                                                 \
            goto err_out;                                             \
        }                                                             \
    }

static inline hid_t nc_type_to_hdf5_type (nc_type nctype) {
    switch (nctype) {
        case NC_INT:
            return H5T_NATIVE_INT;
        case NC_FLOAT:
            return H5T_NATIVE_FLOAT;
        case NC_DOUBLE:
            return H5T_NATIVE_DOUBLE;
        case NC_CHAR:
            return H5T_NATIVE_CHAR;
        default:
            printf ("Error at line %d in %s: Unknown type %d\n", __LINE__, __FILE__, nctype);
            DEBUG_ABORT
    }

    return -1;
}

static inline hid_t mpi_type_to_hdf5_type (MPI_Datatype type) {
    if (type == MPI_DOUBLE){
        return H5T_NATIVE_DOUBLE;
    }
    else if (type == MPI_FLOAT){
        return H5T_NATIVE_FLOAT;
    }
    else if (type == MPI_INT){
        return H5T_NATIVE_INT;
    }
    else if (type == MPI_LONG_LONG){
        return H5T_NATIVE_INT64;
    }
    else if (type == MPI_CHAR){
        return H5T_NATIVE_CHAR;
    }
    else if (type == MPI_BYTE){
        return H5T_NATIVE_UINT8;
    }

    return -1;
}
