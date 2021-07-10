#pragma once

#include <mpi.h>

#ifdef E3SM_IO_DEBUG
#include <stdlib.h>
#include <string.h>
#define DEBUG_ABORT                                        \
    {                                                      \
        char *val = getenv ("E3SM_IO_DEBUG_ABORT_ON_ERR"); \
        if (val && (strcmp (val, "1") == 0)) { abort (); } \
    }
#else
#define DEBUG_ABORT
#endif

#define CHECK_ERR {                                                      \
    if (err < 0) {                                                       \
        printf ("Error in %s line %d function %s\n", __FILE__, __LINE__, \
                __func__);                                               \
        DEBUG_ABORT                                                      \
        goto err_out;                                                    \
    }                                                                    \
}

#define CHECK_MPIERR {                                                       \
    if (err != MPI_SUCCESS) {                                                \
        int el = 256;                                                        \
        char errstr[256];                                                    \
        MPI_Error_string (err, errstr, &el);                                 \
        printf ("Error in %s line %d function %s: %s\n", __FILE__, __LINE__, \
                __func__, errstr);                                           \
        err = -1;                                                            \
        DEBUG_ABORT                                                          \
        goto err_out;                                                        \
    }                                                                        \
}

#define CHECK_PTR(A)                                                              \
    {                                                                             \
        if (A == NULL) {                                                          \
            printf ("Error in %s line %d function %s\n", __FILE__, __LINE__,  \
                __func__);                                                \
            err = -1;                                                             \
            DEBUG_ABORT                                                           \
            goto err_out;                                                         \
        }                                                                         \
    }

#define ERR_OUT(msg) {                                                  \
    printf("Error in %s line %d function %s: %s\n", __FILE__, __LINE__, \
           __func__, msg);                                              \
    err = -1;                                                           \
    DEBUG_ABORT                                                         \
    goto err_out;                                                       \
}

