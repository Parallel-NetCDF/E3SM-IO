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

#define CHECK_NERR                                                             \
    {                                                                         \
        if (nerrs != 0) {                                                        \
            printf ("Stop at line %d in %s due to %d errors\n", __LINE__, __FILE__, nerrs); \
            DEBUG_ABORT                                                       \
            goto err_out;                                                     \
        }                                                                     \
    }

#define CHECK_ERR                                                             \
    {                                                                         \
        if (err != 0) {                                                        \
            printf ("Error at line %d in %s: %d\n", __LINE__, __FILE__, err); \
            nerrs++;                                                          \
            DEBUG_ABORT                                                       \
            goto err_out;                                                     \
        }                                                                     \
    }

#define CHECK_MPIERR                                                             \
    {                                                                            \
        if (err != MPI_SUCCESS) {                                             \
            int el = 256;                                                        \
            char errstr[256];                                                    \
            MPI_Error_string (err, errstr, &el);                              \
            printf ("Error at line %d in %s: %s\n", __LINE__, __FILE__, errstr); \
            nerrs++;                                                             \
            DEBUG_ABORT                                                          \
            goto err_out;                                                        \
        }                                                                        \
    }

#define CHECK_PTR(A)                                                              \
    {                                                                             \
        if (A == NULL) {                                                          \
            printf ("Error at line %d in %s: " #A "=NULL\n", __LINE__, __FILE__); \
            nerrs++;                                                              \
            DEBUG_ABORT                                                           \
            goto err_out;                                                         \
        }                                                                         \
    }

#define ERR_OUT(A)                                                      \
    {                                                                   \
        printf ("Error at line %d in %s: %s\n", __LINE__, __FILE__, A); \
        DEBUG_ABORT                                                     \
        goto err_out;                                                   \
    }

#define RET_ERR(A)  \
    {               \
        nerrs++;    \
        ERR_OUT (A) \
    }
