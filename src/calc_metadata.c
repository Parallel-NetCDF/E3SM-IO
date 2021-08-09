/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program splits the MPI communicator of I/O processes into multiple
 * sub-communicators, such that there is one sub-communicator per compute node
 * and all MPI processes on the same compute node belong to the node's
 * sub-communicator.
 *
 * Based on the scope of sub-communicator, this program also calculate the
 * metadata aggregated across all processes in the sub-communicator, such as
 * the number of requests, number of offset-length pairs, and each process's
 * starting offset array indices for all decompositions.
 *
 *********************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <e3sm_io.h>

#define CHECK_MPI_ERROR(mpi_errorcode, err_msg) {                      \
    if (mpi_errorcode != MPI_SUCCESS) {                                \
        char errorString[MPI_MAX_ERROR_STRING];                        \
        int errorStringLen;                                            \
        MPI_Error_string(mpi_errorcode, errorString, &errorStringLen); \
        printf("%2d: MPI Failure at line %d of %s (%s : %s)\n",        \
               rank, __LINE__, __FILE__, err_msg, errorString);        \
        return -1;                                                     \
    }                                                                  \
}

/*----< split_communicator() >-----------------------------------------------*/
static
int split_communicator(MPI_Comm comm,
                       MPI_Comm *sub_comm)
{
    int err, rank;
    MPI_Comm_rank(comm, &rank);

    /* split communicator to create one sub-communicator per compute node and
     * all processes on the same compute node are on the same sub-communicator.
     */
    err = MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL,
                              sub_comm);
    CHECK_MPI_ERROR(err, "MPI_Comm_split_type")

    return 0;
}

/*----< blob_metadata() >----------------------------------------------------*/
static
int blob_metadata(e3sm_io_config *cfg,
                  e3sm_io_decom  *decom)
{
    int i, j, err, rank, color, int_msg[2], *nnprocs;
    MPI_Comm comm_roots;

    /* split communicator to create one sub-communicator per compute node */
    err = split_communicator(cfg->io_comm, &cfg->sub_comm);
    if (err < 0) return err;

    /* calculate subfile ID and take care of both process rank assignments:
     * block-based (MPICH_RANK_REORDER_METHOD=1) or
     * round-robin (MPICH_RANK_REORDER_METHOD=0)
     */
    rank = cfg->rank;
    MPI_Comm_rank(cfg->sub_comm, &cfg->sub_rank);
    color = (cfg->sub_rank == 0) ? 1 : 0;
    err = MPI_Comm_split(cfg->io_comm, color, cfg->rank, &comm_roots);
    CHECK_MPI_ERROR(err, "MPI_Comm_split")

    MPI_Comm_size(comm_roots, &cfg->num_subfiles);
    MPI_Comm_rank(comm_roots, &cfg->subfile_ID);

    /* print the number of MPI processes per node (subfile) */
    MPI_Comm_size(cfg->sub_comm, &cfg->sub_nprocs);
    nnprocs = (int*) malloc(cfg->num_subfiles * sizeof(int));
    MPI_Gather(&cfg->sub_nprocs, 1, MPI_INT, nnprocs, 1, MPI_INT, 0,
               comm_roots);
    if (rank == 0) {
        char str[64], *msg=cfg->node_info;
        sprintf(msg,"Total number of compute nodes: %d\n", cfg->num_subfiles);
        for (j=0, i=1; i<cfg->num_subfiles; i++) {
            if (nnprocs[i] != nnprocs[j]) {
                if (j+1 == i) sprintf(str,"Node %d runs %d processes; ",j,nnprocs[j]);
                else sprintf(str,"Nodes %d to %d run %d processes each; ",j,i-1,nnprocs[j]);
                j = i;
                strcat(msg, str);
            }
        }
        if (j+1 == i) {
            if (nnprocs[j] == 1)
                sprintf(str,"Node %d runs %d process.",j,nnprocs[j]);
            else
                sprintf(str,"Node %d runs %d processes.",j,nnprocs[j]);
        }
        else sprintf(str,"Nodes %d to %d run %d processes each.",j,i-1,nnprocs[j]);
        strcat(msg, str);
    }
    free(nnprocs);
    MPI_Comm_free(&comm_roots);

    int_msg[0] = cfg->num_subfiles;
    int_msg[1] = cfg->subfile_ID;
    err = MPI_Bcast(int_msg, 2, MPI_INT, 0, cfg->sub_comm);
    CHECK_MPI_ERROR(err, "MPI_Bcast")
    cfg->num_subfiles = int_msg[0];
    cfg->subfile_ID   = int_msg[1];
    if (cfg->verbose && cfg->sub_rank == 0)
        printf("cfg->rank=%5d sub_rank=%5d color=%d subfile_ID=%5d\n",
           cfg->rank, cfg->sub_rank, color, cfg->subfile_ID);

    for (i=0; i<decom->num_decomp; i++) {
        decom->start[i] = 0;
        decom->count[i] = 0;
        /* total request amount per decomposition by this process */
        for (j=0; j<decom->contig_nreqs[i]; j++)
            decom->count[i] += decom->blocklens[i][j];

        /* calculate max no. requests per decomposition by this process */
        decom->max_nreqs[i] = decom->contig_nreqs[i];
    }
    if (cfg->verbose && cfg->rank == 0)
        printf("cfg->rank %d max_nreqs= %d %d %d\n",cfg->rank,
               decom->max_nreqs[0],decom->max_nreqs[1],decom->max_nreqs[2]);

    /* calculate starting offset per decomposition of this process */
    err = MPI_Exscan(decom->count, decom->start, decom->num_decomp, MPI_OFFSET,
                     MPI_SUM, cfg->sub_comm);
    CHECK_MPI_ERROR(err," MPI_Exscan")

    /* decom->nelems is the total number of array elements in this subfile */
    err = MPI_Allreduce(decom->count, decom->nelems, decom->num_decomp,
                        MPI_OFFSET, MPI_SUM, cfg->sub_comm);
    CHECK_MPI_ERROR(err, "MPI_Allreduce")

    /* decom->max_nreqs is the max number of requests in this subfile */
    err = MPI_Allreduce(MPI_IN_PLACE, decom->max_nreqs, decom->num_decomp,
                        MPI_INT, MPI_MAX, cfg->sub_comm);
    CHECK_MPI_ERROR(err, "MPI_Allreduce")

    if (cfg->verbose && cfg->rank == 0)
        printf("cfg->rank %d nelems= %lld %lld %lld start=%lld %lld %lld count=%lld %lld %lld\n",cfg->rank,
               decom->nelems[0],decom->nelems[1],decom->nelems[2],
               decom->start[0],decom->start[1],decom->start[2],
               decom->count[0],decom->count[1],decom->count[2]);

    return 0;
}

/*----< set_starts_counts() >------------------------------------------------*/
static
int set_starts_counts(e3sm_io_decom *dp)
{
    int i, j;

    for (i=0; i<dp->num_decomp; i++) {
        int nreqs = dp->contig_nreqs[i];

        /* construct starts[] and counts[] for iput_varn */
        dp->w_starts[i] = (MPI_Offset**) malloc(nreqs * 4 *
                                                sizeof(MPI_Offset*));
        dp->w_counts[i] = dp->w_starts[i] + nreqs;
        dp->w_startx[i] = dp->w_counts[i] + nreqs;
        dp->w_countx[i] = dp->w_startx[i] + nreqs;

        dp->w_starts[i][0] = (MPI_Offset*) malloc(nreqs * MAX_NDIMS * 4 *
                                                  sizeof(MPI_Offset));
        dp->w_counts[i][0] = dp->w_starts[i][0] + nreqs * MAX_NDIMS;
        dp->w_startx[i][0] = dp->w_counts[i][0] + nreqs * MAX_NDIMS;
        dp->w_countx[i][0] = dp->w_startx[i][0] + nreqs * MAX_NDIMS;
        for (j=1; j<nreqs; j++) {
            dp->w_starts[i][j] = dp->w_starts[i][j-1] + MAX_NDIMS;
            dp->w_counts[i][j] = dp->w_counts[i][j-1] + MAX_NDIMS;
            dp->w_startx[i][j] = dp->w_startx[i][j-1] + MAX_NDIMS;
            dp->w_countx[i][j] = dp->w_countx[i][j-1] + MAX_NDIMS;
        }

        for (j=0; j<nreqs; j++) {
            dp->w_starts[i][j][0] = 0;
            dp->w_counts[i][j][0] = 1;
            if (dp->ndims[i] == 1) { /* decomposition is 1D */
                dp->w_starts[i][j][1] = dp->disps[i][j];
                dp->w_counts[i][j][1] = dp->blocklens[i][j];
                dp->w_startx[i][j][0] = dp->w_starts[i][j][1];
                dp->w_countx[i][j][0] = dp->w_counts[i][j][1];
            }
            else if (dp->ndims[i] == 2) { /* decomposition is 2D */
                dp->w_starts[i][j][1] = dp->disps[i][j] / dp->dims[i][1];
                dp->w_starts[i][j][2] = dp->disps[i][j] % dp->dims[i][1];
                dp->w_counts[i][j][1] = 1;
                dp->w_counts[i][j][2] = dp->blocklens[i][j];
                dp->w_startx[i][j][0] = dp->w_starts[i][j][1];
                dp->w_startx[i][j][1] = dp->w_starts[i][j][2];
                dp->w_countx[i][j][0] = dp->w_counts[i][j][1];
                dp->w_countx[i][j][1] = dp->w_counts[i][j][2];
            }
            else if (dp->ndims[i] == 3) { /* decomposition is 3D */
                int xy = dp->dims[i][2] * dp->dims[i][1];

                dp->w_starts[i][j][1] = dp->disps[i][j] / xy;
                dp->w_starts[i][j][2] = dp->disps[i][j] % xy / dp->dims[i][2];
                dp->w_starts[i][j][3] = dp->disps[i][j] % dp->dims[i][2];
                dp->w_counts[i][j][1] = 1;
                dp->w_counts[i][j][2] = 1;
                dp->w_counts[i][j][3] = dp->blocklens[i][j];

                dp->w_startx[i][j][0] = dp->w_starts[i][j][1];
                dp->w_startx[i][j][1] = dp->w_starts[i][j][2];
                dp->w_startx[i][j][2] = dp->w_starts[i][j][3];
                dp->w_countx[i][j][0] = dp->w_counts[i][j][1];
                dp->w_countx[i][j][1] = dp->w_counts[i][j][2];
                dp->w_countx[i][j][2] = dp->w_counts[i][j][3];
            }
            /* each blocklens[j] is no bigger than last dims[] */
        }
    }

    return 0;
}

/*----< calc_metadata() >----------------------------------------------------*/
int calc_metadata(e3sm_io_config *cfg,
                  e3sm_io_decom  *decom)
{
    int i, j;

    /* Note adios blob I/O also uses canonical metadata */
    if (cfg->strategy == blob && cfg->api != adios)
        return blob_metadata(cfg, decom);

    /* for canonical order I/O */
    for (i=0; i<decom->num_decomp; i++) {
        decom->count[i] = 0;
        /* total request amount per decomposition by this process */
        for (j=0; j<decom->contig_nreqs[i]; j++)
            decom->count[i] += decom->blocklens[i][j];
    }

    if (cfg->api != adios) set_starts_counts(decom);

    return 0;
}

