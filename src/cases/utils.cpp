/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 *********************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <mpi.h>

#include <e3sm_io.h>
#include <e3sm_io_case.hpp>
#include <e3sm_io_driver.hpp>

/*----< check_malloc() >-----------------------------------------------------*/
int check_malloc(e3sm_io_config *cfg,
                 e3sm_io_driver *driver)
{
    int err=0, global_rank;
    MPI_Offset m_alloc, s_alloc, x_alloc;

    if (!cfg->verbose || cfg->api != pnetcdf) return 0;

    MPI_Comm_rank(cfg->io_comm, &global_rank);

    /* check if there is any PnetCDF internal malloc residue */
    err = driver->inq_malloc_size(&m_alloc);
    if (err == NC_NOERR) {
        MPI_Reduce(&m_alloc, &s_alloc, 1, MPI_OFFSET, MPI_SUM, 0, cfg->io_comm);
        if (global_rank == 0 && s_alloc > 0) {
            printf("-------------------------------------------------------\n");
            printf("Residue heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   s_alloc);
        }
    }

    /* find the high water mark among all processes */
    driver->inq_malloc_max_size(&m_alloc);
    MPI_Reduce(&m_alloc, &x_alloc, 1, MPI_OFFSET, MPI_MAX, 0, cfg->io_comm);
    if (global_rank == 0)
        printf("High water mark of heap memory allocated by PnetCDF internally is %.2f MiB\n",
               (float)x_alloc / 1048576);

    return err;
}

/*----< wr_buf_init() >------------------------------------------------------*/
void wr_buf_init(io_buffers &buf,
                 int         gap)
{
    buf.gap             = gap;

    buf.fix_txt_buflen  = 0;
    buf.fix_int_buflen  = 0;
    buf.fix_dbl_buflen  = 0;
    buf.fix_buflen      = 0;
    buf.rec_txt_buflen  = 0;
    buf.rec_int_buflen  = 0;
    buf.rec_dbl_buflen  = 0;
    buf.rec_buflen      = 0;

    buf.fix_txt_buf  = NULL;
    buf.fix_int_buf  = NULL;
    buf.fix_dbl_buf  = NULL;
    buf.fix_buf      = NULL;
    buf.rec_txt_buf  = NULL;
    buf.rec_int_buf  = NULL;
    buf.rec_dbl_buf  = NULL;
    buf.rec_buf      = NULL;
}

/*----< wr_buf_malloc() >----------------------------------------------------*/
int wr_buf_malloc(e3sm_io_config &cfg,
                  int             one_flush,
                  io_buffers     &buf)
{
    int rank, nrecs;
    size_t j;

    MPI_Comm_rank(cfg.io_comm, &rank);

    if (cfg.run_case == F) {
        if (cfg.hist == h0) nrecs = cfg.F_case_h0.nrecs;
        else                nrecs = cfg.F_case_h1.nrecs;
    }
    else if (cfg.run_case == G) nrecs = cfg.G_case.nrecs;
    else if (cfg.run_case == I) {
        if (cfg.hist == h0) nrecs = cfg.I_case_h0.nrecs;
        else                nrecs = cfg.I_case_h1.nrecs;
    }

    if (one_flush && cfg.api == pnetcdf) {
        /* write buffers should not be touched when using PnetCDF iput before
         * ncmpi_wait_all is called. For HDF5 and ADIOS blob I/O, write data
         * will be copied and cached into internally allocated buffers and user
         * buffers can be reused after put call returned.
         */
        buf.rec_dbl_buflen *= nrecs;
        buf.rec_int_buflen *= nrecs;
        buf.rec_txt_buflen *= nrecs;
        buf.rec_buflen     *= nrecs;
    }

    /* allocate and initialize write buffers */
    buf.fix_buf     = (vtype*)  malloc(buf.fix_buflen     * sizeof(vtype));
    buf.fix_txt_buf = (char*)   malloc(buf.fix_txt_buflen * sizeof(char));
    buf.fix_int_buf = (int*)    malloc(buf.fix_int_buflen * sizeof(int));
    buf.fix_dbl_buf = (double*) malloc(buf.fix_dbl_buflen * sizeof(double));
    buf.rec_txt_buf = (char*)   malloc(buf.rec_txt_buflen * sizeof(char));
    buf.rec_int_buf = (int*)    malloc(buf.rec_int_buflen * sizeof(int));
    buf.rec_dbl_buf = (double*) malloc(buf.rec_dbl_buflen * sizeof(double));
    buf.rec_buf     = (vtype*)  malloc(buf.rec_buflen     * sizeof(vtype));

    for (j=0; j<buf.fix_txt_buflen; j++) buf.fix_txt_buf[j] = 'a' + rank;
    for (j=0; j<buf.fix_int_buflen; j++) buf.fix_int_buf[j] = rank;
    for (j=0; j<buf.fix_dbl_buflen; j++) buf.fix_dbl_buf[j] = rank;
    for (j=0; j<buf.fix_buflen;     j++) buf.fix_buf[j]     = rank;
    for (j=0; j<buf.rec_txt_buflen; j++) buf.rec_txt_buf[j] = 'a' + rank;
    for (j=0; j<buf.rec_int_buflen; j++) buf.rec_int_buf[j] = rank;
    for (j=0; j<buf.rec_dbl_buflen; j++) buf.rec_dbl_buf[j] = rank;
    for (j=0; j<buf.rec_buflen;     j++) buf.rec_buf[j]     = rank;

    return 0;
}

/*----< wr_buf_free() >------------------------------------------------------*/
void wr_buf_free(io_buffers &buf)
{
    if (buf.fix_txt_buf != NULL) free(buf.fix_txt_buf);
    if (buf.fix_int_buf != NULL) free(buf.fix_int_buf);
    if (buf.fix_dbl_buf != NULL) free(buf.fix_dbl_buf);
    if (buf.fix_buf     != NULL) free(buf.fix_buf);
    if (buf.rec_txt_buf != NULL) free(buf.rec_txt_buf);
    if (buf.rec_int_buf != NULL) free(buf.rec_int_buf);
    if (buf.rec_dbl_buf != NULL) free(buf.rec_dbl_buf);
    if (buf.rec_buf     != NULL) free(buf.rec_buf);

    buf.fix_txt_buf  = NULL;
    buf.fix_int_buf  = NULL;
    buf.fix_dbl_buf  = NULL;
    buf.fix_buf      = NULL;
    buf.rec_txt_buf  = NULL;
    buf.rec_int_buf  = NULL;
    buf.rec_dbl_buf  = NULL;
    buf.rec_buf      = NULL;

    buf.fix_txt_buflen  = 0;
    buf.fix_int_buflen  = 0;
    buf.fix_dbl_buflen  = 0;
    buf.fix_buflen      = 0;
    buf.rec_txt_buflen  = 0;
    buf.rec_int_buflen  = 0;
    buf.rec_dbl_buflen  = 0;
    buf.rec_buflen      = 0;
}

