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
int e3sm_io_case::check_malloc(e3sm_io_config *cfg,
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
void e3sm_io_case::wr_buf_init(int gap)
{
    wr_buf.gap = gap;

    wr_buf.fix_txt_buflen = 0;
    wr_buf.fix_int_buflen = 0;
    wr_buf.fix_flt_buflen = 0;
    wr_buf.fix_dbl_buflen = 0;
    wr_buf.fix_lld_buflen = 0;
    wr_buf.rec_txt_buflen = 0;
    wr_buf.rec_int_buflen = 0;
    wr_buf.rec_flt_buflen = 0;
    wr_buf.rec_dbl_buflen = 0;
    wr_buf.rec_lld_buflen = 0;

    wr_buf.fix_txt_buf = NULL;
    wr_buf.fix_int_buf = NULL;
    wr_buf.fix_flt_buf = NULL;
    wr_buf.fix_dbl_buf = NULL;
    wr_buf.fix_lld_buf = NULL;
    wr_buf.rec_txt_buf = NULL;
    wr_buf.rec_int_buf = NULL;
    wr_buf.rec_flt_buf = NULL;
    wr_buf.rec_dbl_buf = NULL;
    wr_buf.rec_lld_buf = NULL;
}

/*----< wr_buf_malloc() >----------------------------------------------------*/
int e3sm_io_case::wr_buf_malloc(e3sm_io_config &cfg, int ffreq)
{
    int rank;
    size_t j;

    MPI_Comm_rank(cfg.io_comm, &rank);

    if (cfg.api == adios) {
        wr_buf.fix_txt_buflen += 64;
        wr_buf.fix_int_buflen += 64;
        wr_buf.fix_flt_buflen += 64;
        wr_buf.fix_dbl_buflen += 64;
        wr_buf.fix_lld_buflen += 64;
        wr_buf.rec_txt_buflen += 64;
        wr_buf.rec_int_buflen += 64;
        wr_buf.rec_flt_buflen += 64;
        wr_buf.rec_dbl_buflen += 64;
        wr_buf.rec_lld_buflen += 64;
    }

    if (cfg.api != adios && !(cfg.strategy == blob && cfg.api == hdf5)) {
        /* Note HDF5 and ADIOS blob I/O copy write data into their internal
         * buffers and only flush them out at file close. Thus, write buffers
         * can be reused for these two I/O methods. For others, such as PnetCDF
         * and HDF5 log-based VOL, write buffers should not be touched as they
         * will later be used during the flushing is called.
         */
        wr_buf.rec_txt_buflen *= ffreq;
        wr_buf.rec_int_buflen *= ffreq;
        wr_buf.rec_flt_buflen *= ffreq;
        wr_buf.rec_dbl_buflen *= ffreq;
        wr_buf.rec_lld_buflen *= ffreq;
    }

    /* allocate and initialize write buffers */
    wr_buf.fix_txt_buf = (char*)   malloc(wr_buf.fix_txt_buflen * sizeof(char));
    wr_buf.fix_int_buf = (int*)    malloc(wr_buf.fix_int_buflen * sizeof(int));
    wr_buf.fix_flt_buf = (float*)  malloc(wr_buf.fix_flt_buflen * sizeof(float));
    wr_buf.fix_dbl_buf = (double*) malloc(wr_buf.fix_dbl_buflen * sizeof(double));
    wr_buf.fix_lld_buf = (long long*) malloc(wr_buf.fix_lld_buflen * sizeof(long long));
    wr_buf.rec_txt_buf = (char*)   malloc(wr_buf.rec_txt_buflen * sizeof(char));
    wr_buf.rec_int_buf = (int*)    malloc(wr_buf.rec_int_buflen * sizeof(int));
    wr_buf.rec_flt_buf = (float*)  malloc(wr_buf.rec_flt_buflen * sizeof(float));
    wr_buf.rec_dbl_buf = (double*) malloc(wr_buf.rec_dbl_buflen * sizeof(double));
    wr_buf.rec_lld_buf = (long long*) malloc(wr_buf.rec_lld_buflen * sizeof(long long));

    for (j=0; j<wr_buf.fix_txt_buflen; j++) wr_buf.fix_txt_buf[j] = 'a' + rank;
    for (j=0; j<wr_buf.fix_int_buflen; j++) wr_buf.fix_int_buf[j] = rank;
    for (j=0; j<wr_buf.fix_flt_buflen; j++) wr_buf.fix_flt_buf[j] = rank;
    for (j=0; j<wr_buf.fix_dbl_buflen; j++) wr_buf.fix_dbl_buf[j] = rank;
    for (j=0; j<wr_buf.fix_lld_buflen; j++) wr_buf.fix_lld_buf[j] = rank;
    for (j=0; j<wr_buf.rec_txt_buflen; j++) wr_buf.rec_txt_buf[j] = 'a' + rank;
    for (j=0; j<wr_buf.rec_int_buflen; j++) wr_buf.rec_int_buf[j] = rank;
    for (j=0; j<wr_buf.rec_flt_buflen; j++) wr_buf.rec_flt_buf[j] = rank;
    for (j=0; j<wr_buf.rec_dbl_buflen; j++) wr_buf.rec_dbl_buf[j] = rank;
    for (j=0; j<wr_buf.rec_lld_buflen; j++) wr_buf.rec_lld_buf[j] = rank;

    return 0;
}

/*----< wr_buf_free() >------------------------------------------------------*/
void e3sm_io_case::wr_buf_free(void)
{
    if (wr_buf.fix_txt_buf != NULL) free(wr_buf.fix_txt_buf);
    if (wr_buf.fix_int_buf != NULL) free(wr_buf.fix_int_buf);
    if (wr_buf.fix_flt_buf != NULL) free(wr_buf.fix_flt_buf);
    if (wr_buf.fix_dbl_buf != NULL) free(wr_buf.fix_dbl_buf);
    if (wr_buf.fix_lld_buf != NULL) free(wr_buf.fix_lld_buf);
    if (wr_buf.rec_txt_buf != NULL) free(wr_buf.rec_txt_buf);
    if (wr_buf.rec_int_buf != NULL) free(wr_buf.rec_int_buf);
    if (wr_buf.rec_flt_buf != NULL) free(wr_buf.rec_flt_buf);
    if (wr_buf.rec_dbl_buf != NULL) free(wr_buf.rec_dbl_buf);
    if (wr_buf.rec_lld_buf != NULL) free(wr_buf.rec_lld_buf);

    wr_buf.fix_txt_buf = NULL;
    wr_buf.fix_int_buf = NULL;
    wr_buf.fix_flt_buf = NULL;
    wr_buf.fix_dbl_buf = NULL;
    wr_buf.fix_lld_buf = NULL;
    wr_buf.rec_txt_buf = NULL;
    wr_buf.rec_int_buf = NULL;
    wr_buf.rec_flt_buf = NULL;
    wr_buf.rec_dbl_buf = NULL;
    wr_buf.rec_lld_buf = NULL;

    wr_buf.fix_txt_buflen = 0;
    wr_buf.fix_int_buflen = 0;
    wr_buf.fix_flt_buflen = 0;
    wr_buf.fix_dbl_buflen = 0;
    wr_buf.fix_lld_buflen = 0;
    wr_buf.rec_txt_buflen = 0;
    wr_buf.rec_int_buflen = 0;
    wr_buf.rec_flt_buflen = 0;
    wr_buf.rec_dbl_buflen = 0;
    wr_buf.rec_lld_buflen = 0;
}

