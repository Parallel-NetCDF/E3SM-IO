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
#include <string.h>
#include <assert.h>
#include <unistd.h> /* unlink() */

#include <mpi.h>

#include <e3sm_io.h>
#include <e3sm_io_err.h>
#include <e3sm_io_case.hpp>
#include <e3sm_io_case_scorpio.hpp>
#include <e3sm_io_driver.hpp>

#define VAR_ITYPE REC_ITYPE

#define CHECK_VAR_ERR(varid) {                                            \
    if (err != 0) {                                                       \
        char var_name[64];                                                \
        driver.inq_var_name(ncid, varid.vid, var_name);                   \
        printf("Error in %s:%d: %s() var %s\n"     ,                      \
               __FILE__, __LINE__, __func__, var_name);                   \
        goto err_out;                                                     \
    }                                                                     \
}
#define FILE_CREATE(filename) {                                           \
    err = driver.create(filename, comm, cfg.info, &ncid);                 \
    CHECK_ERR                                                             \
}
#define FILE_CLOSE {                                                      \
    err = driver.close(ncid);                                             \
    CHECK_ERR                                                             \
}
#define ENDDEF {                                                          \
    err = driver.enddef(ncid);                                            \
    CHECK_ERR                                                             \
}
#define INQ_DIM_LEN(name, size) {                                         \
    int dimid;                                                            \
    err = driver.inq_dim(ncid, name, &dimid);                             \
    CHECK_ERR                                                             \
    err = driver.inq_dimlen(ncid, dimid, &size);                          \
    CHECK_ERR                                                             \
}
#define INQ_PUT_SIZE(size) {                                              \
    err = driver.inq_put_size(&size);                                     \
    CHECK_ERR                                                             \
}
#define INQ_FILE_INFO(info) {                                             \
    err = driver.inq_file_info(ncid, &info);                              \
    CHECK_ERR                                                             \
}
#define IPUT_VARA_NOADV(itype, buf) {                                     \
    err = driver.put_vara(ncid, varid, itype, start, count, buf, nb);     \
    CHECK_VAR_ERR(varid)                                                  \
    my_nreqs++;                                                           \
    varid++;                                                              \
}
#define IPUT_VARA(varid, itype, adv, buf) {                               \
    err = e3sm_io_scorpio_write_var(driver, rec_no, ncid, varid, itype, buf, nb); \
    CHECK_VAR_ERR(varid)                                                  \
    buf += (adv);                                                         \
    my_nreqs++;                                                           \
}
#define IPUT_VAR(varid, itype, adv, buf) {                                \
    err = e3sm_io_scorpio_write_var(driver, -1, ncid, varid, itype, buf, nb); \
    CHECK_VAR_ERR(varid)                                                  \
    buf += (adv);                                                         \
    my_nreqs++;                                                           \
}
#define WAIT_ALL_REQS {                                                   \
    err = driver.wait(ncid);                                              \
    CHECK_ERR                                                             \
    nflushes++;                                                           \
}
#define FIX_VAR_IPUT(varid, dp, itype, buf) {                             \
    err = e3sm_io_scorpio_write_var(driver, -1, ncid, varid, itype, buf, nb);     \
    my_nreqs += decom.contig_nreqs[dp];                                   \
    CHECK_VAR_ERR(varid)                                                  \
    buf += decom.raw_nreqs[dp] + gap;                                     \
    nvars_D[dp]++;                                                        \
}

#define REC_VAR_IPUT(varid, dp, itype, buf) {                             \
    err = e3sm_io_scorpio_write_var(driver, rec_no, ncid, varid, itype, buf, nb); \
    my_nreqs += decom.contig_nreqs[dp];                                   \
    CHECK_VAR_ERR(varid)                                                  \
    buf += decom.raw_nreqs[dp] + gap;                                     \
    if (rec_no == 0) nvars_D[dp]++;                                       \
}

/*----< var_io_I_case_scorpio.cpp() >----------------------------------------*/
int var_wr_all_cases_scorpio(e3sm_io_config &cfg,
                             e3sm_io_decom  &decom,
                             e3sm_io_driver &driver,
                             case_meta      *cmeta)
{
    char *fix_txt_buf_ptr, *rec_txt_buf_ptr;
    int i, j, err=0, sub_rank, global_rank, ncid=-1, nflushes=0;
    int rec_no, ffreq, gap=0, my_nreqs, nvars, num_decomp_vars, *nvars_D;
    int *fix_int_buf_ptr, *rec_int_buf_ptr;
    double *fix_dbl_buf_ptr, *rec_dbl_buf_ptr, timing;
    MPI_Offset previous_size, metadata_size, total_size;
    MPI_Comm comm;
    vtype *fix_buf_ptr, *rec_buf_ptr;
    io_buffers wr_buf;

    var_meta_scorpio *vars;
    int scorpiovars[7];

    MPI_Barrier(cfg.io_comm); /*---------------------------------------------*/
    cmeta->end2end_time = timing = MPI_Wtime();

    cmeta->post_time  = 0.0;
    cmeta->flush_time = 0.0;

    if (cfg.strategy == blob && cfg.api != adios)
        comm = cfg.sub_comm;
    else
        comm = cfg.io_comm;

    MPI_Comm_rank(cfg.io_comm,  &global_rank);
    MPI_Comm_rank(comm,         &sub_rank);

    /* I/O amount from previous I/O */
    INQ_PUT_SIZE(previous_size)

    if (cfg.verbose) {
        if (cfg.strategy == blob && sub_rank == 0)
            printf("global_rank=%d sub_rank=%d outfile=%s\n",
                   global_rank,sub_rank,cmeta->outfile);
        else if (global_rank == 0)
            printf("global_rank=%d outfile=%s\n",global_rank,cmeta->outfile);
    }

    /* create the output file */
    FILE_CREATE(cmeta->outfile)

    cmeta->open_time = MPI_Wtime() - timing;

    MPI_Barrier(comm); /*----------------------------------------------------*/
    timing = MPI_Wtime();

    /* there are num_decomp_vars number of decomposition variables */
    if (cfg.strategy == blob && cfg.api != adios)
        num_decomp_vars = decom.num_decomp * NVARS_DECOMP;
    else
        num_decomp_vars = 0;

    nvars = cmeta->nvars + num_decomp_vars;

    /* allocate space to store variable metadata */
    vars = (var_meta_scorpio*) malloc(nvars * sizeof(var_meta_scorpio));

    if (cfg.non_contig_buf) gap = BUF_GAP;

    wr_buf_init(wr_buf, gap);

    if (cfg.run_case == F)
        err = def_F_case_scorpio(driver, cfg, decom, ncid, vars, scorpiovars,
                                 &wr_buf);
    else if (cfg.run_case == G)
        err = def_G_case_scorpio(cfg, decom, driver, ncid, vars, scorpiovars,
                                 &wr_buf);
    else if (cfg.run_case == I)
        err = def_I_case_scorpio(cfg, decom, driver, ncid, vars, scorpiovars, 
                                 &wr_buf);
    CHECK_ERR

    /* exit define mode and enter data mode */
    ENDDEF

    /* I/O amount so far */
    INQ_PUT_SIZE(metadata_size)
    metadata_size -= previous_size;

    INQ_FILE_INFO(cmeta->info_used)

    cmeta->def_time = MPI_Wtime() - timing;

    MPI_Barrier(comm); /*----------------------------------------------------*/
    timing = MPI_Wtime();

    /* flush drequency only affects pnetcdf API.  Note for HDF5 and ADIOS blob
     * I/O, write data is copied into their internal buffers and only flushed
     * at file close. Calling driver.wait() takes no effect.
     */
    ffreq = (cfg.api == pnetcdf) ? cmeta->ffreq : 1;

    /* allocate write buffers */
    wr_buf_malloc(cfg, ffreq, wr_buf);

    nvars_D = cmeta->nvars_D;
    for (j=0; j<decom.num_decomp; j++)
        nvars_D[j] = 0; /* number of variables using decomposition j */

    cmeta->pre_time = MPI_Wtime() - timing;

    MPI_Barrier(comm); /*----------------------------------------------------*/
    timing = MPI_Wtime();

    my_nreqs = 0; /* number of noncontiguous requests written by this rank */

    for (j=0; j<decom.num_decomp; j++) {
        err = driver.put_varl(ncid, scorpiovars[j], MPI_LONG_LONG,
                              decom.raw_offsets[j], nb);
        CHECK_ERR
    }
    // Nproc only written by rank 0
    if (cfg.rank == 0) {
        err = driver.put_varl(ncid, scorpiovars[j], MPI_INT, &cfg.np, nb);
        CHECK_ERR
    }

    fix_txt_buf_ptr = wr_buf.fix_txt_buf;
    fix_int_buf_ptr = wr_buf.fix_int_buf;
    fix_dbl_buf_ptr = wr_buf.fix_dbl_buf;
    fix_buf_ptr     = wr_buf.fix_buf;

    for (rec_no=0; rec_no<cmeta->nrecs; rec_no++) {

        if (rec_no % ffreq == 0) {
            /* reset the pointers to the beginning of the buffers */
            rec_txt_buf_ptr = wr_buf.rec_txt_buf;
            rec_int_buf_ptr = wr_buf.rec_int_buf;
            rec_dbl_buf_ptr = wr_buf.rec_dbl_buf;
            rec_buf_ptr     = wr_buf.rec_buf;
        }

        if (cfg.strategy == canonical) {
            /* set the start index for the next record */
            for (i=0; i<decom.num_decomp; i++) {
                for (j=0; j<decom.contig_nreqs[i]; j++)
                    decom.w_starts[i][j][0] = rec_no;
            }
        }

        /* write all climate variables */
        for (j=num_decomp_vars; j<nvars; j++) {
            var_meta_scorpio varid = vars[j];
            int          dp    = vars[j].decomp_id;
            MPI_Datatype itype = vars[j].itype;
            size_t       adv   = vars[j].vlen + gap;

            if (vars[j].decomp_id >= 0) { /* this variable is partitioned */
                if (vars[j].isRecVar) { /* this is a record variable */
                    REC_VAR_IPUT(varid, dp, VAR_ITYPE, rec_buf_ptr)
                }
                else if (rec_no == 0) { /* this is a fixed-size variable */
                    if (itype == VAR_ITYPE)
                        FIX_VAR_IPUT(varid, dp, VAR_ITYPE, fix_buf_ptr)
                    else if (itype == MPI_INT)
                        FIX_VAR_IPUT(varid, dp, MPI_INT, fix_int_buf_ptr)
                }
            }
            else if (sub_rank == 0) {
                /* not-partitioned variables are written by root only */

                if (vars[j].isRecVar) { /* this is a record variable */
                    if (itype == VAR_ITYPE)
                        IPUT_VARA(varid, itype, adv,     rec_buf_ptr)
                    else if (itype == MPI_INT)
                        IPUT_VARA(varid, itype, adv, rec_int_buf_ptr)
                    else if (itype == MPI_CHAR)
                        IPUT_VARA(varid, itype, adv, rec_txt_buf_ptr)
                    else if (itype == MPI_DOUBLE)
                        IPUT_VARA(varid, itype, adv, rec_dbl_buf_ptr)
                }
                else if (rec_no == 0) { /* this is a fixed-size variable */
                    if (itype == VAR_ITYPE)
                        IPUT_VAR(varid, itype, adv,     fix_buf_ptr)
                    else if (itype == MPI_INT)
                        IPUT_VAR(varid, itype, adv, fix_int_buf_ptr)
                    else if (itype == MPI_CHAR)
                        IPUT_VAR(varid, itype, adv, fix_txt_buf_ptr)
                    else if (itype == MPI_DOUBLE)
                        IPUT_VAR(varid, itype, adv, fix_dbl_buf_ptr)
                }
            }
        }

        /* flush out the pending iput requests */
        if ((rec_no + 1) % ffreq == 0 || (rec_no + 1) == cmeta->nrecs) {
            cmeta->post_time += MPI_Wtime() - timing;

            MPI_Barrier(comm); /*--------------------------------------------*/
            timing = MPI_Wtime();

            /* flush once per time record */
            WAIT_ALL_REQS
            cmeta->flush_time += MPI_Wtime() - timing;

            if ((rec_no + 1) < cmeta->nrecs) timing = MPI_Wtime();
        }
    }

    MPI_Barrier(comm); /*----------------------------------------------------*/
    timing = MPI_Wtime();

    /* free up previously allocated heap memory space */
    wr_buf_free(wr_buf);
    if (vars != NULL) free(vars);

    FILE_CLOSE

    cmeta->close_time = MPI_Wtime() - timing;

    /* obtain the write amount tracked by the driver */
    INQ_PUT_SIZE(total_size)
    total_size -= previous_size;

    cmeta->num_flushes     = nflushes;
    cmeta->num_decomp_vars = num_decomp_vars;
    cmeta->my_nreqs        = my_nreqs;
    cmeta->metadata_WR     = metadata_size;
    cmeta->amount_WR       = total_size;
    cmeta->end2end_time    = MPI_Wtime() - cmeta->end2end_time;

    /* check if there is any PnetCDF internal malloc residue */
    check_malloc(&cfg, &driver);

    /* note inquiring file size may be expensive on some machines */
    if (cfg.verbose && global_rank == 0)
        driver.inq_file_size(cmeta->outfile, &cmeta->file_size);

err_out:
    if (err < 0 && ncid >= 0) driver.close(ncid);
    if (!cfg.keep_outfile && sub_rank == 0) unlink(cmeta->outfile);

    return err;
}

int e3sm_io_scorpio_define_dim(e3sm_io_driver &driver,
                               int fid,
                               std::string name,
                               MPI_Offset size,
                               std::map<int, std::string> &dnames,
                               int *did)
{
    int err;

    err = driver.def_dim (fid, name, size, did);
    CHECK_ERR

    dnames[*did] = name;

err_out:
    return err;
}

int e3sm_io_scorpio_define_var(e3sm_io_driver &driver,
                               e3sm_io_config &cfg,
                               std::map<int, std::string> &dnames,
                               e3sm_io_decom &decom,
                               int decomp_id,
                               int fid,
                               std::string name,
                               nc_type xtype,
                               int ndims,
                               int *dimids,
                               var_meta_scorpio *var)
{
    int i, err, ibuf;
    char cbuf[64];
    std::vector<const char*> dnames_array (ndims);

    var->piodecomid = decomp_id + 512;
    var->ndims = ndims;

    for(i = 0; i < ndims; i++){
        dnames_array[i] = dnames[dimids[i]].c_str();
    }

    // If there is a decomposition map associated with the variable,
    // created 2 associated scalar variables frame_id (timesteps) and decom_id (decomposition map)
    if (decomp_id >= 0) {
        MPI_Offset one = 1;

        err = driver.def_local_var (fid, name, xtype, 1, &(decom.raw_nreqs[decomp_id]), &var->vid);
        CHECK_ERR

        assert(ndims > 0);

        err = driver.def_local_var(fid, "frame_id/" + name, NC_INT, 1,
                                   &one, &var->frame_id);
        CHECK_ERR

        err = driver.def_local_var(fid, "decomp_id/" + name, NC_INT, 1,
                                   &one, &var->decom_id);
        CHECK_ERR

        // Double vars have an additional fillval_id
        if (xtype == NC_DOUBLE) {
            err = driver.def_local_var(fid, "fillval_id/" + name, xtype,
                                       1, &one, &(var->fillval_id));
            CHECK_ERR
        } else {
            var->fillval_id = -1;
        }

        // Attributes
        // err = driver.put_att (fid, var->vid, "_FillValue", var->type, 1, cbuf);
        // CHECK_ERR

        // Scorpio attributes are only written by rank 0
        if (cfg.rank == 0) {
            // Decomposition map
            sprintf(cbuf, "%d", var->piodecomid);
            err = driver.put_att (fid, var->vid, "__pio__/decomp", NC_CHAR, strlen(cbuf), &cbuf);
            CHECK_ERR

            err = driver.put_att (fid, var->vid, "__pio__/dims", NC_STRING, ndims, dnames_array.data());
            CHECK_ERR

            // Type of NetCDF API called
            err = driver.put_att (fid, var->vid, "__pio__/ncop", NC_CHAR, 7, (void *)"darray");
            CHECK_ERR

            // NetCDF type enum
            ibuf = xtype;
            err  = driver.put_att (fid, var->vid, "__pio__/nctype", NC_INT, 1, &ibuf);
            CHECK_ERR

            // Number of dimensions
            err = driver.put_att (fid, var->vid, "__pio__/ndims", NC_INT, 1, &ndims);
            CHECK_ERR
        }
    } else { /* this variable is not partitioned */
        std::vector<MPI_Offset> dsize (ndims);
        MPI_Offset vsize = 1;
        int esize;

        for (i = 0; i < ndims; i++) {
            err = driver.inq_dimlen (fid, dimids[i], &dsize[i]);
            CHECK_ERR

            // Time dim is always 0, but block size should be 1
            if (dsize[i] == 0) dsize[i] = 1;
            var->dims[i] = dsize[i];
            vsize *= dsize[i];
        }

        // flatten into 1 dim only apply to non-scalar variables
        if (ndims){
            // Convert into byte array
            err = e3sm_io_xlen_nc_type(xtype, &esize);
            CHECK_MPIERR
            vsize *= esize;
            vsize += 8 * 2 * ndims; // Include start and count array
            err = driver.def_local_var (fid, name, NC_UBYTE, 1, &vsize, &var->vid);
            CHECK_ERR
        }
        else{
            err = driver.def_local_var (fid, name, xtype, ndims, NULL, &var->vid);
            CHECK_ERR
        }

        // Attributes for non-constant small vars
        // Scorpio attributes are only written by rank 0
        if (cfg.rank == 0) {
            // ADIOS type enum, variables without decomposition map are stored as byte array
            ibuf = (int)e3sm_io_type_nc2adios (xtype);
            err  = driver.put_att (fid, var->vid, "__pio__/adiostype", NC_INT, 1, &ibuf);
            CHECK_ERR

            // Scalar var does not have dims
            if (ndims > 0) {
                err = driver.put_att (fid, var->vid, "__pio__/dims", NC_STRING, ndims, dnames_array.data());
                CHECK_ERR
            }

            // Type of NetCDF API called
            err = driver.put_att (fid, var->vid, "__pio__/ncop", NC_CHAR, 7, (void *)"put_var");
            CHECK_ERR

            // NetCDF type enum
            ibuf = xtype;
            err  = driver.put_att (fid, var->vid, "__pio__/nctype", NC_INT, 1, &ibuf);
            CHECK_ERR

            // Number of dimensions
            err = driver.put_att (fid, var->vid, "__pio__/ndims", NC_INT, 1, &ndims);
            CHECK_ERR
        }

        var->decom_id   = -1;
        var->frame_id   = -1;
        var->fillval_id = -1;
    }
err_out:
    return err;
}

int e3sm_io_scorpio_write_var (e3sm_io_driver &driver,
                                  int frameid,
                                  int fid,
                                  var_meta_scorpio &var,
                                  MPI_Datatype itype,
                                  void *buf,
                                  e3sm_io_op_mode mode) {
    int err, decomid;
    void *wbuf;

    if (var.isRecVar) frameid = -1;

    /* prepend start and count to write buffer for small, not-partitioned,
     * non-scalar variables This must be done by allocating another buffer,
     * add start/count, and copy over user write buffer.
     * These variables are stored as arrays of type byte.
     */
    if (var.decomp_id < 0 && var.ndims) {
        int esize;
        size_t cp_len;

#ifdef DEBUG
        size_t i, vlen=1;
        for (i=0; i<var.ndims; i++) vlen *= var.dims[i];
        assert(vlen == var.vlen);
#endif

        MPI_Type_size(itype, &esize);
        cp_len = var.ndims * sizeof(int64_t);
        wbuf = (void*) malloc(2 * cp_len + var.vlen * esize);
        memset(wbuf, 0, cp_len);
        memcpy((char*)wbuf+cp_len,   var.dims, cp_len);
        memcpy((char*)wbuf+cp_len*2, buf,      var.vlen * esize);
    }
    else
        wbuf = buf;

    if ((var.frame_id < 0) && var.ndims){
        itype = MPI_BYTE;
    }

    err = driver.put_varl (fid, var.vid, itype, wbuf, nbe);
    CHECK_ERR

    if (var.frame_id >= 0) {

        err = driver.put_varl (fid, var.frame_id, MPI_INT, &frameid, nbe);
        CHECK_ERR

        decomid = var.piodecomid ;
        if (var.fillval_id < 0) {
            decomid *= -1;
        }
        err = driver.put_varl (fid, var.decom_id, MPI_INT, &decomid, nbe);
        CHECK_ERR

        if (var.fillval_id >= 0) {
            double fbuf = 1e+20;

            err = driver.put_varl (fid, var.fillval_id, MPI_DOUBLE, &fbuf, nbe);
            CHECK_ERR
        }
    }

err_out:
    if (wbuf != buf) free(wbuf);

    return err;
}

int e3sm_io_scorpio_put_att (e3sm_io_driver &driver,
                                int fid,
                                int vid,
                                std::string name,
                                nc_type xtype,
                                MPI_Offset size,
                                void *buf) {
    return driver.put_att (fid, vid, name, xtype, size, buf);
}

int e3sm_io_scorpio_put_att (e3sm_io_driver &driver,
                                int fid,
                                var_meta_scorpio &var,
                                std::string name,
                                nc_type xtype,
                                MPI_Offset size,
                                void *buf) {
    return driver.put_att (fid, var.vid, name, xtype, size, buf);
}

