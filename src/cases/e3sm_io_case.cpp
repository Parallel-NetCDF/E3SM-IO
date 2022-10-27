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
#include <string.h> /* strcpy(), strcat(), strdup() */
#include <assert.h>

#include <e3sm_io.h>
#include <e3sm_io_err.h>
#include <e3sm_io_case.hpp>

e3sm_io_case::e3sm_io_case () {}
e3sm_io_case::~e3sm_io_case () {}

int e3sm_io_case::wr_test(e3sm_io_config &cfg,
                          e3sm_io_decom  &decom,
                          e3sm_io_driver &driver)
{
    int err=0;
    char *ext;
    case_meta *cmeta;

    /* construct I/O metadata */
    err = calc_metadata(&cfg, &decom);
    CHECK_ERR

    ext = strrchr(cfg.out_path, '.');

    if (cfg.run_case == G) {
        cmeta            = &cfg.G_case;
        cmeta->nrecs     =  cfg.G_case.nrecs;
        cmeta->nvars     = NVARS_G_CASE;
        cmeta->info_used = MPI_INFO_NULL;

        /* set flush frequency, once per ffreq time steps */
        if (cmeta->ffreq == -1 || cmeta->ffreq > cmeta->nrecs)
            cmeta->ffreq = cmeta->nrecs;

        /* construct file name */
        strcpy(cmeta->outfile, cfg.out_path);
        if (cfg.strategy == blob && cfg.api != adios) {
            /* append subfile ID to subfile name */
            char sub_str[8];
            sprintf(sub_str, ".%04d", cfg.subfile_ID);
            strcat(cmeta->outfile, sub_str);
        }

        cfg.nvars = cmeta->nvars;
        err = var_wr_case(cfg, decom, driver, cmeta);
        CHECK_ERR

        goto err_out;
    }

    /* Only F anc I cases write to h0 and h1 files */
    if (cfg.hx == 0 || cfg.hx == -1) {  /* h0 file */
        if (cfg.run_case == F) {
            cmeta        = &cfg.F_case_h0;
            cmeta->nrecs =  cfg.F_case_h0.nrecs;
            cmeta->nvars = NVARS_F_CASE_H0;
        }
        else if (cfg.run_case == I) {
            cmeta        = &cfg.I_case_h0;
            cmeta->nrecs =  cfg.I_case_h0.nrecs;
            cmeta->nvars = NVARS_I_CASE_H0;
        }
        cmeta->info_used = MPI_INFO_NULL;

        /* set flush frequency, once per ffreq time steps */
        if (cmeta->ffreq == -1 || cmeta->ffreq > cmeta->nrecs)
            cmeta->ffreq = cmeta->nrecs;

        /* construct file name */
        if (ext == NULL || (strcmp(ext, ".nc") && strcmp(ext, ".h5") && strcmp(ext, ".nc4") && strcmp(ext, ".bp")))
            sprintf(cmeta->outfile, "%s_h0", cfg.out_path);
        else { /* add "_h0" before file extension */
            strcpy(cmeta->outfile, cfg.out_path);
            sprintf(cmeta->outfile + (ext - cfg.out_path), "_h0%s", ext);
        }
        if (cfg.strategy == blob && cfg.api != adios) {
            /* append subfile ID to subfile name */
            char sub_str[8];
            sprintf(sub_str, ".%04d", cfg.subfile_ID);
            strcat(cmeta->outfile, sub_str);
        }

        cfg.hist = h0;
        cfg.nvars = cmeta->nvars;
        err = var_wr_case(cfg, decom, driver, cmeta);
        CHECK_ERR
    }

    if (cfg.hx == 1 || cfg.hx == -1) {  /* h1 file */
        if (cfg.run_case == F) {
            cmeta        = &cfg.F_case_h1;
            cmeta->nrecs =  cfg.F_case_h1.nrecs;
            cmeta->nvars = NVARS_F_CASE_H1;
        }
        else if (cfg.run_case == I) {
            cmeta        = &cfg.I_case_h1;
            cmeta->nrecs =  cfg.I_case_h1.nrecs;
            cmeta->nvars = NVARS_I_CASE_H1;
        }
        cmeta->info_used = MPI_INFO_NULL;

        /* set flush frequency, once per ffreq time steps */
        if (cmeta->ffreq == -1 || cmeta->ffreq > cmeta->nrecs)
            cmeta->ffreq = cmeta->nrecs;

        /* construct file name */
        if (ext == NULL || (strcmp(ext, ".nc") && strcmp(ext, ".h5") && strcmp(ext, ".nc4") && strcmp(ext, ".bp")))
            sprintf(cmeta->outfile, "%s_h1", cfg.out_path);
        else { /* add "_h1" before file extension */
            strcpy(cmeta->outfile, cfg.out_path);
            sprintf(cmeta->outfile + (ext - cfg.out_path), "_h1%s", ext);
        }
        if (cfg.strategy == blob && cfg.api != adios) {
            /* append subfile ID to subfile name */
            char sub_str[8];
            sprintf(sub_str, ".%04d", cfg.subfile_ID);
            strcat(cmeta->outfile, sub_str);
        }

        cfg.hist = h1;
        cfg.nvars = cmeta->nvars;
        err = var_wr_case(cfg, decom, driver, cmeta);
        CHECK_ERR
    }

err_out:
    if (cfg.sub_comm != MPI_COMM_NULL) {
        MPI_Comm_free(&cfg.sub_comm);
        cfg.sub_comm = MPI_COMM_NULL;
    }

    return err;
}

int e3sm_io_case::rd_test(e3sm_io_config &cfg,
                          e3sm_io_decom &decom,
                          e3sm_io_driver &driver)
{
    int err=0;
    char *ext;
    case_meta *cmeta;

    if (cfg.strategy != canonical)
        ERR_OUT ("Read is only supported in canonical storage layout")

    if (cfg.api == hdf5)
        ERR_OUT ("Read currently does not support HDF5 methods")

    /* construct I/O metadata */
    err = calc_metadata(&cfg, &decom);
    CHECK_ERR

    ext = strrchr(cfg.in_path, '.');

    if (cfg.run_case == G) {
        cmeta            = &cfg.G_case;
        cmeta->nrecs     =  cfg.G_case.nrecs;
        cmeta->nvars     = NVARS_G_CASE;
        cmeta->info_used = MPI_INFO_NULL;

        /* set flush frequency, once per ffreq time steps */
        if (cmeta->ffreq == -1 || cmeta->ffreq > cmeta->nrecs)
            cmeta->ffreq = cmeta->nrecs;

        /* construct file name */
        strcpy(cmeta->outfile, cfg.in_path);

        cfg.nvars = cmeta->nvars;
        err = var_rd_case(cfg, decom, driver, cmeta);
        CHECK_ERR

        goto err_out;
    }

    /* Only F anc I cases write to h0 and h1 files */
    if (cfg.hx == 0 || cfg.hx == -1) {  /* h0 file */
        if (cfg.run_case == F) {
            cmeta        = &cfg.F_case_h0;
            cmeta->nrecs =  cfg.F_case_h0.nrecs;
            cmeta->nvars = NVARS_F_CASE_H0;
        }
        else if (cfg.run_case == I) {
            cmeta        = &cfg.I_case_h0;
            cmeta->nrecs =  cfg.I_case_h0.nrecs;
            cmeta->nvars = NVARS_I_CASE_H0;
        }
        cmeta->info_used = MPI_INFO_NULL;

        /* set flush frequency, once per ffreq time steps */
        if (cmeta->ffreq == -1 || cmeta->ffreq > cmeta->nrecs)
            cmeta->ffreq = cmeta->nrecs;

        /* construct file name */
        if (ext == NULL || (strcmp(ext, ".nc") && strcmp(ext, ".h5") && strcmp(ext, ".nc4")))
            sprintf(cmeta->outfile, "%s_h0", cfg.in_path);
        else { /* add "_h0" before file extension */
            strcpy(cmeta->outfile, cfg.in_path);
            sprintf(cmeta->outfile + (ext - cfg.in_path), "_h0%s", ext);
        }

        cfg.hist = h0;
        cfg.nvars = cmeta->nvars;
        err = var_rd_case(cfg, decom, driver, cmeta);
        CHECK_ERR
    }

    if (cfg.hx == 1 || cfg.hx == -1) {  /* h1 file */
        if (cfg.run_case == F) {
            cmeta        = &cfg.F_case_h1;
            cmeta->nrecs =  cfg.F_case_h1.nrecs;
            cmeta->nvars = NVARS_F_CASE_H1;
        }
        else if (cfg.run_case == I) {
            cmeta        = &cfg.I_case_h1;
            cmeta->nrecs =  cfg.I_case_h1.nrecs;
            cmeta->nvars = NVARS_I_CASE_H1;
        }
        cmeta->info_used = MPI_INFO_NULL;

        /* set flush frequency, once per ffreq time steps */
        if (cmeta->ffreq == -1 || cmeta->ffreq > cmeta->nrecs)
            cmeta->ffreq = cmeta->nrecs;

        /* construct file name */
        if (ext == NULL || (strcmp(ext, ".nc") && strcmp(ext, ".h5") && strcmp(ext, ".nc4")))
            sprintf(cmeta->outfile, "%s_h1", cfg.in_path);
        else { /* add "_h1" before file extension */
            strcpy(cmeta->outfile, cfg.in_path);
            sprintf(cmeta->outfile + (ext - cfg.in_path), "_h1%s", ext);
        }

        cfg.hist = h1;
        cfg.nvars = cmeta->nvars;
        err = var_rd_case(cfg, decom, driver, cmeta);
        CHECK_ERR
    }

err_out:
    if (cfg.sub_comm != MPI_COMM_NULL) {
        MPI_Comm_free(&cfg.sub_comm);
        cfg.sub_comm = MPI_COMM_NULL;
    }

    return err;
}

int e3sm_io_case::def_var(e3sm_io_config            &cfg,
                          e3sm_io_decom             &decom,
                          e3sm_io_driver            &driver,
                          case_meta                 *cmeta,
                          int                        ncid,
                          std::string                name,
                          std::map<int, std::string> &dnames,
                          int                        xtype,
                          int                        nDims,
                          int                        dim_time,
                          int                       *dimids,
                          MPI_Datatype               itype,
                          int                        decomid,
                          var_meta                  *varp)
{
    /* nDims and dimids are canonical dimensions */
    int err=0, j, *_dimids = dimids;
    varp->_name     = strdup(name.c_str());
    varp->ndims     = nDims;   /* number of dimensions */
    varp->iType     = itype;   /* internal data type of write buffer */
    varp->xType     = xtype;   /* external data type of variable in file */
    varp->decomp_id = decomid; /* decomposition map ID */
    varp->isRecVar  = (nDims != 0 && *_dimids == dim_time);
    /* calculate variable size */
    for (varp->vlen=1, j=0; j<nDims; j++) {
        err = driver.inq_dimlen(ncid, _dimids[j], &varp->dims[j]);
        CHECK_ERR
        if (j == 0 && varp->isRecVar) varp->dims[j] = 1;
        varp->vlen *= varp->dims[j];
    }
    /* define a new variable */
    if (cfg.api == adios) {
        err = scorpio_define_var(cfg, decom, driver, dnames, decomid, ncid,
                                 name, xtype, nDims, dimids, varp);
        if (decomid >= 0) varp->vlen = decom.raw_nreqs[decomid];
    } else if (cfg.strategy == blob && decomid >= 0) {
        /* use blob dimensions to define blob variables */
        int ival, _ndims;
        if (varp->isRecVar) {
            _ndims = 2;  /* all blob record variables are 2D */
            _dimids = rec_dimids[decomid];
        } else {
            _ndims = 1;  /* all blob fixed-size variables are 1D */
            _dimids = &fix_dimids[decomid];
        }
        err = driver.def_var(ncid, name, xtype, _ndims, _dimids, &varp->vid);
        /* save the canonical dimensions as attributes */
        ival = decomid + 1;
        PUT_ATTR_INT("decomposition_ID", 1, &ival)
        PUT_ATTR_INT("global_dimids", nDims, dimids)
        varp->vlen = decom.count[decomid];
    } else { /* cfg.strategy == canonical or log or decomid == -1 */
        err = driver.def_var(ncid, name, xtype, nDims, dimids, &varp->vid);
        if (decomid >= 0) varp->vlen = decom.count[decomid];
    }
    if (err != 0) {
        printf("Error in %s line %d: def_var %s\n", __FILE__, __LINE__,varp->_name);
        return err;
    }
    /* increment I/O buffer sizes */
    if (varp->isRecVar) {
        if (varp->iType == MPI_DOUBLE)
            wr_buf.rec_dbl_buflen += varp->vlen + wr_buf.gap;
        else if (varp->iType == MPI_INT)
            wr_buf.rec_int_buflen += varp->vlen + wr_buf.gap;
        else if (varp->iType == MPI_CHAR)
            wr_buf.rec_txt_buflen += varp->vlen + wr_buf.gap;
        else if (varp->iType == MPI_FLOAT)
            wr_buf.rec_flt_buflen += varp->vlen + wr_buf.gap;
        else if (varp->iType == MPI_LONG_LONG)
            wr_buf.rec_lld_buflen += varp->vlen + wr_buf.gap;
        else assert(0);
    } else {
        if (varp->iType == MPI_DOUBLE)
            wr_buf.fix_dbl_buflen += varp->vlen + wr_buf.gap;
        else if (varp->iType == MPI_INT)
            wr_buf.fix_int_buflen += varp->vlen + wr_buf.gap;
        else if (varp->iType == MPI_CHAR)
            wr_buf.fix_txt_buflen += varp->vlen + wr_buf.gap;
        else if (varp->iType == MPI_FLOAT)
            wr_buf.fix_flt_buflen += varp->vlen + wr_buf.gap;
        else if (varp->iType == MPI_LONG_LONG)
            wr_buf.fix_lld_buflen += varp->vlen + wr_buf.gap;
        else assert(0);
    }
err_out:
    return err;
}

int e3sm_io_case::inq_var(e3sm_io_config &cfg,
                          e3sm_io_decom  &decom,
                          e3sm_io_driver &driver,
                          case_meta      *cmeta,
                          int             ncid,
                          const char     *name,
                          int             dim_time,
                          int            *dimids,
                          MPI_Datatype    itype,
                          int             decomid,
                          var_meta       *varp)
{
    /* nDims and dimids are canonical dimensions */
    int err=0, j;

    err = driver.inq_varid(ncid, name, &varp->vid);
if (err != 0) printf("name=%s\n",name);
    CHECK_ERR

    varp->_name = strdup(name);

    err = driver.inq_var(ncid, varp->vid, NULL, &varp->xType,
                         &varp->ndims, dimids, NULL);
    CHECK_ERR

    varp->iType     = itype;   /* internal data type of write buffer */
    varp->decomp_id = decomid; /* decomposition map ID */
    varp->isRecVar  = (varp->ndims != 0 && dimids[0] == dim_time);

    /* calculate variable size (1 record only) */
    for (varp->vlen=1, j=0; j<varp->ndims; j++) {
        err = driver.inq_dimlen(ncid, dimids[j], &varp->dims[j]);
        CHECK_ERR
        if (j == 0 && varp->isRecVar) varp->dims[j] = 1;
        varp->vlen *= varp->dims[j];
    }
    if (cfg.api == adios) {
        if (decomid >= 0) varp->vlen = decom.raw_nreqs[decomid];
    } else if (cfg.strategy == blob && decomid >= 0) {
        /* use blob dimensions to define blob variables */
        if (varp->isRecVar)
            assert(varp->ndims == 2); /* all blob record variables are 2D */
        else
            assert(varp->ndims == 1); /* all blob fixed-size variables are 1D */
        /* decomposition map ID */
        GET_ATTR_INT("decomposition_ID", 1, &decomid)
        decomid--; /* start with 0 */
        varp->decomp_id = decomid;
        /* the canonical dimensions as attributes */
        GET_ATTR_INT("global_dimids", varp->ndims, dimids)
        varp->vlen = decom.count[decomid];
    } else { /* cfg.strategy == canonical or log or decomid == -1 */
        if (decomid >= 0) varp->vlen = decom.count[decomid];
    }
    if (err != 0) {
        printf("Error in %s line %d: def_var %s\n", __FILE__, __LINE__,varp->_name);
        return err;
    }
    /* increment I/O buffer sizes */
    if (varp->isRecVar) {
        if (varp->iType == MPI_DOUBLE)
            wr_buf.rec_dbl_buflen += varp->vlen + wr_buf.gap;
        else if (varp->iType == MPI_INT)
            wr_buf.rec_int_buflen += varp->vlen + wr_buf.gap;
        else if (varp->iType == MPI_CHAR)
            wr_buf.rec_txt_buflen += varp->vlen + wr_buf.gap;
        else if (varp->iType == MPI_FLOAT)
            wr_buf.rec_flt_buflen += varp->vlen + wr_buf.gap;
        else if (varp->iType == MPI_LONG_LONG)
            wr_buf.rec_lld_buflen += varp->vlen + wr_buf.gap;
        else assert(0);
    } else {
        if (varp->iType == MPI_DOUBLE)
            wr_buf.fix_dbl_buflen += varp->vlen + wr_buf.gap;
        else if (varp->iType == MPI_INT)
            wr_buf.fix_int_buflen += varp->vlen + wr_buf.gap;
        else if (varp->iType == MPI_CHAR)
            wr_buf.fix_txt_buflen += varp->vlen + wr_buf.gap;
        else if (varp->iType == MPI_FLOAT)
            wr_buf.fix_flt_buflen += varp->vlen + wr_buf.gap;
        else if (varp->iType == MPI_LONG_LONG)
            wr_buf.fix_lld_buflen += varp->vlen + wr_buf.gap;
        else assert(0);
    }
err_out:
    return err;
}

