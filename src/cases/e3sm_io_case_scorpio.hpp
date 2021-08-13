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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//
#include <cstring>
#include <map>
//
#include <e3sm_io.h>
#include <e3sm_io_err.h>

#include <e3sm_io_case.hpp>
#include <e3sm_io_driver.hpp>
#include <e3sm_io_driver_adios2.hpp>
#include <e3sm_io_driver_pnc.hpp>

typedef struct {
    int vid;        /* variable ID (ADIOS or NetCDF) */

    int frame_id;    /* ADIOS ID */
    int fillval_id;  /* ADIOS ID */
    int decom_id;    /* ADIOS ID of decomposition map variable */
    int piodecomid;  /* map IDs used on Scorpio starting at 512 */
    int64_t dims[3]; /* dimension sizes */
    int ndims;

    int decomp_id;      /* decomposition map ID */
    int isRecVar;       /* whether is a record variable */
    size_t vlen;        /* length to be written by this rank */
    MPI_Datatype itype; /* memory buffer of internal data type */
} var_meta_scorpio;


extern
int e3sm_io_scorpio_define_dim (e3sm_io_driver &driver,
                                   int fid,
                                   std::string name,
                                   MPI_Offset size,
                                   std::map<int, std::string> &dnames,
                                   int *did);

extern
int e3sm_io_scorpio_define_var (e3sm_io_driver &driver,
                                   e3sm_io_config &cfg,
                                   std::map<int, std::string> &dnames,
                                   e3sm_io_decom &decom,
                                   int decomid,
                                   int fid,
                                   std::string name,
                                   nc_type xtype,
                                   int ndims,
                                   int *dimids,
                                   var_meta_scorpio *var);

extern
int e3sm_io_scorpio_write_var (e3sm_io_driver &driver,
                                  int frameid,
                                  int fid,
                                  var_meta_scorpio &var,
                                  MPI_Datatype itype,
                                  void *buf,
                                  e3sm_io_op_mode mode);

extern
int e3sm_io_scorpio_put_att (e3sm_io_driver &driver,
                                int fid,
                                int vid,
                                std::string name,
                                nc_type xtype,
                                MPI_Offset size,
                                void *buf);

extern
int e3sm_io_scorpio_put_att (e3sm_io_driver &driver,
                                int fid,
                                var_meta_scorpio &var,
                                std::string name,
                                nc_type xtype,
                                MPI_Offset size,
                                void *buf);

extern
int def_F_case_scorpio(e3sm_io_driver   &driver,
                       e3sm_io_config   &cfg,
                       e3sm_io_decom    &decom,
                       int               ncid,
                       var_meta_scorpio *vars,
                       int              *scorpiovars,
                       io_buffers       *wr_buf);

extern
int def_G_case_scorpio(e3sm_io_config   &cfg,
                       e3sm_io_decom    &decom,
                       e3sm_io_driver   &driver,
                       int               ncid, /* file ID */
                       var_meta_scorpio *vars, /* variable IDs */
                       int              *scorpiovars,
                       io_buffers       *wr_buf);

extern
int def_I_case_scorpio(e3sm_io_config   &cfg,
                       e3sm_io_decom    &decom,
                       e3sm_io_driver   &driver,
                       int               ncid,    /* file ID */
                       var_meta_scorpio *vars,    /* variable metadata */
                       int              *scorpiovars,
                       io_buffers       *wr_buf);

extern
int var_wr_all_cases_scorpio(e3sm_io_config &cfg,
                             e3sm_io_decom  &decom,
                             e3sm_io_driver &driver,
                             case_meta      *cmeta);
