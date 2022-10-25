#pragma once

#include <string>
//
#include "e3sm_io.h"

#define E3SM_IO_DRIVER_MAX_RANK 32

typedef enum e3sm_io_op_mode { coll, indep, nb, nbe } e3sm_io_op_mode;

class e3sm_io_driver {
   protected:
    MPI_Offset amount_WR;
    MPI_Offset amount_RD;

   public:
    e3sm_io_config *cfg;

    e3sm_io_driver (e3sm_io_config *cfg) : amount_WR(0), amount_RD(0), cfg(cfg) {};
    virtual ~e3sm_io_driver () {};
    virtual int create (std::string path, MPI_Comm comm, MPI_Info info, int *fid) = 0;
    virtual int open (std::string path, MPI_Comm comm, MPI_Info info, int *fid)   = 0;
    virtual int close (int fid)                                                   = 0;
    virtual int inq_file_size (std::string path, MPI_Offset *size)                = 0;
    virtual int inq_file_info (int fid, MPI_Info *info)                           = 0;
    virtual int inq_put_size (MPI_Offset *size)                          = 0;
    virtual int inq_get_size (MPI_Offset *size)                          = 0;
    virtual int inq_malloc_size (MPI_Offset *size)                                = 0;
    virtual int inq_malloc_max_size (MPI_Offset *size)                            = 0;
    virtual int inq_rec_size (int fid, MPI_Offset *size)                          = 0;
    virtual int expand_rec_size (int fid, MPI_Offset size)                        = 0;
    virtual int def_var (
        int fid, std::string name, nc_type xtype, int ndim, int *dimids, int *did) = 0;
    virtual int def_local_var (
        int fid, std::string name, nc_type xtype, int ndim, MPI_Offset *dsize, int *did) = 0;
    virtual int inq_varid(int fid, const char *name, int *did) = 0;
    virtual int inq_var(int fid, int varid, char *name, nc_type *xtypep, int *ndimsp,
                        int *dimids, int *nattsp) = 0;
    virtual int inq_var_name (int fid, int did, char *name)                                  = 0;
    virtual int inq_var_off (int fid, int vid, MPI_Offset *off)                              = 0;
    virtual int def_dim (int fid, std::string name, MPI_Offset size, int *dimid)             = 0;
    virtual int inq_dim (int fid, std::string name, int *dimid)                              = 0;
    virtual int inq_dimlen (int fid, int dimid, MPI_Offset *size)                            = 0;
    virtual int enddef (int fid)                                                             = 0;
    virtual int redef (int fid)                                                              = 0;
    virtual int wait (int fid)                                                               = 0;
    virtual int put_att (
        int fid, int vid, std::string name, nc_type xtype, MPI_Offset size, const void *buf)= 0;
    virtual int get_att (int fid, int vid, std::string name, void *buf)                         = 0;
    virtual int inq_att (int fid, int vid, std::string name, MPI_Offset *size)               = 0;
    virtual int put_varl (int fid, int vid, MPI_Datatype itype, void *buf, e3sm_io_op_mode mode) = 0;
    virtual int put_vara (int fid,
                          int vid,
                          MPI_Datatype itype,
                          MPI_Offset *start,
                          MPI_Offset *count,
                          void *buf,
                          e3sm_io_op_mode mode)                                                 = 0;
    virtual int put_varn (int fid,
                          int vid,
                          MPI_Datatype itype,
                          int nreq,
                          MPI_Offset **starts,
                          MPI_Offset **counts,
                          void *buf,
                          e3sm_io_op_mode mode)                                                 = 0;
    virtual int get_vara (int fid,
                          int vid,
                          MPI_Datatype itype,
                          MPI_Offset *start,
                          MPI_Offset *count,
                          void *buf,
                          e3sm_io_op_mode mode)                                                 = 0;
    virtual int get_varn (int fid,
                          int vid,
                          MPI_Datatype itype,
                          int nreq,
                          MPI_Offset **starts,
                          MPI_Offset **counts,
                          void *buf,
                          e3sm_io_op_mode mode)                                                 = 0;
};

e3sm_io_driver *e3sm_io_get_driver (const char *filename, e3sm_io_config *cfg);
