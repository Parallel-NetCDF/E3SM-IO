#pragma once

#include <string>
//
#include "e3sm_io.h"

#define E3SM_IO_DRIVER_MAX_RANK 32

typedef enum e3sm_io_op_mode { coll, indep, nb, nbe } e3sm_io_op_mode;

class e3sm_io_driver {
   protected:
    e3sm_io_config *cfg;

   public:
    e3sm_io_driver (e3sm_io_config *cfg) : cfg (cfg) {};
    virtual ~e3sm_io_driver (){   };
    virtual int create (std::string path, MPI_Comm comm, MPI_Info info, int *fid) = 0;
    virtual int open (std::string path, MPI_Comm comm, MPI_Info info, int *fid)   = 0;
    virtual int close (int fid)                                                   = 0;
    virtual int inq_file_info (int fid, MPI_Info *info)                           = 0;
    virtual int inq_put_size (int fid, MPI_Offset *size)                          = 0;
    virtual int inq_get_size (int fid, MPI_Offset *size)                          = 0;
    virtual int inq_malloc_size (MPI_Offset *size)                                = 0;
    virtual int inq_malloc_max_size (MPI_Offset *size)                            = 0;
    virtual int inq_rec_size (int fid, MPI_Offset *size)                          = 0;
    virtual int def_var (
        int fid, std::string name, MPI_Datatype type, int ndim, int *dimids, int *did) = 0;
    virtual int def_local_var (
        int fid, std::string name, MPI_Datatype type, int ndim, MPI_Offset *dsize, int *did) = 0;
    virtual int inq_var (int fid, std::string name, int *did)                                = 0;
    virtual int inq_var_off (int fid, int vid, MPI_Offset *off)                              = 0;
    virtual int def_dim (int fid, std::string name, MPI_Offset size, int *dimid)             = 0;
    virtual int inq_dim (int fid, std::string name, int *dimid)                              = 0;
    virtual int inq_dimlen (int fid, int dimid, MPI_Offset *size)                            = 0;
    virtual int enddef (int fid)                                                             = 0;
    virtual int redef (int fid)                                                              = 0;
    virtual int wait (int fid)                                                               = 0;
    virtual int put_att (
        int fid, int vid, std::string name, MPI_Datatype type, MPI_Offset size, void *buf)      = 0;
    virtual int get_att (int fid, int vid, std::string name, void *buf)                         = 0;
    virtual int put_varl (int fid, int vid, MPI_Datatype type, void *buf, e3sm_io_op_mode mode) = 0;
    virtual int put_vara (int fid,
                          int vid,
                          MPI_Datatype type,
                          MPI_Offset *start,
                          MPI_Offset *count,
                          void *buf,
                          e3sm_io_op_mode mode)                                                 = 0;
    virtual int put_vars (int fid,
                          int vid,
                          MPI_Datatype type,
                          MPI_Offset *start,
                          MPI_Offset *count,
                          MPI_Offset *stride,
                          void *buf,
                          e3sm_io_op_mode mode)                                                 = 0;
    virtual int put_varn (int fid,
                          int vid,
                          MPI_Datatype type,
                          int nreq,
                          MPI_Offset **starts,
                          MPI_Offset **counts,
                          void *buf,
                          e3sm_io_op_mode mode)                                                 = 0;
    virtual int put_vard (int fid,
                          int vid,
                          MPI_Datatype type,
                          MPI_Offset nelem,
                          MPI_Datatype ftype,
                          void *buf,
                          e3sm_io_op_mode mode)                                                 = 0;
    virtual int get_vara (int fid,
                          int vid,
                          MPI_Datatype type,
                          MPI_Offset *start,
                          MPI_Offset *count,
                          void *buf,
                          e3sm_io_op_mode mode)                                                 = 0;
    virtual int get_vars (int fid,
                          int vid,
                          MPI_Datatype type,
                          MPI_Offset *start,
                          MPI_Offset *count,
                          MPI_Offset *stride,
                          void *buf,
                          e3sm_io_op_mode mode)                                                 = 0;
    virtual int get_varn (int fid,
                          int vid,
                          MPI_Datatype type,
                          int nreq,
                          MPI_Offset **starts,
                          MPI_Offset **counts,
                          void *buf,
                          e3sm_io_op_mode mode)                                                 = 0;
    virtual int get_vard (int fid,
                          int vid,
                          MPI_Datatype type,
                          MPI_Offset nelem,
                          MPI_Datatype ftype,
                          void *buf,
                          e3sm_io_op_mode mode)                                                 = 0;
};