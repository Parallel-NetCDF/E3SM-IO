#pragma once
#include <config.h>
//
#include <map>
#include <vector>
//
#include <adios2_c.h>
#include <mpi.h>
//
#include "e3sm_io.h"
#include "e3sm_io_driver.hpp"

class e3sm_io_driver_adios2 : public e3sm_io_driver {
    typedef struct adios2_file {
        std::string path;
        adios2_adios *adp;
        adios2_io *iop;
        adios2_engine *ep;
        int read;
        std::vector<adios2_variable *> dids;
        std::vector<int> ndims;
        std::vector<size_t> dsizes;
        MPI_Offset recsize = 0;
    } adios2_file;
    std::vector<adios2_file *> files;

    double tsel, twrite, tread, text;

   public:
    e3sm_io_driver_adios2 ();
    ~e3sm_io_driver_adios2 ();
    int create (std::string path, MPI_Comm comm, MPI_Info info, int *fid);
    int open (std::string path, MPI_Comm comm, MPI_Info info, int *fid);
    int close (int fid);
    int inq_file_info (int fid, MPI_Info *info);
    int inq_put_size (int fid, MPI_Offset *size);
    int inq_get_size (int fid, MPI_Offset *size);
    int inq_malloc_size (MPI_Offset *size);
    int inq_malloc_max_size (MPI_Offset *size);
    int inq_rec_size (int fid, MPI_Offset *size);
    int def_var (int fid, std::string name, MPI_Datatype type, int ndim, int *dimids, int *did);
    int inq_var (int fid, std::string name, int *did);
    int inq_var_off (int fid, int vid, MPI_Offset *off);
    int def_dim (int fid, std::string name, MPI_Offset size, int *dimid);
    int inq_dim (int fid, std::string name, int *dimid);
    int inq_dimlen (int fid, int dimid, MPI_Offset *size);
    int enddef (int fid);
    int redef (int fid);
    int wait (int fid);
    int put_att (int fid, int vid, std::string name, MPI_Datatype type, MPI_Offset size, void *buf);
    int get_att (int fid, int vid, std::string name, void *buf);
    int put_vara (int fid,
                  int vid,
                  MPI_Datatype type,
                  MPI_Offset *start,
                  MPI_Offset *count,
                  void *buf,
                  e3sm_io_op_mode mode);
    int put_vars (int fid,
                  int vid,
                  MPI_Datatype type,
                  MPI_Offset *start,
                  MPI_Offset *count,
                  MPI_Offset *stride,
                  void *buf,
                  e3sm_io_op_mode mode);
    int put_varn (int fid,
                  int vid,
                  MPI_Datatype type,
                  int nreq,
                  MPI_Offset **starts,
                  MPI_Offset **counts,
                  void *buf,
                  e3sm_io_op_mode mode);
    int put_vard (int fid,
                  int vid,
                  MPI_Datatype type,
                  MPI_Offset nelem,
                  MPI_Datatype ftype,
                  void *buf,
                  e3sm_io_op_mode mode);
    int get_vara (int fid,
                  int vid,
                  MPI_Datatype type,
                  MPI_Offset *start,
                  MPI_Offset *count,
                  void *buf,
                  e3sm_io_op_mode mode);
    int get_vars (int fid,
                  int vid,
                  MPI_Datatype type,
                  MPI_Offset *start,
                  MPI_Offset *count,
                  MPI_Offset *stride,
                  void *buf,
                  e3sm_io_op_mode mode);
    int get_varn (int fid,
                  int vid,
                  MPI_Datatype type,
                  int nreq,
                  MPI_Offset **starts,
                  MPI_Offset **counts,
                  void *buf,
                  e3sm_io_op_mode mode);
    int get_vard (int fid,
                  int vid,
                  MPI_Datatype type,
                  MPI_Offset nelem,
                  MPI_Datatype ftype,
                  void *buf,
                  e3sm_io_op_mode mode);
};
