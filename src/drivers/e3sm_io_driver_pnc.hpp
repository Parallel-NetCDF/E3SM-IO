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
#include <map>
#include <vector>
//
#include <mpi.h>
//
#include <e3sm_io.h>

#include <e3sm_io_driver.hpp>

inline nc_type mpitype2nctype (MPI_Datatype type) {
    switch (type) {
        case MPI_INT:
            return NC_INT;
        case MPI_FLOAT:
            return NC_FLOAT;
        case MPI_DOUBLE:
            return NC_DOUBLE;
        case MPI_CHAR:
            return NC_CHAR;
        default:
            throw "Unsupported datatype";
    }
}

class e3sm_io_driver_pnc : public e3sm_io_driver {
    MPI_Comm comm;
    MPI_Info info;
    std::vector<MPI_Offset> var_nelems;
    std::vector<int> var_ndims;
    std::vector<MPI_Offset> dim_lens;
    std::map<int, MPI_Info> file_infos;

   public:
    e3sm_io_driver_pnc (e3sm_io_config *cfg);
    ~e3sm_io_driver_pnc ();
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
    int def_local_var (
        int fid, std::string name, MPI_Datatype type, int ndim, MPI_Offset *dsize, int *did);
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
    int put_varl (int fid, int vid, MPI_Datatype type, void *buf, e3sm_io_op_mode mode);
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