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
#include <hdf5.h>
#include <mpi.h>
//
#include <e3sm_io.h>

#include <e3sm_io_driver.hpp>

class e3sm_io_driver_hdf5_log : public e3sm_io_driver_hdf5 {      
    // Config
    bool use_logvol_varn  = true;
    bool merge_varn       = false;
    int num_subfiles      = 0;

   public:
    e3sm_io_driver_hdf5_log (e3sm_io_config *cfg);
    ~e3sm_io_driver_hdf5_log ();
    int create (std::string path, MPI_Comm comm, MPI_Info info, int *fid) override;
    int open (std::string path, MPI_Comm comm, MPI_Info info, int *fid) override;
    int inq_file_size (std::string path, MPI_Offset *size) override;

    int put_vara (int fid,
                  int vid,
                  MPI_Datatype itype,
                  MPI_Offset *start,
                  MPI_Offset *count,
                  void *buf,
                  e3sm_io_op_mode mode) override;
    int put_varn (int fid,
                  int vid,
                  MPI_Datatype itype,
                  int nreq,
                  MPI_Offset **starts,
                  MPI_Offset **counts,
                  void *buf,
                  e3sm_io_op_mode mode) override;
    int get_vara (int fid,
                  int vid,
                  MPI_Datatype itype,
                  MPI_Offset *start,
                  MPI_Offset *count,
                  void *buf,
                  e3sm_io_op_mode mode) override;
    int get_varn (int fid,
                  int vid,
                  MPI_Datatype itype,
                  int nreq,
                  MPI_Offset **starts,
                  MPI_Offset **counts,
                  void *buf,
                  e3sm_io_op_mode mode) override;
};


