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

class e3sm_io_driver_nc4 : public e3sm_io_driver {
	std::map<std::pair<int, int>, MPI_Offset> var_size;

   public:
	static bool compatible (std::string path);
	e3sm_io_driver_nc4 (e3sm_io_config *cfg);
	~e3sm_io_driver_nc4 ();
	int create (std::string path, MPI_Comm comm, MPI_Info info, int *fid);
	int open (std::string path, MPI_Comm comm, MPI_Info info, int *fid);
	int close (int fid);
	int inq_file_info (int fid, MPI_Info *info);
	int inq_file_size (std::string path, MPI_Offset *size);
	int inq_put_size (MPI_Offset *size);
	int inq_get_size (MPI_Offset *size);
	int inq_malloc_size (MPI_Offset *size);
	int inq_malloc_max_size (MPI_Offset *size);
	int inq_rec_size (int fid, MPI_Offset *size);
	int expand_rec_size (int fid, MPI_Offset size);
	int def_var (int fid, std::string name, nc_type xtype, int ndim, int *dimids, int *did);
	int def_local_var (
		int fid, std::string name, nc_type xtype, int ndim, MPI_Offset *dsize, int *did);
        int inq_varid(int fid, const char *name, int *did);
        int inq_var(int fid, int varid, char *name, nc_type *xtypep, int *ndimsp,
                    int *dimids, int *nattsp);
	int inq_var_name (int ncid, int varid, char *name);
	int inq_var_off (int fid, int vid, MPI_Offset *off);
	int def_dim (int fid, std::string name, MPI_Offset size, int *dimid);
	int inq_dim (int fid, std::string name, int *dimid);
	int inq_dimlen (int fid, int dimid, MPI_Offset *size);
	int enddef (int fid);
	int redef (int fid);
	int wait (int fid);
	int put_att (
		int fid, int vid, std::string name, nc_type xtype, MPI_Offset size, const void *buf);
	int get_att (int fid, int vid, std::string name, void *buf);
	int inq_att (int fid, int vid, std::string name, MPI_Offset *size);
	int put_varl (int fid, int vid, MPI_Datatype itype, void *buf, e3sm_io_op_mode mode);
	int put_vara (int fid,
				  int vid,
				  MPI_Datatype itype,
				  MPI_Offset *start,
				  MPI_Offset *count,
				  void *buf,
				  e3sm_io_op_mode mode);
	int put_varn (int fid,
				  int vid,
				  MPI_Datatype itype,
				  int nreq,
				  MPI_Offset **starts,
				  MPI_Offset **counts,
				  void *buf,
				  e3sm_io_op_mode mode);
	int get_vara (int fid,
				  int vid,
				  MPI_Datatype itype,
				  MPI_Offset *start,
				  MPI_Offset *count,
				  void *buf,
				  e3sm_io_op_mode mode);
	int get_varn (int fid,
				  int vid,
				  MPI_Datatype itype,
				  int nreq,
				  MPI_Offset **starts,
				  MPI_Offset **counts,
				  void *buf,
				  e3sm_io_op_mode mode);
};
