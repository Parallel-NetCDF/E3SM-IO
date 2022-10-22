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
#include <adios2_c.h>
#include <mpi.h>
//
#include <e3sm_io.h>

#include <e3sm_io_driver.hpp>

inline adios2_type mpi_type_to_adios2_type (MPI_Datatype type) {
    if (type == MPI_DOUBLE){
        return adios2_type_double;
    }
    else if (type == MPI_FLOAT){
        return adios2_type_float;
    }
    else if (type == MPI_INT){
        return adios2_type_int32_t;
    }
    else if (type == MPI_LONG_LONG){
        return adios2_type_int64_t;
    }
    else if (type == MPI_CHAR){
        return adios2_type_int8_t;
    }
    else if (type == MPI_WCHAR){
        return adios2_type_string;
    }
    else if (type == MPI_BYTE){
        return adios2_type_uint8_t;
    }

    int name_len;
    char type_name[MPI_MAX_OBJECT_NAME];
    MPI_Type_get_name(type, type_name, &name_len);
    printf ("Error at line %d in %s: Unknown MPI Datatype %s\n", __LINE__, __FILE__, type_name);
    DEBUG_ABORT

    return adios2_type_unknown;
}

/*----< e3sm_io_type_nc2adios() >---------------------------------------------*/
inline adios2_type
e3sm_io_type_nc2adios(nc_type xtype)
{
    switch(xtype) {
        case NC_BYTE:
        case NC_CHAR:   return adios2_type_int8_t;
        case NC_UBYTE:  return adios2_type_uint8_t;  
        case NC_INT:    return adios2_type_int32_t;
        case NC_FLOAT : return adios2_type_float;
        case NC_DOUBLE: return adios2_type_double;
        case NC_INT64:  return adios2_type_int64_t;
        case NC_STRING: return adios2_type_string;
        case NC_SHORT:
        case NC_USHORT:
        case NC_UINT:
        case NC_UINT64:
        default: return adios2_type_unknown;
    }
}

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
        std::vector<adios2_variable *> ddids;
        std::map<std::string, int> dimmap;
        MPI_Offset recsize = 0;
        adios2_operator *op;
        MPI_Offset putsize = 0;
        MPI_Offset getsize = 0;
        bool wr;
        int rank;
    } adios2_file;
    std::vector<adios2_file *> files;

   public:
    static bool compatible (std::string path);
    e3sm_io_driver_adios2 (e3sm_io_config *cfg);
    ~e3sm_io_driver_adios2 ();
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
    int inq_var_name(int ncid, int varid, char *name);
    int inq_var_off (int fid, int vid, MPI_Offset *off);
    int def_dim (int fid, std::string name, MPI_Offset size, int *dimid);
    int inq_dim (int fid, std::string name, int *dimid);
    int inq_dimlen (int fid, int dimid, MPI_Offset *size);
    int enddef (int fid);
    int redef (int fid);
    int wait (int fid);
    int put_att (int fid, int vid, std::string name, nc_type xtype, MPI_Offset size, const void *buf);
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
