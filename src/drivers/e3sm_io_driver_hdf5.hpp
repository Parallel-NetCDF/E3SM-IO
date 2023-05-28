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


#define CHECK_HERR {                                              \
    if (herr < 0) {                                               \
        printf ("Error at line %d in %s:\n", __LINE__, __FILE__); \
        H5Eprint1 (stdout);                                       \
        DEBUG_ABORT;                                              \
        err = -1;                                                 \
        goto err_out;                                             \
    }                                                             \
}

#define CHECK_HID(A) {                                            \
    if (A < 0) {                                                  \
        printf ("Error at line %d in %s:\n", __LINE__, __FILE__); \
        H5Eprint1 (stdout);                                       \
        DEBUG_ABORT;                                              \
        err = -1;                                                 \
        goto err_out;                                             \
    }                                                             \
}

class e3sm_io_driver_hdf5 : public e3sm_io_driver {
   protected:
    class hdf5_file {
       public:
        MPI_Comm comm;
        int rank;

        e3sm_io_driver_hdf5 &driver;
        e3sm_io_config *cfg;
        hid_t id;                    // HDF5 file ID
        std::vector<hid_t> dids;     // HDF5 dataset IDs
        std::vector<hsize_t> dsizes; // Size of dimensions
        MPI_Offset recsize = 0;
        MPI_Offset putsize = 0;
        MPI_Offset getsize = 0;
        MPI_Info info_used = MPI_INFO_NULL;

        std::vector<bool>  dset_isRec;      /* whether a dataset is a record variable */
        std::vector<hid_t> wdset_ids;       /* dataset ID */
        std::vector<hid_t> wmem_type_ids;   /* memory datatype ID */
        std::vector<hid_t> wmem_space_ids;  /* memory selection dataspace ID */
        std::vector<hid_t> wdset_space_ids; /* dataset selection dataspace ID */
        std::vector<void*> wbufs;           /* pointer to data buffer */

        std::vector<hid_t> rdset_ids;       /* dataset ID */
        std::vector<hid_t> rmem_type_ids;   /* memory datatype ID */
        std::vector<hid_t> rmem_space_ids;  /* memory selection dataspace ID */
        std::vector<hid_t> rdset_space_ids; /* dataset selection dataspace ID */
        std::vector<void*> rbufs;           /* pointer to data buffer */

        std::vector<hid_t> memspace_recycle;
        std::vector<hid_t> dataspace_recycle;

        hdf5_file (e3sm_io_driver_hdf5 &x) : driver (x), cfg{x.cfg} {};

        int flush_multidatasets ();
    };

    hid_t log_vlid;

    std::vector<hdf5_file *> files;
    hid_t dxplid_coll;
    hid_t dxplid_indep;
    hid_t dxplid_coll_nb;
    hid_t dxplid_indep_nb;

    hsize_t one[H5S_MAX_RANK];
    int nfixVars; /* number of fixed-size variables posted in varn calls */

   private:
    // Config
    bool use_dwrite_multi = false;

   public:
    e3sm_io_driver_hdf5 (e3sm_io_config *cfg);
    ~e3sm_io_driver_hdf5 ();
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

   private:
    int varn_expand (int fid,
                     int vid,
                     MPI_Datatype itype,
                     int nreq,
                     MPI_Offset **starts,
                     MPI_Offset **counts,
                     void *buf,
                     bool isWrite);
    int post_varn (int fid,
                   int vid,
                   MPI_Datatype itype,
                   int nreq,
                   MPI_Offset **starts,
                   MPI_Offset **counts,
                   void *buf,
                   bool mode);
};

extern hid_t e3sm_io_type_mpi2hdf5(MPI_Datatype itype);
extern hid_t e3sm_io_type_nc2hdf5(nc_type xtype);

