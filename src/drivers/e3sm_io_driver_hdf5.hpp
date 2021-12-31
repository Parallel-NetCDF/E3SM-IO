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

class e3sm_io_driver_hdf5 : public e3sm_io_driver {
    typedef struct Index_order {
        hsize_t index;
        hsize_t coverage;
        char *data;
    } Index_order;

   protected:
    class hdf5_file {
       public:
        MPI_Comm comm;
        e3sm_io_driver_hdf5 &driver;
        e3sm_io_config *cfg;
        hid_t id;
        std::vector<hid_t> dids; // HDF5 dataset IDs
        std::vector<hsize_t> dsizes; // Size of dimensions
        std::map<hid_t, int> inv_dids;  // Inverse of dids
        MPI_Offset recsize = 0;
        MPI_Offset putsize = 0;
        MPI_Offset getsize = 0;
        int rank;

#ifndef HDF5_HAVE_DWRITE_MULTI
        typedef struct H5D_rw_multi_t {
            hid_t dset_id;       /* dataset ID */
            hid_t dset_space_id; /* dataset selection dataspace ID */
            hid_t mem_type_id;   /* memory datatype ID */
            hid_t mem_space_id;  /* memory selection dataspace ID */
            void *buf;           /* pointer to data buffer */
        } H5D_rw_multi_t;
#endif
        std::vector<H5D_rw_multi_t> wreqs;
        std::vector<H5D_rw_multi_t> rreqs;

        std::vector<std::vector<e3sm_io_driver_hdf5::Index_order> > dataset_segments;

        std::vector<hid_t> memspace_recycle;
        std::vector<hid_t> dataspace_recycle;

        hdf5_file (e3sm_io_driver_hdf5 &x) : driver (x), cfg{x.cfg} {};

// #ifdef HDF5_HAVE_DWRITE_MULTI
        herr_t register_multidataset (
            void *buf, hid_t did, hid_t dsid, hid_t msid, hid_t mtype, int write);
        herr_t pull_multidatasets ();
        int flush_multidatasets ();
// #endif
    };

    std::vector<hdf5_file *> files;
    hid_t dxplid_coll;
    hid_t dxplid_indep;
    hid_t dxplid_coll_nb;
    hid_t dxplid_indep_nb;

    hsize_t one[H5S_MAX_RANK];

   private:
    // Config
    bool use_dwrite_multi = false;
    bool merge_varn       = false;

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
    int inq_var (int fid, std::string name, int *did);
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
    int put_varl (int fid, int vid, MPI_Datatype itype, void *buf, e3sm_io_op_mode mode);
    int put_vara (int fid,
                  int vid,
                  MPI_Datatype itype,
                  MPI_Offset *start,
                  MPI_Offset *count,
                  void *buf,
                  e3sm_io_op_mode mode);
    int put_vars (int fid,
                  int vid,
                  MPI_Datatype itype,
                  MPI_Offset *start,
                  MPI_Offset *count,
                  MPI_Offset *stride,
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
    int get_vars (int fid,
                  int vid,
                  MPI_Datatype itype,
                  MPI_Offset *start,
                  MPI_Offset *count,
                  MPI_Offset *stride,
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
    static int index_order_cmp (const void *a, const void *b);
    int pack_data (Index_order *index_order,
                   int *index,
                   char *src,
                   hsize_t esize,
                   int ndim,
                   const hsize_t *dims,
                   const hsize_t *start,
                   const hsize_t *block);
    int copy_index_buf (Index_order *index_order, int total_blocks, char *out_buf);

    int put_varn_expand (int fid,
                         int vid,
                         MPI_Datatype itype,
                         int nreq,
                         MPI_Offset **starts,
                         MPI_Offset **counts,
                         void *buf,
                         e3sm_io_op_mode mode);

    int get_varn_expand (int fid,
                         int vid,
                         MPI_Datatype itype,
                         int nreq,
                         MPI_Offset **starts,
                         MPI_Offset **counts,
                         void *buf,
                         e3sm_io_op_mode mode);
    int put_varn_merge (int fid,
                        int vid,
                        MPI_Datatype itype,
                        int nreq,
                        MPI_Offset **starts,
                        MPI_Offset **counts,
                        void *buf,
                        e3sm_io_op_mode mode);

    int get_varn_merge (int fid,
                        int vid,
                        MPI_Datatype itype,
                        int nreq,
                        MPI_Offset **starts,
                        MPI_Offset **counts,
                        void *buf,
                        e3sm_io_op_mode mode);
};

/*----< e3sm_io_type_nc2hdf5() >---------------------------------------------*/
inline hid_t
e3sm_io_type_nc2hdf5(nc_type xtype)
{
    switch(xtype) {
        case NC_BYTE:   return H5T_NATIVE_UINT8;
        case NC_UBYTE:  return H5T_NATIVE_UCHAR;
        case NC_CHAR:   return H5T_NATIVE_CHAR;
        case NC_SHORT:  return H5T_NATIVE_SHORT;
        case NC_USHORT: return H5T_NATIVE_USHORT;
        case NC_INT:    return H5T_NATIVE_INT;
        case NC_UINT:   return H5T_NATIVE_UINT;
        case NC_FLOAT : return H5T_NATIVE_FLOAT;
        case NC_DOUBLE: return H5T_NATIVE_DOUBLE;
        case NC_INT64:  return H5T_NATIVE_INT64;
        case NC_UINT64: return H5T_NATIVE_UINT64;
        default: return -1;
    }
}

