/*********************************************************************
 *
 * Copyright (C) 2022, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 *********************************************************************/
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstring>
//
#include <string.h> /* strdup() */
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
//
#include "e3sm_io.h"
#include "e3sm_io_driver.hpp"
#include "e3sm_io_err.h"

#ifdef ENABLE_HDF5
#include <hdf5.h>
#include <e3sm_io_driver_hdf5.hpp>
#include <e3sm_io_driver_h5blob.hpp>
#ifdef ENABLE_LOGVOL
#include <e3sm_io_driver_hdf5_log.hpp>
#endif
#endif

#ifdef ENABLE_PNC
#include <e3sm_io_driver_pnc.hpp>
#endif

#ifdef ENABLE_NETCDF4
#include <e3sm_io_driver_nc4.hpp>
#endif

#ifdef ENABLE_ADIOS2
#include <adios2_c.h>

#include <e3sm_io_driver_adios2.hpp>
#endif

/* File system types recognized by ROMIO in MPICH 4.0.0 */
static const char* fstypes[] = {"ufs", "nfs", "xfs", "pvfs2", "gpfs", "panfs", "lustre", "daos", "testfs", "ime", "quobyte", NULL};

/* Return a pointer to filename by removing the file system type prefix name if
 * there is any.  For example, when filename = "lustre:/home/foo/testfile.nc",
 * remove "lustre:" to return a pointer to "/home/foo/testfile.nc", so the name
 * can be used in POSIX open() calls.
 */
char* remove_file_system_type_prefix(const char *filename)
{
    char *prefix, *colon, *ret_filename;

    if (filename == NULL) return NULL;

    ret_filename = (char*)filename;
    prefix = strdup(filename);

    colon = strchr(prefix, ':');
    if (colon != NULL) { /* there is a prefix end with : */
        int i=0;
        *colon = '\0';
        /* check if prefix is one of recognized file system types */
        while (fstypes[i] != NULL) {
            if (!strcmp(prefix, fstypes[i])) { /* found */
                ret_filename += colon - prefix + 1;
                break;
            }
            i++;
        }
    }
    free(prefix);

    return ret_filename;
}

e3sm_io_driver *e3sm_io_get_driver (const char *filename, /* NULL is for read */
                                    e3sm_io_config *cfg)
{
    int err=0, fd;
    const char *cdf_signature="CDF", *path;
    char signature[8];
    ssize_t rlen;
    e3sm_io_driver *driver = NULL;
#ifdef ENABLE_HDF5
    const char *hdf5_signature = "\211HDF\r\n\032\n";
    off_t offset;
#endif

    /* root checks the file format and broadcasts the finding */
    if (filename && cfg->rank == 0) {
        /* remove the file system type prefix name if there is any.  For example,
         * when filename = "lustre:/home/foo/testfile.nc", remove "lustre:" to make
         * path pointing to "/home/foo/testfile.nc", so it can be used in POSIX
         * open() below
         */
        path = remove_file_system_type_prefix(filename);

        /* must include config.h on 32-bit machines, as AC_SYS_LARGEFILE is
         * called at the configure time and it defines _FILE_OFFSET_BITS to 64
         * if large file feature is supported.
         */
        if ((fd = open(path, O_RDONLY, 00400)) == -1) { /* open for read */
            char msg[128];
            sprintf(msg, "Cannot open file \"%s\"", path);
            ERR_OUT(msg)
        }

        /* get first 8 bytes of file */
        rlen = read (fd, signature, 8);
        if (rlen != 8) { ERR_OUT ("Cannot read file.") }

        // PnetCDF ?
        if (memcmp (signature, cdf_signature, 3) == 0) {
            cfg->api = pnetcdf;
            goto done_check;
        }

        // HDF5 ?
#ifdef ENABLE_HDF5
        offset = 512;
        while (rlen == 8 && memcmp (signature, hdf5_signature, 8)) {
            lseek (fd, offset, SEEK_SET);
            offset <<= 1;
            rlen = read (fd, signature, 8);
        }
        if (rlen == 8) {  // It is HDF5, check if it is Log VOL
            hid_t fid = -1;
            htri_t isnc, islog;

            hid_t faplid = H5Pcreate (H5P_FILE_ACCESS);
            H5Pset_fapl_mpio(faplid, MPI_COMM_SELF, MPI_INFO_NULL);
            fid = H5Fopen (path, H5F_ACC_RDONLY, faplid);
            if (fid < 0) { ERR_OUT ("HDF5 header detected, but not a HDF5 file"); }

            // Check for NetCDF4
            isnc = H5Aexists(fid, "_NCProperties");
            
            if(isnc == false){
                islog = H5Lexists(fid, "_LOG", H5P_DEFAULT);
                if (islog == true) {
                    cfg->api = hdf5_log;
                } else {
                    cfg->api = hdf5;
                }
            }
            if (fid >= 0) { H5Fclose (fid); }
            if (faplid >= 0) { H5Pclose(faplid); }

            if (isnc == false) {
                goto done_check;
            }
        }
#endif

        // NetCDF4 ?
        // NC4 will try to open HDF5 files even if it is not nc4 file, check HDF5 files first
#ifdef ENABLE_NETCDF4
        if(e3sm_io_driver_nc4::compatible(std::string(path))) {
            cfg->api = netcdf4;
            goto done_check;
        }
#endif

// ADIOS2?
#ifdef ENABLE_ADIOS2
        if(e3sm_io_driver_adios2::compatible(std::string(path))) {
            cfg->api = adios;
            goto done_check;
        }
#endif


done_check:
        if (cfg->rank == 0 && fd >= 0) { close (fd); }
    }
    if (filename) {
        err = MPI_Bcast (&(cfg->api), 1, MPI_INT, 0, cfg->io_comm);
        CHECK_MPIERR
    }

    /* For read tests, the read driver is determined by the format of input
     * file. However, if the input file is an HDF5 in canonical layout, the
     * read driver will be selected based on command-line option -a, because
     * the file can be read either by the official release of HDF5 library or
     * the develop branch that implements the multi-dataset APIs, not yet in
     * the official release.
     */
    switch (cfg->api) {
        case pnetcdf:
#ifdef ENABLE_PNC
            driver = new e3sm_io_driver_pnc (cfg);
            break;
#else
            ERR_OUT ("PnetCDF support was not enabled in this build")
#endif
        case netcdf4:
#ifdef ENABLE_NETCDF4
            driver = new e3sm_io_driver_nc4 (cfg);
            break;
#else
            ERR_OUT ("NetCDF4 support was not enabled in this build")
#endif
        case hdf5:
#ifdef ENABLE_HDF5
            if (cfg->strategy == blob)
                driver = new e3sm_io_driver_h5blob (cfg);
            else
                driver = new e3sm_io_driver_hdf5 (cfg);
#else
            ERR_OUT ("HDF5 support was not enabled in this build")
#endif
            break;
        case hdf5_md:
#ifdef HDF5_HAVE_MULTI_DATASET_API
            driver = new e3sm_io_driver_hdf5 (cfg);
            break;
#else
            ERR_OUT("HDF5 library is not built with support of multi-dataset APIs")
#endif
        case hdf5_log:
#ifdef ENABLE_LOGVOL
            driver = new e3sm_io_driver_hdf5_log (cfg);
            break;
#else
            ERR_OUT ("Log VOL support was not enabled in this build")
#endif
        case adios:
#ifdef ENABLE_ADIOS2
            driver = new e3sm_io_driver_adios2 (cfg);
            break;
#else
            ERR_OUT ("ADIOS2 support was not enabled in this build")
#endif
        default:
            ERR_OUT ("I/O API is not set")
            break;
    }

err_out:;
    if (err) {
        if (driver) {
            delete driver;
            driver = NULL;
        }
    }
    return driver;
}

/*----< e3sm_io_xlen_nc_type() >---------------------------------------------*/
int e3sm_io_xlen_nc_type (nc_type xtype, int *size) {
	int err	  = 0;
	size_t sz = 0;
	switch (xtype) {
		case NC_NAT:
			sz = 0;
			break;
		case NC_BYTE:
			sz = sizeof (signed char);
			break;
		case NC_CHAR:
			sz = sizeof (char);
			break;
		case NC_SHORT:
			sz = sizeof (short);
			break;
		case NC_INT:
			sz = sizeof (int);
			break;
		case NC_FLOAT:
			sz = sizeof (float);
			break;
		case NC_DOUBLE:
			sz = sizeof (double);
			break;
		case NC_INT64:
			sz = sizeof (signed long long);
			break;
		case NC_UBYTE:
			sz = sizeof (unsigned char);
			break;
		case NC_USHORT:
			sz = sizeof (unsigned short);
			break;
		case NC_UINT:
			sz = sizeof (unsigned int);
			break;
		case NC_UINT64:
			sz = sizeof (unsigned long long);
			break;
#ifdef USE_NETCDF4
		case NC_STRING:
			sz = sizeof (char *);
			break;
#endif
		default:
			ERR_OUT ("Unknown type")
	}
	*size = (int)sz;

err_out:;
	return err;
}
