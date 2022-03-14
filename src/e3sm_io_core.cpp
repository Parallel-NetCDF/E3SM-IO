/*********************************************************************
 *
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 *********************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <e3sm_io.h>
#include <e3sm_io_err.h>
#include <e3sm_io_case.hpp>
#include <e3sm_io_driver.hpp>
#include <e3sm_io_profile.hpp>

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

e3sm_io_driver *e3sm_io_get_driver (const char *filename, e3sm_io_config *cfg) {
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
        /* remove the file system type prefix name if there is any.  For
         * example, when filename = "lustre:/home/foo/testfile.nc", remove
         * "lustre:" to make path = "/home/foo/testfile.nc" in open() below
         */
        path = strchr (filename, ':');
        if (path == NULL)
            path = filename; /* no prefix */
        else
            path++;

        /* must include config.h on 32-bit machines, as AC_SYS_LARGEFILE is
         * called at the configure time and it defines _FILE_OFFSET_BITS to 64
         * if large file feature is supported.
         */
        if ((fd = open(path, O_RDONLY, 00400)) == -1) { /* open for read */
            ERR_OUT ("Cannot open file.")
        }

        /* get first 8 bytes of file */
        rlen = read (fd, signature, 8);
        if (rlen != 8) { ERR_OUT ("Cannot read file.") }

        // PnetCDF ?
        if (memcmp (signature, cdf_signature, 3) == 0) {
            cfg->api = pnetcdf;
            goto done_check;
        }

        // NetCDF4 ?
        // NC4 is HDF5, so NC4 must be checked first
#ifdef ENABLE_NETCDF4
        if(e3sm_io_driver_nc4::compatible(std::string(path))) {
            cfg->api = netcdf4;
            goto done_check;
        }
#endif

        // HDF5 ?
#ifdef ENABLE_HDF5
        offset = 512;
        while (rlen == 8 && memcmp (signature, hdf5_signature, 8)) {
            lseek (fd, offset, SEEK_SET);
            offset <<= 1;
            rlen = read (fd, signature, 8);
        }
        if (rlen == 8) {  // It is HDF5, check if it is Log VOL
            hid_t fid = -1, gid = -1;

            fid = H5Fopen (path, H5F_ACC_RDONLY, H5P_DEFAULT);
            if (fid < 0) { ERR_OUT ("HDF5 header detected, but not a HDF5 file"); }

            gid = H5Gopen (fid, "_LOG", H5P_DEFAULT);
            if (gid >= 0) {
                cfg->api = hdf5_log;
            } else {
                cfg->api = hdf5;
            }
            if (gid >= 0) { H5Gclose (gid); }
            if (fid >= 0) { H5Fclose (fid); }

            goto done_check;
        }
#endif

// ADIOS2?
#ifdef ENABLE_ADIOS
        adios2_error aerr;
        adios2_adios *adp = NULL;
        adios2_io *iop    = NULL;
        adios2_engine *ep = NULL;
        adios2_bool result;

        adp = adios2_init (cfg->io_comm, "");
        CHECK_APTR (adp)

        iop = adios2_declare_io (adp, "e3sm_check");
        CHECK_APTR (iop)

        aerr = adios2_set_engine (fp->iop, "BP3");
        CHECK_AERR

        ep = adios2_open (iop, path, adios2_mode_read);
        if (ep) { cfg->api = adios; }

        adios2_close (ep);
        adios2_remove_io (&result, adp, "e3sm_check");
        adios2_finalize (adp);
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
#ifdef ENABLE_PNC
        case pnetcdf:
            driver = new e3sm_io_driver_pnc (cfg);
            break;
#endif
#ifdef ENABLE_NETCDF4
        case netcdf4:
            driver = new e3sm_io_driver_nc4 (cfg);
            break;
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
            ERR_OUT("HDF5 does not support multi-dataset APIs")
            break;
        case hdf5_log:
#ifdef ENABLE_HDF5
#ifdef ENABLE_LOGVOL
            driver = new e3sm_io_driver_hdf5_log (cfg);
#else
            ERR_OUT ("Log VOL support was not enabled in this build")
#endif
#else
            ERR_OUT ("HDF5 support was not enabled in this build")
#endif
            break;
        case adios:
#ifdef ENABLE_ADIOS2
            driver = new e3sm_io_driver_adios2 (cfg);
#else
            ERR_OUT ("ADIOS2 support was not enabled in this build")
#endif
            break;
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

extern "C" int e3sm_io_core (e3sm_io_config *cfg, e3sm_io_decom *decom) {
    int err=0;
    e3sm_io_case *tcase    = NULL;
    e3sm_io_driver *driver = NULL;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_TOTAL)
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_CORE)

    /* Select test case */
    E3SM_IO_TIMER_START (E3SM_IO_TIMER_INIT_CASE)
    tcase = new e3sm_io_case();
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_INIT_CASE)

    /* perform read */
    if (cfg->rd) {
        e3sm_io_api api_tmp = cfg->api;

        E3SM_IO_TIMER_START (E3SM_IO_TIMER_INIT_DRIVER)
        driver = e3sm_io_get_driver (cfg->in_path, cfg);
        CHECK_PTR (driver)
        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_INIT_DRIVER)

        E3SM_IO_TIMER_START (E3SM_IO_TIMER_RD)
        err = tcase->rd_test (*cfg, *decom, *driver);
        CHECK_ERR
        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_RD)

        cfg->api = api_tmp;  // Restore API to user setting for write test

        if (driver) {
            delete driver;
            driver = NULL;
        }
    }

    /* perform write */
    if (cfg->wr) {
        E3SM_IO_TIMER_START (E3SM_IO_TIMER_INIT_DRIVER)
        driver = e3sm_io_get_driver (NULL, cfg);
        CHECK_PTR (driver)
        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_INIT_DRIVER)

        E3SM_IO_TIMER_START (E3SM_IO_TIMER_WR)
        err = tcase->wr_test (*cfg, *decom, *driver);
        CHECK_ERR
        E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_WR)

        if (driver) {
            delete driver;
            driver = NULL;
        }
    }

err_out:
    if (driver) { delete driver; }
    if (tcase) { delete tcase; }

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_CORE)
    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_TOTAL)

    return err;
}
