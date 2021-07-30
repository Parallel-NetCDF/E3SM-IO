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
//
#include <cstdio>
#include <cstdlib>
#include <cstring>
//
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
//
#include <e3sm_io.h>
#include <e3sm_io_err.h>

#include <e3sm_io_case.hpp>
#include <e3sm_io_driver.hpp>
#include <e3sm_io_driver_pnc.hpp>
#include <e3sm_io_profile.hpp>
#ifdef ENABLE_HDF5
#include <hdf5.h>

#include <e3sm_io_driver_hdf5.hpp>
#include <e3sm_io_driver_h5blob.hpp>
#endif
#ifdef ENABLE_ADIOS2
#include <adios2_c.h>

#include <e3sm_io_driver_adios2.hpp>
#endif

e3sm_io_driver *e3sm_io_get_driver (const char *filename, e3sm_io_config *cfg) {
    int err                    = 0;
    const char *cdf_signature  = "CDF";
    const char *hdf5_signature = "\211HDF\r\n\032\n";
    const char *path;
    char signature[8];
    bool not_done = true;
    int fd;
    ssize_t rlen;
    e3sm_io_driver *driver = NULL;

    if (filename) {  // Reset API by file content
        if (cfg->rank == 0) {
            /* remove the file system type prefix name if there is any.
             * For example, when filename = "lustre:/home/foo/testfile.nc", remove
             * "lustre:" to make path = "/home/foo/testfile.nc" in open() below
             */
            path = strchr (filename, ':');
            if (path == NULL)
                path = filename; /* no prefix */
            else
                path++;

            /* must include config.h on 32-bit machines, as AC_SYS_LARGEFILE is called
             * at the configure time and it defines _FILE_OFFSET_BITS to 64 if large
             * file feature is supported.
             */
            if ((fd = open (path, O_RDONLY, 00400)) == -1) { /* open for read */
                ERR_OUT ("Cannot open file.")
            }

            /* get first 8 bytes of file */
            rlen = read (fd, signature, 8);
            if (rlen != 8) { ERR_OUT ("Cannot read file.") }

            // PnetCDF ?
            if (memcmp (signature, cdf_signature, 3) == 0) {
                cfg->api = pnetcdf;
                not_done = false;
            }

            // HDF5 ?
#ifdef ENABLE_HDF5
            if (not_done) {
                off_t offset = 512;

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
                    } else if (!((cfg->api == hdf5_md) || (cfg->api == hdf5))) {
                        ERR_OUT ("Selected API not compatible with input file format");
                    }

                    not_done = false;

                    if (gid >= 0) { H5Gclose (gid); }
                    if (fid >= 0) { H5Fclose (fid); }
                }
            }
#endif

// ADIOS2?
#ifdef ENABLE_ADIOS
            if (not_done) {
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
            }
#endif
            if (fd >= 0) { close (fd); }
        }

        err = MPI_Bcast (&(cfg->api), 1, MPI_INT, 0, cfg->io_comm);
        CHECK_MPIERR
    }

    // Cannot infer by file content, follow hint
    switch (cfg->api) {
        case pnetcdf:
            driver = new e3sm_io_driver_pnc (cfg);
            break;
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
        case hdf5_log:
#ifdef ENABLE_HDF5
            driver = new e3sm_io_driver_hdf5 (cfg);
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
    int err, nerrs = 0;
    e3sm_io_case *tcase    = NULL;
    e3sm_io_driver *driver = NULL;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_TOTAL)

    if ((cfg->verbose >= 0) && (cfg->rank == 0)) {
        printf ("Total number of MPI processes      = %d\n", cfg->np);
        printf ("Number of IO processes             = %d\n", cfg->num_iotasks);
        printf ("Input decomposition file           = %s\n", cfg->cfg_path);
        printf ("Number of decompositions           = %d\n", decom->num_decomp);
        if (decom->num_decomp == 3)
            printf("==== Benchmarking F case =============================\n");
        else if (decom->num_decomp == 6)
            printf("==== Benchmarking G case =============================\n");
        else if (decom->num_decomp == 5)
            printf("==== Benchmarking I case =============================\n");
    }

    /* Select test case */
    if (decom->num_decomp == 3) {
        /* F case has 3 decompositions */
        cfg->nvars = 414;
        switch (cfg->strategy) {
            case canonical:
            case log:
                tcase = new e3sm_io_case_F ();
                break;
            case blob:
                if (cfg->api == pnetcdf || cfg->api == hdf5)
                    tcase = new e3sm_io_case_F();
#ifdef ENABLE_ADIOS2
                else
                    tcase = new e3sm_io_case_F_scorpio();
#endif
                break;
            default:
                ERR_OUT ("I/O strategy is not set")
                break;
        }
    } else if (decom->num_decomp == 6) {
        /* G case has 6 decompositions */
        cfg->nrecs = 1;  /* only one record for G case */
        cfg->nvars = 52;
        switch (cfg->strategy) {
            case canonical:
            case log:
                tcase = new e3sm_io_case_G ();
                break;
            case blob:
                if (cfg->api == pnetcdf || cfg->api == hdf5)
                    tcase = new e3sm_io_case_G();
#ifdef ENABLE_ADIOS2
                else
                    tcase = new e3sm_io_case_G_scorpio ();
#endif
                break;
            default:
                ERR_OUT ("I/O strategy is not set")
                break;
        }
    } else if (decom->num_decomp == 5) {
        /* I case has 5 decompositions */
        cfg->nvars = 560;
        switch (cfg->strategy) {
            case canonical:
            case log:
                tcase = new e3sm_io_case_I ();
                break;
            case blob:
                if (cfg->api == pnetcdf || cfg->api == hdf5)
                    tcase = new e3sm_io_case_I();
#ifdef ENABLE_ADIOS2
                // else
                    // tcase = new e3sm_io_case_I_scorpio ();
#endif
                break;
            default:
                ERR_OUT ("I/O strategy is not set")
                break;
        }
    } else {
        ERR_OUT ("Unknown decom file")
    }

    /* perform read */
    if (cfg->rd) {
        if ((cfg->verbose >= 0) && (cfg->rank == 0)) {
            printf ("Input f ile/directory              = %s\n", cfg->in_path);
            printf ("Using noncontiguous read buffer    = %s\n", cfg->non_contig_buf ? "yes" : "no");
            if (cfg->two_buf)
                printf("Variable read order: 2D variables then 3D variables\n");
            else
                printf("Variable read order: same as variables are defined\n");
        }

        e3sm_io_api api_tmp = cfg->api;

        driver = e3sm_io_get_driver (cfg->in_path, cfg);
        CHECK_PTR (driver)

        err = tcase->rd_test (*cfg, *decom, *driver);
        CHECK_ERR

        cfg->api = api_tmp;  // Restore API to user setting for write test

        if (driver) {
            delete driver;
            driver = NULL;
        }
    }

    /* perform write */
    if (cfg->wr) {
        if ((cfg->verbose >= 0) && (cfg->rank == 0)) {
            printf ("Output file/directory              = %s\n", cfg->out_path);
            printf ("Using noncontiguous write buffer   = %s\n", cfg->non_contig_buf ? "yes" : "no");
            if (cfg->two_buf)
                printf("Variable write order: 2D variables then 3D variables\n");
            else
                printf("Variable write order: same as variables are defined\n");

            if (cfg->strategy == canonical) {
                if (cfg->api == pnetcdf)
                    printf("==== PnetCDF canonical I/O using varn API ============\n");
                else if (cfg->api == hdf5_md)
                    printf("==== HDF5 canonical I/O using multi-dataset API ======\n");
                else
                    ERR_OUT ("I/O strategy and API used is not supported yet")
            }
            else if (cfg->strategy == log) {
                if (cfg->api == hdf5_log)
                    printf("==== HDF5 using log-based VOL ========================\n");
                else
                    ERR_OUT ("I/O strategy and API used is not supported yet")
            }
            else if (cfg->strategy == blob) {
                if (cfg->api == pnetcdf)
                    printf("==== PnetCDF blob I/O ================================\n");
                else if (cfg->api == hdf5)
                    printf("==== HDF5 blob I/O ===================================\n");
                else if (cfg->api == adios)
                    printf("==== ADIOS blob I/O ==================================\n");
                else
                    ERR_OUT ("I/O strategy and API used is not supported yet")
            }
        }

        driver = e3sm_io_get_driver (NULL, cfg);
        CHECK_PTR (driver)

        err = tcase->wr_test (*cfg, *decom, *driver);
        CHECK_ERR

        if (driver) {
            delete driver;
            driver = NULL;
        }
    }

err_out:
    if (driver) { delete driver; }
    if (tcase) { delete tcase; }

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_TOTAL)

    if (cfg->verbose) e3sm_io_print_profile (cfg);

    return (err != 0) ? err : nerrs;
}
