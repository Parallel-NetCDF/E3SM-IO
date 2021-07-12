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
#include <e3sm_io.h>
#include <e3sm_io_err.h>

#include <e3sm_io_case.hpp>
#include <e3sm_io_driver.hpp>
#include <e3sm_io_driver_pnc.hpp>
#include <e3sm_io_profile.hpp>
#ifdef ENABLE_HDF5
#include <e3sm_io_driver_hdf5.hpp>
#endif
#ifdef ENABLE_ADIOS2
#include <e3sm_io_driver_adios2.hpp>
#endif

extern "C" int e3sm_io_core (e3sm_io_config *cfg, e3sm_io_decom *decom) {
    int err, nerrs = 0;
    e3sm_io_case *tcase    = NULL;
    e3sm_io_driver *driver = NULL;

    E3SM_IO_TIMER_START (E3SM_IO_TIMER_TOTAL)

    switch (cfg->api) {
        case pnetcdf:
            driver = new e3sm_io_driver_pnc (cfg);
            break;
        case hdf5_ra:
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
            ERR_OUT ("Unknown driver")
            break;
    }

    /* F case has 3 decompositions, G case has 6 */
    if (decom->num_decomp == 3) {
        cfg->nvars = 414;
        switch (cfg->strategy) {
            case canonical:
            case log:
                tcase = new e3sm_io_case_F ();
                break;
            case blob:
                if (cfg->api == pnetcdf)
                    tcase = new e3sm_io_case_F ();
#ifdef ENABLE_ADIOS2
                else
                    tcase = new e3sm_io_case_F_scorpio ();
#endif
                break;
        }
    } else if (decom->num_decomp == 6) {
        cfg->nvars = 52;
        switch (cfg->strategy) {
            case canonical:
            case log:
                tcase = new e3sm_io_case_G ();
                break;
            case blob:
                if (cfg->api == pnetcdf)
                    tcase = new e3sm_io_case_G ();
#ifdef ENABLE_ADIOS2
/* TODO:
                else
                    tcase = new e3sm_io_case_G_scorpio ();
*/
#endif
                break;
        }
    } else {
        ERR_OUT ("Unknown decom file")
    }

    if (cfg->wr) {
        err = tcase->wr_test (*cfg, *decom, *driver);
        if (err != 0) goto err_out;
    }
    if (cfg->rd) {
        err = tcase->rd_test (*cfg, *decom, *driver);
        if (err != 0) goto err_out;
    }

err_out:
    if (driver) { delete driver; }
    if (tcase) { delete tcase; }

    E3SM_IO_TIMER_STOP (E3SM_IO_TIMER_TOTAL)

    if (cfg->verbose) e3sm_io_print_profile (cfg);

    return (err != 0) ? err : nerrs;
}
