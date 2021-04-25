/*********************************************************************
 *
 * Copyright (C) 2018, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 *********************************************************************/

#include <e3sm_io.h>
#include <e3sm_io_err.h>
#include <e3sm_io_case.hpp>
#include <e3sm_io_driver.hpp>
#include <e3sm_io_driver_pnc.hpp>

extern "C" int e3sm_io_core (e3sm_io_config *cfg, e3sm_io_decom *decom) {
    int nerrs=0;
    e3sm_io_case *tcase    = NULL;
    e3sm_io_driver *driver = NULL;

    switch (cfg->api) {
        case pnc:
            driver = new e3sm_io_driver_pnc ();
            break;
        default:
            RET_ERR ("Unknown driver")
            break;
    }

    /* F case has 3 decompositions, G case has 6 */
    if (decom->num_decomp == 3) {
        cfg->nvars = 414;
        tcase      = new e3sm_io_case_F ();
    } else if (decom->num_decomp == 6) {
        cfg->nvars = 52;
        tcase      = new e3sm_io_case_G ();
    } else {
        RET_ERR ("Unknown decom file")
    }

    if (cfg->wr) { tcase->wr_test (*cfg, *decom, *driver); }
    if (cfg->rd) { tcase->rd_test (*cfg, *decom, *driver); }

err_out:
    if (driver) { delete driver; }
    if (tcase) { delete tcase; }

    return nerrs;
}