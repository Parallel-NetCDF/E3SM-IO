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

#include <e3sm_io.h>
#include <e3sm_io_err.h>
#include <e3sm_io_case.hpp>
#include <e3sm_io_driver.hpp>
#include <e3sm_io_profile.hpp>

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
        char path[1028], *ext;
        ext = strrchr(cfg->in_path, '.');

        E3SM_IO_TIMER_START (E3SM_IO_TIMER_INIT_DRIVER)
        /* construct file name */
        if (cfg->run_case == F || cfg->run_case == I) {
            if (ext == NULL || (strcmp(ext, ".nc") && strcmp(ext, ".h5") && strcmp(ext, ".nc4")))
                if (cfg->hx == 0 || cfg->hx == -1){
                    sprintf(path, "%s_h0", cfg->in_path);
                }
                else{
                    sprintf(path, "%s_h1", cfg->in_path);
                }
            else { /* add "_h0" before file extension */
                strcpy(path, cfg->in_path);
                if (cfg->hx == 0 || cfg->hx == -1){
                    sprintf(path + (ext - cfg->in_path), "_h0%s", ext);
                }
                else{
                    sprintf(path + (ext - cfg->in_path), "_h1%s", ext);
                }
            }
        }
        else{
            strcpy(path, cfg->in_path);
        }
        if (cfg->strategy == blob && cfg->api != adios) {
            /* append subfile ID to subfile name */
            char sub_str[8];
            sprintf(sub_str, ".%04d", cfg->subfile_ID);
            strcat(path, sub_str);
        }

        driver = e3sm_io_get_driver (path, cfg);
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
