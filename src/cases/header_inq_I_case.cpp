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

#include <string.h>
#include <assert.h>
#include <mpi.h>

#include <e3sm_io.h>
#include <e3sm_io_err.h>
#include <e3sm_io_case.hpp>
#include <e3sm_io_driver.hpp>

/*----< get_gattrs() >-------------------------------------------------------*/
/* add global attributes */
static
int get_gattrs(e3sm_io_config &cfg,
               e3sm_io_decom  &decom,
               e3sm_io_driver &driver,
               case_meta      *cmeta,
               int             ncid)
{
    char txtBuf[1024];
    std::string prefix("");
    int err=0, nprocs, intBuf;

    /* save number of processes as global attributes */
    if (cfg.strategy == blob) {
        if (cfg.api == adios) {
            GET_GATTR_INT("/__pio__/fillmode", &intBuf)
            prefix = "pio_global/";
        }
        else {
            MPI_Comm_size(cfg.io_comm, &nprocs);
            GET_GATTR_INT("global_nprocs", &nprocs)
            GET_GATTR_INT("num_decompositions", &decom.num_decomp)
            GET_GATTR_INT("num_subfiles", &cfg.num_subfiles)
        }
    }

    txtBuf[0] = '\0';
    /* 28 global attributes: */
    GET_GATTR_TXT("title", txtBuf)
    GET_GATTR_TXT("source", txtBuf)
    GET_GATTR_TXT("source_id", txtBuf)
    GET_GATTR_TXT("product", txtBuf)
    GET_GATTR_TXT("realm", txtBuf)
    GET_GATTR_TXT("case", txtBuf)
    GET_GATTR_TXT("username", txtBuf)
    GET_GATTR_TXT("hostname", txtBuf)
    GET_GATTR_TXT("git_version", txtBuf)
    GET_GATTR_TXT("history", txtBuf)
    GET_GATTR_TXT("institution_id", txtBuf)
    GET_GATTR_TXT("institution", txtBuf)
    GET_GATTR_TXT("contact", txtBuf)
    GET_GATTR_TXT("Conventions", txtBuf)
    GET_GATTR_TXT("comment", txtBuf)
    GET_GATTR_TXT("Surface_dataset", txtBuf)
    GET_GATTR_TXT("Initial_conditions_dataset", txtBuf)
    GET_GATTR_TXT("PFT_physiological_constants_dataset", txtBuf)

    GET_GATTR_INT("ltype_vegetated_or_bare_soil", &intBuf)
    GET_GATTR_INT("ltype_crop", &intBuf)
    GET_GATTR_INT("ltype_landice", &intBuf)
    GET_GATTR_INT("ltype_landice_multiple_elevation_classes", &intBuf)
    GET_GATTR_INT("ltype_deep_lake", &intBuf)
    GET_GATTR_INT("ltype_wetland", &intBuf)
    GET_GATTR_INT("ltype_urban_tbd", &intBuf)
    GET_GATTR_INT("ltype_urban_hd", &intBuf)
    GET_GATTR_INT("ltype_urban_md", &intBuf)

    if (cfg.hist == h1) { /* h1 only */
        GET_GATTR_TXT("Time_constant_3Dvars_filenamae", txtBuf)
        GET_GATTR_TXT("Time_constant_3Dvars", txtBuf)
    }

    assert(txtBuf[0] != '\0');
err_out:
    return err;
}

/*----< inq_I_case() >-------------------------------------------------------*/
int e3sm_io_case::inq_I_case(e3sm_io_config   &cfg,
                             e3sm_io_decom    &decom,
                             e3sm_io_driver   &driver,
                             int               ncid)    /* file ID */
{
    /* Total 52 variables */
    char txtBuf[1024];
    int i, err, nprocs, dimids[4], nvars_decomp, intBuf;
    int dim_time, dim_lon, dim_lat, dim_gridcell, dim_topounit, dim_landunit;
    int dim_column, dim_pft, dim_levgrnd, dim_levurb, dim_levlak, dim_numrad;
    int dim_month, dim_levsno, dim_ltype, dim_nvegwcs, dim_natpft, dim_nblobs;
    int dim_string_length, dim_levdcmp, dim_levtrc, dim_hist_interval;
    int dim_nelems[MAX_NUM_DECOMP], dim_max_nreqs[MAX_NUM_DECOMP];
    int g_dimids[MAX_NUM_DECOMP][MAX_NDIMS];
    // MPI_Offset lat, lon, levgrnd, levdcmp, levlak, ltype, natpft;
    // MPI_Offset string_length, hist_interval;
    float fillv = 1.e+36f, missv = 1.e+36f;
    std::map<int, std::string> dnames;
    var_meta *varp;
    case_meta *cmeta;

    if (cfg.run_case == F) {
        if (cfg.hist == h0) cmeta = &cfg.F_case_h0;
        else                cmeta = &cfg.F_case_h1;
    } else if (cfg.run_case == I) {
        if (cfg.hist == h0) cmeta = &cfg.I_case_h0;
        else                cmeta = &cfg.I_case_h1;
    } else if (cfg.run_case == G)
        cmeta = &cfg.G_case;

    /* add global attributes */
    err = get_gattrs(cfg, decom, driver, cmeta, ncid);
    CHECK_ERR

    /* define dimensions */

    /* piodecomp1344tasks_D1_lat_lon.dat
     *
     *  num_decomp = 5 ;
     *  decomp_nprocs = 1344 ;
     *
     *  dimensions:
     *    lon     = 720 ;
     *    lat     = 360 ;
     *    levgrnd =  15 ;
     *    levdcmp =  15 ;
     *    levlak  =  10 ;
     *    ltype   =   9 ;
     *    natpft  =  17 ;
     *                                         total: 560 variables
     *     D1:    360 720                lat x lon    465 variables
     *     D2: 15 360 720      levgrnd x lat x lon     19 variables
     *                         levdcmp x lat x lon     56 variables
     *     D3: 10 360 720      levlak  x lat x lon      4 variables
     *     D4:  9 360 720      ltype   x lat x lon      1 variables
     *     D5: 17 360 720      natpft  x lat x lon      1 variables
     *                                 not decomposed: 14 variables
     */

    /*
    lat           = decom.dims[0][0];
    lon           = decom.dims[0][1];
    levgrnd       = decom.dims[1][0];
    levdcmp       = decom.dims[1][0];
    levlak        = decom.dims[2][0];
    ltype         = decom.dims[3][0];
    natpft        = decom.dims[4][0];
    string_length =  16;
    hist_interval =   2;
    */

    /* define dimensions */
    INQ_DIM("time",           NC_UNLIMITED, &dim_time)
    INQ_DIM("lon",                     lon, &dim_lon)
    INQ_DIM("lat",                     lat, &dim_lat)
    INQ_DIM("gridcell",              62482, &dim_gridcell)
    INQ_DIM("topounit",              62482, &dim_topounit)
    INQ_DIM("landunit",             271014, &dim_landunit)
    INQ_DIM("column",              1020642, &dim_column)
    INQ_DIM("pft",                 2020354, &dim_pft)
    INQ_DIM("levgrnd",             levgrnd, &dim_levgrnd)
    INQ_DIM("levurb",                    5, &dim_levurb)
    INQ_DIM("levlak",               levlak, &dim_levlak)
    INQ_DIM("numrad",                    2, &dim_numrad)
    INQ_DIM("month",                    12, &dim_month)
    INQ_DIM("levsno",                    5, &dim_levsno)
    INQ_DIM("ltype",                 ltype, &dim_ltype)
    INQ_DIM("nvegwcs",                   4, &dim_nvegwcs)
    INQ_DIM("natpft",               natpft, &dim_natpft)
    INQ_DIM("string_length", string_length, &dim_string_length)
    INQ_DIM("levdcmp",             levdcmp, &dim_levdcmp)
    INQ_DIM("levtrc",                   10, &dim_levtrc)
    INQ_DIM("hist_interval", hist_interval, &dim_hist_interval)

    if (cfg.strategy == blob && cfg.api != adios) {
        char name[64];

        /* g_dimids[] are the original dimension IDs */
        g_dimids[0][0] = dim_lat;
        g_dimids[0][1] = dim_lon;
        for (i=1; i<decom.num_decomp; i++) {
            g_dimids[i][1] = dim_lat;
            g_dimids[i][2] = dim_lon;
        }
        g_dimids[1][0] = dim_levgrnd   /* same size as levdcmp */;
        g_dimids[2][0] = dim_levlak;
        g_dimids[3][0] = dim_ltype;
        g_dimids[4][0] = dim_natpft;

        /* additional dimensions to be used by decomposition variables */
        MPI_Comm_size(cfg.sub_comm, &nprocs);
        INQ_DIM("nblobs", (MPI_Offset)nprocs, &dim_nblobs)

        for (i=0; i<decom.num_decomp; i++) {
            sprintf(name, "D%d.nelems", i+1);
            INQ_DIM(name, decom.nelems[i], &dim_nelems[i])
        }

        for (i=0; i<decom.num_decomp; i++) {
            sprintf(name, "D%d.max_nreqs", i+1);
            INQ_DIM(name, (MPI_Offset)decom.max_nreqs[i], &dim_max_nreqs[i])
        }

        for (i=0; i<decom.num_decomp; i++) {
            fix_dimids[i]    = dim_nelems[i];
            rec_dimids[i][0] = dim_time;
            rec_dimids[i][1] = dim_nelems[i];
        }
    }

    /* define variables related to decompositions */
    nvars_decomp = 0;
    if (cfg.strategy == blob) {
        if (cfg.api == adios)
            nvars_decomp = decom.num_decomp + 1;
        else
            nvars_decomp = NVARS_DECOMP * decom.num_decomp;

        err = inq_var_decomp(cfg, decom, driver, cmeta, ncid, dim_time,
                             dim_nblobs, dim_max_nreqs, g_dimids);
        CHECK_ERR
    }

    /* For h0 file, There are 560 climate variables:
     *    0 scalar variables     + 560 array variables
     *   18 fixed-size variables + 542 record variables
     *   14 not partitioned      + 546 partitioned
     *
     * For h1 file, There are 552 climate variables:
     *    0 scalar variables     + 552 array variables
     *   10 fixed-size variables + 542 record variables
     *   14 not partitioned      + 538 partitioned
     *
     * For both h0 and h1 files, those 14 not partitioned variables are written
     * by root only:
     *     5 fixed-size, 9 record variables
     */
    varp = vars + nvars_decomp - 1;

    /* float levgrnd(levgrnd) */
    INQ_VAR("levgrnd", REC_XTYPE, 1, &dim_levgrnd, REC_ITYPE, -1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)

    /* float levlak(levlak) */
    INQ_VAR("levlak", REC_XTYPE, 1, &dim_levlak, REC_ITYPE, -1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)

    /* float levdcmp(levdcmp) */
    INQ_VAR("levdcmp", REC_XTYPE, 1, &dim_levdcmp, REC_ITYPE, -1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)

    /* float time(time) */
    INQ_VAR("time", REC_XTYPE, 1, &dim_time, REC_ITYPE, -1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("calendar", txtBuf)
    GET_ATTR_TXT("bounds", txtBuf)

    /* int mcdate(time) */
    INQ_VAR("mcdate", NC_INT, 1, &dim_time, MPI_INT, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* int mcsec(time) */
    INQ_VAR("mcsec", NC_INT, 1, &dim_time, MPI_INT, -1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)

    /* int mdcur(time) */
    INQ_VAR("mdcur", NC_INT, 1, &dim_time, MPI_INT, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* int mscur(time) */
    INQ_VAR("mscur", NC_INT, 1, &dim_time, MPI_INT, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* int nstep(time) */
    INQ_VAR("nstep", NC_INT, 1, &dim_time, MPI_INT, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* double time_bounds(time, hist_interval) */
    dimids[0] = dim_time;
    dimids[1] = dim_hist_interval;
    INQ_VAR("time_bounds", NC_DOUBLE, 2, dimids, MPI_DOUBLE, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* char date_written(time, string_length) */
    dimids[0] = dim_time;
    dimids[1] = dim_string_length;
    INQ_VAR("date_written", NC_CHAR, 2, dimids, MPI_CHAR, -1)

    /* char time_written(time, string_length) */
    dimids[0] = dim_time;
    dimids[1] = dim_string_length;
    INQ_VAR("time_written", NC_CHAR, 2, dimids, MPI_CHAR, -1)

    /* float lon(lon) */
    INQ_VAR("lon", REC_XTYPE, 1, &dim_lon, REC_ITYPE, -1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float lat(lat) */
    INQ_VAR("lat", REC_XTYPE, 1, &dim_lat, REC_ITYPE, -1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float area(lat, lon) */
    dimids[0] = dim_lat;
    dimids[1] = dim_lon;
    INQ_VAR("area", REC_XTYPE, 2, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float topo(lat, lon) */
    dimids[0] = dim_lat;
    dimids[1] = dim_lon;
    INQ_VAR("topo", REC_XTYPE, 2, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float landfrac(lat, lon) */
    dimids[0] = dim_lat;
    dimids[1] = dim_lon;
    INQ_VAR("landfrac", REC_XTYPE, 2, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* int landmask(lat, lon) */
    dimids[0] = dim_lat;
    dimids[1] = dim_lon;
    INQ_VAR("landmask", NC_INT, 2, dimids, MPI_INT, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_FILL(intBuf)
    GET_ATTR_INT1("missing_value", &intBuf)

    /* int pftmask(lat, lon) */
    dimids[0] = dim_lat;
    dimids[1] = dim_lon;
    INQ_VAR("pftmask", NC_INT, 2, dimids, MPI_INT, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_FILL(intBuf)
    GET_ATTR_INT1("missing_value", &intBuf)

    if (cfg.hist == h0) {  /* h0 only */
        /* float ZSOI(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        INQ_VAR("ZSOI", REC_XTYPE, 3, dimids, REC_ITYPE, 1)
        GET_ATTR_TXT("long_name", txtBuf)
        GET_ATTR_TXT("units", txtBuf)
        GET_ATTR_FILL(fillv)
        GET_ATTR_FLT1("missing_value", &missv)

        /* float DZSOI(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        INQ_VAR("DZSOI", REC_XTYPE, 3, dimids, REC_ITYPE, 1)
        GET_ATTR_TXT("long_name", txtBuf)
        GET_ATTR_TXT("units", txtBuf)
        GET_ATTR_FILL(fillv)
        GET_ATTR_FLT1("missing_value", &missv)

        /* float WATSAT(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        INQ_VAR("WATSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 1)
        GET_ATTR_TXT("long_name", txtBuf)
        GET_ATTR_TXT("units", txtBuf)
        GET_ATTR_FILL(fillv)
        GET_ATTR_FLT1("missing_value", &missv)

        /* float SUCSAT(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        INQ_VAR("SUCSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 1)
        GET_ATTR_TXT("long_name", txtBuf)
        GET_ATTR_TXT("units", txtBuf)
        GET_ATTR_FILL(fillv)
        GET_ATTR_FLT1("missing_value", &missv)

        /* float BSW(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        INQ_VAR("BSW", REC_XTYPE, 3, dimids, REC_ITYPE, 1)
        GET_ATTR_TXT("long_name", txtBuf)
        GET_ATTR_TXT("units", txtBuf)
        GET_ATTR_FILL(fillv)
        GET_ATTR_FLT1("missing_value", &missv)

        /* float HKSAT(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        INQ_VAR("HKSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 1)
        GET_ATTR_TXT("long_name", txtBuf)
        GET_ATTR_TXT("units", txtBuf)
        GET_ATTR_FILL(fillv)
        GET_ATTR_FLT1("missing_value", &missv)

        /* float ZLAKE(levlak, lat, lon) */
        dimids[0] = dim_levlak;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        INQ_VAR("ZLAKE", REC_XTYPE, 3, dimids, REC_ITYPE, 2)
        GET_ATTR_TXT("long_name", txtBuf)
        GET_ATTR_TXT("units", txtBuf)
        GET_ATTR_FILL(fillv)
        GET_ATTR_FLT1("missing_value", &missv)

        /* float DZLAKE(levlak, lat, lon) */
        dimids[0] = dim_levlak;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        INQ_VAR("DZLAKE", REC_XTYPE, 3, dimids, REC_ITYPE, 2)
        GET_ATTR_TXT("long_name", txtBuf)
        GET_ATTR_TXT("units", txtBuf)
        GET_ATTR_FILL(fillv)
        GET_ATTR_FLT1("missing_value", &missv)
    }

    /* float ACTUAL_IMMOB(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ACTUAL_IMMOB", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float ACTUAL_IMMOB_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ACTUAL_IMMOB_P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float ADSORBTION_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ADSORBTION_P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float AGNPP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("AGNPP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float AGWDNPP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("AGWDNPP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float ALT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ALT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float ALTMAX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ALTMAX", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float ALTMAX_LASTYEAR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ALTMAX_LASTYEAR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float AR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("AR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float AVAILC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("AVAILC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float AVAIL_RETRANSP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("AVAIL_RETRANSP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float BAF_CROP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("BAF_CROP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float BAF_PEATF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("BAF_PEATF", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float BCDEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("BCDEP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float BGNPP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("BGNPP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float BIOCHEM_PMIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("BIOCHEM_PMIN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float BIOCHEM_PMIN_TO_PLANT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("BIOCHEM_PMIN_TO_PLANT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float BTRAN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("BTRAN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float BUILDHEAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("BUILDHEAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CH4PROD(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CH4PROD", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CH4_SURF_AERE_SAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CH4_SURF_AERE_SAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CH4_SURF_AERE_UNSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CH4_SURF_AERE_UNSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CH4_SURF_DIFF_SAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CH4_SURF_DIFF_SAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CH4_SURF_DIFF_UNSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CH4_SURF_DIFF_UNSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CH4_SURF_EBUL_SAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CH4_SURF_EBUL_SAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CH4_SURF_EBUL_UNSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CH4_SURF_EBUL_UNSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float COL_PTRUNC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("COL_PTRUNC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CONC_CH4_SAT(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("CONC_CH4_SAT", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CONC_CH4_UNSAT(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("CONC_CH4_UNSAT", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CONC_O2_SAT(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("CONC_O2_SAT", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CONC_O2_UNSAT(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("CONC_O2_UNSAT", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CPOOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CWDC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CWDC_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDC_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CWDC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDC_LOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CWDC_TO_LITR2C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDC_TO_LITR2C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CWDC_TO_LITR3C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDC_TO_LITR3C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CWDC_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("CWDC_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CWDN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CWDN_TO_LITR2N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDN_TO_LITR2N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CWDN_TO_LITR3N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDN_TO_LITR3N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CWDN_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("CWDN_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CWDP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CWDP_TO_LITR2P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDP_TO_LITR2P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CWDP_TO_LITR3P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDP_TO_LITR3P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float CWDP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("CWDP_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DEADCROOTC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DEADCROOTC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DEADCROOTN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DEADCROOTN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DEADCROOTP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DEADCROOTP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DEADSTEMC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DEADSTEMC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DEADSTEMN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DEADSTEMN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DEADSTEMP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DEADSTEMP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DEFICIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DEFICIT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DENIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DENIT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DESORPTION_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DESORPTION_P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DISPVEGC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DISPVEGC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DISPVEGN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DISPVEGN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DISPVEGP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DISPVEGP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DSTDEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DSTDEP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DSTFLXT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DSTFLXT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DWB(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DWB", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DWT_CONV_CFLUX_DRIBBLED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DWT_CONV_CFLUX_DRIBBLED", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DWT_CONV_CFLUX_GRC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DWT_CONV_CFLUX_GRC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DWT_CONV_NFLUX_GRC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DWT_CONV_NFLUX_GRC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DWT_CONV_PFLUX_GRC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DWT_CONV_PFLUX_GRC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DWT_SLASH_CFLUX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DWT_SLASH_CFLUX", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DWT_SLASH_NFLUX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DWT_SLASH_NFLUX", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float DWT_SLASH_PFLUX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DWT_SLASH_PFLUX", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float EFLX_DYNBAL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("EFLX_DYNBAL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float EFLX_GRND_LAKE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("EFLX_GRND_LAKE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float EFLX_LH_TOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("EFLX_LH_TOT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float EFLX_LH_TOT_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("EFLX_LH_TOT_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float EFLX_LH_TOT_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("EFLX_LH_TOT_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float ELAI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ELAI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float ER(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ER", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float ERRH2O(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ERRH2O", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float ERRH2OSNO(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ERRH2OSNO", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float ERRSEB(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ERRSEB", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float ERRSOI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ERRSOI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float ERRSOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ERRSOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float ESAI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ESAI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FAREA_BURNED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FAREA_BURNED", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FCEV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FCEV", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FCH4(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FCH4", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FCH4TOCO2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FCH4TOCO2", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FCH4_DFSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FCH4_DFSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FCOV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FCOV", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FCTR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FCTR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FGEV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FGEV", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FGR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FGR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FGR12(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FGR12", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FGR_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FGR_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FGR_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FGR_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FH2OSFC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FH2OSFC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FINUNDATED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FINUNDATED", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FINUNDATED_LAG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FINUNDATED_LAG", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FIRA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FIRA", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FIRA_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FIRA_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FIRA_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FIRA_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FIRE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FIRE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FIRE_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FIRE_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FIRE_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FIRE_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FLDS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FLDS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FPG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FPG", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FPG_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FPG_P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FPI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FPI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FPI_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FPI_P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FPI_P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("FPI_P_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FPI_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("FPI_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FPSN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FPSN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FPSN_WC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FPSN_WC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FPSN_WJ(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FPSN_WJ", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FPSN_WP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FPSN_WP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FROOTC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FROOTC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FROOTC_ALLOC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FROOTC_ALLOC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FROOTC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FROOTC_LOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FROOTN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FROOTN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FROOTP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FROOTP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FROST_TABLE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FROST_TABLE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSA", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSA_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSA_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSA_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSA_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSDS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSDS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSDSND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSDSND", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSDSNDLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSDSNDLN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSDSNI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSDSNI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSDSVD(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSDSVD", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSDSVDLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSDSVDLN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSDSVI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSDSVI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSDSVILN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSDSVILN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSH(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSH", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSH_G(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSH_G", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSH_NODYNLNDUSE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSH_NODYNLNDUSE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSH_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSH_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSH_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSH_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSH_V(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSH_V", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSM(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSM", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSM_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSM_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSM_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSM_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSNO(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSNO", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSNO_EFF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSNO_EFF", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSRND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSRND", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSRNDLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSRNDLN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSRNI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSRNI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSRVD(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSRVD", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSRVDLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSRVDLN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float FSRVI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSRVI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float F_CO2_SOIL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("F_CO2_SOIL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float F_CO2_SOIL_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("F_CO2_SOIL_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float F_DENIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("F_DENIT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float F_DENIT_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("F_DENIT_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float F_N2O_DENIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("F_N2O_DENIT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float F_N2O_NIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("F_N2O_NIT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float F_NIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("F_NIT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float F_NIT_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("F_NIT_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float GC_HEAT1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("GC_HEAT1", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float GC_ICE1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("GC_ICE1", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float GC_LIQ1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("GC_LIQ1", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float GPP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("GPP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float GR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("GR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float GROSS_NMIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("GROSS_NMIN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float GROSS_PMIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("GROSS_PMIN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float H2OCAN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("H2OCAN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float H2OSFC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("H2OSFC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float H2OSNO(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("H2OSNO", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float H2OSNO_TOP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("H2OSNO_TOP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float H2OSOI(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("H2OSOI", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float HC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("HC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float HCSOI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("HCSOI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float HEAT_FROM_AC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("HEAT_FROM_AC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float HR_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("HR_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float HTOP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("HTOP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float INT_SNOW(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("INT_SNOW", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LABILEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LABILEP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LABILEP_TO_SECONDP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LABILEP_TO_SECONDP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LABILEP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LABILEP_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LAISHA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LAISHA", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LAISUN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LAISUN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LAKEICEFRAC(time, levlak, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levlak;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LAKEICEFRAC", REC_XTYPE, 4, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LAKEICETHICK(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LAKEICETHICK", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LAND_UPTAKE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LAND_UPTAKE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LAND_USE_FLUX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LAND_USE_FLUX", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LEAFC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LEAFC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LEAFC_ALLOC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LEAFC_ALLOC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LEAFC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LEAFC_LOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LEAFC_TO_LITTER(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LEAFC_TO_LITTER", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LEAFN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LEAFN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LEAFP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LEAFP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LEAF_MR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LEAF_MR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LFC2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LFC2", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITFALL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITFALL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITHR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITHR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR1C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR1C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR1C_TO_SOIL1C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR1C_TO_SOIL1C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR1C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR1C_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR1N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR1N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR1N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR1N_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR1N_TO_SOIL1N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR1N_TO_SOIL1N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR1N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR1N_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR1P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR1P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR1P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR1P_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR1P_TO_SOIL1P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR1P_TO_SOIL1P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR1P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR1P_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR1_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR1_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR2C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR2C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR2C_TO_SOIL2C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR2C_TO_SOIL2C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR2C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR2C_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR2N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR2N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR2N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR2N_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR2N_TO_SOIL2N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR2N_TO_SOIL2N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR2N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR2N_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR2P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR2P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR2P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR2P_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR2P_TO_SOIL2P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR2P_TO_SOIL2P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR2P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR2P_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR2_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR2_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR3C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR3C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR3C_TO_SOIL3C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR3C_TO_SOIL3C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR3C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR3C_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR3N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR3N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR3N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR3N_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR3N_TO_SOIL3N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR3N_TO_SOIL3N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR3N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR3N_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR3P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR3P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR3P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR3P_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR3P_TO_SOIL3P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR3P_TO_SOIL3P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR3P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR3P_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITR3_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR3_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITTERC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITTERC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITTERC_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITTERC_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LITTERC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITTERC_LOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LIVECROOTC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LIVECROOTC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LIVECROOTN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LIVECROOTN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LIVECROOTP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LIVECROOTP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LIVESTEMC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LIVESTEMC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LIVESTEMN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LIVESTEMN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float LIVESTEMP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LIVESTEMP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float MR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("MR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float M_LITR1C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("M_LITR1C_TO_LEACHING", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float M_LITR2C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("M_LITR2C_TO_LEACHING", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float M_LITR3C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("M_LITR3C_TO_LEACHING", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float M_SOIL1C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("M_SOIL1C_TO_LEACHING", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float M_SOIL2C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("M_SOIL2C_TO_LEACHING", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float M_SOIL3C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("M_SOIL3C_TO_LEACHING", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float M_SOIL4C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("M_SOIL4C_TO_LEACHING", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float NBP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NBP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float NDEPLOY(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NDEPLOY", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float NDEP_TO_SMINN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NDEP_TO_SMINN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float NEE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NEE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float NEM(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NEM", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float NEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NEP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float NET_NMIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NET_NMIN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float NET_PMIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NET_PMIN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float NFIRE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NFIRE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float NFIX_TO_SMINN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NFIX_TO_SMINN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float NPP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NPP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float OCCLP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("OCCLP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float OCCLP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("OCCLP_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float OCDEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("OCDEP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float O_SCALAR(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("O_SCALAR", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PARVEGLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PARVEGLN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PBOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PBOT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PCH4(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PCH4", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PCO2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PCO2", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PCT_LANDUNIT(time, ltype, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_ltype;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("PCT_LANDUNIT", REC_XTYPE, 4, dimids, REC_ITYPE, 3)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PCT_NAT_PFT(time, natpft, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_natpft;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("PCT_NAT_PFT", REC_XTYPE, 4, dimids, REC_ITYPE, 4)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PDEPLOY(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PDEPLOY", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PDEP_TO_SMINP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PDEP_TO_SMINP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PFT_FIRE_CLOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PFT_FIRE_CLOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PFT_FIRE_NLOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PFT_FIRE_NLOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PLANT_CALLOC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PLANT_CALLOC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PLANT_NDEMAND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PLANT_NDEMAND", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PLANT_NDEMAND_COL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PLANT_NDEMAND_COL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PLANT_PALLOC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PLANT_PALLOC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PLANT_PDEMAND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PLANT_PDEMAND", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PLANT_PDEMAND_COL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PLANT_PDEMAND_COL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float POTENTIAL_IMMOB(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("POTENTIAL_IMMOB", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float POTENTIAL_IMMOB_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("POTENTIAL_IMMOB_P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float POT_F_DENIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("POT_F_DENIT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float POT_F_NIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("POT_F_NIT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PRIMP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PRIMP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PRIMP_TO_LABILEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PRIMP_TO_LABILEP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PRIMP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("PRIMP_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PROD1P_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PROD1P_LOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PSNSHA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PSNSHA", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PSNSHADE_TO_CPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PSNSHADE_TO_CPOOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PSNSUN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PSNSUN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float PSNSUN_TO_CPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PSNSUN_TO_CPOOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float Q2M(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("Q2M", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QBOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QBOT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QCHARGE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QCHARGE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QDRAI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QDRAI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QDRAI_PERCH(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QDRAI_PERCH", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QDRAI_XS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QDRAI_XS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QDRIP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QDRIP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QFLOOD(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QFLOOD", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QFLX_ICE_DYNBAL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QFLX_ICE_DYNBAL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QFLX_LIQ_DYNBAL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QFLX_LIQ_DYNBAL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QH2OSFC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QH2OSFC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QINFL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QINFL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QINTR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QINTR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QIRRIG_GRND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QIRRIG_GRND", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QIRRIG_ORIG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QIRRIG_ORIG", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QIRRIG_REAL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QIRRIG_REAL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QIRRIG_SURF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QIRRIG_SURF", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QIRRIG_WM(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QIRRIG_WM", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QOVER(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QOVER", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QOVER_LAG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QOVER_LAG", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QRGWL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QRGWL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QRUNOFF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QRUNOFF", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QRUNOFF_NODYNLNDUSE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QRUNOFF_NODYNLNDUSE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QRUNOFF_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QRUNOFF_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QRUNOFF_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QRUNOFF_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QSNOMELT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QSNOMELT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QSNWCPICE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QSNWCPICE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QSNWCPICE_NODYNLNDUSE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QSNWCPICE_NODYNLNDUSE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QSOIL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QSOIL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QVEGE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QVEGE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float QVEGT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QVEGT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float RAIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RAIN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float RETRANSN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RETRANSN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float RETRANSN_TO_NPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RETRANSN_TO_NPOOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float RETRANSP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RETRANSP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float RETRANSP_TO_PPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RETRANSP_TO_PPOOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float RH2M(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RH2M", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float RH2M_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RH2M_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float RH2M_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RH2M_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float RR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SABG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SABG", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SABG_PEN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SABG_PEN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SABV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SABV", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SCALARAVG_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SCALARAVG_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SECONDP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SECONDP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SECONDP_TO_LABILEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SECONDP_TO_LABILEP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SECONDP_TO_OCCLP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SECONDP_TO_OCCLP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SECONDP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SECONDP_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SEEDC_GRC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SEEDC_GRC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINN_TO_NPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN_TO_NPOOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINN_TO_PLANT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN_TO_PLANT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINN_TO_SOIL1N_L1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN_TO_SOIL1N_L1", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINN_TO_SOIL2N_L2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN_TO_SOIL2N_L2", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINN_TO_SOIL2N_S1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN_TO_SOIL2N_S1", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINN_TO_SOIL3N_L3(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN_TO_SOIL3N_L3", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINN_TO_SOIL3N_S2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN_TO_SOIL3N_S2", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINN_TO_SOIL4N_S3(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN_TO_SOIL4N_S3", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINP_LEACHED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_LEACHED", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINP_TO_PLANT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_TO_PLANT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINP_TO_PPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_TO_PPOOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINP_TO_SOIL1P_L1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_TO_SOIL1P_L1", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINP_TO_SOIL2P_L2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_TO_SOIL2P_L2", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINP_TO_SOIL2P_S1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_TO_SOIL2P_S1", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINP_TO_SOIL3P_L3(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_TO_SOIL3P_L3", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINP_TO_SOIL3P_S2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_TO_SOIL3P_S2", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINP_TO_SOIL4P_S3(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_TO_SOIL4P_S3", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMINP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SMINP_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMIN_NH4(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMIN_NH4", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMIN_NH4_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SMIN_NH4_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMIN_NO3(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMIN_NO3", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMIN_NO3_LEACHED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMIN_NO3_LEACHED", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMIN_NO3_RUNOFF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMIN_NO3_RUNOFF", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SMIN_NO3_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SMIN_NO3_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SNOBCMCL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOBCMCL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SNOBCMSL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOBCMSL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SNODSTMCL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNODSTMCL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SNODSTMSL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNODSTMSL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SNOINTABS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOINTABS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SNOOCMCL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOOCMCL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SNOOCMSL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOOCMSL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SNOW(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOW", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SNOWDP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOWDP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SNOWICE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOWICE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SNOWLIQ(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOWLIQ", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SNOW_DEPTH(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOW_DEPTH", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SNOW_SINKS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOW_SINKS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SNOW_SOURCES(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOW_SOURCES", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL1C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL1C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL1C_TO_SOIL2C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL1C_TO_SOIL2C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL1C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL1C_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL1N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL1N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL1N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL1N_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL1N_TO_SOIL2N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL1N_TO_SOIL2N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL1N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL1N_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL1P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL1P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL1P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL1P_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL1P_TO_SOIL2P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL1P_TO_SOIL2P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL1P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL1P_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL1_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL1_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL2C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL2C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL2C_TO_SOIL3C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL2C_TO_SOIL3C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL2C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL2C_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL2N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL2N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL2N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL2N_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL2N_TO_SOIL3N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL2N_TO_SOIL3N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL2N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL2N_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL2P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL2P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL2P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL2P_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL2P_TO_SOIL3P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL2P_TO_SOIL3P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL2P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL2P_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL2_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL2_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL3C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL3C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL3C_TO_SOIL4C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL3C_TO_SOIL4C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL3C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL3C_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL3N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL3N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL3N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL3N_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL3N_TO_SOIL4N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL3N_TO_SOIL4N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL3N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL3N_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL3P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL3P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL3P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL3P_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL3P_TO_SOIL4P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL3P_TO_SOIL4P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL3P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL3P_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL3_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL3_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL4C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL4C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL4C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL4C_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL4N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL4N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL4N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL4N_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL4N_TO_SMINN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL4N_TO_SMINN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL4N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL4N_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL4P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL4P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL4P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL4P_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL4P_TO_SMINP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL4P_TO_SMINP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL4P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL4P_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOIL4_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL4_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOILC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOILC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOILC_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOILC_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOILC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOILC_LOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOILICE(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOILICE", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOILICE_ICE(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOILICE_ICE", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOILLIQ(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOILLIQ", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOILLIQ_ICE(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOILLIQ_ICE", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOILPSI(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOILPSI", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOILWATER_10CM(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOILWATER_10CM", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOLUTIONP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOLUTIONP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOLUTIONP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOLUTIONP_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOMHR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOMHR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SOM_C_LEACHED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOM_C_LEACHED", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float STORVEGC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("STORVEGC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float STORVEGN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("STORVEGN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float STORVEGP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("STORVEGP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SUPPLEMENT_TO_SMINN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SUPPLEMENT_TO_SMINN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SUPPLEMENT_TO_SMINP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SUPPLEMENT_TO_SMINP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SUPPLY(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SUPPLY", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SoilAlpha(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SoilAlpha", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float SoilAlpha_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SoilAlpha_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TAUX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TAUX", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TAUY(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TAUY", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TBOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TBOT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TBUILD(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TBUILD", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TCS_MONTH_BEGIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TCS_MONTH_BEGIN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TCS_MONTH_END(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TCS_MONTH_END", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TG", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TG_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TG_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TG_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TG_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TH2OSFC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TH2OSFC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float THBOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("THBOT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TKE1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TKE1", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TLAI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TLAI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TLAKE(time, levlak, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levlak;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("TLAKE", REC_XTYPE, 4, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTCOLC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTCOLC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTCOLCH4(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTCOLCH4", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTCOLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTCOLN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTCOLP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTCOLP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTECOSYSC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTECOSYSC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTECOSYSN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTECOSYSN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTECOSYSP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTECOSYSP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTLITC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTLITC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTLITC_1m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTLITC_1m", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTLITN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTLITN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTLITP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTLITP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTLITP_1m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTLITP_1m", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTPFTC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTPFTC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTPFTN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTPFTN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTPFTP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTPFTP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTSOMC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTSOMC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTSOMC_1m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTSOMC_1m", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTSOMN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTSOMN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTSOMP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTSOMP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTSOMP_1m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTSOMP_1m", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTVEGC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTVEGC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTVEGC_ABG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTVEGC_ABG", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTVEGN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTVEGN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TOTVEGP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTVEGP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TREFMNAV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TREFMNAV", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TREFMNAV_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TREFMNAV_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TREFMNAV_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TREFMNAV_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TREFMXAV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TREFMXAV", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TREFMXAV_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TREFMXAV_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TREFMXAV_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TREFMXAV_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TSA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TSA", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TSAI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TSAI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TSA_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TSA_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TSA_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TSA_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TSOI(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("TSOI", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TSOI_10CM(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TSOI_10CM", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TSOI_ICE(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("TSOI_ICE", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TV", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TWS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TWS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TWS_MONTH_BEGIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TWS_MONTH_BEGIN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float TWS_MONTH_END(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TWS_MONTH_END", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float T_SCALAR(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("T_SCALAR", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float U10(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("U10", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float URBAN_AC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("URBAN_AC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float URBAN_HEAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("URBAN_HEAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float VOLR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("VOLR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float VOLRMCH(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("VOLRMCH", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float WA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WA", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float WASTEHEAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WASTEHEAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float WF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WF", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float WIND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WIND", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float WOODC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WOODC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float WOODC_ALLOC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WOODC_ALLOC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float WOODC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WOODC_LOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float WOOD_HARVESTC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WOOD_HARVESTC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float WOOD_HARVESTN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WOOD_HARVESTN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float WTGQ(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WTGQ", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float W_SCALAR(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("W_SCALAR", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float XR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("XR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float XSMRPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("XSMRPOOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float ZBOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ZBOT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float ZWT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ZWT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float ZWT_CH4_UNSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ZWT_CH4_UNSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float ZWT_PERCH(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ZWT_PERCH", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float cn_scalar(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("cn_scalar", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float cp_scalar(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("cp_scalar", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float leaf_npimbalance(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("leaf_npimbalance", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float nlim_m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("nlim_m", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float o2_decomp_depth_unsat(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("o2_decomp_depth_unsat", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    /* float plim_m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("plim_m", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)

    assert(varp - vars + 1 == cfg.nvars + nvars_decomp);

err_out:
    return err;
}

