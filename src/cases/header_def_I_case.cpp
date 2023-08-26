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

/*----< add_gattrs() >-------------------------------------------------------*/
/* add global attributes */
static
int add_gattrs(e3sm_io_config &cfg,
               e3sm_io_decom  &decom,
               e3sm_io_driver &driver,
               case_meta      *cmeta,
               int             ncid)
{
    std::string prefix("");
    int err=0, nprocs;

    /* save number of processes as global attributes */
    if (cfg.strategy == blob) {
        if (cfg.api == adios) {
            PUT_GATTR_INT("/__pio__/fillmode", 256)
            prefix = "/__pio__/global/";
        }
        else {
            MPI_Comm_size(cfg.io_comm, &nprocs);
            PUT_GATTR_INT("global_nprocs", nprocs)
            PUT_GATTR_INT("num_decompositions", decom.num_decomp)
            PUT_GATTR_INT("num_subfiles", cfg.num_subfiles)
        }
    }

    /* 28 global attributes: */
    PUT_GATTR_TXT("title", "ELM History file information")
    PUT_GATTR_TXT("source", "E3SM Land Model")
    PUT_GATTR_TXT("source_id", "6eb829238a")
    PUT_GATTR_TXT("product", "model-output")
    PUT_GATTR_TXT("realm", "land")
    PUT_GATTR_TXT("case", "I1850GSWCNPRDCTCBC_hcru_hcru")
    PUT_GATTR_TXT("username", "E3SM")
    PUT_GATTR_TXT("hostname", "cori-knl")
    PUT_GATTR_TXT("git_version", "6eb829238a")
    PUT_GATTR_TXT("history", "created on 07/13/21 20:37:15")
    PUT_GATTR_TXT("institution_id", "E3SM-Project")
    PUT_GATTR_TXT("institution", "LLNL (Lawrence Livermore National Laboratory, Livermore, CA 94550, USA); ANL (Argonne National Laboratory, Argonne, IL 60439, USA); BNL (Brookhaven National Laboratory, Upton, NY 11973, USA); LANL (Los Alamos National Laboratory, Los Alamos, NM 87545, USA); LBNL (Lawrence Berkeley National Laboratory, Berkeley, CA 94720, USA); ORNL (Oak Ridge National Laboratory, Oak Ridge, TN 37831, USA); PNNL (Pacific Northwest National Laboratory, Richland, WA 99352, USA); SNL (Sandia National Laboratories, Albuquerque, NM 87185, USA). Mailing address: LLNL Climate Program, c/o David C. Bader, Principal Investigator, L-103, 7000 East Avenue, Livermore, CA 94550, USA")
    PUT_GATTR_TXT("contact", "e3sm-data-support@listserv.llnl.gov")
    PUT_GATTR_TXT("Conventions", "CF-1.7")
    PUT_GATTR_TXT("comment", "NOTE: None of the variables are weighted by land fraction!")
    PUT_GATTR_TXT("Surface_dataset", "surfdata_360x720cru_simyr1850_c180216.nc")
    PUT_GATTR_TXT("Initial_conditions_dataset", "arbitrary initialization")
    PUT_GATTR_TXT("PFT_physiological_constants_dataset", "clm_params_c180524.nc")

    PUT_GATTR_INT("ltype_vegetated_or_bare_soil", 1)
    PUT_GATTR_INT("ltype_crop", 2)
    PUT_GATTR_INT("ltype_landice", 3)
    PUT_GATTR_INT("ltype_landice_multiple_elevation_classes", 4)
    PUT_GATTR_INT("ltype_deep_lake", 5)
    PUT_GATTR_INT("ltype_wetland", 6)
    PUT_GATTR_INT("ltype_urban_tbd", 7)
    PUT_GATTR_INT("ltype_urban_hd", 8)
    PUT_GATTR_INT("ltype_urban_md", 9)

    if (cfg.hist == h1) { /* h1 only */
        PUT_GATTR_TXT("Time_constant_3Dvars_filenamae", "./I1850GSWCNPRDCTCBC_hcru_hcru.elm.h0.0001-01-01-00000.nc")
        PUT_GATTR_TXT("Time_constant_3Dvars", "ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE")
    }

err_out:
    return err;
}

/*----< def_I_case() >-------------------------------------------------------*/
int e3sm_io_case::def_I_case(e3sm_io_config   &cfg,
                             e3sm_io_decom    &decom,
                             e3sm_io_driver   &driver,
                             int               ncid)    /* file ID */
{
    /* Total 52 variables */
    int i, err, nprocs, dimids[4], nvars_decomp;
    int dim_time, dim_lon, dim_lat, dim_gridcell, dim_topounit, dim_landunit;
    int dim_column, dim_pft, dim_levgrnd, dim_levurb, dim_levlak, dim_numrad;
    int dim_month, dim_levsno, dim_ltype, dim_nvegwcs, dim_natpft, dim_nblobs;
    int dim_string_length, dim_levdcmp, dim_levtrc, dim_hist_interval;
    int dim_nelems[MAX_NUM_DECOMP], dim_max_nreqs[MAX_NUM_DECOMP];
    int g_dimids[MAX_NUM_DECOMP][MAX_NDIMS];
    MPI_Offset lat, lon, levgrnd, levdcmp, levlak, ltype, natpft;
    MPI_Offset string_length, hist_interval;
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
    err = add_gattrs(cfg, decom, driver, cmeta, ncid);
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

    lat           = decom.dims[0][0];
    lon           = decom.dims[0][1];
    levgrnd       = decom.dims[1][0];
    levdcmp       = decom.dims[1][0];
    levlak        = decom.dims[2][0];
    ltype         = decom.dims[3][0];
    natpft        = decom.dims[4][0];
    string_length =  16;
    hist_interval =   2;

    /* define dimensions */
    DEF_DIM("time",           NC_UNLIMITED, &dim_time)
    DEF_DIM("lon",                     lon, &dim_lon)
    DEF_DIM("lat",                     lat, &dim_lat)
    DEF_DIM("gridcell",              62482, &dim_gridcell)
    DEF_DIM("topounit",              62482, &dim_topounit)
    DEF_DIM("landunit",             271014, &dim_landunit)
    DEF_DIM("column",              1020642, &dim_column)
    DEF_DIM("pft",                 2020354, &dim_pft)
    DEF_DIM("levgrnd",             levgrnd, &dim_levgrnd)
    DEF_DIM("levurb",                    5, &dim_levurb)
    DEF_DIM("levlak",               levlak, &dim_levlak)
    DEF_DIM("numrad",                    2, &dim_numrad)
    DEF_DIM("month",                    12, &dim_month)
    DEF_DIM("levsno",                    5, &dim_levsno)
    DEF_DIM("ltype",                 ltype, &dim_ltype)
    DEF_DIM("nvegwcs",                   4, &dim_nvegwcs)
    DEF_DIM("natpft",               natpft, &dim_natpft)
    DEF_DIM("string_length", string_length, &dim_string_length)
    DEF_DIM("levdcmp",             levdcmp, &dim_levdcmp)
    DEF_DIM("levtrc",                   10, &dim_levtrc)
    DEF_DIM("hist_interval", hist_interval, &dim_hist_interval)

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
        DEF_DIM("nblobs", (MPI_Offset)nprocs, &dim_nblobs)

        for (i=0; i<decom.num_decomp; i++) {
            sprintf(name, "D%d.nelems", i+1);
            DEF_DIM(name, decom.nelems[i], &dim_nelems[i])
        }

        for (i=0; i<decom.num_decomp; i++) {
            sprintf(name, "D%d.max_nreqs", i+1);
            DEF_DIM(name, (MPI_Offset)decom.max_nreqs[i], &dim_max_nreqs[i])
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
            nvars_decomp = (2 * decom.num_decomp) + 3;
        else
            nvars_decomp = NVARS_DECOMP * decom.num_decomp;

        err = def_var_decomp(cfg, decom, driver, cmeta, ncid, dim_time,
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

    /* define fixed-size variables */

    /* int landmask(lat, lon) */
    dimids[0] = dim_lat;
    dimids[1] = dim_lon;
    DEF_VAR("landmask", NC_INT, 2, dimids, MPI_INT, 0)
    PUT_ATTR_TXT("long_name", "land/ocean mask (0.=ocean and 1.=land)")
    PUT_ATTR_FILL(-9999)
    PUT_ATTR_INT1("missing_value", -9999)

    /* int pftmask(lat, lon) */
    dimids[0] = dim_lat;
    dimids[1] = dim_lon;
    DEF_VAR("pftmask", NC_INT, 2, dimids, MPI_INT, 0)
    PUT_ATTR_TXT("long_name", "pft real/fake mask (0.=fake and 1.=real)")
    PUT_ATTR_FILL(-9999)
    PUT_ATTR_INT1("missing_value", -9999)

    /* float levgrnd(levgrnd) */
    DEF_VAR("levgrnd", REC_XTYPE, 1, &dim_levgrnd, REC_ITYPE, -1)
    PUT_ATTR_TXT("long_name", "coordinate soil levels")
    PUT_ATTR_TXT("units", "m")

    /* float levlak(levlak) */
    DEF_VAR("levlak", REC_XTYPE, 1, &dim_levlak, REC_ITYPE, -1)
    PUT_ATTR_TXT("long_name", "coordinate lake levels")
    PUT_ATTR_TXT("units", "m")

    /* float levdcmp(levdcmp) */
    DEF_VAR("levdcmp", REC_XTYPE, 1, &dim_levdcmp, REC_ITYPE, -1)
    PUT_ATTR_TXT("long_name", "coordinate soil levels")
    PUT_ATTR_TXT("units", "m")

    /* float lon(lon) */
    DEF_VAR("lon", REC_XTYPE, 1, &dim_lon, REC_ITYPE, -1)
    PUT_ATTR_TXT("long_name", "coordinate longitude")
    PUT_ATTR_TXT("units", "degrees_east")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float lat(lat) */
    DEF_VAR("lat", REC_XTYPE, 1, &dim_lat, REC_ITYPE, -1)
    PUT_ATTR_TXT("long_name", "coordinate latitude")
    PUT_ATTR_TXT("units", "degrees_north")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float area(lat, lon) */
    dimids[0] = dim_lat;
    dimids[1] = dim_lon;
    DEF_VAR("area", REC_XTYPE, 2, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "grid cell areas")
    PUT_ATTR_TXT("units", "km^2")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float topo(lat, lon) */
    dimids[0] = dim_lat;
    dimids[1] = dim_lon;
    DEF_VAR("topo", REC_XTYPE, 2, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "grid cell topography")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float landfrac(lat, lon) */
    dimids[0] = dim_lat;
    dimids[1] = dim_lon;
    DEF_VAR("landfrac", REC_XTYPE, 2, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "land fraction")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    if (cfg.hist == h0) {  /* h0 only */
        /* float ZSOI(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        DEF_VAR("ZSOI", REC_XTYPE, 3, dimids, REC_ITYPE, 1)
        PUT_ATTR_TXT("long_name", "soil depth")
        PUT_ATTR_TXT("units", "m")
        PUT_ATTR_FILL(fillv)
        PUT_ATTR_FLT1("missing_value", missv)

        /* float DZSOI(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        DEF_VAR("DZSOI", REC_XTYPE, 3, dimids, REC_ITYPE, 1)
        PUT_ATTR_TXT("long_name", "soil thickness")
        PUT_ATTR_TXT("units", "m")
        PUT_ATTR_FILL(fillv)
        PUT_ATTR_FLT1("missing_value", missv)

        /* float WATSAT(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        DEF_VAR("WATSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 1)
        PUT_ATTR_TXT("long_name", "saturated soil water content (porosity)")
        PUT_ATTR_TXT("units", "mm3/mm3")
        PUT_ATTR_FILL(fillv)
        PUT_ATTR_FLT1("missing_value", missv)

        /* float SUCSAT(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        DEF_VAR("SUCSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 1)
        PUT_ATTR_TXT("long_name", "saturated soil matric potential")
        PUT_ATTR_TXT("units", "mm")
        PUT_ATTR_FILL(fillv)
        PUT_ATTR_FLT1("missing_value", missv)

        /* float BSW(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        DEF_VAR("BSW", REC_XTYPE, 3, dimids, REC_ITYPE, 1)
        PUT_ATTR_TXT("long_name", "slope of soil water retention curve")
        PUT_ATTR_TXT("units", "1")
        PUT_ATTR_FILL(fillv)
        PUT_ATTR_FLT1("missing_value", missv)

        /* float HKSAT(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        DEF_VAR("HKSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 1)
        PUT_ATTR_TXT("long_name", "saturated hydraulic conductivity")
        PUT_ATTR_TXT("units", "1")
        PUT_ATTR_FILL(fillv)
        PUT_ATTR_FLT1("missing_value", missv)

        /* float ZLAKE(levlak, lat, lon) */
        dimids[0] = dim_levlak;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        DEF_VAR("ZLAKE", REC_XTYPE, 3, dimids, REC_ITYPE, 2)
        PUT_ATTR_TXT("long_name", "lake layer node depth")
        PUT_ATTR_TXT("units", "m")
        PUT_ATTR_FILL(fillv)
        PUT_ATTR_FLT1("missing_value", missv)

        /* float DZLAKE(levlak, lat, lon) */
        dimids[0] = dim_levlak;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        DEF_VAR("DZLAKE", REC_XTYPE, 3, dimids, REC_ITYPE, 2)
        PUT_ATTR_TXT("long_name", "lake layer thickness")
        PUT_ATTR_TXT("units", "m")
        PUT_ATTR_FILL(fillv)
        PUT_ATTR_FLT1("missing_value", missv)
    }

    /* define record variables */

    /* char date_written(time, string_length) */
    dimids[0] = dim_time;
    dimids[1] = dim_string_length;
    DEF_VAR("date_written", NC_CHAR, 2, dimids, MPI_CHAR, -1)

    /* char time_written(time, string_length) */
    dimids[0] = dim_time;
    dimids[1] = dim_string_length;
    DEF_VAR("time_written", NC_CHAR, 2, dimids, MPI_CHAR, -1)

    /* int mcdate(time) */
    DEF_VAR("mcdate", NC_INT, 1, &dim_time, MPI_INT, -1)
    PUT_ATTR_TXT("long_name", "current date (YYYYMMDD)")

    /* int mcsec(time) */
    DEF_VAR("mcsec", NC_INT, 1, &dim_time, MPI_INT, -1)
    PUT_ATTR_TXT("long_name", "current seconds of current date")
    PUT_ATTR_TXT("units", "s")

    /* int mdcur(time) */
    DEF_VAR("mdcur", NC_INT, 1, &dim_time, MPI_INT, -1)
    PUT_ATTR_TXT("long_name", "current day (from base day)")

    /* int mscur(time) */
    DEF_VAR("mscur", NC_INT, 1, &dim_time, MPI_INT, -1)
    PUT_ATTR_TXT("long_name", "current seconds of current day")

    /* int nstep(time) */
    DEF_VAR("nstep", NC_INT, 1, &dim_time, MPI_INT, -1)
    PUT_ATTR_TXT("long_name", "time step")

    /* double time_bounds(time, hist_interval) */
    dimids[0] = dim_time;
    dimids[1] = dim_hist_interval;
    DEF_VAR("time_bounds", NC_DOUBLE, 2, dimids, MPI_DOUBLE, -1)
    PUT_ATTR_TXT("long_name", "history time interval endpoints")

    /* float time(time) */
    DEF_VAR("time", REC_XTYPE, 1, &dim_time, REC_ITYPE, -1)
    PUT_ATTR_TXT("long_name", "time")
    PUT_ATTR_TXT("units", "days since 0001-01-01 00:00:00")
    PUT_ATTR_TXT("calendar", "noleap")
    PUT_ATTR_TXT("bounds", "time_bounds")

    /* float ACTUAL_IMMOB(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ACTUAL_IMMOB", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "actual N immobilization")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float ACTUAL_IMMOB_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ACTUAL_IMMOB_P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "actual P immobilization")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float ADSORBTION_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ADSORBTION_P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "adsorb P flux")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float AGNPP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("AGNPP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "aboveground NPP")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float AGWDNPP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("AGWDNPP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "aboveground wood NPP")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float ALT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ALT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "current active layer thickness")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float ALTMAX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ALTMAX", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "maximum annual active layer thickness")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float ALTMAX_LASTYEAR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ALTMAX_LASTYEAR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "maximum prior year active layer thickness")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float AR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("AR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "autotrophic respiration (MR + GR)")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float AVAILC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("AVAILC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "C flux available for allocation")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float AVAIL_RETRANSP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("AVAIL_RETRANSP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "P flux available from retranslocation pool")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float BAF_CROP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("BAF_CROP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "fractional area burned for crop")
    PUT_ATTR_TXT("units", "proportion/sec")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float BAF_PEATF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("BAF_PEATF", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "fractional area burned in peatland")
    PUT_ATTR_TXT("units", "proportion/sec")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float BCDEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("BCDEP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total BC deposition (dry+wet) from atmosphere")
    PUT_ATTR_TXT("units", "kg/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float BGNPP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("BGNPP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "belowground NPP")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float BIOCHEM_PMIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("BIOCHEM_PMIN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "biochemical rate of P mineralization")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float BIOCHEM_PMIN_TO_PLANT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("BIOCHEM_PMIN_TO_PLANT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "plant uptake of biochemical P mineralization")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float BTRAN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("BTRAN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "transpiration beta factor")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float BUILDHEAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("BUILDHEAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "heat flux from urban building interior to walls and roof")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CH4PROD(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CH4PROD", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Gridcell total production of CH4")
    PUT_ATTR_TXT("units", "gC/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CH4_SURF_AERE_SAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CH4_SURF_AERE_SAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "aerenchyma surface CH4 flux for inundated area; (+ to atm)")
    PUT_ATTR_TXT("units", "mol/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CH4_SURF_AERE_UNSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CH4_SURF_AERE_UNSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "aerenchyma surface CH4 flux for non-inundated area; (+ to atm)")
    PUT_ATTR_TXT("units", "mol/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CH4_SURF_DIFF_SAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CH4_SURF_DIFF_SAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "diffusive surface CH4 flux for inundated / lake area; (+ to atm)")
    PUT_ATTR_TXT("units", "mol/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CH4_SURF_DIFF_UNSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CH4_SURF_DIFF_UNSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "diffusive surface CH4 flux for non-inundated area; (+ to atm)")
    PUT_ATTR_TXT("units", "mol/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CH4_SURF_EBUL_SAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CH4_SURF_EBUL_SAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "ebullition surface CH4 flux for inundated / lake area; (+ to atm)")
    PUT_ATTR_TXT("units", "mol/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CH4_SURF_EBUL_UNSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CH4_SURF_EBUL_UNSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "ebullition surface CH4 flux for non-inundated area; (+ to atm)")
    PUT_ATTR_TXT("units", "mol/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float COL_PTRUNC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("COL_PTRUNC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "column-level sink for P truncation")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CONC_CH4_SAT(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("CONC_CH4_SAT", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "CH4 soil Concentration for inundated / lake area")
    PUT_ATTR_TXT("units", "mol/m3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CONC_CH4_UNSAT(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("CONC_CH4_UNSAT", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "CH4 soil Concentration for non-inundated area")
    PUT_ATTR_TXT("units", "mol/m3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CONC_O2_SAT(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("CONC_O2_SAT", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "O2 soil Concentration for inundated / lake area")
    PUT_ATTR_TXT("units", "mol/m3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CONC_O2_UNSAT(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("CONC_O2_UNSAT", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "O2 soil Concentration for non-inundated area")
    PUT_ATTR_TXT("units", "mol/m3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CPOOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "temporary photosynthate C pool")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CWDC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CWDC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "CWD C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CWDC_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CWDC_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "coarse woody debris C heterotrophic respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CWDC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CWDC_LOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "coarse woody debris C loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CWDC_TO_LITR2C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CWDC_TO_LITR2C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of coarse woody debris C to litter 2 C")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CWDC_TO_LITR3C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CWDC_TO_LITR3C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of coarse woody debris C to litter 3 C")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CWDC_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("CWDC_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "CWD C (vertically resolved)")
    PUT_ATTR_TXT("units", "gC/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CWDN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CWDN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "CWD N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CWDN_TO_LITR2N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CWDN_TO_LITR2N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of coarse woody debris N to litter 2 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CWDN_TO_LITR3N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CWDN_TO_LITR3N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of coarse woody debris N to litter 3 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CWDN_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("CWDN_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "CWD N (vertically resolved)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CWDP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CWDP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "CWD P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CWDP_TO_LITR2P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CWDP_TO_LITR2P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of coarse woody debris P to litter 2 N")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CWDP_TO_LITR3P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("CWDP_TO_LITR3P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of coarse woody debris P to litter 3 N")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float CWDP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("CWDP_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "CWD P (vertically resolved)")
    PUT_ATTR_TXT("units", "gP/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DEADCROOTC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DEADCROOTC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "dead coarse root C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DEADCROOTN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DEADCROOTN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "dead coarse root N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DEADCROOTP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DEADCROOTP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "dead coarse root P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DEADSTEMC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DEADSTEMC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "dead stem C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DEADSTEMN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DEADSTEMN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "dead stem N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DEADSTEMP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DEADSTEMP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "dead stem P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DEFICIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DEFICIT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "runoff supply deficit")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DENIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DENIT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total rate of denitrification")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DESORPTION_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DESORPTION_P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "desorp P flux")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DISPVEGC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DISPVEGC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "displayed veg carbon, excluding storage and cpool")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DISPVEGN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DISPVEGN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "displayed vegetation nitrogen")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DISPVEGP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DISPVEGP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "displayed vegetation phosphorus")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DSTDEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DSTDEP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total dust deposition (dry+wet) from atmosphere")
    PUT_ATTR_TXT("units", "kg/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DSTFLXT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DSTFLXT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total surface dust emission")
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DWB(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DWB", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "net change in total water mass")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DWT_CONV_CFLUX_DRIBBLED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DWT_CONV_CFLUX_DRIBBLED", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "conversion C flux (immediate loss to atm), dribbled throughout the year")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DWT_CONV_CFLUX_GRC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DWT_CONV_CFLUX_GRC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DWT_CONV_NFLUX_GRC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DWT_CONV_NFLUX_GRC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DWT_CONV_PFLUX_GRC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DWT_CONV_PFLUX_GRC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DWT_SLASH_CFLUX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DWT_SLASH_CFLUX", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "slash C flux to litter and CWD due to land use")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DWT_SLASH_NFLUX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DWT_SLASH_NFLUX", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "slash N flux to litter and CWD due to land use")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float DWT_SLASH_PFLUX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("DWT_SLASH_PFLUX", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "slash P flux to litter and CWD due to land use")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float EFLX_DYNBAL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("EFLX_DYNBAL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "dynamic land cover change conversion energy flux")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float EFLX_GRND_LAKE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("EFLX_GRND_LAKE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "net heat flux into lake/snow surface, excluding light transmission")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float EFLX_LH_TOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("EFLX_LH_TOT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total latent heat flux [+ to atm]")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float EFLX_LH_TOT_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("EFLX_LH_TOT_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Rural total evaporation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float EFLX_LH_TOT_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("EFLX_LH_TOT_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Urban total evaporation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float ELAI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ELAI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "exposed one-sided leaf area index")
    PUT_ATTR_TXT("units", "m^2/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float ER(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ER", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total ecosystem respiration, autotrophic + heterotrophic")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float ERRH2O(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ERRH2O", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total water conservation error")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float ERRH2OSNO(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ERRH2OSNO", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "imbalance in snow depth (liquid water)")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float ERRSEB(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ERRSEB", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "surface energy conservation error")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float ERRSOI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ERRSOI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil/lake energy conservation error")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float ERRSOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ERRSOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "solar radiation conservation error")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float ESAI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ESAI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "exposed one-sided stem area index")
    PUT_ATTR_TXT("units", "m^2/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FAREA_BURNED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FAREA_BURNED", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "timestep fractional area burned")
    PUT_ATTR_TXT("units", "proportion")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FCEV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FCEV", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "canopy evaporation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FCH4(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FCH4", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Gridcell surface CH4 flux to atmosphere (+ to atm)")
    PUT_ATTR_TXT("units", "kgC/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FCH4TOCO2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FCH4TOCO2", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Gridcell oxidation of CH4 to CO2")
    PUT_ATTR_TXT("units", "gC/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FCH4_DFSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FCH4_DFSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "CH4 additional flux due to changing fsat, vegetated landunits only")
    PUT_ATTR_TXT("units", "kgC/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FCOV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FCOV", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "fractional impermeable area")
    PUT_ATTR_TXT("units", "unitless")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FCTR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FCTR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "canopy transpiration")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FGEV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FGEV", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "ground evaporation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FGR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FGR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "heat flux into soil/snow including snow melt and lake / snow light transmission")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FGR12(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FGR12", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "heat flux between soil layers 1 and 2")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FGR_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FGR_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Rural heat flux into soil/snow including snow melt and snow light transmission")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FGR_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FGR_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Urban heat flux into soil/snow including snow melt")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FH2OSFC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FH2OSFC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "fraction of ground covered by surface water")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FINUNDATED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FINUNDATED", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "fractional inundated area of vegetated columns")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FINUNDATED_LAG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FINUNDATED_LAG", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "time-lagged inundated fraction of vegetated columns")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FIRA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FIRA", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "net infrared (longwave) radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FIRA_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FIRA_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Rural net infrared (longwave) radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FIRA_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FIRA_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Urban net infrared (longwave) radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FIRE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FIRE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "emitted infrared (longwave) radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FIRE_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FIRE_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Rural emitted infrared (longwave) radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FIRE_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FIRE_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Urban emitted infrared (longwave) radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FLDS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FLDS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "atmospheric longwave radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FPG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FPG", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "fraction of potential gpp due to N limitation")
    PUT_ATTR_TXT("units", "proportion")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FPG_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FPG_P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "fraction of potential gpp due to P limitation")
    PUT_ATTR_TXT("units", "proportion")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FPI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FPI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "fraction of potential immobilization of nitrogen")
    PUT_ATTR_TXT("units", "proportion")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FPI_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FPI_P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "fraction of potential immobilization of phosphorus")
    PUT_ATTR_TXT("units", "proportion")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FPI_P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("FPI_P_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "fraction of potential immobilization of phosphorus")
    PUT_ATTR_TXT("units", "proportion")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FPI_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("FPI_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "fraction of potential immobilization of nitrogen")
    PUT_ATTR_TXT("units", "proportion")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FPSN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FPSN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "photosynthesis")
    PUT_ATTR_TXT("units", "umol/m2s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FPSN_WC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FPSN_WC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Rubisco-limited photosynthesis")
    PUT_ATTR_TXT("units", "umol/m2s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FPSN_WJ(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FPSN_WJ", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "RuBP-limited photosynthesis")
    PUT_ATTR_TXT("units", "umol/m2s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FPSN_WP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FPSN_WP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Product-limited photosynthesis")
    PUT_ATTR_TXT("units", "umol/m2s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FROOTC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FROOTC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "fine root C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FROOTC_ALLOC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FROOTC_ALLOC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "fine root C allocation")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FROOTC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FROOTC_LOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "fine root C loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FROOTN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FROOTN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "fine root N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FROOTP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FROOTP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "fine root P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FROST_TABLE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FROST_TABLE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "frost table depth (vegetated landunits only)")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSA", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "absorbed solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "fractional area with water table at surface")
    PUT_ATTR_TXT("units", "unitless")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSA_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSA_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Rural absorbed solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSA_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSA_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Urban absorbed solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSDS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSDS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "atmospheric incident solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSDSND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSDSND", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "direct nir incident solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSDSNDLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSDSNDLN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "direct nir incident solar radiation at local noon")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSDSNI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSDSNI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "diffuse nir incident solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSDSVD(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSDSVD", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "direct vis incident solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSDSVDLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSDSVDLN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "direct vis incident solar radiation at local noon")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSDSVI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSDSVI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "diffuse vis incident solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSDSVILN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSDSVILN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "diffuse vis incident solar radiation at local noon")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSH(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSH", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "sensible heat")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSH_G(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSH_G", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "sensible heat from ground")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSH_NODYNLNDUSE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSH_NODYNLNDUSE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "sensible heat not including correction for land use change")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSH_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSH_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Rural sensible heat")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSH_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSH_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Urban sensible heat")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSH_V(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSH_V", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "sensible heat from veg")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSM(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSM", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "snow melt heat flux")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSM_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSM_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Rural snow melt heat flux")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSM_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSM_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Urban snow melt heat flux")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSNO(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSNO", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "fraction of ground covered by snow")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSNO_EFF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSNO_EFF", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "effective fraction of ground covered by snow")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "reflected solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSRND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSRND", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "direct nir reflected solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSRNDLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSRNDLN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "direct nir reflected solar radiation at local noon")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSRNI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSRNI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "diffuse nir reflected solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSRVD(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSRVD", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "direct vis reflected solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSRVDLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSRVDLN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "direct vis reflected solar radiation at local noon")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float FSRVI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("FSRVI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "diffuse vis reflected solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float F_CO2_SOIL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("F_CO2_SOIL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total soil-atm. CO2 exchange")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float F_CO2_SOIL_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("F_CO2_SOIL_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "total vertically resolved soil-atm. CO2 exchange")
    PUT_ATTR_TXT("units", "gC/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float F_DENIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("F_DENIT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "denitrification flux")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float F_DENIT_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("F_DENIT_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "denitrification flux")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float F_N2O_DENIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("F_N2O_DENIT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "denitrification N2O flux")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float F_N2O_NIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("F_N2O_NIT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "nitrification N2O flux")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float F_NIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("F_NIT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "nitrification flux")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float F_NIT_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("F_NIT_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "nitrification flux")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float GC_HEAT1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("GC_HEAT1", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "initial gridcell total heat content")
    PUT_ATTR_TXT("units", "J/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float GC_ICE1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("GC_ICE1", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "initial gridcell total ice content")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float GC_LIQ1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("GC_LIQ1", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "initial gridcell total liq content")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float GPP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("GPP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "gross primary production")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float GR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("GR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total growth respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float GROSS_NMIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("GROSS_NMIN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "gross rate of N mineralization")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float GROSS_PMIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("GROSS_PMIN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "gross rate of P mineralization")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float H2OCAN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("H2OCAN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "intercepted water")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float H2OSFC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("H2OSFC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "surface water depth")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float H2OSNO(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("H2OSNO", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "snow depth (liquid water)")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float H2OSNO_TOP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("H2OSNO_TOP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mass of snow in top snow layer")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float H2OSOI(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("H2OSOI", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "volumetric soil water (vegetated landunits only)")
    PUT_ATTR_TXT("units", "mm3/mm3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float HC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("HC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "heat content of soil/snow/lake")
    PUT_ATTR_TXT("units", "MJ/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float HCSOI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("HCSOI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil heat content")
    PUT_ATTR_TXT("units", "MJ/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float HEAT_FROM_AC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("HEAT_FROM_AC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "sensible heat flux put into canyon due to heat removed from air conditioning")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total heterotrophic respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float HR_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("HR_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "total vertically resolved heterotrophic respiration")
    PUT_ATTR_TXT("units", "gC/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float HTOP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("HTOP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "canopy top")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float INT_SNOW(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("INT_SNOW", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "accumulated swe (vegetated landunits only)")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LABILEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LABILEP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil Labile P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LABILEP_TO_SECONDP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LABILEP_TO_SECONDP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "LABILE P TO SECONDARY MINERAL P")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LABILEP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("LABILEP_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil labile P (vert. res.)")
    PUT_ATTR_TXT("units", "gp/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LAISHA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LAISHA", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "shaded projected leaf area index")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LAISUN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LAISUN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "sunlit projected leaf area index")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LAKEICEFRAC(time, levlak, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levlak;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("LAKEICEFRAC", REC_XTYPE, 4, dimids, REC_ITYPE, 2)
    PUT_ATTR_TXT("long_name", "lake layer ice mass fraction")
    PUT_ATTR_TXT("units", "unitless")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LAKEICETHICK(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LAKEICETHICK", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "thickness of lake ice (including physical expansion on freezing)")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LAND_UPTAKE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LAND_UPTAKE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "NEE minus LAND_USE_FLUX, negative for update")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LAND_USE_FLUX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LAND_USE_FLUX", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total C emitted from land cover conversion and wood product pools")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LEAFC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LEAFC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "leaf C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LEAFC_ALLOC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LEAFC_ALLOC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "leaf C allocation")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LEAFC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LEAFC_LOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "leaf C loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LEAFC_TO_LITTER(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LEAFC_TO_LITTER", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "leaf C litterfall")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LEAFN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LEAFN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "leaf N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LEAFP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LEAFP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "leaf P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LEAF_MR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LEAF_MR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "leaf maintenance respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LFC2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LFC2", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "conversion area fraction of BET and BDT that burned")
    PUT_ATTR_TXT("units", "per sec")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITFALL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITFALL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "litterfall (leaves and fine roots)")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITHR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITHR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "litter heterotrophic respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR1C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR1C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "LITR1 C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR1C_TO_SOIL1C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR1C_TO_SOIL1C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of litter 1 C to soil 1 C")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR1C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("LITR1C_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "LITR1 C (vertically resolved)")
    PUT_ATTR_TXT("units", "gC/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR1N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR1N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "LITR1 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR1N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("LITR1N_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "litter 1 N tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR1N_TO_SOIL1N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR1N_TO_SOIL1N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of litter 1 N to soil 1 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR1N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("LITR1N_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "LITR1 N (vertically resolved)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR1P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR1P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "LITR1 P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR1P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("LITR1P_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "litter 1 P tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gP/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR1P_TO_SOIL1P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR1P_TO_SOIL1P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of litter 1 P to soil 1 N")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR1P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("LITR1P_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "LITR1 P (vertically resolved)")
    PUT_ATTR_TXT("units", "gP/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR1_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR1_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Het. Resp. from litter 1")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR2C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR2C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "LITR2 C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR2C_TO_SOIL2C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR2C_TO_SOIL2C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of litter 2 C to soil 2 C")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR2C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("LITR2C_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "LITR2 C (vertically resolved)")
    PUT_ATTR_TXT("units", "gC/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR2N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR2N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "LITR2 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR2N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("LITR2N_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "litter 2 N tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR2N_TO_SOIL2N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR2N_TO_SOIL2N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of litter 2 N to soil 2 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR2N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("LITR2N_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "LITR2 N (vertically resolved)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR2P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR2P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "LITR2 P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR2P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("LITR2P_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "litter 2 P tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gP/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR2P_TO_SOIL2P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR2P_TO_SOIL2P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of litter 2 P to soil 2 N")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR2P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("LITR2P_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "LITR2 P (vertically resolved)")
    PUT_ATTR_TXT("units", "gP/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR2_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR2_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Het. Resp. from litter 2")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR3C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR3C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "LITR3 C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR3C_TO_SOIL3C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR3C_TO_SOIL3C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of litter 3 C to soil 3 C")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR3C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("LITR3C_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "LITR3 C (vertically resolved)")
    PUT_ATTR_TXT("units", "gC/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR3N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR3N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "LITR3 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR3N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("LITR3N_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "litter 3 N tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR3N_TO_SOIL3N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR3N_TO_SOIL3N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of litter 3 N to soil 3 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR3N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("LITR3N_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "LITR3 N (vertically resolved)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR3P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR3P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "LITR3 P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR3P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("LITR3P_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "litter 3 P tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gP/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR3P_TO_SOIL3P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR3P_TO_SOIL3P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of litter 3 P to soil 3 N")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR3P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("LITR3P_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "LITR3 P (vertically resolved)")
    PUT_ATTR_TXT("units", "gP/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITR3_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITR3_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Het. Resp. from litter 3")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITTERC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITTERC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "litter C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITTERC_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITTERC_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "litter C heterotrophic respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LITTERC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LITTERC_LOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "litter C loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LIVECROOTC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LIVECROOTC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "live coarse root C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LIVECROOTN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LIVECROOTN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "live coarse root N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LIVECROOTP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LIVECROOTP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "live coarse root P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LIVESTEMC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LIVESTEMC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "live stem C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LIVESTEMN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LIVESTEMN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "live stem N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float LIVESTEMP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("LIVESTEMP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "live stem P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float MR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("MR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "maintenance respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float M_LITR1C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("M_LITR1C_TO_LEACHING", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "litter 1 C leaching loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float M_LITR2C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("M_LITR2C_TO_LEACHING", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "litter 2 C leaching loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float M_LITR3C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("M_LITR3C_TO_LEACHING", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "litter 3 C leaching loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float M_SOIL1C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("M_SOIL1C_TO_LEACHING", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil 1 C leaching loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float M_SOIL2C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("M_SOIL2C_TO_LEACHING", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil 2 C leaching loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float M_SOIL3C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("M_SOIL3C_TO_LEACHING", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil 3 C leaching loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float M_SOIL4C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("M_SOIL4C_TO_LEACHING", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil 4 C leaching loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float NBP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("NBP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "net biome production, includes fire, landuse, and harvest flux, positive for sink")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float NDEPLOY(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("NDEPLOY", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total N deployed in new growth")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float NDEP_TO_SMINN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("NDEP_TO_SMINN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "atmospheric N deposition to soil mineral N")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float NEE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("NEE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "net ecosystem exchange of carbon, includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float NEM(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("NEM", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Gridcell net adjustment to NEE passed to atm. for methane production")
    PUT_ATTR_TXT("units", "gC/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float NEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("NEP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "net ecosystem production, excludes fire, landuse, and harvest flux, positive for sink")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float NET_NMIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("NET_NMIN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "net rate of N mineralization")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float NET_PMIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("NET_PMIN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "net rate of P mineralization")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float NFIRE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("NFIRE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "fire counts valid only in Reg.C")
    PUT_ATTR_TXT("units", "counts/km2/sec")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float NFIX_TO_SMINN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("NFIX_TO_SMINN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "symbiotic/asymbiotic N fixation to soil mineral N")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float NPP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("NPP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "net primary production")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float OCCLP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("OCCLP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil occluded P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float OCCLP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("OCCLP_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil occluded P (vert. res.)")
    PUT_ATTR_TXT("units", "gp/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float OCDEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("OCDEP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total OC deposition (dry+wet) from atmosphere")
    PUT_ATTR_TXT("units", "kg/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float O_SCALAR(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("O_SCALAR", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "fraction by which decomposition is reduced due to anoxia")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PARVEGLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PARVEGLN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "absorbed par by vegetation at local noon")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PBOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PBOT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "atmospheric pressure")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PCH4(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PCH4", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "atmospheric partial pressure of CH4")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PCO2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PCO2", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "atmospheric partial pressure of CO2")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PCT_LANDUNIT(time, ltype, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_ltype;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("PCT_LANDUNIT", REC_XTYPE, 4, dimids, REC_ITYPE, 3)
    PUT_ATTR_TXT("long_name", "% of each landunit on grid cell")
    PUT_ATTR_TXT("units", "%")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PCT_NAT_PFT(time, natpft, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_natpft;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("PCT_NAT_PFT", REC_XTYPE, 4, dimids, REC_ITYPE, 4)
    PUT_ATTR_TXT("long_name", "% of each PFT on the natural vegetation (i.e., soil) landunit")
    PUT_ATTR_TXT("units", "%")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PDEPLOY(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PDEPLOY", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total P deployed in new growth")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PDEP_TO_SMINP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PDEP_TO_SMINP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "atmospheric P deposition to soil mineral P")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PFT_FIRE_CLOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PFT_FIRE_CLOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total patch-level fire C loss for non-peat fires outside land-type converted region")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PFT_FIRE_NLOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PFT_FIRE_NLOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total pft-level fire N loss")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PLANT_CALLOC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PLANT_CALLOC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total allocated C flux")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PLANT_NDEMAND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PLANT_NDEMAND", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "N flux required to support initial GPP")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PLANT_NDEMAND_COL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PLANT_NDEMAND_COL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "N flux required to support initial GPP")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PLANT_PALLOC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PLANT_PALLOC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total allocated P flux")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PLANT_PDEMAND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PLANT_PDEMAND", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "P flux required to support initial GPP")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PLANT_PDEMAND_COL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PLANT_PDEMAND_COL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "P flux required to support initial GPP")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float POTENTIAL_IMMOB(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("POTENTIAL_IMMOB", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "potential N immobilization")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float POTENTIAL_IMMOB_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("POTENTIAL_IMMOB_P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "potential P immobilization")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float POT_F_DENIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("POT_F_DENIT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "potential denitrification flux")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float POT_F_NIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("POT_F_NIT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "potential nitrification flux")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PRIMP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PRIMP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil primary P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PRIMP_TO_LABILEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PRIMP_TO_LABILEP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "PRIMARY MINERAL P TO LABILE P")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PRIMP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("PRIMP_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil primary P (vert. res.)")
    PUT_ATTR_TXT("units", "gp/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PROD1P_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PROD1P_LOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "loss from 1-yr crop product pool")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PSNSHA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PSNSHA", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "shaded leaf photosynthesis")
    PUT_ATTR_TXT("units", "umolCO2/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PSNSHADE_TO_CPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PSNSHADE_TO_CPOOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "C fixation from shaded canopy")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PSNSUN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PSNSUN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "sunlit leaf photosynthesis")
    PUT_ATTR_TXT("units", "umolCO2/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float PSNSUN_TO_CPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("PSNSUN_TO_CPOOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "C fixation from sunlit canopy")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float Q2M(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("Q2M", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "2m specific humidity")
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QBOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QBOT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "atmospheric specific humidity")
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QCHARGE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QCHARGE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "aquifer recharge rate (vegetated landunits only)")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QDRAI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QDRAI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "sub-surface drainage")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QDRAI_PERCH(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QDRAI_PERCH", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "perched wt drainage")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QDRAI_XS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QDRAI_XS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "saturation excess drainage")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QDRIP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QDRIP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "throughfall")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QFLOOD(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QFLOOD", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "runoff from river flooding")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QFLX_ICE_DYNBAL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QFLX_ICE_DYNBAL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "ice dynamic land cover change conversion runoff flux")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QFLX_LIQ_DYNBAL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QFLX_LIQ_DYNBAL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "liq dynamic land cover change conversion runoff flux")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QH2OSFC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QH2OSFC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "surface water runoff")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QINFL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QINFL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "infiltration")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QINTR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QINTR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "interception")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QIRRIG_GRND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QIRRIG_GRND", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Groundwater irrigation")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QIRRIG_ORIG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QIRRIG_ORIG", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Original total irrigation water demand (surface + ground)")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QIRRIG_REAL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QIRRIG_REAL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "actual water added through irrigation (surface + ground)")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QIRRIG_SURF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QIRRIG_SURF", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Surface water irrigation")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QIRRIG_WM(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QIRRIG_WM", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Surface water irrigation demand sent to MOSART/WM")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QOVER(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QOVER", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "surface runoff")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QOVER_LAG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QOVER_LAG", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "time-lagged surface runoff for soil columns")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QRGWL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QRGWL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "surface runoff at glaciers (liquid only), wetlands, lakes")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QRUNOFF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QRUNOFF", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total liquid runoff (does not include QSNWCPICE)")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QRUNOFF_NODYNLNDUSE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QRUNOFF_NODYNLNDUSE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total liquid runoff (does not include QSNWCPICE) not including correction for land use change")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QRUNOFF_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QRUNOFF_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Rural total runoff")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QRUNOFF_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QRUNOFF_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Urban total runoff")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QSNOMELT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QSNOMELT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "snow melt")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QSNWCPICE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QSNWCPICE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "excess snowfall due to snow capping")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QSNWCPICE_NODYNLNDUSE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QSNWCPICE_NODYNLNDUSE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "excess snowfall due to snow capping not including correction for land use change")
    PUT_ATTR_TXT("units", "mm H2O/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QSOIL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QSOIL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Ground evaporation (soil/snow evaporation + soil/snow sublimation - dew)")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QVEGE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QVEGE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "canopy evaporation")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float QVEGT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("QVEGT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "canopy transpiration")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float RAIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("RAIN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "atmospheric rain")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float RETRANSN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("RETRANSN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "plant pool of retranslocated N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float RETRANSN_TO_NPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("RETRANSN_TO_NPOOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "deployment of retranslocated N")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float RETRANSP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("RETRANSP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "plant pool of retranslocated P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float RETRANSP_TO_PPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("RETRANSP_TO_PPOOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "deployment of retranslocated P")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float RH2M(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("RH2M", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "2m relative humidity")
    PUT_ATTR_TXT("units", "%")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float RH2M_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("RH2M_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Rural 2m specific humidity")
    PUT_ATTR_TXT("units", "%")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float RH2M_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("RH2M_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Urban 2m relative humidity")
    PUT_ATTR_TXT("units", "%")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float RR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("RR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "root respiration (fine root MR + total root GR)")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SABG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SABG", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "solar rad absorbed by ground")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SABG_PEN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SABG_PEN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Rural solar rad penetrating top soil or snow layer")
    PUT_ATTR_TXT("units", "watt/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SABV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SABV", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "solar rad absorbed by veg")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SCALARAVG_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SCALARAVG_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "average of decomposition scalar")
    PUT_ATTR_TXT("units", "fraction")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SECONDP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SECONDP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil secondary P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SECONDP_TO_LABILEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SECONDP_TO_LABILEP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "SECONDARY MINERAL P TO LABILE P")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SECONDP_TO_OCCLP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SECONDP_TO_OCCLP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "SECONDARY MINERAL P TO OCCLUDED P")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SECONDP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SECONDP_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil secondary P (vert. res.)")
    PUT_ATTR_TXT("units", "gp/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SEEDC_GRC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SEEDC_GRC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "pool for seeding new PFTs via dynamic landcover")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil mineral N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINN_TO_NPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINN_TO_NPOOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "deployment of soil mineral N uptake")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINN_TO_PLANT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINN_TO_PLANT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "plant uptake of soil mineral N")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINN_TO_SOIL1N_L1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINN_TO_SOIL1N_L1", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mineral N flux for decomp. of LITR1to SOIL1")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINN_TO_SOIL2N_L2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINN_TO_SOIL2N_L2", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mineral N flux for decomp. of LITR2to SOIL2")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINN_TO_SOIL2N_S1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINN_TO_SOIL2N_S1", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mineral N flux for decomp. of SOIL1to SOIL2")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINN_TO_SOIL3N_L3(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINN_TO_SOIL3N_L3", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mineral N flux for decomp. of LITR3to SOIL3")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINN_TO_SOIL3N_S2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINN_TO_SOIL3N_S2", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mineral N flux for decomp. of SOIL2to SOIL3")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINN_TO_SOIL4N_S3(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINN_TO_SOIL4N_S3", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mineral N flux for decomp. of SOIL3to SOIL4")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil mineral P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINP_LEACHED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINP_LEACHED", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil mineral P pool loss to leaching")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINP_TO_PLANT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINP_TO_PLANT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "plant uptake of soil mineral P")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINP_TO_PPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINP_TO_PPOOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "deployment of soil mineral P uptake")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINP_TO_SOIL1P_L1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINP_TO_SOIL1P_L1", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mineral P flux for decomp. of LITR1to SOIL1")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINP_TO_SOIL2P_L2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINP_TO_SOIL2P_L2", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mineral P flux for decomp. of LITR2to SOIL2")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINP_TO_SOIL2P_S1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINP_TO_SOIL2P_S1", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mineral P flux for decomp. of SOIL1to SOIL2")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINP_TO_SOIL3P_L3(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINP_TO_SOIL3P_L3", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mineral P flux for decomp. of LITR3to SOIL3")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINP_TO_SOIL3P_S2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINP_TO_SOIL3P_S2", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mineral P flux for decomp. of SOIL2to SOIL3")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINP_TO_SOIL4P_S3(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMINP_TO_SOIL4P_S3", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mineral P flux for decomp. of SOIL3to SOIL4")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMINP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SMINP_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil mineral P (vert. res.)")
    PUT_ATTR_TXT("units", "gp/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMIN_NH4(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMIN_NH4", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil mineral NH4")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMIN_NH4_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SMIN_NH4_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil mineral NH4 (vert. res.)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMIN_NO3(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMIN_NO3", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil mineral NO3")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMIN_NO3_LEACHED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMIN_NO3_LEACHED", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil NO3 pool loss to leaching")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMIN_NO3_RUNOFF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SMIN_NO3_RUNOFF", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil NO3 pool loss to runoff")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SMIN_NO3_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SMIN_NO3_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil mineral NO3 (vert. res.)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SNOBCMCL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SNOBCMCL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mass of BC in snow column")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SNOBCMSL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SNOBCMSL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mass of BC in top snow layer")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SNODSTMCL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SNODSTMCL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mass of dust in snow column")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SNODSTMSL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SNODSTMSL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mass of dust in top snow layer")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SNOINTABS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SNOINTABS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Percent of incoming solar absorbed by lower snow layers")
    PUT_ATTR_TXT("units", "%")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SNOOCMCL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SNOOCMCL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mass of OC in snow column")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SNOOCMSL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SNOOCMSL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mass of OC in top snow layer")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SNOW(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SNOW", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "atmospheric snow")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SNOWDP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SNOWDP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "gridcell mean snow height")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SNOWICE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SNOWICE", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "snow ice")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SNOWLIQ(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SNOWLIQ", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "snow liquid water")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SNOW_DEPTH(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SNOW_DEPTH", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "snow height of snow covered area")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SNOW_SINKS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SNOW_SINKS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "snow sinks (liquid water)")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SNOW_SOURCES(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SNOW_SOURCES", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "snow sources (liquid water)")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL1C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL1C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "SOIL1 C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL1C_TO_SOIL2C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL1C_TO_SOIL2C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of soil 1 C to soil 2 C")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL1C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL1C_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "SOIL1 C (vertically resolved)")
    PUT_ATTR_TXT("units", "gC/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL1N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL1N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "SOIL1 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL1N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL1N_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil 1 N tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL1N_TO_SOIL2N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL1N_TO_SOIL2N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of soil 1 N to soil 2 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL1N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL1N_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "SOIL1 N (vertically resolved)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL1P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL1P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "SOIL1 P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL1P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL1P_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil 1 P tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gP/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL1P_TO_SOIL2P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL1P_TO_SOIL2P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of soil 1 P to soil 2 N")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL1P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL1P_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "SOIL1 P (vertically resolved)")
    PUT_ATTR_TXT("units", "gP/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL1_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL1_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Het. Resp. from soil 1")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL2C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL2C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "SOIL2 C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL2C_TO_SOIL3C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL2C_TO_SOIL3C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of soil 2 C to soil 3 C")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL2C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL2C_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "SOIL2 C (vertically resolved)")
    PUT_ATTR_TXT("units", "gC/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL2N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL2N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "SOIL2 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL2N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL2N_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil 2 N tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL2N_TO_SOIL3N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL2N_TO_SOIL3N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of soil 2 N to soil 3 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL2N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL2N_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "SOIL2 N (vertically resolved)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL2P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL2P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "SOIL2 P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL2P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL2P_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil 2 P tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gP/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL2P_TO_SOIL3P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL2P_TO_SOIL3P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of soil 2 P to soil 3 N")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL2P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL2P_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "SOIL2 P (vertically resolved)")
    PUT_ATTR_TXT("units", "gP/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL2_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL2_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Het. Resp. from soil 2")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL3C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL3C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "SOIL3 C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL3C_TO_SOIL4C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL3C_TO_SOIL4C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of soil 3 C to soil 4 C")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL3C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL3C_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "SOIL3 C (vertically resolved)")
    PUT_ATTR_TXT("units", "gC/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL3N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL3N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "SOIL3 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL3N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL3N_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil 3 N tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL3N_TO_SOIL4N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL3N_TO_SOIL4N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of soil 3 N to soil 4 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL3N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL3N_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "SOIL3 N (vertically resolved)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL3P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL3P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "SOIL3 P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL3P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL3P_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil 3 P tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gP/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL3P_TO_SOIL4P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL3P_TO_SOIL4P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "decomp. of soil 3 P to soil 4 N")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL3P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL3P_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "SOIL3 P (vertically resolved)")
    PUT_ATTR_TXT("units", "gP/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL3_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL3_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Het. Resp. from soil 3")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL4C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL4C", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "SOIL4 C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL4C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL4C_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "SOIL4 C (vertically resolved)")
    PUT_ATTR_TXT("units", "gC/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL4N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL4N", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "SOIL4 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL4N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL4N_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil 4 N tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL4N_TO_SMINN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL4N_TO_SMINN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mineral N flux for decomp. of SOIL4")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL4N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL4N_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "SOIL4 N (vertically resolved)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL4P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL4P", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "SOIL4 P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL4P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL4P_TNDNCY_VERT_TRANS", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil 4 P tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gP/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL4P_TO_SMINP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL4P_TO_SMINP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "mineral P flux for decomp. of SOIL4")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL4P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOIL4P_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "SOIL4 P (vertically resolved)")
    PUT_ATTR_TXT("units", "gP/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOIL4_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOIL4_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Het. Resp. from soil 4")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOILC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOILC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOILC_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOILC_HR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil C heterotrophic respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOILC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOILC_LOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil C loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOILICE(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOILICE", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil ice (vegetated landunits only)")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOILICE_ICE(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOILICE_ICE", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil ice (ice landunits only)")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOILLIQ(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOILLIQ", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil liquid water (vegetated landunits only)")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOILLIQ_ICE(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOILLIQ_ICE", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil liquid water (ice landunits only)")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOILPSI(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOILPSI", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil water potential in each soil layer")
    PUT_ATTR_TXT("units", "MPa")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOILWATER_10CM(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOILWATER_10CM", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil liquid water + ice in top 10cm of soil (veg landunits only)")
    PUT_ATTR_TXT("standard_name", "mass_content_of_water_in_soil_layer")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOLUTIONP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOLUTIONP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil solution P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOLUTIONP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("SOLUTIONP_vr", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil solution P (vert. res.)")
    PUT_ATTR_TXT("units", "gp/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOMHR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOMHR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil organic matter heterotrophic respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SOM_C_LEACHED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SOM_C_LEACHED", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total flux of C from SOM pools due to leaching")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total soil respiration (HR + root resp)")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float STORVEGC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("STORVEGC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "stored vegetation carbon, excluding cpool")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float STORVEGN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("STORVEGN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "stored vegetation nitrogen")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float STORVEGP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("STORVEGP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "stored vegetation phosphorus")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SUPPLEMENT_TO_SMINN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SUPPLEMENT_TO_SMINN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "supplemental N supply")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SUPPLEMENT_TO_SMINP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SUPPLEMENT_TO_SMINP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "supplemental P supply")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SUPPLY(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SUPPLY", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "runoff supply for land use")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SoilAlpha(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SoilAlpha", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "factor limiting ground evap")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float SoilAlpha_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("SoilAlpha_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "urban factor limiting ground evap")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TAUX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TAUX", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "zonal surface stress")
    PUT_ATTR_TXT("units", "kg/m/s^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TAUY(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TAUY", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "meridional surface stress")
    PUT_ATTR_TXT("units", "kg/m/s^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TBOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TBOT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "atmospheric air temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TBUILD(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TBUILD", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "internal urban building temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TCS_MONTH_BEGIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TCS_MONTH_BEGIN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total carbon storage at the beginning of a month")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TCS_MONTH_END(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TCS_MONTH_END", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total carbon storage at the end of a month")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TG", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "ground temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TG_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TG_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Rural ground temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TG_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TG_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Urban ground temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TH2OSFC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TH2OSFC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "surface water temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float THBOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("THBOT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "atmospheric air potential temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TKE1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TKE1", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "top lake level eddy thermal conductivity")
    PUT_ATTR_TXT("units", "W/(mK)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TLAI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TLAI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total projected leaf area index")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TLAKE(time, levlak, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levlak;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("TLAKE", REC_XTYPE, 4, dimids, REC_ITYPE, 2)
    PUT_ATTR_TXT("long_name", "lake temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTCOLC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTCOLC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total column carbon, incl veg and cpool but excl product pools")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTCOLCH4(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTCOLCH4", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total belowground CH4, (0 for non-lake special landunits)")
    PUT_ATTR_TXT("units", "gC/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTCOLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTCOLN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total column-level N but excl product pools")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTCOLP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTCOLP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total column-level P but excl product pools")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTECOSYSC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTECOSYSC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total ecosystem carbon, incl veg but excl cpool but excl product pools")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTECOSYSN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTECOSYSN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total ecosystem N but excl product pools")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTECOSYSP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTECOSYSP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total ecosystem P but excl product pools")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTLITC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTLITC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total litter carbon")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTLITC_1m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTLITC_1m", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total litter carbon to 1 meter depth")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTLITN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTLITN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total litter N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTLITP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTLITP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total litter P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTLITP_1m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTLITP_1m", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total litter P to 1 meter")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTPFTC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTPFTC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total patch-level carbon, including cpool")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTPFTN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTPFTN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total PFT-level nitrogen")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTPFTP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTPFTP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total PFT-level phosphorus")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTSOMC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTSOMC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total soil organic matter carbon")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTSOMC_1m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTSOMC_1m", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total soil organic matter carbon to 1 meter depth")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTSOMN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTSOMN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total soil organic matter N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTSOMP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTSOMP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total soil organic matter P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTSOMP_1m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTSOMP_1m", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total soil organic matter P to 1 meter")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTVEGC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTVEGC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total vegetation carbon, excluding cpool")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTVEGC_ABG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTVEGC_ABG", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total aboveground vegetation carbon, excluding cpool")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTVEGN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTVEGN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total vegetation nitrogen")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TOTVEGP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TOTVEGP", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total vegetation phosphorus")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TREFMNAV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TREFMNAV", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "daily minimum of average 2-m temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TREFMNAV_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TREFMNAV_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Rural daily minimum of average 2-m temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TREFMNAV_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TREFMNAV_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Urban daily minimum of average 2-m temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TREFMXAV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TREFMXAV", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "daily maximum of average 2-m temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TREFMXAV_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TREFMXAV_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Rural daily maximum of average 2-m temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TREFMXAV_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TREFMXAV_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Urban daily maximum of average 2-m temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TSA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TSA", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "2m air temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TSAI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TSAI", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total projected stem area index")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TSA_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TSA_R", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Rural 2m air temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TSA_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TSA_U", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "Urban 2m air temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TSOI(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("TSOI", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil temperature (vegetated landunits only)")
    PUT_ATTR_TXT("standard_name", "soil_temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TSOI_10CM(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TSOI_10CM", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil temperature in top 10cm of soil")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TSOI_ICE(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("TSOI_ICE", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "soil temperature (ice landunits only)")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TV", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "vegetation temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TWS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TWS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total water storage")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TWS_MONTH_BEGIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TWS_MONTH_BEGIN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total water storage at the beginning of a month")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float TWS_MONTH_END(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("TWS_MONTH_END", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total water storage at the end of a month")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float T_SCALAR(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("T_SCALAR", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "temperature inhibition of decomposition")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float U10(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("U10", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "10-m wind")
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float URBAN_AC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("URBAN_AC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "urban air conditioning flux")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float URBAN_HEAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("URBAN_HEAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "urban heating flux")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float VOLR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("VOLR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "river channel total water storage")
    PUT_ATTR_TXT("units", "m3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float VOLRMCH(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("VOLRMCH", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "river channel main channel water storage")
    PUT_ATTR_TXT("units", "m3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float WA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("WA", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "water in the unconfined aquifer (vegetated landunits only)")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float WASTEHEAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("WASTEHEAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "sensible heat flux from heating/cooling sources of urban waste heat")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float WF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("WF", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "soil water as frac. of whc for top 0.05 m")
    PUT_ATTR_TXT("units", "proportion")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float WIND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("WIND", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "atmospheric wind velocity magnitude")
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float WOODC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("WOODC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "wood C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float WOODC_ALLOC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("WOODC_ALLOC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "wood C eallocation")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float WOODC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("WOODC_LOSS", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "wood C loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float WOOD_HARVESTC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("WOOD_HARVESTC", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "wood harvest carbon (to product pools)")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float WOOD_HARVESTN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("WOOD_HARVESTN", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "wood harvest N (to product pools)")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float WTGQ(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("WTGQ", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "surface tracer conductance")
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float W_SCALAR(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("W_SCALAR", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "Moisture (dryness) inhibition of decomposition")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float XR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("XR", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "total excess respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float XSMRPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("XSMRPOOL", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "temporary photosynthate C pool")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float ZBOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ZBOT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "atmospheric reference height")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float ZWT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ZWT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "water table depth (vegetated landunits only)")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float ZWT_CH4_UNSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ZWT_CH4_UNSAT", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "depth of water table for methane production used in non-inundated area")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float ZWT_PERCH(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("ZWT_PERCH", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "perched water table depth (vegetated landunits only)")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float cn_scalar(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("cn_scalar", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "N limitation factor")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float cp_scalar(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("cp_scalar", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "P limitation factor")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float leaf_npimbalance(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("leaf_npimbalance", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "leaf np imbalance partial C partial P/partial C partial N")
    PUT_ATTR_TXT("units", "gN/gP")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float nlim_m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("nlim_m", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "runmean N limitation factor")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float o2_decomp_depth_unsat(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    DEF_VAR("o2_decomp_depth_unsat", REC_XTYPE, 4, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("long_name", "o2_decomp_depth_unsat")
    PUT_ATTR_TXT("units", "mol/m3/2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    /* float plim_m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    DEF_VAR("plim_m", REC_XTYPE, 3, dimids, REC_ITYPE, 0)
    PUT_ATTR_TXT("long_name", "runmean P limitation factor")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)

    assert(varp - vars + 1 == cfg.nvars + nvars_decomp);

err_out:
    return err;
}

