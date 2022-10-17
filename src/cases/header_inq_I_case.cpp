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
    std::string prefix("");
    int err=0, nprocs;

    /* save number of processes as global attributes */
    if (cfg.strategy == blob) {
        if (cfg.api == adios) {
            GET_GATTR_INT("/__pio__/fillmode", 256)
            prefix = "pio_global/";
        }
        else {
            MPI_Comm_size(cfg.io_comm, &nprocs);
            GET_GATTR_INT("global_nprocs", nprocs)
            GET_GATTR_INT("num_decompositions", decom.num_decomp)
            GET_GATTR_INT("num_subfiles", cfg.num_subfiles)
        }
    }

    /* 28 global attributes: */
    GET_GATTR_TXT("title", "ELM History file information")
    GET_GATTR_TXT("source", "E3SM Land Model")
    GET_GATTR_TXT("source_id", "6eb829238a")
    GET_GATTR_TXT("product", "model-output")
    GET_GATTR_TXT("realm", "land")
    GET_GATTR_TXT("case", "I1850GSWCNPRDCTCBC_hcru_hcru")
    GET_GATTR_TXT("username", "E3SM")
    GET_GATTR_TXT("hostname", "cori-knl")
    GET_GATTR_TXT("git_version", "6eb829238a")
    GET_GATTR_TXT("history", "created on 07/13/21 20:37:15")
    GET_GATTR_TXT("institution_id", "E3SM-Project")
    GET_GATTR_TXT("institution", "LLNL (Lawrence Livermore National Laboratory, Livermore, CA 94550, USA); ANL (Argonne National Laboratory, Argonne, IL 60439, USA); BNL (Brookhaven National Laboratory, Upton, NY 11973, USA); LANL (Los Alamos National Laboratory, Los Alamos, NM 87545, USA); LBNL (Lawrence Berkeley National Laboratory, Berkeley, CA 94720, USA); ORNL (Oak Ridge National Laboratory, Oak Ridge, TN 37831, USA); PNNL (Pacific Northwest National Laboratory, Richland, WA 99352, USA); SNL (Sandia National Laboratories, Albuquerque, NM 87185, USA). Mailing address: LLNL Climate Program, c/o David C. Bader, Principal Investigator, L-103, 7000 East Avenue, Livermore, CA 94550, USA")
    GET_GATTR_TXT("contact", "e3sm-data-support@listserv.llnl.gov")
    GET_GATTR_TXT("Conventions", "CF-1.7")
    GET_GATTR_TXT("comment", "NOTE: None of the variables are weighted by land fraction!")
    GET_GATTR_TXT("Surface_dataset", "surfdata_360x720cru_simyr1850_c180216.nc")
    GET_GATTR_TXT("Initial_conditions_dataset", "arbitrary initialization")
    GET_GATTR_TXT("PFT_physiological_constants_dataset", "clm_params_c180524.nc")

    GET_GATTR_INT("ltype_vegetated_or_bare_soil", 1)
    GET_GATTR_INT("ltype_crop", 2)
    GET_GATTR_INT("ltype_landice", 3)
    GET_GATTR_INT("ltype_landice_multiple_elevation_classes", 4)
    GET_GATTR_INT("ltype_deep_lake", 5)
    GET_GATTR_INT("ltype_wetland", 6)
    GET_GATTR_INT("ltype_urban_tbd", 7)
    GET_GATTR_INT("ltype_urban_hd", 8)
    GET_GATTR_INT("ltype_urban_md", 9)

    if (cfg.hist == h1) { /* h1 only */
        GET_GATTR_TXT("Time_constant_3Dvars_filenamae", "./I1850GSWCNPRDCTCBC_hcru_hcru.elm.h0.0001-01-01-00000.nc")
        GET_GATTR_TXT("Time_constant_3Dvars", "ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE")
    }

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
    int i, err, nprocs, dimids[4], nvars_decomp;
    int dim_time, dim_lon, dim_lat, dim_gridcell, dim_topounit, dim_landunit;
    int dim_column, dim_pft, dim_levgrnd, dim_levurb, dim_levlak, dim_numrad;
    int dim_month, dim_levsno, dim_ltype, dim_nvegwcs, dim_natpft, dim_nblobs;
    int dim_string_length, dim_levdcmp, dim_levtrc, dim_hist_interval;
    int dim_nelems[MAX_NUM_DECOMP], dim_max_nreqs[MAX_NUM_DECOMP];
    int g_dimids[MAX_NUM_DECOMP][4];
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
    INQ_VAR("levgrnd", NC_FLOAT, 1, &dim_levgrnd, REC_ITYPE, -1)
    GET_ATTR_TXT("long_name", "coordinate soil levels")
    GET_ATTR_TXT("units", "m")

    /* float levlak(levlak) */
    INQ_VAR("levlak", NC_FLOAT, 1, &dim_levlak, REC_ITYPE, -1)
    GET_ATTR_TXT("long_name", "coordinate lake levels")
    GET_ATTR_TXT("units", "m")

    /* float levdcmp(levdcmp) */
    INQ_VAR("levdcmp", NC_FLOAT, 1, &dim_levdcmp, REC_ITYPE, -1)
    GET_ATTR_TXT("long_name", "coordinate soil levels")
    GET_ATTR_TXT("units", "m")

    /* float time(time) */
    INQ_VAR("time", NC_FLOAT, 1, &dim_time, REC_ITYPE, -1)
    GET_ATTR_TXT("long_name", "time")
    GET_ATTR_TXT("units", "days since 0001-01-01 00:00:00")
    GET_ATTR_TXT("calendar", "noleap")
    GET_ATTR_TXT("bounds", "time_bounds")

    /* int mcdate(time) */
    INQ_VAR("mcdate", NC_INT, 1, &dim_time, MPI_INT, -1)
    GET_ATTR_TXT("long_name", "current date (YYYYMMDD)")

    /* int mcsec(time) */
    INQ_VAR("mcsec", NC_INT, 1, &dim_time, MPI_INT, -1)
    GET_ATTR_TXT("long_name", "current seconds of current date")
    GET_ATTR_TXT("units", "s")

    /* int mdcur(time) */
    INQ_VAR("mdcur", NC_INT, 1, &dim_time, MPI_INT, -1)
    GET_ATTR_TXT("long_name", "current day (from base day)")

    /* int mscur(time) */
    INQ_VAR("mscur", NC_INT, 1, &dim_time, MPI_INT, -1)
    GET_ATTR_TXT("long_name", "current seconds of current day")

    /* int nstep(time) */
    INQ_VAR("nstep", NC_INT, 1, &dim_time, MPI_INT, -1)
    GET_ATTR_TXT("long_name", "time step")

    /* double time_bounds(time, hist_interval) */
    dimids[0] = dim_time;
    dimids[1] = dim_hist_interval;
    INQ_VAR("time_bounds", NC_DOUBLE, 2, dimids, MPI_DOUBLE, -1)
    GET_ATTR_TXT("long_name", "history time interval endpoints")

    /* char date_written(time, string_length) */
    dimids[0] = dim_time;
    dimids[1] = dim_string_length;
    INQ_VAR("date_written", NC_CHAR, 2, dimids, MPI_CHAR, -1)

    /* char time_written(time, string_length) */
    dimids[0] = dim_time;
    dimids[1] = dim_string_length;
    INQ_VAR("time_written", NC_CHAR, 2, dimids, MPI_CHAR, -1)

    /* float lon(lon) */
    INQ_VAR("lon", NC_FLOAT, 1, &dim_lon, REC_ITYPE, -1)
    GET_ATTR_TXT("long_name", "coordinate longitude")
    GET_ATTR_TXT("units", "degrees_east")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float lat(lat) */
    INQ_VAR("lat", NC_FLOAT, 1, &dim_lat, REC_ITYPE, -1)
    GET_ATTR_TXT("long_name", "coordinate latitude")
    GET_ATTR_TXT("units", "degrees_north")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float area(lat, lon) */
    dimids[0] = dim_lat;
    dimids[1] = dim_lon;
    INQ_VAR("area", NC_FLOAT, 2, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "grid cell areas")
    GET_ATTR_TXT("units", "km^2")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float topo(lat, lon) */
    dimids[0] = dim_lat;
    dimids[1] = dim_lon;
    INQ_VAR("topo", NC_FLOAT, 2, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "grid cell topography")
    GET_ATTR_TXT("units", "m")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float landfrac(lat, lon) */
    dimids[0] = dim_lat;
    dimids[1] = dim_lon;
    INQ_VAR("landfrac", NC_FLOAT, 2, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "land fraction")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* int landmask(lat, lon) */
    dimids[0] = dim_lat;
    dimids[1] = dim_lon;
    INQ_VAR("landmask", NC_INT, 2, dimids, MPI_INT, 0)
    GET_ATTR_TXT("long_name", "land/ocean mask (0.=ocean and 1.=land)")
    GET_ATTR_FILL(-9999)
    GET_ATTR_INT1("missing_value", -9999)

    /* int pftmask(lat, lon) */
    dimids[0] = dim_lat;
    dimids[1] = dim_lon;
    INQ_VAR("pftmask", NC_INT, 2, dimids, MPI_INT, 0)
    GET_ATTR_TXT("long_name", "pft real/fake mask (0.=fake and 1.=real)")
    GET_ATTR_FILL(-9999)
    GET_ATTR_INT1("missing_value", -9999)

    if (cfg.hist == h0) {  /* h0 only */
        /* float ZSOI(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        INQ_VAR("ZSOI", NC_FLOAT, 3, dimids, REC_ITYPE, 1)
        GET_ATTR_TXT("long_name", "soil depth")
        GET_ATTR_TXT("units", "m")
        GET_ATTR_FILL(fillv)
        GET_ATTR_FLT1("missing_value", missv)

        /* float DZSOI(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        INQ_VAR("DZSOI", NC_FLOAT, 3, dimids, REC_ITYPE, 1)
        GET_ATTR_TXT("long_name", "soil thickness")
        GET_ATTR_TXT("units", "m")
        GET_ATTR_FILL(fillv)
        GET_ATTR_FLT1("missing_value", missv)

        /* float WATSAT(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        INQ_VAR("WATSAT", NC_FLOAT, 3, dimids, REC_ITYPE, 1)
        GET_ATTR_TXT("long_name", "saturated soil water content (porosity)")
        GET_ATTR_TXT("units", "mm3/mm3")
        GET_ATTR_FILL(fillv)
        GET_ATTR_FLT1("missing_value", missv)

        /* float SUCSAT(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        INQ_VAR("SUCSAT", NC_FLOAT, 3, dimids, REC_ITYPE, 1)
        GET_ATTR_TXT("long_name", "saturated soil matric potential")
        GET_ATTR_TXT("units", "mm")
        GET_ATTR_FILL(fillv)
        GET_ATTR_FLT1("missing_value", missv)

        /* float BSW(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        INQ_VAR("BSW", NC_FLOAT, 3, dimids, REC_ITYPE, 1)
        GET_ATTR_TXT("long_name", "slope of soil water retention curve")
        GET_ATTR_TXT("units", "1")
        GET_ATTR_FILL(fillv)
        GET_ATTR_FLT1("missing_value", missv)

        /* float HKSAT(levgrnd, lat, lon) */
        dimids[0] = dim_levgrnd;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        INQ_VAR("HKSAT", NC_FLOAT, 3, dimids, REC_ITYPE, 1)
        GET_ATTR_TXT("long_name", "saturated hydraulic conductivity")
        GET_ATTR_TXT("units", "1")
        GET_ATTR_FILL(fillv)
        GET_ATTR_FLT1("missing_value", missv)

        /* float ZLAKE(levlak, lat, lon) */
        dimids[0] = dim_levlak;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        INQ_VAR("ZLAKE", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
        GET_ATTR_TXT("long_name", "lake layer node depth")
        GET_ATTR_TXT("units", "m")
        GET_ATTR_FILL(fillv)
        GET_ATTR_FLT1("missing_value", missv)

        /* float DZLAKE(levlak, lat, lon) */
        dimids[0] = dim_levlak;
        dimids[1] = dim_lat;
        dimids[2] = dim_lon;
        INQ_VAR("DZLAKE", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
        GET_ATTR_TXT("long_name", "lake layer thickness")
        GET_ATTR_TXT("units", "m")
        GET_ATTR_FILL(fillv)
        GET_ATTR_FLT1("missing_value", missv)
    }

    /* float ACTUAL_IMMOB(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ACTUAL_IMMOB", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "actual N immobilization")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float ACTUAL_IMMOB_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ACTUAL_IMMOB_P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "actual P immobilization")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float ADSORBTION_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ADSORBTION_P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "adsorb P flux")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float AGNPP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("AGNPP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "aboveground NPP")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float AGWDNPP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("AGWDNPP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "aboveground wood NPP")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float ALT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ALT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "current active layer thickness")
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float ALTMAX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ALTMAX", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "maximum annual active layer thickness")
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float ALTMAX_LASTYEAR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ALTMAX_LASTYEAR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "maximum prior year active layer thickness")
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float AR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("AR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "autotrophic respiration (MR + GR)")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float AVAILC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("AVAILC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "C flux available for allocation")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float AVAIL_RETRANSP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("AVAIL_RETRANSP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "P flux available from retranslocation pool")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float BAF_CROP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("BAF_CROP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "fractional area burned for crop")
    GET_ATTR_TXT("units", "proportion/sec")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float BAF_PEATF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("BAF_PEATF", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "fractional area burned in peatland")
    GET_ATTR_TXT("units", "proportion/sec")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float BCDEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("BCDEP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total BC deposition (dry+wet) from atmosphere")
    GET_ATTR_TXT("units", "kg/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float BGNPP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("BGNPP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "belowground NPP")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float BIOCHEM_PMIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("BIOCHEM_PMIN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "biochemical rate of P mineralization")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float BIOCHEM_PMIN_TO_PLANT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("BIOCHEM_PMIN_TO_PLANT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "plant uptake of biochemical P mineralization")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float BTRAN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("BTRAN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "transpiration beta factor")
    GET_ATTR_TXT("units", "1")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float BUILDHEAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("BUILDHEAT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "heat flux from urban building interior to walls and roof")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CH4PROD(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CH4PROD", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Gridcell total production of CH4")
    GET_ATTR_TXT("units", "gC/m2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CH4_SURF_AERE_SAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CH4_SURF_AERE_SAT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "aerenchyma surface CH4 flux for inundated area; (+ to atm)")
    GET_ATTR_TXT("units", "mol/m2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CH4_SURF_AERE_UNSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CH4_SURF_AERE_UNSAT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "aerenchyma surface CH4 flux for non-inundated area; (+ to atm)")
    GET_ATTR_TXT("units", "mol/m2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CH4_SURF_DIFF_SAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CH4_SURF_DIFF_SAT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "diffusive surface CH4 flux for inundated / lake area; (+ to atm)")
    GET_ATTR_TXT("units", "mol/m2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CH4_SURF_DIFF_UNSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CH4_SURF_DIFF_UNSAT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "diffusive surface CH4 flux for non-inundated area; (+ to atm)")
    GET_ATTR_TXT("units", "mol/m2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CH4_SURF_EBUL_SAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CH4_SURF_EBUL_SAT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "ebullition surface CH4 flux for inundated / lake area; (+ to atm)")
    GET_ATTR_TXT("units", "mol/m2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CH4_SURF_EBUL_UNSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CH4_SURF_EBUL_UNSAT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "ebullition surface CH4 flux for non-inundated area; (+ to atm)")
    GET_ATTR_TXT("units", "mol/m2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float COL_PTRUNC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("COL_PTRUNC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "column-level sink for P truncation")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CONC_CH4_SAT(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("CONC_CH4_SAT", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "CH4 soil Concentration for inundated / lake area")
    GET_ATTR_TXT("units", "mol/m3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CONC_CH4_UNSAT(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("CONC_CH4_UNSAT", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "CH4 soil Concentration for non-inundated area")
    GET_ATTR_TXT("units", "mol/m3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CONC_O2_SAT(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("CONC_O2_SAT", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "O2 soil Concentration for inundated / lake area")
    GET_ATTR_TXT("units", "mol/m3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CONC_O2_UNSAT(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("CONC_O2_UNSAT", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "O2 soil Concentration for non-inundated area")
    GET_ATTR_TXT("units", "mol/m3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CPOOL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "temporary photosynthate C pool")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CWDC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "CWD C")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CWDC_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDC_HR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "coarse woody debris C heterotrophic respiration")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CWDC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDC_LOSS", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "coarse woody debris C loss")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CWDC_TO_LITR2C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDC_TO_LITR2C", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of coarse woody debris C to litter 2 C")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CWDC_TO_LITR3C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDC_TO_LITR3C", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of coarse woody debris C to litter 3 C")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CWDC_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("CWDC_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "CWD C (vertically resolved)")
    GET_ATTR_TXT("units", "gC/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CWDN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "CWD N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CWDN_TO_LITR2N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDN_TO_LITR2N", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of coarse woody debris N to litter 2 N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CWDN_TO_LITR3N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDN_TO_LITR3N", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of coarse woody debris N to litter 3 N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CWDN_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("CWDN_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "CWD N (vertically resolved)")
    GET_ATTR_TXT("units", "gN/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CWDP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "CWD P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CWDP_TO_LITR2P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDP_TO_LITR2P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of coarse woody debris P to litter 2 N")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CWDP_TO_LITR3P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("CWDP_TO_LITR3P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of coarse woody debris P to litter 3 N")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float CWDP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("CWDP_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "CWD P (vertically resolved)")
    GET_ATTR_TXT("units", "gP/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DEADCROOTC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DEADCROOTC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "dead coarse root C")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DEADCROOTN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DEADCROOTN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "dead coarse root N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DEADCROOTP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DEADCROOTP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "dead coarse root P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DEADSTEMC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DEADSTEMC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "dead stem C")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DEADSTEMN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DEADSTEMN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "dead stem N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DEADSTEMP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DEADSTEMP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "dead stem P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float INQICIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("INQICIT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "runoff supply deficit")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DENIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DENIT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total rate of denitrification")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DESORPTION_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DESORPTION_P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "desorp P flux")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DISPVEGC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DISPVEGC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "displayed veg carbon, excluding storage and cpool")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DISPVEGN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DISPVEGN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "displayed vegetation nitrogen")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DISPVEGP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DISPVEGP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "displayed vegetation phosphorus")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DSTDEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DSTDEP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total dust deposition (dry+wet) from atmosphere")
    GET_ATTR_TXT("units", "kg/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DSTFLXT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DSTFLXT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total surface dust emission")
    GET_ATTR_TXT("units", "kg/m2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DWB(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DWB", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "net change in total water mass")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DWT_CONV_CFLUX_DRIBBLED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DWT_CONV_CFLUX_DRIBBLED", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "conversion C flux (immediate loss to atm), dribbled throughout the year")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DWT_CONV_CFLUX_GRC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DWT_CONV_CFLUX_GRC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DWT_CONV_NFLUX_GRC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DWT_CONV_NFLUX_GRC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DWT_CONV_PFLUX_GRC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DWT_CONV_PFLUX_GRC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DWT_SLASH_CFLUX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DWT_SLASH_CFLUX", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "slash C flux to litter and CWD due to land use")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DWT_SLASH_NFLUX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DWT_SLASH_NFLUX", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "slash N flux to litter and CWD due to land use")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float DWT_SLASH_PFLUX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("DWT_SLASH_PFLUX", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "slash P flux to litter and CWD due to land use")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float EFLX_DYNBAL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("EFLX_DYNBAL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "dynamic land cover change conversion energy flux")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float EFLX_GRND_LAKE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("EFLX_GRND_LAKE", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "net heat flux into lake/snow surface, excluding light transmission")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float EFLX_LH_TOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("EFLX_LH_TOT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total latent heat flux [+ to atm]")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float EFLX_LH_TOT_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("EFLX_LH_TOT_R", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Rural total evaporation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float EFLX_LH_TOT_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("EFLX_LH_TOT_U", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Urban total evaporation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float ELAI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ELAI", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "exposed one-sided leaf area index")
    GET_ATTR_TXT("units", "m^2/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float ER(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ER", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total ecosystem respiration, autotrophic + heterotrophic")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float ERRH2O(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ERRH2O", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total water conservation error")
    GET_ATTR_TXT("units", "mm")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float ERRH2OSNO(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ERRH2OSNO", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "imbalance in snow depth (liquid water)")
    GET_ATTR_TXT("units", "mm")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float ERRSEB(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ERRSEB", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "surface energy conservation error")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float ERRSOI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ERRSOI", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil/lake energy conservation error")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float ERRSOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ERRSOL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "solar radiation conservation error")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float ESAI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ESAI", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "exposed one-sided stem area index")
    GET_ATTR_TXT("units", "m^2/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FAREA_BURNED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FAREA_BURNED", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "timestep fractional area burned")
    GET_ATTR_TXT("units", "proportion")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FCEV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FCEV", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "canopy evaporation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FCH4(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FCH4", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Gridcell surface CH4 flux to atmosphere (+ to atm)")
    GET_ATTR_TXT("units", "kgC/m2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FCH4TOCO2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FCH4TOCO2", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Gridcell oxidation of CH4 to CO2")
    GET_ATTR_TXT("units", "gC/m2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FCH4_DFSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FCH4_DFSAT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "CH4 additional flux due to changing fsat, vegetated landunits only")
    GET_ATTR_TXT("units", "kgC/m2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FCOV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FCOV", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "fractional impermeable area")
    GET_ATTR_TXT("units", "unitless")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FCTR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FCTR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "canopy transpiration")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FGEV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FGEV", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "ground evaporation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FGR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FGR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "heat flux into soil/snow including snow melt and lake / snow light transmission")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FGR12(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FGR12", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "heat flux between soil layers 1 and 2")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FGR_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FGR_R", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Rural heat flux into soil/snow including snow melt and snow light transmission")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FGR_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FGR_U", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Urban heat flux into soil/snow including snow melt")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FH2OSFC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FH2OSFC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "fraction of ground covered by surface water")
    GET_ATTR_TXT("units", "1")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FINUNDATED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FINUNDATED", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "fractional inundated area of vegetated columns")
    GET_ATTR_TXT("units", "1")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FINUNDATED_LAG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FINUNDATED_LAG", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "time-lagged inundated fraction of vegetated columns")
    GET_ATTR_TXT("units", "1")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FIRA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FIRA", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "net infrared (longwave) radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FIRA_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FIRA_R", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Rural net infrared (longwave) radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FIRA_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FIRA_U", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Urban net infrared (longwave) radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FIRE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FIRE", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "emitted infrared (longwave) radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FIRE_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FIRE_R", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Rural emitted infrared (longwave) radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FIRE_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FIRE_U", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Urban emitted infrared (longwave) radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FLDS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FLDS", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "atmospheric longwave radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FPG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FPG", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "fraction of potential gpp due to N limitation")
    GET_ATTR_TXT("units", "proportion")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FPG_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FPG_P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "fraction of potential gpp due to P limitation")
    GET_ATTR_TXT("units", "proportion")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FPI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FPI", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "fraction of potential immobilization of nitrogen")
    GET_ATTR_TXT("units", "proportion")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FPI_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FPI_P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "fraction of potential immobilization of phosphorus")
    GET_ATTR_TXT("units", "proportion")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FPI_P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("FPI_P_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "fraction of potential immobilization of phosphorus")
    GET_ATTR_TXT("units", "proportion")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FPI_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("FPI_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "fraction of potential immobilization of nitrogen")
    GET_ATTR_TXT("units", "proportion")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FPSN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FPSN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "photosynthesis")
    GET_ATTR_TXT("units", "umol/m2s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FPSN_WC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FPSN_WC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Rubisco-limited photosynthesis")
    GET_ATTR_TXT("units", "umol/m2s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FPSN_WJ(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FPSN_WJ", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "RuBP-limited photosynthesis")
    GET_ATTR_TXT("units", "umol/m2s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FPSN_WP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FPSN_WP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Product-limited photosynthesis")
    GET_ATTR_TXT("units", "umol/m2s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FROOTC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FROOTC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "fine root C")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FROOTC_ALLOC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FROOTC_ALLOC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "fine root C allocation")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FROOTC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FROOTC_LOSS", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "fine root C loss")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FROOTN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FROOTN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "fine root N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FROOTP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FROOTP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "fine root P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FROST_TABLE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FROST_TABLE", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "frost table depth (vegetated landunits only)")
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSA", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "absorbed solar radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSAT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "fractional area with water table at surface")
    GET_ATTR_TXT("units", "unitless")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSA_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSA_R", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Rural absorbed solar radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSA_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSA_U", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Urban absorbed solar radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSDS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSDS", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "atmospheric incident solar radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSDSND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSDSND", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "direct nir incident solar radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSDSNDLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSDSNDLN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "direct nir incident solar radiation at local noon")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSDSNI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSDSNI", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "diffuse nir incident solar radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSDSVD(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSDSVD", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "direct vis incident solar radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSDSVDLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSDSVDLN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "direct vis incident solar radiation at local noon")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSDSVI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSDSVI", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "diffuse vis incident solar radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSDSVILN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSDSVILN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "diffuse vis incident solar radiation at local noon")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSH(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSH", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "sensible heat")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSH_G(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSH_G", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "sensible heat from ground")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSH_NODYNLNDUSE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSH_NODYNLNDUSE", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "sensible heat not including correction for land use change")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSH_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSH_R", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Rural sensible heat")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSH_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSH_U", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Urban sensible heat")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSH_V(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSH_V", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "sensible heat from veg")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSM(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSM", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "snow melt heat flux")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSM_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSM_R", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Rural snow melt heat flux")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSM_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSM_U", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Urban snow melt heat flux")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSNO(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSNO", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "fraction of ground covered by snow")
    GET_ATTR_TXT("units", "1")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSNO_EFF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSNO_EFF", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "effective fraction of ground covered by snow")
    GET_ATTR_TXT("units", "1")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "reflected solar radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSRND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSRND", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "direct nir reflected solar radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSRNDLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSRNDLN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "direct nir reflected solar radiation at local noon")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSRNI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSRNI", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "diffuse nir reflected solar radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSRVD(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSRVD", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "direct vis reflected solar radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSRVDLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSRVDLN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "direct vis reflected solar radiation at local noon")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float FSRVI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("FSRVI", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "diffuse vis reflected solar radiation")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float F_CO2_SOIL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("F_CO2_SOIL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total soil-atm. CO2 exchange")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float F_CO2_SOIL_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("F_CO2_SOIL_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "total vertically resolved soil-atm. CO2 exchange")
    GET_ATTR_TXT("units", "gC/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float F_DENIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("F_DENIT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "denitrification flux")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float F_DENIT_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("F_DENIT_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "denitrification flux")
    GET_ATTR_TXT("units", "gN/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float F_N2O_DENIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("F_N2O_DENIT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "denitrification N2O flux")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float F_N2O_NIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("F_N2O_NIT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "nitrification N2O flux")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float F_NIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("F_NIT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "nitrification flux")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float F_NIT_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("F_NIT_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "nitrification flux")
    GET_ATTR_TXT("units", "gN/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float GC_HEAT1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("GC_HEAT1", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "initial gridcell total heat content")
    GET_ATTR_TXT("units", "J/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float GC_ICE1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("GC_ICE1", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "initial gridcell total ice content")
    GET_ATTR_TXT("units", "mm")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float GC_LIQ1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("GC_LIQ1", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "initial gridcell total liq content")
    GET_ATTR_TXT("units", "mm")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float GPP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("GPP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "gross primary production")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float GR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("GR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total growth respiration")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float GROSS_NMIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("GROSS_NMIN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "gross rate of N mineralization")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float GROSS_PMIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("GROSS_PMIN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "gross rate of P mineralization")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float H2OCAN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("H2OCAN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "intercepted water")
    GET_ATTR_TXT("units", "mm")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float H2OSFC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("H2OSFC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "surface water depth")
    GET_ATTR_TXT("units", "mm")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float H2OSNO(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("H2OSNO", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "snow depth (liquid water)")
    GET_ATTR_TXT("units", "mm")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float H2OSNO_TOP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("H2OSNO_TOP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mass of snow in top snow layer")
    GET_ATTR_TXT("units", "kg/m2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float H2OSOI(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("H2OSOI", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "volumetric soil water (vegetated landunits only)")
    GET_ATTR_TXT("units", "mm3/mm3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float HC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("HC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "heat content of soil/snow/lake")
    GET_ATTR_TXT("units", "MJ/m2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float HCSOI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("HCSOI", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil heat content")
    GET_ATTR_TXT("units", "MJ/m2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float HEAT_FROM_AC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("HEAT_FROM_AC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "sensible heat flux put into canyon due to heat removed from air conditioning")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("HR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total heterotrophic respiration")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float HR_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("HR_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "total vertically resolved heterotrophic respiration")
    GET_ATTR_TXT("units", "gC/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float HTOP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("HTOP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "canopy top")
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float INT_SNOW(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("INT_SNOW", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "accumulated swe (vegetated landunits only)")
    GET_ATTR_TXT("units", "mm")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LABILEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LABILEP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil Labile P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LABILEP_TO_SECONDP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LABILEP_TO_SECONDP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "LABILE P TO SECONDARY MINERAL P")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LABILEP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LABILEP_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil labile P (vert. res.)")
    GET_ATTR_TXT("units", "gp/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LAISHA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LAISHA", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "shaded projected leaf area index")
    GET_ATTR_TXT("units", "1")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LAISUN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LAISUN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "sunlit projected leaf area index")
    GET_ATTR_TXT("units", "1")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LAKEICEFRAC(time, levlak, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levlak;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LAKEICEFRAC", NC_FLOAT, 4, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("long_name", "lake layer ice mass fraction")
    GET_ATTR_TXT("units", "unitless")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LAKEICETHICK(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LAKEICETHICK", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "thickness of lake ice (including physical expansion on freezing)")
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LAND_UPTAKE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LAND_UPTAKE", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "NEE minus LAND_USE_FLUX, negative for update")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LAND_USE_FLUX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LAND_USE_FLUX", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total C emitted from land cover conversion and wood product pools")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LEAFC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LEAFC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "leaf C")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LEAFC_ALLOC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LEAFC_ALLOC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "leaf C allocation")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LEAFC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LEAFC_LOSS", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "leaf C loss")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LEAFC_TO_LITTER(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LEAFC_TO_LITTER", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "leaf C litterfall")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LEAFN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LEAFN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "leaf N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LEAFP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LEAFP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "leaf P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LEAF_MR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LEAF_MR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "leaf maintenance respiration")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LFC2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LFC2", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "conversion area fraction of BET and BDT that burned")
    GET_ATTR_TXT("units", "per sec")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITFALL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITFALL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "litterfall (leaves and fine roots)")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITHR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITHR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "litter heterotrophic respiration")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR1C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR1C", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "LITR1 C")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR1C_TO_SOIL1C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR1C_TO_SOIL1C", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of litter 1 C to soil 1 C")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR1C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR1C_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "LITR1 C (vertically resolved)")
    GET_ATTR_TXT("units", "gC/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR1N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR1N", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "LITR1 N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR1N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR1N_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "litter 1 N tendency due to vertical transport")
    GET_ATTR_TXT("units", "gN/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR1N_TO_SOIL1N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR1N_TO_SOIL1N", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of litter 1 N to soil 1 N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR1N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR1N_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "LITR1 N (vertically resolved)")
    GET_ATTR_TXT("units", "gN/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR1P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR1P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "LITR1 P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR1P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR1P_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "litter 1 P tendency due to vertical transport")
    GET_ATTR_TXT("units", "gP/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR1P_TO_SOIL1P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR1P_TO_SOIL1P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of litter 1 P to soil 1 N")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR1P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR1P_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "LITR1 P (vertically resolved)")
    GET_ATTR_TXT("units", "gP/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR1_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR1_HR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Het. Resp. from litter 1")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR2C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR2C", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "LITR2 C")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR2C_TO_SOIL2C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR2C_TO_SOIL2C", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of litter 2 C to soil 2 C")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR2C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR2C_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "LITR2 C (vertically resolved)")
    GET_ATTR_TXT("units", "gC/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR2N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR2N", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "LITR2 N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR2N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR2N_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "litter 2 N tendency due to vertical transport")
    GET_ATTR_TXT("units", "gN/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR2N_TO_SOIL2N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR2N_TO_SOIL2N", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of litter 2 N to soil 2 N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR2N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR2N_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "LITR2 N (vertically resolved)")
    GET_ATTR_TXT("units", "gN/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR2P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR2P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "LITR2 P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR2P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR2P_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "litter 2 P tendency due to vertical transport")
    GET_ATTR_TXT("units", "gP/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR2P_TO_SOIL2P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR2P_TO_SOIL2P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of litter 2 P to soil 2 N")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR2P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR2P_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "LITR2 P (vertically resolved)")
    GET_ATTR_TXT("units", "gP/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR2_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR2_HR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Het. Resp. from litter 2")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR3C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR3C", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "LITR3 C")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR3C_TO_SOIL3C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR3C_TO_SOIL3C", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of litter 3 C to soil 3 C")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR3C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR3C_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "LITR3 C (vertically resolved)")
    GET_ATTR_TXT("units", "gC/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR3N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR3N", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "LITR3 N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR3N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR3N_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "litter 3 N tendency due to vertical transport")
    GET_ATTR_TXT("units", "gN/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR3N_TO_SOIL3N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR3N_TO_SOIL3N", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of litter 3 N to soil 3 N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR3N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR3N_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "LITR3 N (vertically resolved)")
    GET_ATTR_TXT("units", "gN/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR3P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR3P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "LITR3 P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR3P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR3P_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "litter 3 P tendency due to vertical transport")
    GET_ATTR_TXT("units", "gP/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR3P_TO_SOIL3P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR3P_TO_SOIL3P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of litter 3 P to soil 3 N")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR3P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("LITR3P_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "LITR3 P (vertically resolved)")
    GET_ATTR_TXT("units", "gP/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITR3_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITR3_HR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Het. Resp. from litter 3")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITTERC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITTERC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "litter C")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITTERC_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITTERC_HR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "litter C heterotrophic respiration")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LITTERC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LITTERC_LOSS", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "litter C loss")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LIVECROOTC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LIVECROOTC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "live coarse root C")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LIVECROOTN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LIVECROOTN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "live coarse root N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LIVECROOTP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LIVECROOTP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "live coarse root P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LIVESTEMC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LIVESTEMC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "live stem C")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LIVESTEMN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LIVESTEMN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "live stem N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float LIVESTEMP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("LIVESTEMP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "live stem P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float MR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("MR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "maintenance respiration")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float M_LITR1C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("M_LITR1C_TO_LEACHING", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "litter 1 C leaching loss")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float M_LITR2C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("M_LITR2C_TO_LEACHING", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "litter 2 C leaching loss")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float M_LITR3C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("M_LITR3C_TO_LEACHING", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "litter 3 C leaching loss")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float M_SOIL1C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("M_SOIL1C_TO_LEACHING", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil 1 C leaching loss")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float M_SOIL2C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("M_SOIL2C_TO_LEACHING", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil 2 C leaching loss")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float M_SOIL3C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("M_SOIL3C_TO_LEACHING", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil 3 C leaching loss")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float M_SOIL4C_TO_LEACHING(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("M_SOIL4C_TO_LEACHING", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil 4 C leaching loss")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float NBP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NBP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "net biome production, includes fire, landuse, and harvest flux, positive for sink")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float NDEPLOY(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NDEPLOY", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total N deployed in new growth")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float NDEP_TO_SMINN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NDEP_TO_SMINN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "atmospheric N deposition to soil mineral N")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float NEE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NEE", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "net ecosystem exchange of carbon, includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float NEM(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NEM", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Gridcell net adjustment to NEE passed to atm. for methane production")
    GET_ATTR_TXT("units", "gC/m2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float NEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NEP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "net ecosystem production, excludes fire, landuse, and harvest flux, positive for sink")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float NET_NMIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NET_NMIN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "net rate of N mineralization")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float NET_PMIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NET_PMIN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "net rate of P mineralization")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float NFIRE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NFIRE", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "fire counts valid only in Reg.C")
    GET_ATTR_TXT("units", "counts/km2/sec")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float NFIX_TO_SMINN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NFIX_TO_SMINN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "symbiotic/asymbiotic N fixation to soil mineral N")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float NPP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("NPP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "net primary production")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float OCCLP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("OCCLP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil occluded P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float OCCLP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("OCCLP_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil occluded P (vert. res.)")
    GET_ATTR_TXT("units", "gp/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float OCDEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("OCDEP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total OC deposition (dry+wet) from atmosphere")
    GET_ATTR_TXT("units", "kg/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float O_SCALAR(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("O_SCALAR", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "fraction by which decomposition is reduced due to anoxia")
    GET_ATTR_TXT("units", "1")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PARVEGLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PARVEGLN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "absorbed par by vegetation at local noon")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PBOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PBOT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "atmospheric pressure")
    GET_ATTR_TXT("units", "Pa")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PCH4(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PCH4", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "atmospheric partial pressure of CH4")
    GET_ATTR_TXT("units", "Pa")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PCO2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PCO2", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "atmospheric partial pressure of CO2")
    GET_ATTR_TXT("units", "Pa")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PCT_LANDUNIT(time, ltype, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_ltype;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("PCT_LANDUNIT", NC_FLOAT, 4, dimids, REC_ITYPE, 3)
    GET_ATTR_TXT("long_name", "% of each landunit on grid cell")
    GET_ATTR_TXT("units", "%")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PCT_NAT_PFT(time, natpft, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_natpft;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("PCT_NAT_PFT", NC_FLOAT, 4, dimids, REC_ITYPE, 4)
    GET_ATTR_TXT("long_name", "% of each PFT on the natural vegetation (i.e., soil) landunit")
    GET_ATTR_TXT("units", "%")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PDEPLOY(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PDEPLOY", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total P deployed in new growth")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PDEP_TO_SMINP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PDEP_TO_SMINP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "atmospheric P deposition to soil mineral P")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PFT_FIRE_CLOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PFT_FIRE_CLOSS", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total patch-level fire C loss for non-peat fires outside land-type converted region")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PFT_FIRE_NLOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PFT_FIRE_NLOSS", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total pft-level fire N loss")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PLANT_CALLOC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PLANT_CALLOC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total allocated C flux")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PLANT_NDEMAND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PLANT_NDEMAND", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "N flux required to support initial GPP")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PLANT_NDEMAND_COL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PLANT_NDEMAND_COL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "N flux required to support initial GPP")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PLANT_PALLOC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PLANT_PALLOC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total allocated P flux")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PLANT_PDEMAND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PLANT_PDEMAND", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "P flux required to support initial GPP")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PLANT_PDEMAND_COL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PLANT_PDEMAND_COL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "P flux required to support initial GPP")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float POTENTIAL_IMMOB(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("POTENTIAL_IMMOB", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "potential N immobilization")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float POTENTIAL_IMMOB_P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("POTENTIAL_IMMOB_P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "potential P immobilization")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float POT_F_DENIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("POT_F_DENIT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "potential denitrification flux")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float POT_F_NIT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("POT_F_NIT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "potential nitrification flux")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PRIMP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PRIMP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil primary P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PRIMP_TO_LABILEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PRIMP_TO_LABILEP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "PRIMARY MINERAL P TO LABILE P")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PRIMP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("PRIMP_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil primary P (vert. res.)")
    GET_ATTR_TXT("units", "gp/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PROD1P_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PROD1P_LOSS", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "loss from 1-yr crop product pool")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PSNSHA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PSNSHA", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "shaded leaf photosynthesis")
    GET_ATTR_TXT("units", "umolCO2/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PSNSHADE_TO_CPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PSNSHADE_TO_CPOOL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "C fixation from shaded canopy")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PSNSUN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PSNSUN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "sunlit leaf photosynthesis")
    GET_ATTR_TXT("units", "umolCO2/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float PSNSUN_TO_CPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("PSNSUN_TO_CPOOL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "C fixation from sunlit canopy")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float Q2M(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("Q2M", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "2m specific humidity")
    GET_ATTR_TXT("units", "kg/kg")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QBOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QBOT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "atmospheric specific humidity")
    GET_ATTR_TXT("units", "kg/kg")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QCHARGE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QCHARGE", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "aquifer recharge rate (vegetated landunits only)")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QDRAI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QDRAI", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "sub-surface drainage")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QDRAI_PERCH(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QDRAI_PERCH", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "perched wt drainage")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QDRAI_XS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QDRAI_XS", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "saturation excess drainage")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QDRIP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QDRIP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "throughfall")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QFLOOD(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QFLOOD", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "runoff from river flooding")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QFLX_ICE_DYNBAL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QFLX_ICE_DYNBAL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "ice dynamic land cover change conversion runoff flux")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QFLX_LIQ_DYNBAL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QFLX_LIQ_DYNBAL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "liq dynamic land cover change conversion runoff flux")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QH2OSFC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QH2OSFC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "surface water runoff")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QINFL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QINFL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "infiltration")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QINTR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QINTR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "interception")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QIRRIG_GRND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QIRRIG_GRND", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Groundwater irrigation")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QIRRIG_ORIG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QIRRIG_ORIG", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Original total irrigation water demand (surface + ground)")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QIRRIG_REAL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QIRRIG_REAL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "actual water added through irrigation (surface + ground)")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QIRRIG_SURF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QIRRIG_SURF", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Surface water irrigation")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QIRRIG_WM(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QIRRIG_WM", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Surface water irrigation demand sent to MOSART/WM")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QOVER(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QOVER", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "surface runoff")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QOVER_LAG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QOVER_LAG", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "time-lagged surface runoff for soil columns")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QRGWL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QRGWL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "surface runoff at glaciers (liquid only), wetlands, lakes")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QRUNOFF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QRUNOFF", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total liquid runoff (does not include QSNWCPICE)")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QRUNOFF_NODYNLNDUSE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QRUNOFF_NODYNLNDUSE", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total liquid runoff (does not include QSNWCPICE) not including correction for land use change")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QRUNOFF_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QRUNOFF_R", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Rural total runoff")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QRUNOFF_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QRUNOFF_U", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Urban total runoff")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QSNOMELT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QSNOMELT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "snow melt")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QSNWCPICE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QSNWCPICE", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "excess snowfall due to snow capping")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QSNWCPICE_NODYNLNDUSE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QSNWCPICE_NODYNLNDUSE", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "excess snowfall due to snow capping not including correction for land use change")
    GET_ATTR_TXT("units", "mm H2O/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QSOIL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QSOIL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Ground evaporation (soil/snow evaporation + soil/snow sublimation - dew)")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QVEGE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QVEGE", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "canopy evaporation")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float QVEGT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("QVEGT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "canopy transpiration")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float RAIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RAIN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "atmospheric rain")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float RETRANSN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RETRANSN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "plant pool of retranslocated N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float RETRANSN_TO_NPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RETRANSN_TO_NPOOL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "deployment of retranslocated N")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float RETRANSP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RETRANSP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "plant pool of retranslocated P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float RETRANSP_TO_PPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RETRANSP_TO_PPOOL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "deployment of retranslocated P")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float RH2M(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RH2M", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "2m relative humidity")
    GET_ATTR_TXT("units", "%")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float RH2M_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RH2M_R", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Rural 2m specific humidity")
    GET_ATTR_TXT("units", "%")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float RH2M_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RH2M_U", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Urban 2m relative humidity")
    GET_ATTR_TXT("units", "%")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float RR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("RR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "root respiration (fine root MR + total root GR)")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SABG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SABG", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "solar rad absorbed by ground")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SABG_PEN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SABG_PEN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Rural solar rad penetrating top soil or snow layer")
    GET_ATTR_TXT("units", "watt/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SABV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SABV", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "solar rad absorbed by veg")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SCALARAVG_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SCALARAVG_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "average of decomposition scalar")
    GET_ATTR_TXT("units", "fraction")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SECONDP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SECONDP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil secondary P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SECONDP_TO_LABILEP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SECONDP_TO_LABILEP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "SECONDARY MINERAL P TO LABILE P")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SECONDP_TO_OCCLP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SECONDP_TO_OCCLP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "SECONDARY MINERAL P TO OCCLUDED P")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SECONDP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SECONDP_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil secondary P (vert. res.)")
    GET_ATTR_TXT("units", "gp/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SEEDC_GRC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SEEDC_GRC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "pool for seeding new PFTs via dynamic landcover")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil mineral N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINN_TO_NPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN_TO_NPOOL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "deployment of soil mineral N uptake")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINN_TO_PLANT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN_TO_PLANT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "plant uptake of soil mineral N")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINN_TO_SOIL1N_L1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN_TO_SOIL1N_L1", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mineral N flux for decomp. of LITR1to SOIL1")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINN_TO_SOIL2N_L2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN_TO_SOIL2N_L2", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mineral N flux for decomp. of LITR2to SOIL2")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINN_TO_SOIL2N_S1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN_TO_SOIL2N_S1", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mineral N flux for decomp. of SOIL1to SOIL2")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINN_TO_SOIL3N_L3(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN_TO_SOIL3N_L3", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mineral N flux for decomp. of LITR3to SOIL3")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINN_TO_SOIL3N_S2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN_TO_SOIL3N_S2", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mineral N flux for decomp. of SOIL2to SOIL3")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINN_TO_SOIL4N_S3(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINN_TO_SOIL4N_S3", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mineral N flux for decomp. of SOIL3to SOIL4")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil mineral P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINP_LEACHED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_LEACHED", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil mineral P pool loss to leaching")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINP_TO_PLANT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_TO_PLANT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "plant uptake of soil mineral P")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINP_TO_PPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_TO_PPOOL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "deployment of soil mineral P uptake")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINP_TO_SOIL1P_L1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_TO_SOIL1P_L1", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mineral P flux for decomp. of LITR1to SOIL1")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINP_TO_SOIL2P_L2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_TO_SOIL2P_L2", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mineral P flux for decomp. of LITR2to SOIL2")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINP_TO_SOIL2P_S1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_TO_SOIL2P_S1", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mineral P flux for decomp. of SOIL1to SOIL2")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINP_TO_SOIL3P_L3(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_TO_SOIL3P_L3", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mineral P flux for decomp. of LITR3to SOIL3")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINP_TO_SOIL3P_S2(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_TO_SOIL3P_S2", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mineral P flux for decomp. of SOIL2to SOIL3")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINP_TO_SOIL4P_S3(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMINP_TO_SOIL4P_S3", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mineral P flux for decomp. of SOIL3to SOIL4")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMINP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SMINP_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil mineral P (vert. res.)")
    GET_ATTR_TXT("units", "gp/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMIN_NH4(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMIN_NH4", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil mineral NH4")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMIN_NH4_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SMIN_NH4_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil mineral NH4 (vert. res.)")
    GET_ATTR_TXT("units", "gN/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMIN_NO3(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMIN_NO3", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil mineral NO3")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMIN_NO3_LEACHED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMIN_NO3_LEACHED", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil NO3 pool loss to leaching")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMIN_NO3_RUNOFF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SMIN_NO3_RUNOFF", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil NO3 pool loss to runoff")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SMIN_NO3_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SMIN_NO3_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil mineral NO3 (vert. res.)")
    GET_ATTR_TXT("units", "gN/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SNOBCMCL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOBCMCL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mass of BC in snow column")
    GET_ATTR_TXT("units", "kg/m2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SNOBCMSL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOBCMSL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mass of BC in top snow layer")
    GET_ATTR_TXT("units", "kg/m2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SNODSTMCL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNODSTMCL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mass of dust in snow column")
    GET_ATTR_TXT("units", "kg/m2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SNODSTMSL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNODSTMSL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mass of dust in top snow layer")
    GET_ATTR_TXT("units", "kg/m2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SNOINTABS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOINTABS", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Percent of incoming solar absorbed by lower snow layers")
    GET_ATTR_TXT("units", "%")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SNOOCMCL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOOCMCL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mass of OC in snow column")
    GET_ATTR_TXT("units", "kg/m2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SNOOCMSL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOOCMSL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mass of OC in top snow layer")
    GET_ATTR_TXT("units", "kg/m2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SNOW(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOW", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "atmospheric snow")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SNOWDP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOWDP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "gridcell mean snow height")
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SNOWICE(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOWICE", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "snow ice")
    GET_ATTR_TXT("units", "kg/m2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SNOWLIQ(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOWLIQ", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "snow liquid water")
    GET_ATTR_TXT("units", "kg/m2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SNOW_DEPTH(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOW_DEPTH", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "snow height of snow covered area")
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SNOW_SINKS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOW_SINKS", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "snow sinks (liquid water)")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SNOW_SOURCES(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SNOW_SOURCES", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "snow sources (liquid water)")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL1C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL1C", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "SOIL1 C")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL1C_TO_SOIL2C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL1C_TO_SOIL2C", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of soil 1 C to soil 2 C")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL1C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL1C_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "SOIL1 C (vertically resolved)")
    GET_ATTR_TXT("units", "gC/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL1N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL1N", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "SOIL1 N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL1N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL1N_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil 1 N tendency due to vertical transport")
    GET_ATTR_TXT("units", "gN/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL1N_TO_SOIL2N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL1N_TO_SOIL2N", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of soil 1 N to soil 2 N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL1N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL1N_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "SOIL1 N (vertically resolved)")
    GET_ATTR_TXT("units", "gN/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL1P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL1P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "SOIL1 P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL1P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL1P_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil 1 P tendency due to vertical transport")
    GET_ATTR_TXT("units", "gP/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL1P_TO_SOIL2P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL1P_TO_SOIL2P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of soil 1 P to soil 2 N")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL1P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL1P_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "SOIL1 P (vertically resolved)")
    GET_ATTR_TXT("units", "gP/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL1_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL1_HR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Het. Resp. from soil 1")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL2C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL2C", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "SOIL2 C")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL2C_TO_SOIL3C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL2C_TO_SOIL3C", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of soil 2 C to soil 3 C")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL2C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL2C_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "SOIL2 C (vertically resolved)")
    GET_ATTR_TXT("units", "gC/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL2N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL2N", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "SOIL2 N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL2N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL2N_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil 2 N tendency due to vertical transport")
    GET_ATTR_TXT("units", "gN/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL2N_TO_SOIL3N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL2N_TO_SOIL3N", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of soil 2 N to soil 3 N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL2N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL2N_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "SOIL2 N (vertically resolved)")
    GET_ATTR_TXT("units", "gN/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL2P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL2P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "SOIL2 P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL2P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL2P_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil 2 P tendency due to vertical transport")
    GET_ATTR_TXT("units", "gP/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL2P_TO_SOIL3P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL2P_TO_SOIL3P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of soil 2 P to soil 3 N")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL2P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL2P_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "SOIL2 P (vertically resolved)")
    GET_ATTR_TXT("units", "gP/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL2_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL2_HR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Het. Resp. from soil 2")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL3C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL3C", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "SOIL3 C")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL3C_TO_SOIL4C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL3C_TO_SOIL4C", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of soil 3 C to soil 4 C")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL3C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL3C_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "SOIL3 C (vertically resolved)")
    GET_ATTR_TXT("units", "gC/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL3N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL3N", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "SOIL3 N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL3N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL3N_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil 3 N tendency due to vertical transport")
    GET_ATTR_TXT("units", "gN/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL3N_TO_SOIL4N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL3N_TO_SOIL4N", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of soil 3 N to soil 4 N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL3N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL3N_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "SOIL3 N (vertically resolved)")
    GET_ATTR_TXT("units", "gN/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL3P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL3P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "SOIL3 P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL3P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL3P_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil 3 P tendency due to vertical transport")
    GET_ATTR_TXT("units", "gP/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL3P_TO_SOIL4P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL3P_TO_SOIL4P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "decomp. of soil 3 P to soil 4 N")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL3P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL3P_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "SOIL3 P (vertically resolved)")
    GET_ATTR_TXT("units", "gP/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL3_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL3_HR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Het. Resp. from soil 3")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL4C(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL4C", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "SOIL4 C")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL4C_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL4C_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "SOIL4 C (vertically resolved)")
    GET_ATTR_TXT("units", "gC/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL4N(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL4N", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "SOIL4 N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL4N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL4N_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil 4 N tendency due to vertical transport")
    GET_ATTR_TXT("units", "gN/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL4N_TO_SMINN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL4N_TO_SMINN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mineral N flux for decomp. of SOIL4")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL4N_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL4N_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "SOIL4 N (vertically resolved)")
    GET_ATTR_TXT("units", "gN/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL4P(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL4P", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "SOIL4 P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL4P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL4P_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil 4 P tendency due to vertical transport")
    GET_ATTR_TXT("units", "gP/m^3/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL4P_TO_SMINP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL4P_TO_SMINP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "mineral P flux for decomp. of SOIL4")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL4P_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOIL4P_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "SOIL4 P (vertically resolved)")
    GET_ATTR_TXT("units", "gP/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOIL4_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOIL4_HR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Het. Resp. from soil 4")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOILC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOILC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil C")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOILC_HR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOILC_HR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil C heterotrophic respiration")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOILC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOILC_LOSS", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil C loss")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOILICE(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOILICE", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil ice (vegetated landunits only)")
    GET_ATTR_TXT("units", "kg/m2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOILICE_ICE(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOILICE_ICE", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil ice (ice landunits only)")
    GET_ATTR_TXT("units", "kg/m2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOILLIQ(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOILLIQ", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil liquid water (vegetated landunits only)")
    GET_ATTR_TXT("units", "kg/m2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOILLIQ_ICE(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOILLIQ_ICE", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil liquid water (ice landunits only)")
    GET_ATTR_TXT("units", "kg/m2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOILPSI(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOILPSI", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil water potential in each soil layer")
    GET_ATTR_TXT("units", "MPa")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOILWATER_10CM(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOILWATER_10CM", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil liquid water + ice in top 10cm of soil (veg landunits only)")
    GET_ATTR_TXT("standard_name", "mass_content_of_water_in_soil_layer")
    GET_ATTR_TXT("units", "kg/m2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOLUTIONP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOLUTIONP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil solution P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOLUTIONP_vr(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("SOLUTIONP_vr", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil solution P (vert. res.)")
    GET_ATTR_TXT("units", "gp/m^3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOMHR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOMHR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil organic matter heterotrophic respiration")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SOM_C_LEACHED(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SOM_C_LEACHED", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total flux of C from SOM pools due to leaching")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total soil respiration (HR + root resp)")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float STORVEGC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("STORVEGC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "stored vegetation carbon, excluding cpool")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float STORVEGN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("STORVEGN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "stored vegetation nitrogen")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float STORVEGP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("STORVEGP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "stored vegetation phosphorus")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SUPPLEMENT_TO_SMINN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SUPPLEMENT_TO_SMINN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "supplemental N supply")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SUPPLEMENT_TO_SMINP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SUPPLEMENT_TO_SMINP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "supplemental P supply")
    GET_ATTR_TXT("units", "gP/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SUPPLY(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SUPPLY", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "runoff supply for land use")
    GET_ATTR_TXT("units", "mm/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SoilAlpha(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SoilAlpha", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "factor limiting ground evap")
    GET_ATTR_TXT("units", "1")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float SoilAlpha_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("SoilAlpha_U", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "urban factor limiting ground evap")
    GET_ATTR_TXT("units", "1")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TAUX(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TAUX", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "zonal surface stress")
    GET_ATTR_TXT("units", "kg/m/s^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TAUY(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TAUY", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "meridional surface stress")
    GET_ATTR_TXT("units", "kg/m/s^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TBOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TBOT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "atmospheric air temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TBUILD(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TBUILD", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "internal urban building temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TCS_MONTH_BEGIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TCS_MONTH_BEGIN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total carbon storage at the beginning of a month")
    GET_ATTR_TXT("units", "mm")
    GET_ATTR_TXT("cell_methods", "time: point")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TCS_MONTH_END(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TCS_MONTH_END", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total carbon storage at the end of a month")
    GET_ATTR_TXT("units", "mm")
    GET_ATTR_TXT("cell_methods", "time: point")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TG", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "ground temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TG_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TG_R", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Rural ground temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TG_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TG_U", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Urban ground temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TH2OSFC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TH2OSFC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "surface water temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float THBOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("THBOT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "atmospheric air potential temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TKE1(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TKE1", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "top lake level eddy thermal conductivity")
    GET_ATTR_TXT("units", "W/(mK)")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TLAI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TLAI", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total projected leaf area index")
    GET_ATTR_TXT("units", "1")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TLAKE(time, levlak, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levlak;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("TLAKE", NC_FLOAT, 4, dimids, REC_ITYPE, 2)
    GET_ATTR_TXT("long_name", "lake temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTCOLC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTCOLC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total column carbon, incl veg and cpool but excl product pools")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTCOLCH4(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTCOLCH4", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total belowground CH4, (0 for non-lake special landunits)")
    GET_ATTR_TXT("units", "gC/m2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTCOLN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTCOLN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total column-level N but excl product pools")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTCOLP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTCOLP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total column-level P but excl product pools")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTECOSYSC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTECOSYSC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total ecosystem carbon, incl veg but excl cpool but excl product pools")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTECOSYSN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTECOSYSN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total ecosystem N but excl product pools")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTECOSYSP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTECOSYSP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total ecosystem P but excl product pools")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTLITC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTLITC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total litter carbon")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTLITC_1m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTLITC_1m", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total litter carbon to 1 meter depth")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTLITN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTLITN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total litter N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTLITP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTLITP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total litter P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTLITP_1m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTLITP_1m", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total litter P to 1 meter")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTPFTC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTPFTC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total patch-level carbon, including cpool")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTPFTN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTPFTN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total PFT-level nitrogen")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTPFTP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTPFTP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total PFT-level phosphorus")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTSOMC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTSOMC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total soil organic matter carbon")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTSOMC_1m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTSOMC_1m", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total soil organic matter carbon to 1 meter depth")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTSOMN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTSOMN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total soil organic matter N")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTSOMP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTSOMP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total soil organic matter P")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTSOMP_1m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTSOMP_1m", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total soil organic matter P to 1 meter")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTVEGC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTVEGC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total vegetation carbon, excluding cpool")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTVEGC_ABG(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTVEGC_ABG", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total aboveground vegetation carbon, excluding cpool")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTVEGN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTVEGN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total vegetation nitrogen")
    GET_ATTR_TXT("units", "gN/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TOTVEGP(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TOTVEGP", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total vegetation phosphorus")
    GET_ATTR_TXT("units", "gP/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TREFMNAV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TREFMNAV", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "daily minimum of average 2-m temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TREFMNAV_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TREFMNAV_R", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Rural daily minimum of average 2-m temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TREFMNAV_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TREFMNAV_U", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Urban daily minimum of average 2-m temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TREFMXAV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TREFMXAV", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "daily maximum of average 2-m temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TREFMXAV_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TREFMXAV_R", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Rural daily maximum of average 2-m temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TREFMXAV_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TREFMXAV_U", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Urban daily maximum of average 2-m temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TSA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TSA", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "2m air temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TSAI(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TSAI", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total projected stem area index")
    GET_ATTR_TXT("units", "1")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TSA_R(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TSA_R", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Rural 2m air temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TSA_U(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TSA_U", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "Urban 2m air temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TSOI(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("TSOI", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil temperature (vegetated landunits only)")
    GET_ATTR_TXT("standard_name", "soil_temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TSOI_10CM(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TSOI_10CM", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil temperature in top 10cm of soil")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TSOI_ICE(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("TSOI_ICE", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "soil temperature (ice landunits only)")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TV(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TV", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "vegetation temperature")
    GET_ATTR_TXT("units", "K")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TWS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TWS", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total water storage")
    GET_ATTR_TXT("units", "mm")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TWS_MONTH_BEGIN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TWS_MONTH_BEGIN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total water storage at the beginning of a month")
    GET_ATTR_TXT("units", "mm")
    GET_ATTR_TXT("cell_methods", "time: point")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float TWS_MONTH_END(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("TWS_MONTH_END", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total water storage at the end of a month")
    GET_ATTR_TXT("units", "mm")
    GET_ATTR_TXT("cell_methods", "time: point")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float T_SCALAR(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("T_SCALAR", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "temperature inhibition of decomposition")
    GET_ATTR_TXT("units", "1")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float U10(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("U10", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "10-m wind")
    GET_ATTR_TXT("units", "m/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float URBAN_AC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("URBAN_AC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "urban air conditioning flux")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float URBAN_HEAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("URBAN_HEAT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "urban heating flux")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float VOLR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("VOLR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "river channel total water storage")
    GET_ATTR_TXT("units", "m3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float VOLRMCH(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("VOLRMCH", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "river channel main channel water storage")
    GET_ATTR_TXT("units", "m3")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float WA(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WA", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "water in the unconfined aquifer (vegetated landunits only)")
    GET_ATTR_TXT("units", "mm")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float WASTEHEAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WASTEHEAT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "sensible heat flux from heating/cooling sources of urban waste heat")
    GET_ATTR_TXT("units", "W/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float WF(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WF", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "soil water as frac. of whc for top 0.05 m")
    GET_ATTR_TXT("units", "proportion")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float WIND(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WIND", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "atmospheric wind velocity magnitude")
    GET_ATTR_TXT("units", "m/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float WOODC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WOODC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "wood C")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float WOODC_ALLOC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WOODC_ALLOC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "wood C eallocation")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float WOODC_LOSS(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WOODC_LOSS", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "wood C loss")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float WOOD_HARVESTC(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WOOD_HARVESTC", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "wood harvest carbon (to product pools)")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float WOOD_HARVESTN(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WOOD_HARVESTN", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "wood harvest N (to product pools)")
    GET_ATTR_TXT("units", "gN/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float WTGQ(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("WTGQ", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "surface tracer conductance")
    GET_ATTR_TXT("units", "m/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float W_SCALAR(time, levdcmp, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levdcmp;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("W_SCALAR", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "Moisture (dryness) inhibition of decomposition")
    GET_ATTR_TXT("units", "1")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float XR(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("XR", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "total excess respiration")
    GET_ATTR_TXT("units", "gC/m^2/s")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float XSMRPOOL(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("XSMRPOOL", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "temporary photosynthate C pool")
    GET_ATTR_TXT("units", "gC/m^2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float ZBOT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ZBOT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "atmospheric reference height")
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float ZWT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ZWT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "water table depth (vegetated landunits only)")
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float ZWT_CH4_UNSAT(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ZWT_CH4_UNSAT", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "depth of water table for methane production used in non-inundated area")
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float ZWT_PERCH(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("ZWT_PERCH", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "perched water table depth (vegetated landunits only)")
    GET_ATTR_TXT("units", "m")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float cn_scalar(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("cn_scalar", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "N limitation factor")
    GET_ATTR_TXT("units", "")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float cp_scalar(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("cp_scalar", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "P limitation factor")
    GET_ATTR_TXT("units", "")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float leaf_npimbalance(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("leaf_npimbalance", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "leaf np imbalance partial C partial P/partial C partial N")
    GET_ATTR_TXT("units", "gN/gP")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float nlim_m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("nlim_m", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "runmean N limitation factor")
    GET_ATTR_TXT("units", "")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float o2_decomp_depth_unsat(time, levgrnd, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_levgrnd;
    dimids[2] = dim_lat;
    dimids[3] = dim_lon;
    INQ_VAR("o2_decomp_depth_unsat", NC_FLOAT, 4, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("long_name", "o2_decomp_depth_unsat")
    GET_ATTR_TXT("units", "mol/m3/2")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    /* float plim_m(time, lat, lon) */
    dimids[0] = dim_time;
    dimids[1] = dim_lat;
    dimids[2] = dim_lon;
    INQ_VAR("plim_m", NC_FLOAT, 3, dimids, REC_ITYPE, 0)
    GET_ATTR_TXT("long_name", "runmean P limitation factor")
    GET_ATTR_TXT("units", "")
    GET_ATTR_TXT("cell_methods", "time: mean")
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", missv)

    assert(varp - vars + 1 == cfg.nvars + nvars_decomp);

err_out:
    return err;
}

