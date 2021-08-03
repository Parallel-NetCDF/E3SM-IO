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
#include <e3sm_io_case_scorpio.hpp>
#include <e3sm_io_driver.hpp>

#define CHECK_VAR_ERR(varid)                                                                  \
    {                                                                                         \
        if (err != 0) {                                                                       \
            char var_name[64];                                                                \
            driver.inq_var_name (ncid, (varid.data), var_name);                               \
            printf ("Error in %s:%d: %s() var %s\n", __FILE__, __LINE__, __func__, var_name); \
            goto err_out;                                                                     \
        }                                                                                     \
    }
#define DEF_DIM(name, num, dimid)                                                  \
    {                                                                              \
        err = e3sm_io_scorpio_define_dim (driver, ncid, name, num, dnames, dimid); \
        CHECK_ERR                                                                  \
    }
#define DEF_VAR(name, type, ndims, dimids)                                                       \
    {                                                                                            \
        if (op) {                                                                          \
        err = e3sm_io_scorpio_define_var (                                                       \
            driver, cfg, dnames, decom, vars[varid].decomp_id,                                    \
            vars[varid].decomp_id, ncid, name, \
            type, ndims, dimids, &var);                                                          \
        if (err != 0) {                                                                          \
            printf ("Error in %s line %d: def_var %s\n", __FILE__, __LINE__, name);              \
            goto err_out;                                                                        \
        }                                                                                        \
        } \
    }

#define PUT_GATTR_TXT(name, buf)                                                             \
    {                                                                                        \
        err = e3sm_io_scorpio_put_att (driver, ncid, NC_GLOBAL, "pio_global/" name, NC_CHAR, \
                                       strlen (buf), (void *)buf);                           \
        CHECK_ERR                                                                            \
    }

#define PUT_PIOATTR_INT1(name, val)                                                               \
    {                                                                                             \
        int buf = val;                                                                            \
        err = e3sm_io_scorpio_put_att (driver, ncid, NC_GLOBAL, name, NC_INT, 1, (void *)(&buf)); \
        CHECK_ERR                                                                                 \
    }
#define PUT_GATTR_INT1(name, val) PUT_PIOATTR_INT1 (("pio_global/" name), val)

#define PUT_GATTR_DBL1(name, val)                                                                 \
    {                                                                                             \
        double buf = val;                                                                         \
        err = e3sm_io_scorpio_put_att (driver, ncid, NC_GLOBAL, "pio_global/" name, NC_DOUBLE, 1, \
                                       (void *)(&buf));                                           \
        CHECK_ERR                                                                                 \
    }

#define PUT_ATTR_TXT(name, buf)                                                            \
    {                                                                                      \
        if (op) {                                                                          \
            err = e3sm_io_scorpio_put_att (driver, ncid, var, name, NC_CHAR, strlen (buf), \
                                           (void *)buf);                                   \
            CHECK_VAR_ERR (var)                                                            \
        }                                                                                  \
    }

#define PUT_ATTR_INT1(name, val)                                                             \
    {                                                                                        \
        if (op) {                                                                            \
            int buf = val;                                                           \
            err = e3sm_io_scorpio_put_att (driver, ncid, var, name, NC_INT, 1, (void*)(&buf)); \
            CHECK_VAR_ERR (var)                                                              \
        }                                                                                    \
    }
#define PUT_ATTR_INT(name, num, buf)                                                           \
    {                                                                                          \
        if (op) {                                                                              \
            err = e3sm_io_scorpio_put_att (driver, ncid, var, name, NC_INT, num, (void *)buf); \
            CHECK_VAR_ERR (var)                                                                \
        }                                                                                      \
    }
#define PUT_ATTR_FLOAT(name, num, buf)                                                            \
    {                                                                                             \
        if (op) {                                                                                 \
            float flt = buf;                                                                      \
            err =                                                                                 \
                e3sm_io_scorpio_put_att (driver, ncid, var, name, NC_FLOAT, num, (void *)(&flt)); \
            CHECK_VAR_ERR (var)                                                                   \
        }                                                                                         \
    }
#define PUT_ATTR_INT64(name, num, buf)                                                         \
    {                                                                                          \
        if (op) {                                                                              \
            err =                                                                       \
                e3sm_io_scorpio_put_att (driver, ncid, var, name, NC_INT64, num, (void *)buf); \
            CHECK_VAR_ERR (var)                                                                \
        }                                                                                      \
    }

#define PUT_ATTR_DECOMP(D, ndims, dimids) \
    {}

#define SET_VAR_META(dtype, id, rec, buflen, varlen)     \
    {                                                    \
        if (op) {                                        \
            vars[varid].vid = var;                       \
            wr_buf->buflen += vars[varid].vlen + wr_buf->gap;      \
        } else {                                         \
            vars[varid].itype     = dtype;               \
            vars[varid].decomp_id = id;                  \
            vars[varid].isRecVar  = rec;                 \
            if (id >= 0) { \
            vars[varid].vlen      = decom.raw_nreqs[id]; \
            } \
            else { \
            vars[varid].vlen      = varlen; \
            } \
        }                                                \
        varid++;                                         \
    }

/*----< add_gattrs() >-------------------------------------------------------*/
static
int add_gattrs(e3sm_io_config &cfg,
               e3sm_io_decom  &decom,
               e3sm_io_driver &driver,
               int             ncid)
{
    int err=0;

    // PIO attributes
    PUT_PIOATTR_INT1 ("/__pio__/fillmode", 256);

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

    PUT_GATTR_INT1("ltype_vegetated_or_bare_soil", 1)
    PUT_GATTR_INT1("ltype_crop", 2)
    PUT_GATTR_INT1("ltype_landice", 3)
    PUT_GATTR_INT1("ltype_landice_multiple_elevation_classes", 4)
    PUT_GATTR_INT1("ltype_deep_lake", 5)
    PUT_GATTR_INT1("ltype_wetland", 6)
    PUT_GATTR_INT1("ltype_urban_tbd", 7)
    PUT_GATTR_INT1("ltype_urban_hd", 8)
    PUT_GATTR_INT1("ltype_urban_md", 9)

    if (cfg.hist == h1) { /* h1 only */
        PUT_GATTR_TXT("Time_constant_3Dvars_filenamae", "./I1850GSWCNPRDCTCBC_hcru_hcru.elm.h0.0001-01-01-00000.nc")
        PUT_GATTR_TXT("Time_constant_3Dvars", "ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE")
    }

err_out:
    return err;
}

/*----< def_I_case() >-------------------------------------------------------*/
int def_I_case_scorpio(e3sm_io_config   &cfg,
               e3sm_io_decom    &decom,
               e3sm_io_driver   &driver,
               int               ncid,    /* file ID */
               var_meta_scorpio         *vars,    /* variable metadata */
               io_buffers       *wr_buf,
               int *scorpiovars,
               bool op
               )
{
    /* Total 52 variables */
    char name[64];
    int i, j, k, err = 0, ndims, dimids[2], varid = 0;
    int fix_dimids[MAX_NUM_DECOMP][4], rec_dimids[MAX_NUM_DECOMP][4];
    int dim_time, dim_lon, dim_lat, dim_gridcell, dim_topounit, dim_landunit;
    int dim_column, dim_pft, dim_levgrnd, dim_levurb, dim_levlak, dim_numrad;
    int dim_month, dim_levsno, dim_ltype, dim_nvegwcs, dim_natpft;
    int dim_string_length, dim_levdcmp, dim_levtrc, dim_hist_interval;
    MPI_Offset lat, lon, levgrnd, levdcmp, levlak, ltype, natpft;
    MPI_Offset string_length, hist_interval;
    std::map<int, std::string> dnames;
    e3sm_io_scorpio_var var;

    if (op) {
        err = add_gattrs(cfg, decom, driver, ncid);
        CHECK_ERR
    }

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

    lat           = 360;
    lon           = 720;
    levgrnd       =  15;
    levdcmp       =  15;
    levlak        =  10;
    ltype         =   9;
    natpft        =  17;
    string_length =  16;
    hist_interval =   2;

    /* define dimensions */
    if(op){
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
    }

    /* canonical or ADIOS blob */
    fix_dimids[0][0] = dim_lat;
    fix_dimids[0][1] = dim_lon;
    rec_dimids[0][0] = dim_time;
    rec_dimids[0][1] = dim_lat;
    rec_dimids[0][2] = dim_lon;

    for (i=1; i<6; i++) {
        fix_dimids[i][1] = dim_lat;
        fix_dimids[i][2] = dim_lon;
        rec_dimids[i][0] = dim_time;
        rec_dimids[i][2] = dim_lat;
        rec_dimids[i][3] = dim_lon;
    }

    fix_dimids[1][0] = dim_levgrnd;   /* same size as levdcmp */
    fix_dimids[2][0] = dim_levlak;
    fix_dimids[3][0] = dim_ltype;
    fix_dimids[4][0] = dim_natpft;
    fix_dimids[5][0] = dim_levdcmp;   /* same size as levgrnd */

    rec_dimids[1][1] = dim_levgrnd;   /* same size as levdcmp */
    rec_dimids[2][1] = dim_levlak;
    rec_dimids[3][1] = dim_ltype;
    rec_dimids[4][1] = dim_natpft;
    rec_dimids[5][1] = dim_levdcmp;   /* same size as levgrnd */

    /* define scorpio decom map variables */
    if (op){
        for (j = 0; j < decom.num_decomp; j++) {
            int piodims[5];

            sprintf (name, "/__pio__/decomp/%d", (j + 512));
            err = driver.def_local_var (ncid, name, NC_INT64, 1, decom.raw_nreqs + j,
                                        scorpiovars + j);
            CHECK_ERR

            for (i = 0; i < decom.ndims[j]; i++) {
                piodims[i] = (int)decom.dims[j][i];
            }
            err = driver.put_att (ncid, scorpiovars[j], "dimlen", NC_INT, decom.ndims[j],
                                piodims);
            CHECK_ERR

            err =
                driver.put_att (ncid, scorpiovars[j], "ndims", NC_INT, 1, decom.ndims + j);
            CHECK_ERR

            k   = 6;
            err = driver.put_att (ncid, scorpiovars[j], "piotype", NC_INT, 1, &k);
            CHECK_ERR
        }
        err = driver.def_local_var (ncid, "/__pio__/info/nproc", NC_INT, 0, NULL, scorpiovars + (j++));
        CHECK_ERR
    }

    /* define 560 climate variables
     *     18 are fixed-size,      542 are record variables
     *     14 are not partitioned, 546 are partitioned.
     * Thus, 14 not partitioned variable (also small) are written by root only:
     *     5 fixed-size, 9 record variables
     */

    /* float levgrnd(levgrnd) */
    DEF_VAR("levgrnd", NC_FLOAT, 1, &dim_levgrnd)
    PUT_ATTR_TXT("long_name", "coordinate soil levels")
    PUT_ATTR_TXT("units", "m")
    SET_VAR_META(REC_ITYPE, -1, 0, fix_buflen, levgrnd)

    /* float levlak(levlak) */
    DEF_VAR("levlak", NC_FLOAT, 1, &dim_levlak)
    PUT_ATTR_TXT("long_name", "coordinate lake levels")
    PUT_ATTR_TXT("units", "m")
    SET_VAR_META(REC_ITYPE, -1, 0, fix_buflen, levlak)

    /* float levdcmp(levdcmp) */
    DEF_VAR("levdcmp", NC_FLOAT, 1, &dim_levdcmp)
    PUT_ATTR_TXT("long_name", "coordinate soil levels")
    PUT_ATTR_TXT("units", "m")
    SET_VAR_META(REC_ITYPE, -1, 0, fix_buflen, levdcmp)

    /* float time(time) */
    DEF_VAR("time", NC_FLOAT, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "time")
    PUT_ATTR_TXT("units", "days since 0001-01-01 00:00:00")
    PUT_ATTR_TXT("calendar", "noleap")
    PUT_ATTR_TXT("bounds", "time_bounds")
    SET_VAR_META(REC_ITYPE, -1, 1, rec_buflen, 1)

    /* int mcdate(time) */
    DEF_VAR("mcdate", NC_INT, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "current date (YYYYMMDD)")
    SET_VAR_META(MPI_INT, -1, 1, rec_int_buflen, 1)

    /* int mcsec(time) */
    DEF_VAR("mcsec", NC_INT, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "current seconds of current date")
    PUT_ATTR_TXT("units", "s")
    SET_VAR_META(MPI_INT, -1, 1, rec_int_buflen, 1)

    /* int mdcur(time) */
    DEF_VAR("mdcur", NC_INT, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "current day (from base day)")
    SET_VAR_META(MPI_INT, -1, 1, rec_int_buflen, 1)

    /* int mscur(time) */
    DEF_VAR("mscur", NC_INT, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "current seconds of current day")
    SET_VAR_META(MPI_INT, -1, 1, rec_int_buflen, 1)

    /* int nstep(time) */
    DEF_VAR("nstep", NC_INT, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "time step")
    SET_VAR_META(MPI_INT, -1, 1, rec_int_buflen, 1)

    /* double time_bounds(time, hist_interval) */
    dimids[0] = dim_time; dimids[1] = dim_hist_interval;
    DEF_VAR("time_bounds", NC_DOUBLE, 2, dimids)
    PUT_ATTR_TXT("long_name", "history time interval endpoints")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, -1, 1, rec_buflen, hist_interval)
#else
    SET_VAR_META(MPI_DOUBLE, -1, 1, rec_dbl_buflen, hist_interval)
#endif

    /* char date_written(time, string_length) */
    dimids[0] = dim_time; dimids[1] = dim_string_length;
    DEF_VAR("date_written", NC_CHAR, 2, dimids)
    SET_VAR_META(MPI_CHAR, -1, 1, rec_txt_buflen, string_length)

    /* char time_written(time, string_length) */
    dimids[0] = dim_time; dimids[1] = dim_string_length;
    DEF_VAR("time_written", NC_CHAR, 2, dimids)
    SET_VAR_META(MPI_CHAR, -1, 1, rec_txt_buflen, string_length)

    /* float lon(lon) */
    DEF_VAR("lon", NC_FLOAT, 1, &dim_lon)
    PUT_ATTR_TXT("long_name", "coordinate longitude")
    PUT_ATTR_TXT("units", "degrees_east")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    SET_VAR_META(REC_ITYPE, -1, 0, fix_buflen, lon)

    /* float lat(lat) */
    DEF_VAR("lat", NC_FLOAT, 1, &dim_lat)
    PUT_ATTR_TXT("long_name", "coordinate latitude")
    PUT_ATTR_TXT("units", "degrees_north")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    SET_VAR_META(REC_ITYPE, -1, 0, fix_buflen, lat)

    /* float area(lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 1 : 2;
    DEF_VAR("area", NC_FLOAT, ndims, fix_dimids[0])
    PUT_ATTR_TXT("long_name", "grid cell areas")
    PUT_ATTR_TXT("units", "km^2")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 2, orig_dimids[0]+1)
    SET_VAR_META(REC_ITYPE, 0, 0, fix_buflen, decom.count[0])

    /* float topo(lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 1 : 2;
    DEF_VAR("topo", NC_FLOAT, ndims, fix_dimids[0])
    PUT_ATTR_TXT("long_name", "grid cell topography")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 2, orig_dimids[0]+1)
    SET_VAR_META(REC_ITYPE, 0, 0, fix_buflen, decom.count[0])

    /* float landfrac(lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 1 : 2;
    DEF_VAR("landfrac", NC_FLOAT, ndims, fix_dimids[0])
    PUT_ATTR_TXT("long_name", "land fraction")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 2, orig_dimids[0]+1)
    SET_VAR_META(REC_ITYPE, 0, 0, fix_buflen, decom.count[0])

    /* int landmask(lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 1 : 2;
    DEF_VAR("landmask", NC_INT, ndims, fix_dimids[0])
    PUT_ATTR_TXT("long_name", "land/ocean mask (0.=ocean and 1.=land)")
    PUT_ATTR_INT1(_FillValue, -9999)
    PUT_ATTR_INT1("missing_value", -9999)
    PUT_ATTR_DECOMP(one, 2, orig_dimids[0]+1)
    SET_VAR_META(MPI_INT, 0, 0, fix_int_buflen, decom.count[0])

    /* int pftmask(lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 1 : 2;
    DEF_VAR("pftmask", NC_INT, ndims, fix_dimids[0])
    PUT_ATTR_TXT("long_name", "pft real/fake mask (0.=fake and 1.=real)")
    PUT_ATTR_INT1(_FillValue, -9999)
    PUT_ATTR_INT1("missing_value", -9999)
    PUT_ATTR_DECOMP(one, 2, orig_dimids[0]+1)
    SET_VAR_META(MPI_INT, 0, 0, fix_int_buflen, decom.count[0])

    if (cfg.hist == h0) {  /* h0 only */
        /* float ZSOI(levgrnd, lat, lon) */
        ndims = (false && (cfg.strategy == blob)) ? 1 : 3;
        DEF_VAR("ZSOI", NC_FLOAT, ndims, fix_dimids[1])
        PUT_ATTR_TXT("long_name", "soil depth")
        PUT_ATTR_TXT("units", "m")
        PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
        PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
        PUT_ATTR_DECOMP(two, 3, orig_dimids[1]+1)
        SET_VAR_META(REC_ITYPE, 1, 0, fix_buflen, decom.count[1])

        /* float DZSOI(levgrnd, lat, lon) */
        ndims = (false && (cfg.strategy == blob)) ? 1 : 3;
        DEF_VAR("DZSOI", NC_FLOAT, ndims, fix_dimids[1])
        PUT_ATTR_TXT("long_name", "soil thickness")
        PUT_ATTR_TXT("units", "m")
        PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
        PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
        PUT_ATTR_DECOMP(two, 3, orig_dimids[1]+1)
        SET_VAR_META(REC_ITYPE, 1, 0, fix_buflen, decom.count[1])

        /* float WATSAT(levgrnd, lat, lon) */
        ndims = (false && (cfg.strategy == blob)) ? 1 : 3;
        DEF_VAR("WATSAT", NC_FLOAT, ndims, fix_dimids[1])
        PUT_ATTR_TXT("long_name", "saturated soil water content (porosity)")
        PUT_ATTR_TXT("units", "mm3/mm3")
        PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
        PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
        PUT_ATTR_DECOMP(two, 3, orig_dimids[1]+1)
        SET_VAR_META(REC_ITYPE, 1, 0, fix_buflen, decom.count[1])

        /* float SUCSAT(levgrnd, lat, lon) */
        ndims = (false && (cfg.strategy == blob)) ? 1 : 3;
        DEF_VAR("SUCSAT", NC_FLOAT, ndims, fix_dimids[1])
        PUT_ATTR_TXT("long_name", "saturated soil matric potential")
        PUT_ATTR_TXT("units", "mm")
        PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
        PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
        PUT_ATTR_DECOMP(two, 3, orig_dimids[1]+1)
        SET_VAR_META(REC_ITYPE, 1, 0, fix_buflen, decom.count[1])

        /* float BSW(levgrnd, lat, lon) */
        ndims = (false && (cfg.strategy == blob)) ? 1 : 3;
        DEF_VAR("BSW", NC_FLOAT, ndims, fix_dimids[1])
        PUT_ATTR_TXT("long_name", "slope of soil water retention curve")
        PUT_ATTR_TXT("units", "1")
        PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
        PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
        PUT_ATTR_DECOMP(two, 3, orig_dimids[1]+1)
        SET_VAR_META(REC_ITYPE, 1, 0, fix_buflen, decom.count[1])

        /* float HKSAT(levgrnd, lat, lon) */
        ndims = (false && (cfg.strategy == blob)) ? 1 : 3;
        DEF_VAR("HKSAT", NC_FLOAT, ndims, fix_dimids[1])
        PUT_ATTR_TXT("long_name", "saturated hydraulic conductivity")
        PUT_ATTR_TXT("units", "1")
        PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
        PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
        PUT_ATTR_DECOMP(two, 3, orig_dimids[1]+1)
        SET_VAR_META(REC_ITYPE, 1, 0, fix_buflen, decom.count[1])

        /* float ZLAKE(levlak, lat, lon) */
        ndims = (false && (cfg.strategy == blob)) ? 1 : 3;
        DEF_VAR("ZLAKE", NC_FLOAT, ndims, fix_dimids[2])
        PUT_ATTR_TXT("long_name", "lake layer node depth")
        PUT_ATTR_TXT("units", "m")
        PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
        PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
        PUT_ATTR_DECOMP(three, 3, orig_dimids[2]+1)
        SET_VAR_META(REC_ITYPE, 2, 0, fix_buflen, decom.count[2])

        /* float DZLAKE(levlak, lat, lon) */
        ndims = (false && (cfg.strategy == blob)) ? 1 : 3;
        DEF_VAR("DZLAKE", NC_FLOAT, ndims, fix_dimids[2])
        PUT_ATTR_TXT("long_name", "lake layer thickness")
        PUT_ATTR_TXT("units", "m")
        PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
        PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
        PUT_ATTR_DECOMP(three, 3, orig_dimids[2]+1)
        SET_VAR_META(REC_ITYPE, 2, 0, fix_buflen, decom.count[2])
    }

    /* float ACTUAL_IMMOB(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ACTUAL_IMMOB", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "actual N immobilization")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float ACTUAL_IMMOB_P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ACTUAL_IMMOB_P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "actual P immobilization")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float ADSORBTION_P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ADSORBTION_P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "adsorb P flux")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float AGNPP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("AGNPP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "aboveground NPP")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float AGWDNPP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("AGWDNPP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "aboveground wood NPP")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float ALT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ALT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "current active layer thickness")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float ALTMAX(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ALTMAX", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "maximum annual active layer thickness")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float ALTMAX_LASTYEAR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ALTMAX_LASTYEAR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "maximum prior year active layer thickness")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float AR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("AR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "autotrophic respiration (MR + GR)")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float AVAILC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("AVAILC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "C flux available for allocation")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float AVAIL_RETRANSP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("AVAIL_RETRANSP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "P flux available from retranslocation pool")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float BAF_CROP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("BAF_CROP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "fractional area burned for crop")
    PUT_ATTR_TXT("units", "proportion/sec")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float BAF_PEATF(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("BAF_PEATF", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "fractional area burned in peatland")
    PUT_ATTR_TXT("units", "proportion/sec")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float BCDEP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("BCDEP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total BC deposition (dry+wet) from atmosphere")
    PUT_ATTR_TXT("units", "kg/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float BGNPP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("BGNPP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "belowground NPP")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float BIOCHEM_PMIN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("BIOCHEM_PMIN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "biochemical rate of P mineralization")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float BIOCHEM_PMIN_TO_PLANT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("BIOCHEM_PMIN_TO_PLANT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "plant uptake of biochemical P mineralization")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float BTRAN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("BTRAN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "transpiration beta factor")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float BUILDHEAT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("BUILDHEAT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "heat flux from urban building interior to walls and roof")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CH4PROD(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CH4PROD", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Gridcell total production of CH4")
    PUT_ATTR_TXT("units", "gC/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CH4_SURF_AERE_SAT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CH4_SURF_AERE_SAT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "aerenchyma surface CH4 flux for inundated area; (+ to atm)")
    PUT_ATTR_TXT("units", "mol/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CH4_SURF_AERE_UNSAT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CH4_SURF_AERE_UNSAT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "aerenchyma surface CH4 flux for non-inundated area; (+ to atm)")
    PUT_ATTR_TXT("units", "mol/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CH4_SURF_DIFF_SAT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CH4_SURF_DIFF_SAT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "diffusive surface CH4 flux for inundated / lake area; (+ to atm)")
    PUT_ATTR_TXT("units", "mol/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CH4_SURF_DIFF_UNSAT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CH4_SURF_DIFF_UNSAT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "diffusive surface CH4 flux for non-inundated area; (+ to atm)")
    PUT_ATTR_TXT("units", "mol/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CH4_SURF_EBUL_SAT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CH4_SURF_EBUL_SAT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "ebullition surface CH4 flux for inundated / lake area; (+ to atm)")
    PUT_ATTR_TXT("units", "mol/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CH4_SURF_EBUL_UNSAT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CH4_SURF_EBUL_UNSAT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "ebullition surface CH4 flux for non-inundated area; (+ to atm)")
    PUT_ATTR_TXT("units", "mol/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float COL_PTRUNC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("COL_PTRUNC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "column-level sink for P truncation")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CONC_CH4_SAT(time, levgrnd, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("CONC_CH4_SAT", NC_FLOAT, ndims, rec_dimids[1])
    PUT_ATTR_TXT("long_name", "CH4 soil Concentration for inundated / lake area")
    PUT_ATTR_TXT("units", "mol/m3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[1])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float CONC_CH4_UNSAT(time, levgrnd, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("CONC_CH4_UNSAT", NC_FLOAT, ndims, rec_dimids[1])
    PUT_ATTR_TXT("long_name", "CH4 soil Concentration for non-inundated area")
    PUT_ATTR_TXT("units", "mol/m3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[1])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float CONC_O2_SAT(time, levgrnd, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("CONC_O2_SAT", NC_FLOAT, ndims, rec_dimids[1])
    PUT_ATTR_TXT("long_name", "O2 soil Concentration for inundated / lake area")
    PUT_ATTR_TXT("units", "mol/m3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[1])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float CONC_O2_UNSAT(time, levgrnd, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("CONC_O2_UNSAT", NC_FLOAT, ndims, rec_dimids[1])
    PUT_ATTR_TXT("long_name", "O2 soil Concentration for non-inundated area")
    PUT_ATTR_TXT("units", "mol/m3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[1])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float CPOOL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CPOOL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "temporary photosynthate C pool")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CWDC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CWDC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "CWD C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CWDC_HR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CWDC_HR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "coarse woody debris C heterotrophic respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CWDC_LOSS(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CWDC_LOSS", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "coarse woody debris C loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CWDC_TO_LITR2C(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CWDC_TO_LITR2C", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of coarse woody debris C to litter 2 C")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CWDC_TO_LITR3C(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CWDC_TO_LITR3C", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of coarse woody debris C to litter 3 C")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CWDC_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("CWDC_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "CWD C (vertically resolved)")
    PUT_ATTR_TXT("units", "gC/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float CWDN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CWDN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "CWD N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CWDN_TO_LITR2N(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CWDN_TO_LITR2N", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of coarse woody debris N to litter 2 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CWDN_TO_LITR3N(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CWDN_TO_LITR3N", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of coarse woody debris N to litter 3 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CWDN_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("CWDN_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "CWD N (vertically resolved)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float CWDP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CWDP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "CWD P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CWDP_TO_LITR2P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CWDP_TO_LITR2P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of coarse woody debris P to litter 2 N")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CWDP_TO_LITR3P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("CWDP_TO_LITR3P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of coarse woody debris P to litter 3 N")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float CWDP_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("CWDP_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "CWD P (vertically resolved)")
    PUT_ATTR_TXT("units", "gP/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float DEADCROOTC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DEADCROOTC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "dead coarse root C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DEADCROOTN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DEADCROOTN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "dead coarse root N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DEADCROOTP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DEADCROOTP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "dead coarse root P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DEADSTEMC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DEADSTEMC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "dead stem C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DEADSTEMN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DEADSTEMN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "dead stem N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DEADSTEMP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DEADSTEMP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "dead stem P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DEFICIT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DEFICIT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "runoff supply deficit")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DENIT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DENIT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total rate of denitrification")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DESORPTION_P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DESORPTION_P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "desorp P flux")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DISPVEGC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DISPVEGC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "displayed veg carbon, excluding storage and cpool")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DISPVEGN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DISPVEGN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "displayed vegetation nitrogen")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DISPVEGP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DISPVEGP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "displayed vegetation phosphorus")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DSTDEP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DSTDEP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total dust deposition (dry+wet) from atmosphere")
    PUT_ATTR_TXT("units", "kg/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DSTFLXT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DSTFLXT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total surface dust emission")
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DWB(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DWB", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "net change in total water mass")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DWT_CONV_CFLUX_DRIBBLED(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DWT_CONV_CFLUX_DRIBBLED", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "conversion C flux (immediate loss to atm), dribbled throughout the year")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DWT_CONV_CFLUX_GRC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DWT_CONV_CFLUX_GRC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DWT_CONV_NFLUX_GRC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DWT_CONV_NFLUX_GRC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DWT_CONV_PFLUX_GRC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DWT_CONV_PFLUX_GRC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DWT_SLASH_CFLUX(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DWT_SLASH_CFLUX", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "slash C flux to litter and CWD due to land use")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DWT_SLASH_NFLUX(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DWT_SLASH_NFLUX", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "slash N flux to litter and CWD due to land use")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float DWT_SLASH_PFLUX(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("DWT_SLASH_PFLUX", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "slash P flux to litter and CWD due to land use")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float EFLX_DYNBAL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("EFLX_DYNBAL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "dynamic land cover change conversion energy flux")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float EFLX_GRND_LAKE(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("EFLX_GRND_LAKE", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "net heat flux into lake/snow surface, excluding light transmission")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float EFLX_LH_TOT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("EFLX_LH_TOT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total latent heat flux [+ to atm]")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float EFLX_LH_TOT_R(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("EFLX_LH_TOT_R", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Rural total evaporation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float EFLX_LH_TOT_U(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("EFLX_LH_TOT_U", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Urban total evaporation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float ELAI(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ELAI", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "exposed one-sided leaf area index")
    PUT_ATTR_TXT("units", "m^2/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float ER(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ER", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total ecosystem respiration, autotrophic + heterotrophic")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float ERRH2O(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ERRH2O", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total water conservation error")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float ERRH2OSNO(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ERRH2OSNO", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "imbalance in snow depth (liquid water)")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float ERRSEB(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ERRSEB", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "surface energy conservation error")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float ERRSOI(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ERRSOI", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil/lake energy conservation error")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float ERRSOL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ERRSOL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "solar radiation conservation error")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float ESAI(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ESAI", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "exposed one-sided stem area index")
    PUT_ATTR_TXT("units", "m^2/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FAREA_BURNED(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FAREA_BURNED", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "timestep fractional area burned")
    PUT_ATTR_TXT("units", "proportion")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FCEV(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FCEV", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "canopy evaporation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FCH4(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FCH4", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Gridcell surface CH4 flux to atmosphere (+ to atm)")
    PUT_ATTR_TXT("units", "kgC/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FCH4TOCO2(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FCH4TOCO2", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Gridcell oxidation of CH4 to CO2")
    PUT_ATTR_TXT("units", "gC/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FCH4_DFSAT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FCH4_DFSAT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "CH4 additional flux due to changing fsat, vegetated landunits only")
    PUT_ATTR_TXT("units", "kgC/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FCOV(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FCOV", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "fractional impermeable area")
    PUT_ATTR_TXT("units", "unitless")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FCTR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FCTR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "canopy transpiration")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FGEV(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FGEV", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "ground evaporation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FGR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FGR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "heat flux into soil/snow including snow melt and lake / snow light transmission")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FGR12(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FGR12", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "heat flux between soil layers 1 and 2")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FGR_R(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FGR_R", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Rural heat flux into soil/snow including snow melt and snow light transmission")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FGR_U(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FGR_U", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Urban heat flux into soil/snow including snow melt")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FH2OSFC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FH2OSFC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "fraction of ground covered by surface water")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FINUNDATED(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FINUNDATED", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "fractional inundated area of vegetated columns")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FINUNDATED_LAG(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FINUNDATED_LAG", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "time-lagged inundated fraction of vegetated columns")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FIRA(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FIRA", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "net infrared (longwave) radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FIRA_R(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FIRA_R", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Rural net infrared (longwave) radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FIRA_U(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FIRA_U", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Urban net infrared (longwave) radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FIRE(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FIRE", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "emitted infrared (longwave) radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FIRE_R(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FIRE_R", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Rural emitted infrared (longwave) radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FIRE_U(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FIRE_U", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Urban emitted infrared (longwave) radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FLDS(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FLDS", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "atmospheric longwave radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FPG(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FPG", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "fraction of potential gpp due to N limitation")
    PUT_ATTR_TXT("units", "proportion")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FPG_P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FPG_P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "fraction of potential gpp due to P limitation")
    PUT_ATTR_TXT("units", "proportion")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FPI(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FPI", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "fraction of potential immobilization of nitrogen")
    PUT_ATTR_TXT("units", "proportion")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FPI_P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FPI_P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "fraction of potential immobilization of phosphorus")
    PUT_ATTR_TXT("units", "proportion")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FPI_P_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("FPI_P_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "fraction of potential immobilization of phosphorus")
    PUT_ATTR_TXT("units", "proportion")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float FPI_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("FPI_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "fraction of potential immobilization of nitrogen")
    PUT_ATTR_TXT("units", "proportion")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float FPSN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FPSN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "photosynthesis")
    PUT_ATTR_TXT("units", "umol/m2s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FPSN_WC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FPSN_WC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Rubisco-limited photosynthesis")
    PUT_ATTR_TXT("units", "umol/m2s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FPSN_WJ(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FPSN_WJ", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "RuBP-limited photosynthesis")
    PUT_ATTR_TXT("units", "umol/m2s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FPSN_WP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FPSN_WP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Product-limited photosynthesis")
    PUT_ATTR_TXT("units", "umol/m2s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FROOTC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FROOTC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "fine root C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FROOTC_ALLOC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FROOTC_ALLOC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "fine root C allocation")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FROOTC_LOSS(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FROOTC_LOSS", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "fine root C loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FROOTN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FROOTN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "fine root N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FROOTP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FROOTP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "fine root P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FROST_TABLE(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FROST_TABLE", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "frost table depth (vegetated landunits only)")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSA(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSA", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "absorbed solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSAT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSAT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "fractional area with water table at surface")
    PUT_ATTR_TXT("units", "unitless")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSA_R(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSA_R", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Rural absorbed solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSA_U(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSA_U", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Urban absorbed solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSDS(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSDS", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "atmospheric incident solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSDSND(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSDSND", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "direct nir incident solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSDSNDLN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSDSNDLN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "direct nir incident solar radiation at local noon")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSDSNI(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSDSNI", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "diffuse nir incident solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSDSVD(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSDSVD", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "direct vis incident solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSDSVDLN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSDSVDLN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "direct vis incident solar radiation at local noon")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSDSVI(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSDSVI", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "diffuse vis incident solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSDSVILN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSDSVILN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "diffuse vis incident solar radiation at local noon")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSH(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSH", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "sensible heat")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSH_G(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSH_G", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "sensible heat from ground")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSH_NODYNLNDUSE(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSH_NODYNLNDUSE", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "sensible heat not including correction for land use change")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSH_R(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSH_R", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Rural sensible heat")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSH_U(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSH_U", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Urban sensible heat")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSH_V(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSH_V", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "sensible heat from veg")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSM(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSM", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "snow melt heat flux")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSM_R(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSM_R", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Rural snow melt heat flux")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSM_U(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSM_U", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Urban snow melt heat flux")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSNO(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSNO", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "fraction of ground covered by snow")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSNO_EFF(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSNO_EFF", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "effective fraction of ground covered by snow")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "reflected solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSRND(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSRND", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "direct nir reflected solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSRNDLN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSRNDLN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "direct nir reflected solar radiation at local noon")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSRNI(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSRNI", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "diffuse nir reflected solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSRVD(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSRVD", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "direct vis reflected solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSRVDLN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSRVDLN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "direct vis reflected solar radiation at local noon")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float FSRVI(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("FSRVI", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "diffuse vis reflected solar radiation")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float F_CO2_SOIL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("F_CO2_SOIL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total soil-atm. CO2 exchange")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float F_CO2_SOIL_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("F_CO2_SOIL_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "total vertically resolved soil-atm. CO2 exchange")
    PUT_ATTR_TXT("units", "gC/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float F_DENIT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("F_DENIT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "denitrification flux")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float F_DENIT_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("F_DENIT_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "denitrification flux")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float F_N2O_DENIT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("F_N2O_DENIT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "denitrification N2O flux")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float F_N2O_NIT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("F_N2O_NIT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "nitrification N2O flux")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float F_NIT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("F_NIT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "nitrification flux")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float F_NIT_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("F_NIT_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "nitrification flux")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float GC_HEAT1(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("GC_HEAT1", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "initial gridcell total heat content")
    PUT_ATTR_TXT("units", "J/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float GC_ICE1(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("GC_ICE1", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "initial gridcell total ice content")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float GC_LIQ1(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("GC_LIQ1", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "initial gridcell total liq content")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float GPP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("GPP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "gross primary production")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float GR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("GR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total growth respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float GROSS_NMIN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("GROSS_NMIN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "gross rate of N mineralization")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float GROSS_PMIN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("GROSS_PMIN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "gross rate of P mineralization")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float H2OCAN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("H2OCAN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "intercepted water")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float H2OSFC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("H2OSFC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "surface water depth")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float H2OSNO(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("H2OSNO", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "snow depth (liquid water)")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float H2OSNO_TOP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("H2OSNO_TOP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mass of snow in top snow layer")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float H2OSOI(time, levgrnd, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("H2OSOI", NC_FLOAT, ndims, rec_dimids[1])
    PUT_ATTR_TXT("long_name", "volumetric soil water (vegetated landunits only)")
    PUT_ATTR_TXT("units", "mm3/mm3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[1])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float HC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("HC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "heat content of soil/snow/lake")
    PUT_ATTR_TXT("units", "MJ/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float HCSOI(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("HCSOI", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil heat content")
    PUT_ATTR_TXT("units", "MJ/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float HEAT_FROM_AC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("HEAT_FROM_AC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "sensible heat flux put into canyon due to heat removed from air conditioning")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float HR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("HR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total heterotrophic respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float HR_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("HR_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "total vertically resolved heterotrophic respiration")
    PUT_ATTR_TXT("units", "gC/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float HTOP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("HTOP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "canopy top")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float INT_SNOW(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("INT_SNOW", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "accumulated swe (vegetated landunits only)")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LABILEP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LABILEP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil Labile P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LABILEP_TO_SECONDP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LABILEP_TO_SECONDP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "LABILE P TO SECONDARY MINERAL P")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LABILEP_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("LABILEP_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "soil labile P (vert. res.)")
    PUT_ATTR_TXT("units", "gp/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float LAISHA(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LAISHA", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "shaded projected leaf area index")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LAISUN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LAISUN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "sunlit projected leaf area index")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LAKEICEFRAC(time, levlak, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("LAKEICEFRAC", NC_FLOAT, ndims, rec_dimids[2])
    PUT_ATTR_TXT("long_name", "lake layer ice mass fraction")
    PUT_ATTR_TXT("units", "unitless")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(three, 4, orig_dimids[2])
    SET_VAR_META(REC_ITYPE, 2, 1, rec_buflen, decom.count[2])

    /* float LAKEICETHICK(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LAKEICETHICK", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "thickness of lake ice (including physical expansion on freezing)")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LAND_UPTAKE(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LAND_UPTAKE", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "NEE minus LAND_USE_FLUX, negative for update")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LAND_USE_FLUX(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LAND_USE_FLUX", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total C emitted from land cover conversion and wood product pools")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LEAFC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LEAFC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "leaf C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LEAFC_ALLOC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LEAFC_ALLOC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "leaf C allocation")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LEAFC_LOSS(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LEAFC_LOSS", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "leaf C loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LEAFC_TO_LITTER(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LEAFC_TO_LITTER", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "leaf C litterfall")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LEAFN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LEAFN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "leaf N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LEAFP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LEAFP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "leaf P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LEAF_MR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LEAF_MR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "leaf maintenance respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LFC2(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LFC2", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "conversion area fraction of BET and BDT that burned")
    PUT_ATTR_TXT("units", "per sec")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITFALL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITFALL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "litterfall (leaves and fine roots)")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITHR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITHR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "litter heterotrophic respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR1C(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR1C", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "LITR1 C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR1C_TO_SOIL1C(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR1C_TO_SOIL1C", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of litter 1 C to soil 1 C")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR1C_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("LITR1C_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "LITR1 C (vertically resolved)")
    PUT_ATTR_TXT("units", "gC/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float LITR1N(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR1N", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "LITR1 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR1N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("LITR1N_TNDNCY_VERT_TRANS", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "litter 1 N tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float LITR1N_TO_SOIL1N(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR1N_TO_SOIL1N", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of litter 1 N to soil 1 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR1N_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("LITR1N_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "LITR1 N (vertically resolved)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float LITR1P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR1P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "LITR1 P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR1P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("LITR1P_TNDNCY_VERT_TRANS", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "litter 1 P tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gP/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float LITR1P_TO_SOIL1P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR1P_TO_SOIL1P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of litter 1 P to soil 1 N")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR1P_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("LITR1P_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "LITR1 P (vertically resolved)")
    PUT_ATTR_TXT("units", "gP/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float LITR1_HR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR1_HR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Het. Resp. from litter 1")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR2C(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR2C", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "LITR2 C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR2C_TO_SOIL2C(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR2C_TO_SOIL2C", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of litter 2 C to soil 2 C")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR2C_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("LITR2C_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "LITR2 C (vertically resolved)")
    PUT_ATTR_TXT("units", "gC/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float LITR2N(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR2N", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "LITR2 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR2N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("LITR2N_TNDNCY_VERT_TRANS", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "litter 2 N tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float LITR2N_TO_SOIL2N(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR2N_TO_SOIL2N", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of litter 2 N to soil 2 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR2N_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("LITR2N_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "LITR2 N (vertically resolved)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float LITR2P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR2P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "LITR2 P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR2P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("LITR2P_TNDNCY_VERT_TRANS", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "litter 2 P tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gP/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float LITR2P_TO_SOIL2P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR2P_TO_SOIL2P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of litter 2 P to soil 2 N")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR2P_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("LITR2P_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "LITR2 P (vertically resolved)")
    PUT_ATTR_TXT("units", "gP/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float LITR2_HR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR2_HR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Het. Resp. from litter 2")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR3C(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR3C", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "LITR3 C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR3C_TO_SOIL3C(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR3C_TO_SOIL3C", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of litter 3 C to soil 3 C")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR3C_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("LITR3C_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "LITR3 C (vertically resolved)")
    PUT_ATTR_TXT("units", "gC/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float LITR3N(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR3N", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "LITR3 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR3N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("LITR3N_TNDNCY_VERT_TRANS", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "litter 3 N tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float LITR3N_TO_SOIL3N(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR3N_TO_SOIL3N", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of litter 3 N to soil 3 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR3N_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("LITR3N_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "LITR3 N (vertically resolved)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float LITR3P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR3P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "LITR3 P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR3P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("LITR3P_TNDNCY_VERT_TRANS", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "litter 3 P tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gP/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float LITR3P_TO_SOIL3P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR3P_TO_SOIL3P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of litter 3 P to soil 3 N")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITR3P_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("LITR3P_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "LITR3 P (vertically resolved)")
    PUT_ATTR_TXT("units", "gP/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float LITR3_HR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITR3_HR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Het. Resp. from litter 3")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITTERC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITTERC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "litter C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITTERC_HR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITTERC_HR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "litter C heterotrophic respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LITTERC_LOSS(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LITTERC_LOSS", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "litter C loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LIVECROOTC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LIVECROOTC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "live coarse root C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LIVECROOTN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LIVECROOTN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "live coarse root N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LIVECROOTP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LIVECROOTP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "live coarse root P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LIVESTEMC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LIVESTEMC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "live stem C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LIVESTEMN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LIVESTEMN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "live stem N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float LIVESTEMP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("LIVESTEMP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "live stem P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float MR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("MR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "maintenance respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float M_LITR1C_TO_LEACHING(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("M_LITR1C_TO_LEACHING", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "litter 1 C leaching loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float M_LITR2C_TO_LEACHING(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("M_LITR2C_TO_LEACHING", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "litter 2 C leaching loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float M_LITR3C_TO_LEACHING(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("M_LITR3C_TO_LEACHING", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "litter 3 C leaching loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float M_SOIL1C_TO_LEACHING(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("M_SOIL1C_TO_LEACHING", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil 1 C leaching loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float M_SOIL2C_TO_LEACHING(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("M_SOIL2C_TO_LEACHING", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil 2 C leaching loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float M_SOIL3C_TO_LEACHING(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("M_SOIL3C_TO_LEACHING", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil 3 C leaching loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float M_SOIL4C_TO_LEACHING(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("M_SOIL4C_TO_LEACHING", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil 4 C leaching loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float NBP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("NBP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "net biome production, includes fire, landuse, and harvest flux, positive for sink")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float NDEPLOY(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("NDEPLOY", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total N deployed in new growth")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float NDEP_TO_SMINN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("NDEP_TO_SMINN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "atmospheric N deposition to soil mineral N")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float NEE(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("NEE", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "net ecosystem exchange of carbon, includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float NEM(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("NEM", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Gridcell net adjustment to NEE passed to atm. for methane production")
    PUT_ATTR_TXT("units", "gC/m2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float NEP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("NEP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "net ecosystem production, excludes fire, landuse, and harvest flux, positive for sink")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float NET_NMIN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("NET_NMIN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "net rate of N mineralization")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float NET_PMIN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("NET_PMIN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "net rate of P mineralization")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float NFIRE(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("NFIRE", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "fire counts valid only in Reg.C")
    PUT_ATTR_TXT("units", "counts/km2/sec")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float NFIX_TO_SMINN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("NFIX_TO_SMINN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "symbiotic/asymbiotic N fixation to soil mineral N")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float NPP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("NPP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "net primary production")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float OCCLP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("OCCLP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil occluded P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float OCCLP_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("OCCLP_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "soil occluded P (vert. res.)")
    PUT_ATTR_TXT("units", "gp/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float OCDEP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("OCDEP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total OC deposition (dry+wet) from atmosphere")
    PUT_ATTR_TXT("units", "kg/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float O_SCALAR(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("O_SCALAR", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "fraction by which decomposition is reduced due to anoxia")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float PARVEGLN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PARVEGLN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "absorbed par by vegetation at local noon")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PBOT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PBOT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "atmospheric pressure")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PCH4(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PCH4", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "atmospheric partial pressure of CH4")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PCO2(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PCO2", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "atmospheric partial pressure of CO2")
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PCT_LANDUNIT(time, ltype, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("PCT_LANDUNIT", NC_FLOAT, ndims, rec_dimids[3])
    PUT_ATTR_TXT("long_name", "% of each landunit on grid cell")
    PUT_ATTR_TXT("units", "%")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(four, 4, orig_dimids[3])
    SET_VAR_META(REC_ITYPE, 3, 1, rec_buflen, decom.count[3])

    /* float PCT_NAT_PFT(time, natpft, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("PCT_NAT_PFT", NC_FLOAT, ndims, rec_dimids[4])
    PUT_ATTR_TXT("long_name", "% of each PFT on the natural vegetation (i.e., soil) landunit")
    PUT_ATTR_TXT("units", "%")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(five, 4, orig_dimids[4])
    SET_VAR_META(REC_ITYPE, 4, 1, rec_buflen, decom.count[4])

    /* float PDEPLOY(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PDEPLOY", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total P deployed in new growth")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PDEP_TO_SMINP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PDEP_TO_SMINP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "atmospheric P deposition to soil mineral P")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PFT_FIRE_CLOSS(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PFT_FIRE_CLOSS", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total patch-level fire C loss for non-peat fires outside land-type converted region")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PFT_FIRE_NLOSS(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PFT_FIRE_NLOSS", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total pft-level fire N loss")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PLANT_CALLOC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PLANT_CALLOC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total allocated C flux")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PLANT_NDEMAND(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PLANT_NDEMAND", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "N flux required to support initial GPP")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PLANT_NDEMAND_COL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PLANT_NDEMAND_COL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "N flux required to support initial GPP")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PLANT_PALLOC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PLANT_PALLOC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total allocated P flux")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PLANT_PDEMAND(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PLANT_PDEMAND", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "P flux required to support initial GPP")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PLANT_PDEMAND_COL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PLANT_PDEMAND_COL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "P flux required to support initial GPP")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float POTENTIAL_IMMOB(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("POTENTIAL_IMMOB", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "potential N immobilization")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float POTENTIAL_IMMOB_P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("POTENTIAL_IMMOB_P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "potential P immobilization")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float POT_F_DENIT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("POT_F_DENIT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "potential denitrification flux")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float POT_F_NIT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("POT_F_NIT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "potential nitrification flux")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PRIMP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PRIMP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil primary P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PRIMP_TO_LABILEP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PRIMP_TO_LABILEP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "PRIMARY MINERAL P TO LABILE P")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PRIMP_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("PRIMP_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "soil primary P (vert. res.)")
    PUT_ATTR_TXT("units", "gp/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float PROD1P_LOSS(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PROD1P_LOSS", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "loss from 1-yr crop product pool")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PSNSHA(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PSNSHA", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "shaded leaf photosynthesis")
    PUT_ATTR_TXT("units", "umolCO2/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PSNSHADE_TO_CPOOL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PSNSHADE_TO_CPOOL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "C fixation from shaded canopy")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PSNSUN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PSNSUN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "sunlit leaf photosynthesis")
    PUT_ATTR_TXT("units", "umolCO2/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float PSNSUN_TO_CPOOL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("PSNSUN_TO_CPOOL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "C fixation from sunlit canopy")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float Q2M(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("Q2M", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "2m specific humidity")
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QBOT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QBOT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "atmospheric specific humidity")
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QCHARGE(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QCHARGE", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "aquifer recharge rate (vegetated landunits only)")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QDRAI(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QDRAI", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "sub-surface drainage")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QDRAI_PERCH(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QDRAI_PERCH", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "perched wt drainage")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QDRAI_XS(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QDRAI_XS", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "saturation excess drainage")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QDRIP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QDRIP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "throughfall")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QFLOOD(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QFLOOD", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "runoff from river flooding")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QFLX_ICE_DYNBAL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QFLX_ICE_DYNBAL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "ice dynamic land cover change conversion runoff flux")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QFLX_LIQ_DYNBAL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QFLX_LIQ_DYNBAL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "liq dynamic land cover change conversion runoff flux")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QH2OSFC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QH2OSFC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "surface water runoff")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QINFL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QINFL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "infiltration")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QINTR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QINTR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "interception")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QIRRIG_GRND(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QIRRIG_GRND", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Groundwater irrigation")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QIRRIG_ORIG(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QIRRIG_ORIG", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Original total irrigation water demand (surface + ground)")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QIRRIG_REAL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QIRRIG_REAL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "actual water added through irrigation (surface + ground)")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QIRRIG_SURF(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QIRRIG_SURF", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Surface water irrigation")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QIRRIG_WM(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QIRRIG_WM", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Surface water irrigation demand sent to MOSART/WM")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QOVER(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QOVER", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "surface runoff")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QOVER_LAG(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QOVER_LAG", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "time-lagged surface runoff for soil columns")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QRGWL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QRGWL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "surface runoff at glaciers (liquid only), wetlands, lakes")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QRUNOFF(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QRUNOFF", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total liquid runoff (does not include QSNWCPICE)")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QRUNOFF_NODYNLNDUSE(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QRUNOFF_NODYNLNDUSE", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total liquid runoff (does not include QSNWCPICE) not including correction for land use change")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QRUNOFF_R(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QRUNOFF_R", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Rural total runoff")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QRUNOFF_U(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QRUNOFF_U", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Urban total runoff")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QSNOMELT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QSNOMELT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "snow melt")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QSNWCPICE(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QSNWCPICE", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "excess snowfall due to snow capping")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QSNWCPICE_NODYNLNDUSE(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QSNWCPICE_NODYNLNDUSE", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "excess snowfall due to snow capping not including correction for land use change")
    PUT_ATTR_TXT("units", "mm H2O/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QSOIL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QSOIL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Ground evaporation (soil/snow evaporation + soil/snow sublimation - dew)")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QVEGE(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QVEGE", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "canopy evaporation")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float QVEGT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("QVEGT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "canopy transpiration")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float RAIN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("RAIN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "atmospheric rain")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float RETRANSN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("RETRANSN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "plant pool of retranslocated N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float RETRANSN_TO_NPOOL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("RETRANSN_TO_NPOOL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "deployment of retranslocated N")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float RETRANSP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("RETRANSP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "plant pool of retranslocated P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float RETRANSP_TO_PPOOL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("RETRANSP_TO_PPOOL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "deployment of retranslocated P")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float RH2M(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("RH2M", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "2m relative humidity")
    PUT_ATTR_TXT("units", "%")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float RH2M_R(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("RH2M_R", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Rural 2m specific humidity")
    PUT_ATTR_TXT("units", "%")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float RH2M_U(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("RH2M_U", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Urban 2m relative humidity")
    PUT_ATTR_TXT("units", "%")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float RR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("RR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "root respiration (fine root MR + total root GR)")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SABG(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SABG", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "solar rad absorbed by ground")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SABG_PEN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SABG_PEN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Rural solar rad penetrating top soil or snow layer")
    PUT_ATTR_TXT("units", "watt/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SABV(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SABV", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "solar rad absorbed by veg")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SCALARAVG_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SCALARAVG_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "average of decomposition scalar")
    PUT_ATTR_TXT("units", "fraction")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SECONDP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SECONDP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil secondary P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SECONDP_TO_LABILEP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SECONDP_TO_LABILEP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "SECONDARY MINERAL P TO LABILE P")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SECONDP_TO_OCCLP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SECONDP_TO_OCCLP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "SECONDARY MINERAL P TO OCCLUDED P")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SECONDP_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SECONDP_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "soil secondary P (vert. res.)")
    PUT_ATTR_TXT("units", "gp/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SEEDC_GRC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SEEDC_GRC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "pool for seeding new PFTs via dynamic landcover")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil mineral N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINN_TO_NPOOL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINN_TO_NPOOL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "deployment of soil mineral N uptake")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINN_TO_PLANT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINN_TO_PLANT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "plant uptake of soil mineral N")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINN_TO_SOIL1N_L1(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINN_TO_SOIL1N_L1", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mineral N flux for decomp. of LITR1to SOIL1")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINN_TO_SOIL2N_L2(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINN_TO_SOIL2N_L2", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mineral N flux for decomp. of LITR2to SOIL2")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINN_TO_SOIL2N_S1(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINN_TO_SOIL2N_S1", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mineral N flux for decomp. of SOIL1to SOIL2")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINN_TO_SOIL3N_L3(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINN_TO_SOIL3N_L3", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mineral N flux for decomp. of LITR3to SOIL3")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINN_TO_SOIL3N_S2(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINN_TO_SOIL3N_S2", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mineral N flux for decomp. of SOIL2to SOIL3")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINN_TO_SOIL4N_S3(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINN_TO_SOIL4N_S3", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mineral N flux for decomp. of SOIL3to SOIL4")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil mineral P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINP_LEACHED(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINP_LEACHED", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil mineral P pool loss to leaching")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINP_TO_PLANT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINP_TO_PLANT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "plant uptake of soil mineral P")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINP_TO_PPOOL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINP_TO_PPOOL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "deployment of soil mineral P uptake")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINP_TO_SOIL1P_L1(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINP_TO_SOIL1P_L1", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mineral P flux for decomp. of LITR1to SOIL1")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINP_TO_SOIL2P_L2(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINP_TO_SOIL2P_L2", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mineral P flux for decomp. of LITR2to SOIL2")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINP_TO_SOIL2P_S1(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINP_TO_SOIL2P_S1", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mineral P flux for decomp. of SOIL1to SOIL2")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINP_TO_SOIL3P_L3(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINP_TO_SOIL3P_L3", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mineral P flux for decomp. of LITR3to SOIL3")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINP_TO_SOIL3P_S2(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINP_TO_SOIL3P_S2", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mineral P flux for decomp. of SOIL2to SOIL3")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINP_TO_SOIL4P_S3(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMINP_TO_SOIL4P_S3", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mineral P flux for decomp. of SOIL3to SOIL4")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMINP_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SMINP_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "soil mineral P (vert. res.)")
    PUT_ATTR_TXT("units", "gp/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SMIN_NH4(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMIN_NH4", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil mineral NH4")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMIN_NH4_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SMIN_NH4_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "soil mineral NH4 (vert. res.)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SMIN_NO3(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMIN_NO3", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil mineral NO3")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMIN_NO3_LEACHED(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMIN_NO3_LEACHED", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil NO3 pool loss to leaching")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMIN_NO3_RUNOFF(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SMIN_NO3_RUNOFF", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil NO3 pool loss to runoff")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SMIN_NO3_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SMIN_NO3_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "soil mineral NO3 (vert. res.)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SNOBCMCL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SNOBCMCL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mass of BC in snow column")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SNOBCMSL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SNOBCMSL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mass of BC in top snow layer")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SNODSTMCL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SNODSTMCL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mass of dust in snow column")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SNODSTMSL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SNODSTMSL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mass of dust in top snow layer")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SNOINTABS(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SNOINTABS", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Percent of incoming solar absorbed by lower snow layers")
    PUT_ATTR_TXT("units", "%")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SNOOCMCL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SNOOCMCL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mass of OC in snow column")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SNOOCMSL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SNOOCMSL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mass of OC in top snow layer")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SNOW(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SNOW", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "atmospheric snow")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SNOWDP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SNOWDP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "gridcell mean snow height")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SNOWICE(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SNOWICE", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "snow ice")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SNOWLIQ(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SNOWLIQ", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "snow liquid water")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SNOW_DEPTH(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SNOW_DEPTH", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "snow height of snow covered area")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SNOW_SINKS(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SNOW_SINKS", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "snow sinks (liquid water)")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SNOW_SOURCES(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SNOW_SOURCES", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "snow sources (liquid water)")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL1C(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL1C", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "SOIL1 C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL1C_TO_SOIL2C(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL1C_TO_SOIL2C", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of soil 1 C to soil 2 C")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL1C_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL1C_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "SOIL1 C (vertically resolved)")
    PUT_ATTR_TXT("units", "gC/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL1N(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL1N", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "SOIL1 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL1N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL1N_TNDNCY_VERT_TRANS", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "soil 1 N tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL1N_TO_SOIL2N(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL1N_TO_SOIL2N", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of soil 1 N to soil 2 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL1N_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL1N_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "SOIL1 N (vertically resolved)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL1P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL1P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "SOIL1 P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL1P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL1P_TNDNCY_VERT_TRANS", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "soil 1 P tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gP/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL1P_TO_SOIL2P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL1P_TO_SOIL2P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of soil 1 P to soil 2 N")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL1P_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL1P_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "SOIL1 P (vertically resolved)")
    PUT_ATTR_TXT("units", "gP/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL1_HR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL1_HR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Het. Resp. from soil 1")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL2C(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL2C", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "SOIL2 C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL2C_TO_SOIL3C(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL2C_TO_SOIL3C", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of soil 2 C to soil 3 C")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL2C_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL2C_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "SOIL2 C (vertically resolved)")
    PUT_ATTR_TXT("units", "gC/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL2N(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL2N", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "SOIL2 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL2N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL2N_TNDNCY_VERT_TRANS", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "soil 2 N tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL2N_TO_SOIL3N(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL2N_TO_SOIL3N", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of soil 2 N to soil 3 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL2N_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL2N_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "SOIL2 N (vertically resolved)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL2P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL2P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "SOIL2 P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL2P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL2P_TNDNCY_VERT_TRANS", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "soil 2 P tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gP/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL2P_TO_SOIL3P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL2P_TO_SOIL3P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of soil 2 P to soil 3 N")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL2P_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL2P_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "SOIL2 P (vertically resolved)")
    PUT_ATTR_TXT("units", "gP/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL2_HR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL2_HR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Het. Resp. from soil 2")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL3C(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL3C", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "SOIL3 C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL3C_TO_SOIL4C(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL3C_TO_SOIL4C", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of soil 3 C to soil 4 C")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL3C_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL3C_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "SOIL3 C (vertically resolved)")
    PUT_ATTR_TXT("units", "gC/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL3N(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL3N", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "SOIL3 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL3N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL3N_TNDNCY_VERT_TRANS", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "soil 3 N tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL3N_TO_SOIL4N(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL3N_TO_SOIL4N", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of soil 3 N to soil 4 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL3N_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL3N_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "SOIL3 N (vertically resolved)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL3P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL3P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "SOIL3 P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL3P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL3P_TNDNCY_VERT_TRANS", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "soil 3 P tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gP/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL3P_TO_SOIL4P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL3P_TO_SOIL4P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "decomp. of soil 3 P to soil 4 N")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL3P_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL3P_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "SOIL3 P (vertically resolved)")
    PUT_ATTR_TXT("units", "gP/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL3_HR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL3_HR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Het. Resp. from soil 3")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL4C(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL4C", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "SOIL4 C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL4C_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL4C_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "SOIL4 C (vertically resolved)")
    PUT_ATTR_TXT("units", "gC/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL4N(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL4N", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "SOIL4 N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL4N_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL4N_TNDNCY_VERT_TRANS", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "soil 4 N tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gN/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL4N_TO_SMINN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL4N_TO_SMINN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mineral N flux for decomp. of SOIL4")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL4N_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL4N_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "SOIL4 N (vertically resolved)")
    PUT_ATTR_TXT("units", "gN/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL4P(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL4P", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "SOIL4 P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL4P_TNDNCY_VERT_TRANS(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL4P_TNDNCY_VERT_TRANS", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "soil 4 P tendency due to vertical transport")
    PUT_ATTR_TXT("units", "gP/m^3/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL4P_TO_SMINP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL4P_TO_SMINP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "mineral P flux for decomp. of SOIL4")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOIL4P_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOIL4P_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "SOIL4 P (vertically resolved)")
    PUT_ATTR_TXT("units", "gP/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOIL4_HR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOIL4_HR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Het. Resp. from soil 4")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOILC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOILC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOILC_HR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOILC_HR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil C heterotrophic respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOILC_LOSS(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOILC_LOSS", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil C loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOILICE(time, levgrnd, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOILICE", NC_FLOAT, ndims, rec_dimids[1])
    PUT_ATTR_TXT("long_name", "soil ice (vegetated landunits only)")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[1])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOILICE_ICE(time, levgrnd, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOILICE_ICE", NC_FLOAT, ndims, rec_dimids[1])
    PUT_ATTR_TXT("long_name", "soil ice (ice landunits only)")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[1])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOILLIQ(time, levgrnd, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOILLIQ", NC_FLOAT, ndims, rec_dimids[1])
    PUT_ATTR_TXT("long_name", "soil liquid water (vegetated landunits only)")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[1])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOILLIQ_ICE(time, levgrnd, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOILLIQ_ICE", NC_FLOAT, ndims, rec_dimids[1])
    PUT_ATTR_TXT("long_name", "soil liquid water (ice landunits only)")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[1])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOILPSI(time, levgrnd, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOILPSI", NC_FLOAT, ndims, rec_dimids[1])
    PUT_ATTR_TXT("long_name", "soil water potential in each soil layer")
    PUT_ATTR_TXT("units", "MPa")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[1])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOILWATER_10CM(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOILWATER_10CM", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil liquid water + ice in top 10cm of soil (veg landunits only)")
    PUT_ATTR_TXT("standard_name", "mass_content_of_water_in_soil_layer")
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOLUTIONP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOLUTIONP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil solution P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOLUTIONP_vr(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("SOLUTIONP_vr", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "soil solution P (vert. res.)")
    PUT_ATTR_TXT("units", "gp/m^3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float SOMHR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOMHR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil organic matter heterotrophic respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SOM_C_LEACHED(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SOM_C_LEACHED", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total flux of C from SOM pools due to leaching")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total soil respiration (HR + root resp)")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float STORVEGC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("STORVEGC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "stored vegetation carbon, excluding cpool")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float STORVEGN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("STORVEGN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "stored vegetation nitrogen")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float STORVEGP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("STORVEGP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "stored vegetation phosphorus")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SUPPLEMENT_TO_SMINN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SUPPLEMENT_TO_SMINN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "supplemental N supply")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SUPPLEMENT_TO_SMINP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SUPPLEMENT_TO_SMINP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "supplemental P supply")
    PUT_ATTR_TXT("units", "gP/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SUPPLY(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SUPPLY", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "runoff supply for land use")
    PUT_ATTR_TXT("units", "mm/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SoilAlpha(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SoilAlpha", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "factor limiting ground evap")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float SoilAlpha_U(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("SoilAlpha_U", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "urban factor limiting ground evap")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TAUX(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TAUX", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "zonal surface stress")
    PUT_ATTR_TXT("units", "kg/m/s^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TAUY(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TAUY", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "meridional surface stress")
    PUT_ATTR_TXT("units", "kg/m/s^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TBOT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TBOT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "atmospheric air temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TBUILD(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TBUILD", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "internal urban building temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TCS_MONTH_BEGIN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TCS_MONTH_BEGIN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total carbon storage at the beginning of a month")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TCS_MONTH_END(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TCS_MONTH_END", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total carbon storage at the end of a month")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TG(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TG", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "ground temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TG_R(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TG_R", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Rural ground temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TG_U(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TG_U", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Urban ground temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TH2OSFC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TH2OSFC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "surface water temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float THBOT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("THBOT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "atmospheric air potential temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TKE1(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TKE1", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "top lake level eddy thermal conductivity")
    PUT_ATTR_TXT("units", "W/(mK)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TLAI(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TLAI", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total projected leaf area index")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TLAKE(time, levlak, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("TLAKE", NC_FLOAT, ndims, rec_dimids[2])
    PUT_ATTR_TXT("long_name", "lake temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(three, 4, orig_dimids[2])
    SET_VAR_META(REC_ITYPE, 2, 1, rec_buflen, decom.count[2])

    /* float TOTCOLC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTCOLC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total column carbon, incl veg and cpool but excl product pools")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTCOLCH4(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTCOLCH4", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total belowground CH4, (0 for non-lake special landunits)")
    PUT_ATTR_TXT("units", "gC/m2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTCOLN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTCOLN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total column-level N but excl product pools")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTCOLP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTCOLP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total column-level P but excl product pools")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTECOSYSC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTECOSYSC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total ecosystem carbon, incl veg but excl cpool but excl product pools")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTECOSYSN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTECOSYSN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total ecosystem N but excl product pools")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTECOSYSP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTECOSYSP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total ecosystem P but excl product pools")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTLITC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTLITC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total litter carbon")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTLITC_1m(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTLITC_1m", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total litter carbon to 1 meter depth")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTLITN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTLITN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total litter N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTLITP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTLITP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total litter P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTLITP_1m(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTLITP_1m", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total litter P to 1 meter")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTPFTC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTPFTC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total patch-level carbon, including cpool")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTPFTN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTPFTN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total PFT-level nitrogen")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTPFTP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTPFTP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total PFT-level phosphorus")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTSOMC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTSOMC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total soil organic matter carbon")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTSOMC_1m(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTSOMC_1m", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total soil organic matter carbon to 1 meter depth")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTSOMN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTSOMN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total soil organic matter N")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTSOMP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTSOMP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total soil organic matter P")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTSOMP_1m(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTSOMP_1m", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total soil organic matter P to 1 meter")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTVEGC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTVEGC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total vegetation carbon, excluding cpool")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTVEGC_ABG(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTVEGC_ABG", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total aboveground vegetation carbon, excluding cpool")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTVEGN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTVEGN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total vegetation nitrogen")
    PUT_ATTR_TXT("units", "gN/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TOTVEGP(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TOTVEGP", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total vegetation phosphorus")
    PUT_ATTR_TXT("units", "gP/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TREFMNAV(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TREFMNAV", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "daily minimum of average 2-m temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TREFMNAV_R(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TREFMNAV_R", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Rural daily minimum of average 2-m temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TREFMNAV_U(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TREFMNAV_U", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Urban daily minimum of average 2-m temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TREFMXAV(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TREFMXAV", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "daily maximum of average 2-m temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TREFMXAV_R(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TREFMXAV_R", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Rural daily maximum of average 2-m temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TREFMXAV_U(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TREFMXAV_U", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Urban daily maximum of average 2-m temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TSA(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TSA", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "2m air temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TSAI(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TSAI", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total projected stem area index")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TSA_R(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TSA_R", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Rural 2m air temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TSA_U(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TSA_U", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "Urban 2m air temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TSOI(time, levgrnd, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("TSOI", NC_FLOAT, ndims, rec_dimids[1])
    PUT_ATTR_TXT("long_name", "soil temperature (vegetated landunits only)")
    PUT_ATTR_TXT("standard_name", "soil_temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[1])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float TSOI_10CM(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TSOI_10CM", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil temperature in top 10cm of soil")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TSOI_ICE(time, levgrnd, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("TSOI_ICE", NC_FLOAT, ndims, rec_dimids[1])
    PUT_ATTR_TXT("long_name", "soil temperature (ice landunits only)")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[1])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float TV(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TV", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "vegetation temperature")
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TWS(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TWS", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total water storage")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TWS_MONTH_BEGIN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TWS_MONTH_BEGIN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total water storage at the beginning of a month")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float TWS_MONTH_END(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("TWS_MONTH_END", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total water storage at the end of a month")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float T_SCALAR(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("T_SCALAR", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "temperature inhibition of decomposition")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float U10(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("U10", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "10-m wind")
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float URBAN_AC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("URBAN_AC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "urban air conditioning flux")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float URBAN_HEAT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("URBAN_HEAT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "urban heating flux")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float VOLR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("VOLR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "river channel total water storage")
    PUT_ATTR_TXT("units", "m3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float VOLRMCH(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("VOLRMCH", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "river channel main channel water storage")
    PUT_ATTR_TXT("units", "m3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float WA(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("WA", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "water in the unconfined aquifer (vegetated landunits only)")
    PUT_ATTR_TXT("units", "mm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float WASTEHEAT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("WASTEHEAT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "sensible heat flux from heating/cooling sources of urban waste heat")
    PUT_ATTR_TXT("units", "W/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float WF(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("WF", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "soil water as frac. of whc for top 0.05 m")
    PUT_ATTR_TXT("units", "proportion")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float WIND(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("WIND", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "atmospheric wind velocity magnitude")
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float WOODC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("WOODC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "wood C")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float WOODC_ALLOC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("WOODC_ALLOC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "wood C eallocation")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float WOODC_LOSS(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("WOODC_LOSS", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "wood C loss")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float WOOD_HARVESTC(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("WOOD_HARVESTC", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "wood harvest carbon (to product pools)")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float WOOD_HARVESTN(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("WOOD_HARVESTN", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "wood harvest N (to product pools)")
    PUT_ATTR_TXT("units", "gN/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float WTGQ(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("WTGQ", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "surface tracer conductance")
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float W_SCALAR(time, levdcmp, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("W_SCALAR", NC_FLOAT, ndims, rec_dimids[5])
    PUT_ATTR_TXT("long_name", "Moisture (dryness) inhibition of decomposition")
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[5])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float XR(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("XR", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "total excess respiration")
    PUT_ATTR_TXT("units", "gC/m^2/s")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float XSMRPOOL(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("XSMRPOOL", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "temporary photosynthate C pool")
    PUT_ATTR_TXT("units", "gC/m^2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float ZBOT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ZBOT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "atmospheric reference height")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float ZWT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ZWT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "water table depth (vegetated landunits only)")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float ZWT_CH4_UNSAT(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ZWT_CH4_UNSAT", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "depth of water table for methane production used in non-inundated area")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float ZWT_PERCH(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("ZWT_PERCH", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "perched water table depth (vegetated landunits only)")
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float cn_scalar(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("cn_scalar", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "N limitation factor")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float cp_scalar(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("cp_scalar", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "P limitation factor")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float leaf_npimbalance(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("leaf_npimbalance", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "leaf np imbalance partial C partial P/partial C partial N")
    PUT_ATTR_TXT("units", "gN/gP")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float nlim_m(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("nlim_m", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "runmean N limitation factor")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    /* float o2_decomp_depth_unsat(time, levgrnd, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 4;
    DEF_VAR("o2_decomp_depth_unsat", NC_FLOAT, ndims, rec_dimids[1])
    PUT_ATTR_TXT("long_name", "o2_decomp_depth_unsat")
    PUT_ATTR_TXT("units", "mol/m3/2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(two, 4, orig_dimids[1])
    SET_VAR_META(REC_ITYPE, 1, 1, rec_buflen, decom.count[1])

    /* float plim_m(time, lat, lon) */
    ndims = (false && (cfg.strategy == blob)) ? 2 : 3;
    DEF_VAR("plim_m", NC_FLOAT, ndims, rec_dimids[0])
    PUT_ATTR_TXT("long_name", "runmean P limitation factor")
    PUT_ATTR_TXT("units", "")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_FLOAT(_FillValue, 1, 1.e+36f)
    PUT_ATTR_FLOAT("missing_value", 1, 1.e+36f)
    PUT_ATTR_DECOMP(one, 3, orig_dimids[0])
    SET_VAR_META(REC_ITYPE, 0, 1, rec_buflen, decom.count[0])

    if (cfg.strategy == blob && (cfg.api == pnetcdf || cfg.api == hdf5))
        assert(varid == cfg.nvars + NVARS_DECOMP*decom.num_decomp);
    else
        assert(varid == cfg.nvars);

err_out:
    return err;
}

