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

#include <string.h>  /* strlen() */
#include <assert.h>
#include <mpi.h>

#include <e3sm_io.h>
#include <e3sm_io_driver.hpp>
#include <e3sm_io_err.h>

#define CHECK_VAR_ERR(varid) {                                \
    if (err != 0) {                                           \
        char var_name[64];                                    \
        driver.inq_var_name(ncid, varid, var_name);           \
        printf("Error in %s:%d: %s() var %s\n",               \
               __FILE__, __LINE__, __func__, var_name);       \
        goto err_out;                                         \
    }                                                         \
}
#define DEF_DIM(name, num, dimid) {                                      \
    err = driver.def_dim(ncid, name, num, dimid);                        \
    CHECK_ERR                                                            \
}
#define DEF_VAR(name, type, ndims, dimids) {                             \
    varid++;                                                             \
    err = driver.def_var(ncid, name, type, ndims, dimids, varid);        \
    if (err != 0) {                                                      \
        printf("Error in %s line %d: def_var %s\n", __FILE__, __LINE__,  \
               name);                                                    \
        goto err_out;                                                    \
    }                                                                    \
}
#define PUT_GATTR_TXT(name, buf) {                                       \
    err = driver.put_att(ncid, NC_GLOBAL, name, NC_CHAR, strlen(buf), buf); \
    CHECK_ERR                                                            \
}
#define PUT_GATTR_INT(name, num, val) {                                  \
    int buf = val;                                                       \
    err = driver.put_att(ncid, NC_GLOBAL, name, NC_INT, num, &buf); \
    CHECK_ERR                                                            \
}
#define PUT_ATTR_TXT(name, buf) {                                        \
    err = driver.put_att(ncid, *varid, name, NC_CHAR, strlen(buf), buf);\
    CHECK_VAR_ERR(*varid)                                                \
}
#define PUT_ATTR_INT(name, num, buf) {                                   \
    err = driver.put_att(ncid, *varid, name, NC_INT, num, buf);         \
    CHECK_VAR_ERR(*varid)                                                \
}
#define PUT_ATTR_FLOAT(name, num, buf) {                                 \
    err = driver.put_att(ncid, *varid, name, NC_FLOAT, num, buf);       \
    CHECK_VAR_ERR(*varid)                                                \
}
#define PUT_ATTR_INT64(name, num, buf) {                                 \
    err = driver.put_att(ncid, *varid, name, NC_INT64, num, buf);   \
    CHECK_VAR_ERR(*varid)                                                \
}
#define PUT_ATTR_DECOMP(D, ndims, dimids) {                                    \
    if (cfg.strategy == blob && (cfg.api == pnetcdf || cfg.api == hdf5)) {     \
        err = driver.put_att(ncid,*varid,"decomposition_ID",NC_INT,1,&D);     \
        CHECK_VAR_ERR(*varid)                                                  \
        err = driver.put_att(ncid,*varid,"global_dimids",NC_INT,ndims,dimids);\
        CHECK_VAR_ERR(*varid)                                                  \
    }                                                                          \
}

/*----< add_gattrs() >-------------------------------------------------------*/
static
int add_gattrs(e3sm_io_config &cfg,
               e3sm_io_decom  &decom,
               e3sm_io_driver &driver,
               int             ncid)
{
    int err, nprocs;

    /* save number of processes as global attributes */
    if (cfg.strategy == blob && (cfg.api == pnetcdf || cfg.api == hdf5)) {
        MPI_Comm_size(cfg.io_comm, &nprocs);
        PUT_GATTR_INT("global_nprocs", 1, nprocs)
        PUT_GATTR_INT("num_decompositions", 1, decom.num_decomp)
        PUT_GATTR_INT("num_subfiles", 1, cfg.num_subfiles)
    }

    PUT_GATTR_INT("ne", 1, 120)
    PUT_GATTR_INT("np", 1, 21600)
    PUT_GATTR_TXT("title", "EAM History file information")
    PUT_GATTR_TXT("source", "E3SM Atmosphere Model")
    PUT_GATTR_TXT("source_id", "025f820fce")
    PUT_GATTR_TXT("product", "model-output")
    PUT_GATTR_TXT("realm", "atmos")
    PUT_GATTR_TXT("case", "FC5AV1C-H01A_ne120_oRRS18v3")
    PUT_GATTR_TXT("username", "E3SM")
    PUT_GATTR_TXT("hostname", "cori-knl")
    PUT_GATTR_TXT("git_version", "025f820fce")
    PUT_GATTR_TXT("history", "created on 06/10/21 11:59:48")
    PUT_GATTR_TXT("Conventions", "CF-1.7")
    PUT_GATTR_TXT("institution_id", "E3SM-Project")
    PUT_GATTR_TXT("institution", "LLNL (Lawrence Livermore National Laboratory, Livermore, CA 94550, USA); ANL (Argonne National Laboratory, Argonne, IL 60439, USA); BNL (Brookhaven National Laboratory, Upton, NY 11973, USA); LANL (Los Alamos National Laboratory, Los Alamos, NM 87545, USA); LBNL (Lawrence Berkeley National Laboratory, Berkeley, CA 94720, USA); ORNL (Oak Ridge National Laboratory, Oak Ridge, TN 37831, USA); PNNL (Pacific Northwest National Laboratory, Richland, WA 99352, USA); SNL (Sandia National Laboratories, Albuquerque,")
    PUT_GATTR_TXT("contact", "e3sm-data-support@listserv.llnl.gov")
    PUT_GATTR_TXT("initial_file", "/global/cfs/cdirs/e3sm/inputdata/atm/cam/inic/homme/cami_mam3_Linoz_0000-01-ne120np4_L72_c160318.nc")
    PUT_GATTR_TXT("topography_file", "/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_ne120np4_16xdel2-PFC-consistentSGH.nc")

err_out:
    return err;
}

/*----< def_dims() >---------------------------------------------------------*/
static
int def_dims(e3sm_io_config &cfg,
             e3sm_io_decom  &decom,
             e3sm_io_driver &driver,
             int             ncid,
             int            *dim_ncol,
             int            *dim_time,
             int            *dim_nbnd,
             int            *dim_chars,
             int            *dim_lev,
             int            *dim_ilev,
             int            *nblobs_ID,
             int             nelems_D[3],
             int             max_nreqs_D[3],
             int             dimids_D[3][3],
             int             fix_D[3][3],
             int             rec_D[3][3])
{
    int i, err=NC_NOERR, nprocs;
    MPI_Offset *dims=decom.dims[2];

    /*
       nbnd  =   2 ;
       chars =   8 ;
       lev   =  72 ;
       ilev  =  73 ;
       ncol  = 866 ;

      num_decomp = 3 ;
      D1 : ncol
      D2 : ncol
      D3 : lev x ncol
    */

    /* define dimensions */
    DEF_DIM("time",  NC_UNLIMITED, dim_time)
    DEF_DIM("nbnd",  2,            dim_nbnd)
    DEF_DIM("chars", 8,            dim_chars)
    DEF_DIM("lev",   dims[0],      dim_lev)
    DEF_DIM("ilev",  dims[0] + 1,  dim_ilev)
    DEF_DIM("ncol",  dims[1],      dim_ncol)

    dimids_D[0][0] = *dim_time;
    dimids_D[1][0] = *dim_time;
    dimids_D[2][0] = *dim_time;
    dimids_D[0][1] = *dim_ncol;
    dimids_D[1][1] = *dim_ncol;
    dimids_D[2][1] = *dim_lev;  dimids_D[2][2] = *dim_ncol;

    if (cfg.strategy == blob && (cfg.api == pnetcdf || cfg.api == hdf5)) {
        MPI_Comm_size(cfg.sub_comm, &nprocs);
        DEF_DIM("nblobs", (MPI_Offset)nprocs, nblobs_ID)

        DEF_DIM("D1.nelems", decom.nelems[0], &nelems_D[0])
        DEF_DIM("D2.nelems", decom.nelems[1], &nelems_D[1])
        DEF_DIM("D3.nelems", decom.nelems[2], &nelems_D[2])

        DEF_DIM("D1.max_nreqs", (MPI_Offset)decom.max_nreqs[0], &max_nreqs_D[0])
        DEF_DIM("D2.max_nreqs", (MPI_Offset)decom.max_nreqs[1], &max_nreqs_D[1])
        DEF_DIM("D3.max_nreqs", (MPI_Offset)decom.max_nreqs[2], &max_nreqs_D[2])

        for (i=0; i<decom.num_decomp; i++) {
            fix_D[i][0] = nelems_D[i];
            rec_D[i][1] = nelems_D[i];
        }
    }
    else { /* canonical or ADIOS blob */
        fix_D[0][0] = *dim_ncol;
        fix_D[1][0] = *dim_ncol;
        fix_D[2][0] = *dim_lev;  fix_D[2][1] = *dim_ncol;
        rec_D[0][1] = *dim_ncol;
        rec_D[1][1] = *dim_ncol;
        rec_D[2][1] = *dim_lev;  rec_D[2][2] = *dim_ncol;
    }

    for (i=0; i<decom.num_decomp; i++)
        rec_D[i][0] = *dim_time;

err_out:
    return err;
}

/*----< def_var_decomp() >---------------------------------------------------*/
static
int def_var_decomp(e3sm_io_config &cfg,
                   e3sm_io_decom  &decom,
                   e3sm_io_driver &driver,
                   int             ncid,
                   int             nblobs_ID,
                   int             max_nreqs_D[3],
                   int             dimids_D[3][3],
                   int            *varids,
                   int            *nvars_decomp)
{
    int err, *varid, dimids[3];

    varid = varids - 1;

    /* define decomposition variables */

    dimids[0] = nblobs_ID;
    dimids[1] = max_nreqs_D[0];
    DEF_VAR("D1.nreqs", NC_INT, 1, dimids)
    PUT_ATTR_TXT("description", "Number of noncontiguous requests per blob")
    PUT_ATTR_INT("global_dimids", 1, dimids_D[0]+1)
    DEF_VAR("D1.blob_start", NC_INT64, 1, dimids)
    PUT_ATTR_TXT("description", "Starting variable array index per blob")
    DEF_VAR("D1.blob_count", NC_INT64, 1, dimids)
    PUT_ATTR_TXT("description", "Number of variable array elements per blob")
    DEF_VAR("D1.offsets", NC_INT, 2, dimids)
    PUT_ATTR_TXT("description", "Flattened starting indices of noncontiguous requests")
    DEF_VAR("D1.lengths", NC_INT, 2, dimids)
    PUT_ATTR_TXT("description", "Lengths of noncontiguous requests")

    dimids[1] = max_nreqs_D[1];
    DEF_VAR("D2.nreqs", NC_INT, 1, dimids)
    PUT_ATTR_TXT("description", "Number of noncontiguous requests per blob")
    PUT_ATTR_INT("global_dimids", 1, dimids_D[1]+1)
    DEF_VAR("D2.blob_start", NC_INT64, 1, dimids)
    PUT_ATTR_TXT("description", "Starting variable array index per blob")
    DEF_VAR("D2.blob_count", NC_INT64, 1, dimids)
    PUT_ATTR_TXT("description", "Number of variable array elements per blob")
    DEF_VAR("D2.offsets", NC_INT, 2, dimids)
    PUT_ATTR_TXT("description", "Flattened starting indices of noncontiguous requests")
    DEF_VAR("D2.lengths", NC_INT, 2, dimids)
    PUT_ATTR_TXT("description", "Lengths of noncontiguous requests")

    dimids[1] = max_nreqs_D[2];
    DEF_VAR("D3.nreqs", NC_INT, 1, dimids)
    PUT_ATTR_TXT("description", "Number of noncontiguous requests per blob")
    PUT_ATTR_INT("global_dimids", 2, dimids_D[2]+1)
    DEF_VAR("D3.blob_start", NC_INT64, 1, dimids)
    PUT_ATTR_TXT("description", "Starting variable array index per blob")
    DEF_VAR("D3.blob_count", NC_INT64, 1, dimids)
    PUT_ATTR_TXT("description", "Number of variable array elements per blob")
    DEF_VAR("D3.offsets", NC_INT, 2, dimids)
    PUT_ATTR_TXT("description", "Flattened starting indices of noncontiguous requests")
    DEF_VAR("D3.lengths", NC_INT, 2, dimids)
    PUT_ATTR_TXT("description", "Lengths of noncontiguous requests")

    *nvars_decomp = varid - varids + 1;

    assert(*nvars_decomp == NVARS_DECOMP * decom.num_decomp);

err_out:
    return err;
}

/*----< def_F_case_h0() >----------------------------------------------------*/
int def_F_case_h0(e3sm_io_config &cfg,
                  e3sm_io_decom  &decom,
                  e3sm_io_driver &driver,
                  int             ncid,     /* file ID */
                  int            *varids)   /* variable IDs */
{
    int one=1, two=2, three=3, *varid;
    int err, ndims, dimids[3], mdims=1, nvars_decomp=0;
    int dim_ncol, dim_time, dim_nbnd, dim_chars, dim_lev, dim_ilev;
    int nblobs_ID, nelems_D[3], max_nreqs_D[3], dimids_D[3][3];
    int fix_D[3][3], rec_D[3][3];
    float fillv = 1.e+20f, missv = 1.e+20f;

    /* add global attributes: */
    err = add_gattrs(cfg, decom, driver, ncid);
    CHECK_ERR

    PUT_GATTR_TXT("time_period_freq", "day_5")

    /* define dimensions */
    err = def_dims(cfg, decom, driver, ncid, &dim_ncol, &dim_time, &dim_nbnd,
                   &dim_chars, &dim_lev, &dim_ilev, &nblobs_ID, nelems_D,
                   max_nreqs_D, dimids_D, fix_D, rec_D);
    CHECK_ERR

    /* define variables related to decompositions */
    if (cfg.strategy == blob && (cfg.api == pnetcdf || cfg.api == hdf5)) {
        err = def_var_decomp(cfg, decom, driver, ncid, nblobs_ID, max_nreqs_D,
                             dimids_D, varids, &nvars_decomp);
        CHECK_ERR

        varid = varids + nvars_decomp - 1;
    }
    else
        varid = varids - 1;

    /* Total 414 climate variables (15 fixed-size and 399 record variables) */

    /* below 10 are fixed-size climate variables ---------------------------*/

    /* double lat(ncol) */
    DEF_VAR("lat", NC_DOUBLE, 1, fix_D[1])
    PUT_ATTR_TXT("long_name", "latitude")
    PUT_ATTR_TXT("units", "degrees_north")
    PUT_ATTR_DECOMP(two, 1, dimids_D[1]+1)

    /* double lon(ncol) */
    DEF_VAR("lon", NC_DOUBLE, 1, fix_D[1])
    PUT_ATTR_TXT("long_name", "longitude")
    PUT_ATTR_TXT("units", "degrees_east")
    PUT_ATTR_DECOMP(two, 1, dimids_D[1]+1)

    /* double area(ncol) */
    DEF_VAR("area", NC_DOUBLE, 1, fix_D[0])
    PUT_ATTR_TXT("long_name", "gll grid areas")
    PUT_ATTR_DECOMP(one, 1, dimids_D[0]+1)

    /* double lev(lev) */
    DEF_VAR("lev", NC_DOUBLE, 1, &dim_lev)
    PUT_ATTR_TXT("long_name", "hybrid level at midpoints (1000*(A+B))")
    PUT_ATTR_TXT("units", "hPa")
    PUT_ATTR_TXT("positive", "down")
    PUT_ATTR_TXT("standard_name", "atmosphere_hybrid_sigma_pressure_coordinate")
    PUT_ATTR_TXT("formula_terms", "a: hyam b: hybm p0: P0 ps: PS")

    /* double hyam(lev) */
    DEF_VAR("hyam", NC_DOUBLE, 1, &dim_lev)
    PUT_ATTR_TXT("long_name", "hybrid A coefficient at layer midpoints")

    /* double hybm(lev) */
    DEF_VAR("hybm", NC_DOUBLE, 1, &dim_lev)
    PUT_ATTR_TXT("long_name", "hybrid B coefficient at layer midpoints")

    /* double P0 */
    DEF_VAR("P0", NC_DOUBLE, 0, NULL)
    PUT_ATTR_TXT("long_name", "reference pressure")
    PUT_ATTR_TXT("units", "Pa")

    /* double ilev(ilev) */
    DEF_VAR("ilev", NC_DOUBLE, 1, &dim_ilev)
    PUT_ATTR_TXT("long_name", "hybrid level at interfaces (1000*(A+B))")
    PUT_ATTR_TXT("units", "hPa")
    PUT_ATTR_TXT("positive", "down")
    PUT_ATTR_TXT("standard_name", "atmosphere_hybrid_sigma_pressure_coordinate")
    PUT_ATTR_TXT("formula_terms", "a: hyai b: hybi p0: P0 ps: PS")

    /* double hyai(ilev) */
    DEF_VAR("hyai", NC_DOUBLE, 1, &dim_ilev)
    PUT_ATTR_TXT("long_name", "hybrid A coefficient at layer interfaces")

    /* double hybi(ilev) */
    DEF_VAR("hybi", NC_DOUBLE, 1, &dim_ilev)
    PUT_ATTR_TXT("long_name", "hybrid B coefficient at layer interfaces")

    /* below 6 are record climate variables --------------------------------*/

    /* double time(time) */
    DEF_VAR("time", NC_DOUBLE, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "time")
    PUT_ATTR_TXT("units", "days since 0001-01-01 00:00:00")
    PUT_ATTR_TXT("calendar", "noleap")
    PUT_ATTR_TXT("bounds", "time_bnds")

    /* int date(time) */
    DEF_VAR("date", NC_INT, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "current date (YYYYMMDD)")

    /* int datesec(time) */
    DEF_VAR("datesec", NC_INT, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "current seconds of current date")

    /* double time_bnds(time, nbnd) */
    dimids[0] = dim_time;
    dimids[1] = dim_nbnd;
    DEF_VAR("time_bnds", NC_DOUBLE, 2, dimids)
    PUT_ATTR_TXT("long_name", "time interval endpoints")

    /* char date_written(time, chars) */
    dimids[0] = dim_time;
    dimids[1] = dim_chars;
    DEF_VAR("date_written", NC_CHAR, 2, dimids)

    /* char time_written(time, chars) */
    dimids[0] = dim_time;
    dimids[1] = dim_chars;
    DEF_VAR("time_written", NC_CHAR, 2, dimids)

    /* below 5 are fixed-size climate variables ----------------------------*/

    /* int ndbase */
    DEF_VAR("ndbase", NC_INT, 0, NULL)
    PUT_ATTR_TXT("long_name", "base day")

    /* int nsbase */
    DEF_VAR("nsbase", NC_INT, 0, NULL)
    PUT_ATTR_TXT("long_name", "seconds of base day")

    /* int nbdate */
    DEF_VAR("nbdate", NC_INT, 0, NULL)
    PUT_ATTR_TXT("long_name", "base date (YYYYMMDD)")

    /* int nbsec */
    DEF_VAR("nbsec", NC_INT, 0, NULL)
    PUT_ATTR_TXT("long_name", "seconds of base date")

    /* int mdt */
    DEF_VAR("mdt", NC_INT, 0, NULL)
    PUT_ATTR_TXT("long_name", "timestep")
    PUT_ATTR_TXT("units", "s")

    /* below 393 are record climate variables ------------------------------*/

    /* int ndcur(time) */
    DEF_VAR("ndcur", NC_INT, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "current day (from base day)")

    /* int nscur(time) */
    DEF_VAR("nscur", NC_INT, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "current seconds of current day")

    /* double co2vmr(time) */
    DEF_VAR("co2vmr", NC_DOUBLE, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "co2 volume mixing ratio")

    /* double ch4vmr(time) */
    DEF_VAR("ch4vmr", NC_DOUBLE, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "ch4 volume mixing ratio")

    /* double n2ovmr(time) */
    DEF_VAR("n2ovmr", NC_DOUBLE, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "n2o volume mixing ratio")

    /* double f11vmr(time) */
    DEF_VAR("f11vmr", NC_DOUBLE, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "f11 volume mixing ratio")

    /* double f12vmr(time) */
    DEF_VAR("f12vmr", NC_DOUBLE, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "f12 volume mixing ratio")

    /* double sol_tsi(time) */
    DEF_VAR("sol_tsi", NC_DOUBLE, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "total solar irradiance")
    PUT_ATTR_TXT("units", "W/m2")

    /* int nsteph(time) */
    DEF_VAR("nsteph", NC_INT, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "current timestep")

    /* float AEROD_v(time, ncol) */
    DEF_VAR("AEROD_v", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Total Aerosol Optical Depth in visible band")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ANRAIN(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("ANRAIN", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m-3")
    PUT_ATTR_TXT("long_name", "Average rain number conc")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float ANSNOW(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("ANSNOW", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m-3")
    PUT_ATTR_TXT("long_name", "Average snow number conc")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float AODABS(time, ncol) */
    DEF_VAR("AODABS", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol absorption optical depth 550 nm")
    PUT_ATTR_TXT("standard_name", "atmosphere_absorption_optical_thickness_due_to_ambient_aerosol_particles")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODABSBC(time, ncol) */
    DEF_VAR("AODABSBC", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol absorption optical depth 550 nm from BC")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODALL(time, ncol) */
    DEF_VAR("AODALL", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "AOD 550 nm for all time and species")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODBC(time, ncol) */
    DEF_VAR("AODBC", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from BC")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODDUST(time, ncol) */
    DEF_VAR("AODDUST", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from dust")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODDUST1(time, ncol) */
    DEF_VAR("AODDUST1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm model 1 from dust")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODDUST3(time, ncol) */
    DEF_VAR("AODDUST3", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm model 3 from dust")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODDUST4(time, ncol) */
    DEF_VAR("AODDUST4", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm model 4 from dust")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODMODE1(time, ncol) */
    DEF_VAR("AODMODE1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm mode 1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODMODE2(time, ncol) */
    DEF_VAR("AODMODE2", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm mode 2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODMODE3(time, ncol) */
    DEF_VAR("AODMODE3", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm mode 3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODMODE4(time, ncol) */
    DEF_VAR("AODMODE4", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm mode 4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODNIR(time, ncol) */
    DEF_VAR("AODNIR", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 850 nm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODPOM(time, ncol) */
    DEF_VAR("AODPOM", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from POM")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODSO4(time, ncol) */
    DEF_VAR("AODSO4", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from SO4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODSOA(time, ncol) */
    DEF_VAR("AODSOA", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from SOA")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODSS(time, ncol) */
    DEF_VAR("AODSS", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from seasalt")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODUV(time, ncol) */
    DEF_VAR("AODUV", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 350 nm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AODVIS(time, ncol) */
    DEF_VAR("AODVIS", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm")
    PUT_ATTR_TXT("standard_name", "atmosphere_optical_thickness_due_to_ambient_aerosol_particles")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AQRAIN(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("AQRAIN", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "Average rain mixing ratio")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float AQSNOW(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("AQSNOW", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "Average snow mixing ratio")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float AQ_DMS(time, ncol) */
    DEF_VAR("AQ_DMS", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "DMS aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AQ_H2O2(time, ncol) */
    DEF_VAR("AQ_H2O2", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "H2O2 aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AQ_H2SO4(time, ncol) */
    DEF_VAR("AQ_H2SO4", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "H2SO4 aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AQ_O3(time, ncol) */
    DEF_VAR("AQ_O3", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "O3 aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AQ_SO2(time, ncol) */
    DEF_VAR("AQ_SO2", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "SO2 aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AQ_SOAG(time, ncol) */
    DEF_VAR("AQ_SOAG", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "SOAG aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float AREI(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("AREI", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "Micron")
    PUT_ATTR_TXT("long_name", "Average ice effective radius")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float AREL(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("AREL", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "Micron")
    PUT_ATTR_TXT("long_name", "Average droplet effective radius")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float AWNC(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("AWNC", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m-3")
    PUT_ATTR_TXT("long_name", "Average cloud water number conc")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float AWNI(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("AWNI", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m-3")
    PUT_ATTR_TXT("long_name", "Average cloud ice number conc")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float BURDEN1(time, ncol) */
    DEF_VAR("BURDEN1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Aerosol burden mode 1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float BURDEN2(time, ncol) */
    DEF_VAR("BURDEN2", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Aerosol burden mode 2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float BURDEN3(time, ncol) */
    DEF_VAR("BURDEN3", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Aerosol burden mode 3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float BURDEN4(time, ncol) */
    DEF_VAR("BURDEN4", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Aerosol burden mode 4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float CCN3(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("CCN3", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/cm3")
    PUT_ATTR_TXT("long_name", "CCN concentration at S=0.1%")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float CDNUMC(time, ncol) */
    DEF_VAR("CDNUMC", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "1/m2")
    PUT_ATTR_TXT("long_name", "Vertically-integrated droplet concentration")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float CLDHGH(time, ncol) */
    DEF_VAR("CLDHGH", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Vertically-integrated high cloud")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float CLDICE(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("CLDICE", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged cloud ice amount")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float CLDLIQ(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("CLDLIQ", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged cloud liquid amount")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float CLDLOW(time, ncol) */
    DEF_VAR("CLDLOW", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Vertically-integrated low cloud")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float CLDMED(time, ncol) */
    DEF_VAR("CLDMED", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Vertically-integrated mid-level cloud")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float CLDTOT(time, ncol) */
    DEF_VAR("CLDTOT", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Vertically-integrated total cloud")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float CLOUD(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("CLOUD", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Cloud fraction")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float CLOUDFRAC_CLUBB(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("CLOUDFRAC_CLUBB", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Cloud Fraction")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float CONCLD(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("CONCLD", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "fraction")
    PUT_ATTR_TXT("long_name", "Convective cloud cover")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float DCQ(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("DCQ", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg/s")
    PUT_ATTR_TXT("long_name", "Q tendency due to moist processes")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float DF_DMS(time, ncol) */
    DEF_VAR("DF_DMS", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float DF_H2O2(time, ncol) */
    DEF_VAR("DF_H2O2", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float DF_H2SO4(time, ncol) */
    DEF_VAR("DF_H2SO4", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float DF_O3(time, ncol) */
    DEF_VAR("DF_O3", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float DF_SO2(time, ncol) */
    DEF_VAR("DF_SO2", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float DF_SOAG(time, ncol) */
    DEF_VAR("DF_SOAG", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float DMS_SRF(time, ncol) */
    DEF_VAR("DMS_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "DMS in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float DP_KCLDBASE(time, ncol) */
    DEF_VAR("DP_KCLDBASE", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Deep conv. cloudbase level index")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float DP_MFUP_MAX(time, ncol) */
    DEF_VAR("DP_MFUP_MAX", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Deep conv. column-max updraft mass flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float DP_WCLDBASE(time, ncol) */
    DEF_VAR("DP_WCLDBASE", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Deep conv. cloudbase vertical velocity")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float DSTSFMBL(time, ncol) */
    DEF_VAR("DSTSFMBL", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Mobilization flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float DTCOND(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("DTCOND", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "K/s")
    PUT_ATTR_TXT("long_name", "T tendency - moist processes")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float DTENDTH(time, ncol) */
    DEF_VAR("DTENDTH", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Dynamic Tendency of Total (vertically integrated) moist static energy")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float DTENDTQ(time, ncol) */
    DEF_VAR("DTENDTQ", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Dynamic Tendency of Total (vertically integrated) specific humidity")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float EXTINCT(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("EXTINCT", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "/m")
    PUT_ATTR_TXT("long_name", "Aerosol extinction")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float FICE(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("FICE", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fractional ice content within cloud")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float FLDS(time, ncol) */
    DEF_VAR("FLDS", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Downwelling longwave flux at surface")
    PUT_ATTR_TXT("standard_name", "surface_downwelling_longwave_flux_in_air")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float FLNS(time, ncol) */
    DEF_VAR("FLNS", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Net longwave flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float FLNSC(time, ncol) */
    DEF_VAR("FLNSC", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky net longwave flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float FLNT(time, ncol) */
    DEF_VAR("FLNT", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Net longwave flux at top of model")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float FLNTC(time, ncol) */
    DEF_VAR("FLNTC", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky net longwave flux at top of model")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float FLUT(time, ncol) */
    DEF_VAR("FLUT", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Upwelling longwave flux at top of model")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float FLUTC(time, ncol) */
    DEF_VAR("FLUTC", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky upwelling longwave flux at top of model")
    PUT_ATTR_TXT("standard_name", "toa_outgoing_longwave_flux_assuming_clear_sky")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float FREQI(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("FREQI", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fractional occurrence of ice")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float FREQL(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("FREQL", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fractional occurrence of liquid")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float FREQR(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("FREQR", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fractional occurrence of rain")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float FREQS(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("FREQS", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fractional occurrence of snow")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float FSDS(time, ncol) */
    DEF_VAR("FSDS", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Downwelling solar flux at surface")
    PUT_ATTR_TXT("standard_name", "surface_downwelling_shortwave_flux_in_air")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float FSDSC(time, ncol) */
    DEF_VAR("FSDSC", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky downwelling solar flux at surface")
    PUT_ATTR_TXT("standard_name", "surface_downwelling_shortwave_flux_in_air_assuming_clear_sky")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float FSNS(time, ncol) */
    DEF_VAR("FSNS", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Net solar flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float FSNSC(time, ncol) */
    DEF_VAR("FSNSC", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky net solar flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float FSNT(time, ncol) */
    DEF_VAR("FSNT", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Net solar flux at top of model")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float FSNTC(time, ncol) */
    DEF_VAR("FSNTC", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky net solar flux at top of model")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float FSNTOA(time, ncol) */
    DEF_VAR("FSNTOA", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Net solar flux at top of atmosphere")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float FSNTOAC(time, ncol) */
    DEF_VAR("FSNTOAC", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky net solar flux at top of atmosphere")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float FSUTOA(time, ncol) */
    DEF_VAR("FSUTOA", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Upwelling solar flux at top of atmosphere")
    PUT_ATTR_TXT("standard_name", "toa_outgoing_shortwave_flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float FSUTOAC(time, ncol) */
    DEF_VAR("FSUTOAC", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky upwelling solar flux at top of atmosphere")
    PUT_ATTR_TXT("standard_name", "toa_outgoing_shortwave_flux_assuming_clear_sky")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float F_eff(time, ncol) */
    DEF_VAR("F_eff", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Effective enrichment factor of marine organic matter")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float H2O2_SRF(time, ncol) */
    DEF_VAR("H2O2_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "H2O2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float H2SO4_SRF(time, ncol) */
    DEF_VAR("H2SO4_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "H2SO4 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float H2SO4_sfgaex1(time, ncol) */
    DEF_VAR("H2SO4_sfgaex1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "H2SO4 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ICEFRAC(time, ncol) */
    DEF_VAR("ICEFRAC", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fraction of sfc area covered by sea-ice")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ICIMR(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("ICIMR", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "Prognostic in-cloud ice mixing ratio")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float ICWMR(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("ICWMR", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "Prognostic in-cloud water mixing ratio")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float IWC(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("IWC", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/m3")
    PUT_ATTR_TXT("long_name", "Grid box average ice water content")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float LANDFRAC(time, ncol) */
    DEF_VAR("LANDFRAC", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fraction of sfc area covered by land")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float LHFLX(time, ncol) */
    DEF_VAR("LHFLX", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Surface latent heat flux")
    PUT_ATTR_TXT("standard_name", "surface_upward_latent_heat_flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float LINOZ_DO3(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("LINOZ_DO3", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/s")
    PUT_ATTR_TXT("long_name", "ozone vmr tendency by linearized ozone chemistry")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float LINOZ_DO3_PSC(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("LINOZ_DO3_PSC", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/s")
    PUT_ATTR_TXT("long_name", "ozone vmr loss by PSCs using Carille et al. (1990)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float LINOZ_O3CLIM(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("LINOZ_O3CLIM", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "climatology of ozone in LINOZ")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float LINOZ_O3COL(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("LINOZ_O3COL", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "DU")
    PUT_ATTR_TXT("long_name", "ozone column above")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float LINOZ_SFCSINK(time, ncol) */
    DEF_VAR("LINOZ_SFCSINK", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "Tg/yr/m2")
    PUT_ATTR_TXT("long_name", "surface o3 sink in LINOZ with an e-fold to a fixed concentration")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float LINOZ_SSO3(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("LINOZ_SSO3", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg")
    PUT_ATTR_TXT("long_name", "steady state ozone in LINOZ")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float LINOZ_SZA(time, ncol) */
    DEF_VAR("LINOZ_SZA", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "degrees")
    PUT_ATTR_TXT("long_name", "solar zenith angle in LINOZ")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float LND_MBL(time, ncol) */
    DEF_VAR("LND_MBL", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Soil erodibility factor")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float LWCF(time, ncol) */
    DEF_VAR("LWCF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Longwave cloud forcing")
    PUT_ATTR_TXT("standard_name", "toa_longwave_cloud_radiative_effect")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float Mass_bc(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("Mass_bc", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of bc mass concentration bc_a1+bc_c1+bc_a3+bc_c3+bc_a4+bc_c4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float Mass_dst(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("Mass_dst", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of dst mass concentration dst_a1+dst_c1+dst_a3+dst_c3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float Mass_mom(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("Mass_mom", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of mom mass concentration mom_a1+mom_c1+mom_a2+mom_c2+mom_a3+mom_c3+mom_a4+mom_c4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float Mass_ncl(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("Mass_ncl", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of ncl mass concentration ncl_a1+ncl_c1+ncl_a2+ncl_c2+ncl_a3+ncl_c3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float Mass_pom(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("Mass_pom", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of pom mass concentration pom_a1+pom_c1+pom_a3+pom_c3+pom_a4+pom_c4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float Mass_so4(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("Mass_so4", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of so4 mass concentration so4_a1+so4_c1+so4_a2+so4_c2+so4_a3+so4_c3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float Mass_soa(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("Mass_soa", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of soa mass concentration soa_a1+soa_c1+soa_a2+soa_c2+soa_a3+soa_c3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float NUMICE(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("NUMICE", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged cloud ice number")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float NUMLIQ(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("NUMLIQ", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged cloud liquid number")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float NUMRAI(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("NUMRAI", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged rain number")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float NUMSNO(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("NUMSNO", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged snow number")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float O3(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("O3", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("mixing_ratio", "dry")
    PUT_ATTR_TXT("long_name", "O3 concentration")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float O3_SRF(time, ncol) */
    DEF_VAR("O3_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "O3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float OCNFRAC(time, ncol) */
    DEF_VAR("OCNFRAC", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fraction of sfc area covered by ocean")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float OMEGA(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("OMEGA", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "Pa/s")
    PUT_ATTR_TXT("long_name", "Vertical velocity (pressure)")
    PUT_ATTR_TXT("standard_name", "lagrangian_tendency_of_air_pressure")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float OMEGA500(time, ncol) */
    DEF_VAR("OMEGA500", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "Pa/s")
    PUT_ATTR_TXT("long_name", "Vertical velocity at 500 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float OMEGAT(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("OMEGAT", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "K Pa/s")
    PUT_ATTR_TXT("long_name", "Vertical heat flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float PBLH(time, ncol) */
    DEF_VAR("PBLH", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "PBL height")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float PHIS(time, ncol) */
    DEF_VAR("PHIS", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m2/s2")
    PUT_ATTR_TXT("long_name", "Surface geopotential")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float PRECC(time, ncol) */
    DEF_VAR("PRECC", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Convective precipitation rate (liq + ice)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float PRECL(time, ncol) */
    DEF_VAR("PRECL", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Large-scale (stable) precipitation rate (liq + ice)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float PRECSC(time, ncol) */
    DEF_VAR("PRECSC", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Convective snow rate (water equivalent)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float PRECSL(time, ncol) */
    DEF_VAR("PRECSL", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Large-scale (stable) snow rate (water equivalent)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float PS(time, ncol) */
    DEF_VAR("PS", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("long_name", "Surface pressure")
    PUT_ATTR_TXT("standard_name", "surface_air_pressure")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float PSL(time, ncol) */
    DEF_VAR("PSL", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("long_name", "Sea level pressure")
    PUT_ATTR_TXT("standard_name", "air_pressure_at_mean_sea_level")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float Q(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("Q", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Specific humidity")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float QFLX(time, ncol) */
    DEF_VAR("QFLX", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Surface water flux")
    PUT_ATTR_TXT("standard_name", "water_evapotranspiration_flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float QREFHT(time, ncol) */
    DEF_VAR("QREFHT", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "Reference height humidity")
    PUT_ATTR_TXT("standard_name", "specific_humidity")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float QRL(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("QRL", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "K/s")
    PUT_ATTR_TXT("long_name", "Longwave heating rate")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float QRS(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("QRS", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "K/s")
    PUT_ATTR_TXT("long_name", "Solar heating rate")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float RAINQM(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("RAINQM", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged rain amount")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float RAM1(time, ncol) */
    DEF_VAR("RAM1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "frac")
    PUT_ATTR_TXT("long_name", "RAM1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float RELHUM(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("RELHUM", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "percent")
    PUT_ATTR_TXT("long_name", "Relative humidity")
    PUT_ATTR_TXT("standard_name", "relative_humidity")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float SFDMS(time, ncol) */
    DEF_VAR("SFDMS", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "DMS surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFH2O2(time, ncol) */
    DEF_VAR("SFH2O2", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "H2O2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFH2SO4(time, ncol) */
    DEF_VAR("SFH2SO4", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "H2SO4 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFO3(time, ncol) */
    DEF_VAR("SFO3", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "O3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFSO2(time, ncol) */
    DEF_VAR("SFSO2", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "SO2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFSOAG(time, ncol) */
    DEF_VAR("SFSOAG", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "SOAG surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFbc_a1(time, ncol) */
    DEF_VAR("SFbc_a1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFbc_a3(time, ncol) */
    DEF_VAR("SFbc_a3", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFbc_a4(time, ncol) */
    DEF_VAR("SFbc_a4", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a4 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFdst_a1(time, ncol) */
    DEF_VAR("SFdst_a1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFdst_a3(time, ncol) */
    DEF_VAR("SFdst_a3", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFmom_a1(time, ncol) */
    DEF_VAR("SFmom_a1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFmom_a2(time, ncol) */
    DEF_VAR("SFmom_a2", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFmom_a3(time, ncol) */
    DEF_VAR("SFmom_a3", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFmom_a4(time, ncol) */
    DEF_VAR("SFmom_a4", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a4 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFncl_a1(time, ncol) */
    DEF_VAR("SFncl_a1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFncl_a2(time, ncol) */
    DEF_VAR("SFncl_a2", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFncl_a3(time, ncol) */
    DEF_VAR("SFncl_a3", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFnum_a1(time, ncol) */
    DEF_VAR("SFnum_a1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFnum_a2(time, ncol) */
    DEF_VAR("SFnum_a2", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFnum_a3(time, ncol) */
    DEF_VAR("SFnum_a3", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFnum_a4(time, ncol) */
    DEF_VAR("SFnum_a4", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a4 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFpom_a1(time, ncol) */
    DEF_VAR("SFpom_a1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFpom_a3(time, ncol) */
    DEF_VAR("SFpom_a3", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFpom_a4(time, ncol) */
    DEF_VAR("SFpom_a4", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a4 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFso4_a1(time, ncol) */
    DEF_VAR("SFso4_a1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFso4_a2(time, ncol) */
    DEF_VAR("SFso4_a2", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFso4_a3(time, ncol) */
    DEF_VAR("SFso4_a3", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFsoa_a1(time, ncol) */
    DEF_VAR("SFsoa_a1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFsoa_a2(time, ncol) */
    DEF_VAR("SFsoa_a2", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SFsoa_a3(time, ncol) */
    DEF_VAR("SFsoa_a3", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SHFLX(time, ncol) */
    DEF_VAR("SHFLX", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Surface sensible heat flux")
    PUT_ATTR_TXT("standard_name", "surface_upward_sensible_heat_flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SH_KCLDBASE(time, ncol) */
    DEF_VAR("SH_KCLDBASE", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Shallow conv. cloudbase level index")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SH_MFUP_MAX(time, ncol) */
    DEF_VAR("SH_MFUP_MAX", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Shallow conv. column-max updraft mass flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SH_WCLDBASE(time, ncol) */
    DEF_VAR("SH_WCLDBASE", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Shallow conv. cloudbase vertical velocity")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SNOWHICE(time, ncol) */
    DEF_VAR("SNOWHICE", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "Snow depth over ice")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SNOWHLND(time, ncol) */
    DEF_VAR("SNOWHLND", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "Water equivalent snow depth")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SNOWQM(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("SNOWQM", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged snow amount")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float SO2(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("SO2", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("mixing_ratio", "dry")
    PUT_ATTR_TXT("long_name", "SO2 concentration")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float SO2_CLXF(time, ncol) */
    DEF_VAR("SO2_CLXF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for SO2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SO2_SRF(time, ncol) */
    DEF_VAR("SO2_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "SO2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SOAG_CLXF(time, ncol) */
    DEF_VAR("SOAG_CLXF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for SOAG")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SOAG_SRF(time, ncol) */
    DEF_VAR("SOAG_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "SOAG in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SOAG_sfgaex1(time, ncol) */
    DEF_VAR("SOAG_sfgaex1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "SOAG gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SOLIN(time, ncol) */
    DEF_VAR("SOLIN", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Solar insolation")
    PUT_ATTR_TXT("standard_name", "toa_incoming_shortwave_flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SSAVIS(time, ncol) */
    DEF_VAR("SSAVIS", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("long_name", "Aerosol singel-scatter albedo")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SSTSFMBL(time, ncol) */
    DEF_VAR("SSTSFMBL", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Mobilization flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SSTSFMBL_OM(time, ncol) */
    DEF_VAR("SSTSFMBL_OM", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Mobilization flux of marine organic matter at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SWCF(time, ncol) */
    DEF_VAR("SWCF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Shortwave cloud forcing")
    PUT_ATTR_TXT("standard_name", "toa_shortwave_cloud_radiative_effect")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float T(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("T", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Temperature")
    PUT_ATTR_TXT("standard_name", "air_temperature")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float TAUGWX(time, ncol) */
    DEF_VAR("TAUGWX", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "N/m2")
    PUT_ATTR_TXT("long_name", "Zonal gravity wave surface stress")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TAUGWY(time, ncol) */
    DEF_VAR("TAUGWY", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "N/m2")
    PUT_ATTR_TXT("long_name", "Meridional gravity wave surface stress")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TAUX(time, ncol) */
    DEF_VAR("TAUX", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "N/m2")
    PUT_ATTR_TXT("long_name", "Zonal surface stress")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TAUY(time, ncol) */
    DEF_VAR("TAUY", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "N/m2")
    PUT_ATTR_TXT("long_name", "Meridional surface stress")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TGCLDCWP(time, ncol) */
    DEF_VAR("TGCLDCWP", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Total grid-box cloud water path (liquid and ice)")
    PUT_ATTR_TXT("standard_name", "atmosphere_mass_content_of_cloud_condensed_water")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TGCLDIWP(time, ncol) */
    DEF_VAR("TGCLDIWP", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Total grid-box cloud ice water path")
    PUT_ATTR_TXT("standard_name", "atmosphere_mass_content_of_cloud_ice")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TGCLDLWP(time, ncol) */
    DEF_VAR("TGCLDLWP", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Total grid-box cloud liquid water path")
    PUT_ATTR_TXT("standard_name", "atmosphere_mass_content_of_cloud_liquid_water")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TH7001000(time, ncol) */
    DEF_VAR("TH7001000", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Theta difference 700 mb - 1000 mb")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TMQ(time, ncol) */
    DEF_VAR("TMQ", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Total (vertically integrated) precipitable water")
    PUT_ATTR_TXT("standard_name", "atmosphere_mass_content_of_water_vapor")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TREFHT(time, ncol) */
    DEF_VAR("TREFHT", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Reference height temperature")
    PUT_ATTR_TXT("standard_name", "air_temperature")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TROP_P(time, ncol) */
    DEF_VAR("TROP_P", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("long_name", "Tropopause Pressure")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TROP_T(time, ncol) */
    DEF_VAR("TROP_T", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Tropopause Temperature")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TS(time, ncol) */
    DEF_VAR("TS", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Surface temperature (radiative)")
    PUT_ATTR_TXT("standard_name", "surface_temperature")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TSMN(time, ncol) */
    DEF_VAR("TSMN", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Minimum surface temperature over output period")
    PUT_ATTR_TXT("cell_methods", "time: minimum")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TSMX(time, ncol) */
    DEF_VAR("TSMX", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Maximum surface temperature over output period")
    PUT_ATTR_TXT("cell_methods", "time: maximum")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TUH(time, ncol) */
    DEF_VAR("TUH", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "W/m")
    PUT_ATTR_TXT("long_name", "Total (vertically integrated) zonal MSE flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TUQ(time, ncol) */
    DEF_VAR("TUQ", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m/s")
    PUT_ATTR_TXT("long_name", "Total (vertically integrated) zonal water flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TVH(time, ncol) */
    DEF_VAR("TVH", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "W/m")
    PUT_ATTR_TXT("long_name", "Total (vertically integrated) meridional MSE flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TVQ(time, ncol) */
    DEF_VAR("TVQ", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m/s")
    PUT_ATTR_TXT("long_name", "Total (vertically integrated) meridional water flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float U(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("U", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Zonal wind")
    PUT_ATTR_TXT("standard_name", "eastward_wind")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float U10(time, ncol) */
    DEF_VAR("U10", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "10m wind speed")
    PUT_ATTR_TXT("standard_name", "wind_speed")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float UU(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("UU", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m2/s2")
    PUT_ATTR_TXT("long_name", "Zonal velocity squared")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float V(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("V", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Meridional wind")
    PUT_ATTR_TXT("standard_name", "northward_wind")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float VQ(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("VQ", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m/s kg/kg")
    PUT_ATTR_TXT("long_name", "Meridional water transport")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float VT(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("VT", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "K m/s")
    PUT_ATTR_TXT("long_name", "Meridional heat transport")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float VU(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("VU", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m2/s2")
    PUT_ATTR_TXT("long_name", "Meridional flux of zonal momentum")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float VV(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("VV", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m2/s2")
    PUT_ATTR_TXT("long_name", "Meridional velocity squared")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float WD_H2O2(time, ncol) */
    DEF_VAR("WD_H2O2", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/s")
    PUT_ATTR_TXT("long_name", "H2O2             wet deposition")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float WD_H2SO4(time, ncol) */
    DEF_VAR("WD_H2SO4", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/s")
    PUT_ATTR_TXT("long_name", "H2SO4            wet deposition")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float WD_SO2(time, ncol) */
    DEF_VAR("WD_SO2", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/s")
    PUT_ATTR_TXT("long_name", "SO2              wet deposition")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float WSUB(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("WSUB", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Diagnostic sub-grid vertical velocity")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float Z3(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("Z3", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "Geopotential Height (above sea level)")
    PUT_ATTR_TXT("standard_name", "geopotential_height")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float aero_water(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("aero_water", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "sum of aerosol water of interstitial modes wat_a1+wat_a2+wat_a3+wat_a4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float airFV(time, ncol) */
    DEF_VAR("airFV", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "frac")
    PUT_ATTR_TXT("long_name", "FV")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_a1DDF(time, ncol) */
    DEF_VAR("bc_a1DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_a1SFWET(time, ncol) */
    DEF_VAR("bc_a1SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_a1_SRF(time, ncol) */
    DEF_VAR("bc_a1_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "bc_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_a1_sfgaex1(time, ncol) */
    DEF_VAR("bc_a1_sfgaex1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_a3DDF(time, ncol) */
    DEF_VAR("bc_a3DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_a3SFWET(time, ncol) */
    DEF_VAR("bc_a3SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_a3_SRF(time, ncol) */
    DEF_VAR("bc_a3_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "bc_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_a4DDF(time, ncol) */
    DEF_VAR("bc_a4DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_a4SFWET(time, ncol) */
    DEF_VAR("bc_a4SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_a4_CLXF(time, ncol) */
    DEF_VAR("bc_a4_CLXF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for bc_a4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_a4_SRF(time, ncol) */
    DEF_VAR("bc_a4_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "bc_a4 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_a4_sfgaex1(time, ncol) */
    DEF_VAR("bc_a4_sfgaex1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a4 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_c1DDF(time, ncol) */
    DEF_VAR("bc_c1DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_c1SFWET(time, ncol) */
    DEF_VAR("bc_c1SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_c3DDF(time, ncol) */
    DEF_VAR("bc_c3DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_c3SFWET(time, ncol) */
    DEF_VAR("bc_c3SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_c4DDF(time, ncol) */
    DEF_VAR("bc_c4DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float bc_c4SFWET(time, ncol) */
    DEF_VAR("bc_c4SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c4 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float chla(time, ncol) */
    DEF_VAR("chla", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "mg L-1")
    PUT_ATTR_TXT("long_name", "ocean input data: chla")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float dst_a1DDF(time, ncol) */
    DEF_VAR("dst_a1DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float dst_a1SF(time, ncol) */
    DEF_VAR("dst_a1SF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a1 dust surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float dst_a1SFWET(time, ncol) */
    DEF_VAR("dst_a1SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float dst_a1_SRF(time, ncol) */
    DEF_VAR("dst_a1_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "dst_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float dst_a3DDF(time, ncol) */
    DEF_VAR("dst_a3DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float dst_a3SF(time, ncol) */
    DEF_VAR("dst_a3SF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a3 dust surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float dst_a3SFWET(time, ncol) */
    DEF_VAR("dst_a3SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float dst_a3_SRF(time, ncol) */
    DEF_VAR("dst_a3_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "dst_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float dst_c1DDF(time, ncol) */
    DEF_VAR("dst_c1DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float dst_c1SFWET(time, ncol) */
    DEF_VAR("dst_c1SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float dst_c3DDF(time, ncol) */
    DEF_VAR("dst_c3DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float dst_c3SFWET(time, ncol) */
    DEF_VAR("dst_c3SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float hstobie_linoz(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("hstobie_linoz", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "fraction of model time")
    PUT_ATTR_TXT("long_name", "Lowest possible Linoz level")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float mlip(time, ncol) */
    DEF_VAR("mlip", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "uM C")
    PUT_ATTR_TXT("long_name", "ocean input data: mlip")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_a1DDF(time, ncol) */
    DEF_VAR("mom_a1DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_a1SF(time, ncol) */
    DEF_VAR("mom_a1SF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a1 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_a1SFWET(time, ncol) */
    DEF_VAR("mom_a1SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_a1_SRF(time, ncol) */
    DEF_VAR("mom_a1_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "mom_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_a1_sfgaex1(time, ncol) */
    DEF_VAR("mom_a1_sfgaex1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_a2DDF(time, ncol) */
    DEF_VAR("mom_a2DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_a2SF(time, ncol) */
    DEF_VAR("mom_a2SF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a2 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_a2SFWET(time, ncol) */
    DEF_VAR("mom_a2SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_a2_SRF(time, ncol) */
    DEF_VAR("mom_a2_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "mom_a2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_a3DDF(time, ncol) */
    DEF_VAR("mom_a3DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_a3SFWET(time, ncol) */
    DEF_VAR("mom_a3SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_a3_SRF(time, ncol) */
    DEF_VAR("mom_a3_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "mom_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_a4DDF(time, ncol) */
    DEF_VAR("mom_a4DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_a4SF(time, ncol) */
    DEF_VAR("mom_a4SF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a4 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_a4SFWET(time, ncol) */
    DEF_VAR("mom_a4SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_a4_SRF(time, ncol) */
    DEF_VAR("mom_a4_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "mom_a4 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_a4_sfgaex1(time, ncol) */
    DEF_VAR("mom_a4_sfgaex1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a4 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_c1DDF(time, ncol) */
    DEF_VAR("mom_c1DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_c1SFWET(time, ncol) */
    DEF_VAR("mom_c1SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_c2DDF(time, ncol) */
    DEF_VAR("mom_c2DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_c2SFWET(time, ncol) */
    DEF_VAR("mom_c2SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c2 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_c3DDF(time, ncol) */
    DEF_VAR("mom_c3DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_c3SFWET(time, ncol) */
    DEF_VAR("mom_c3SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_c4DDF(time, ncol) */
    DEF_VAR("mom_c4DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mom_c4SFWET(time, ncol) */
    DEF_VAR("mom_c4SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c4 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mpoly(time, ncol) */
    DEF_VAR("mpoly", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "uM C")
    PUT_ATTR_TXT("long_name", "ocean input data: mpoly")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float mprot(time, ncol) */
    DEF_VAR("mprot", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "uM C")
    PUT_ATTR_TXT("long_name", "ocean input data: mprot")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_a1DDF(time, ncol) */
    DEF_VAR("ncl_a1DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_a1SF(time, ncol) */
    DEF_VAR("ncl_a1SF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a1 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_a1SFWET(time, ncol) */
    DEF_VAR("ncl_a1SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_a1_SRF(time, ncol) */
    DEF_VAR("ncl_a1_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "ncl_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_a2DDF(time, ncol) */
    DEF_VAR("ncl_a2DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_a2SF(time, ncol) */
    DEF_VAR("ncl_a2SF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a2 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_a2SFWET(time, ncol) */
    DEF_VAR("ncl_a2SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_a2_SRF(time, ncol) */
    DEF_VAR("ncl_a2_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "ncl_a2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_a3DDF(time, ncol) */
    DEF_VAR("ncl_a3DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_a3SF(time, ncol) */
    DEF_VAR("ncl_a3SF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a3 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_a3SFWET(time, ncol) */
    DEF_VAR("ncl_a3SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_a3_SRF(time, ncol) */
    DEF_VAR("ncl_a3_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "ncl_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_c1DDF(time, ncol) */
    DEF_VAR("ncl_c1DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_c1SFWET(time, ncol) */
    DEF_VAR("ncl_c1SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_c2DDF(time, ncol) */
    DEF_VAR("ncl_c2DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_c2SFWET(time, ncol) */
    DEF_VAR("ncl_c2SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c2 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_c3DDF(time, ncol) */
    DEF_VAR("ncl_c3DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float ncl_c3SFWET(time, ncol) */
    DEF_VAR("ncl_c3SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a1DDF(time, ncol) */
    DEF_VAR("num_a1DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a1SF(time, ncol) */
    DEF_VAR("num_a1SF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "num_a1 dust surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a1SFWET(time, ncol) */
    DEF_VAR("num_a1SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a1_CLXF(time, ncol) */
    DEF_VAR("num_a1_CLXF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for num_a1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a1_SRF(time, ncol) */
    DEF_VAR("num_a1_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/kg")
    PUT_ATTR_TXT("long_name", "num_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a1_sfgaex1(time, ncol) */
    DEF_VAR("num_a1_sfgaex1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "num_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a2DDF(time, ncol) */
    DEF_VAR("num_a2DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a2SFWET(time, ncol) */
    DEF_VAR("num_a2SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a2_CLXF(time, ncol) */
    DEF_VAR("num_a2_CLXF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for num_a2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a2_SRF(time, ncol) */
    DEF_VAR("num_a2_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/kg")
    PUT_ATTR_TXT("long_name", "num_a2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a3DDF(time, ncol) */
    DEF_VAR("num_a3DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a3SF(time, ncol) */
    DEF_VAR("num_a3SF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "num_a3 dust surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a3SFWET(time, ncol) */
    DEF_VAR("num_a3SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a3_SRF(time, ncol) */
    DEF_VAR("num_a3_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/kg")
    PUT_ATTR_TXT("long_name", "num_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a4DDF(time, ncol) */
    DEF_VAR("num_a4DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a4SFWET(time, ncol) */
    DEF_VAR("num_a4SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a4_CLXF(time, ncol) */
    DEF_VAR("num_a4_CLXF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for num_a4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a4_SRF(time, ncol) */
    DEF_VAR("num_a4_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/kg")
    PUT_ATTR_TXT("long_name", "num_a4 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_a4_sfgaex1(time, ncol) */
    DEF_VAR("num_a4_sfgaex1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "num_a4 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_c1DDF(time, ncol) */
    DEF_VAR("num_c1DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_c1SFWET(time, ncol) */
    DEF_VAR("num_c1SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_c2DDF(time, ncol) */
    DEF_VAR("num_c2DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_c2SFWET(time, ncol) */
    DEF_VAR("num_c2SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c2 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_c3DDF(time, ncol) */
    DEF_VAR("num_c3DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_c3SFWET(time, ncol) */
    DEF_VAR("num_c3SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_c4DDF(time, ncol) */
    DEF_VAR("num_c4DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float num_c4SFWET(time, ncol) */
    DEF_VAR("num_c4SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c4 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_a1DDF(time, ncol) */
    DEF_VAR("pom_a1DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_a1SFWET(time, ncol) */
    DEF_VAR("pom_a1SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_a1_SRF(time, ncol) */
    DEF_VAR("pom_a1_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "pom_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_a1_sfgaex1(time, ncol) */
    DEF_VAR("pom_a1_sfgaex1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_a3DDF(time, ncol) */
    DEF_VAR("pom_a3DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_a3SFWET(time, ncol) */
    DEF_VAR("pom_a3SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_a3_SRF(time, ncol) */
    DEF_VAR("pom_a3_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "pom_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_a4DDF(time, ncol) */
    DEF_VAR("pom_a4DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_a4SFWET(time, ncol) */
    DEF_VAR("pom_a4SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_a4_CLXF(time, ncol) */
    DEF_VAR("pom_a4_CLXF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for pom_a4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_a4_SRF(time, ncol) */
    DEF_VAR("pom_a4_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "pom_a4 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_a4_sfgaex1(time, ncol) */
    DEF_VAR("pom_a4_sfgaex1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a4 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_c1DDF(time, ncol) */
    DEF_VAR("pom_c1DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_c1SFWET(time, ncol) */
    DEF_VAR("pom_c1SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_c3DDF(time, ncol) */
    DEF_VAR("pom_c3DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_c3SFWET(time, ncol) */
    DEF_VAR("pom_c3SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_c4DDF(time, ncol) */
    DEF_VAR("pom_c4DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float pom_c4SFWET(time, ncol) */
    DEF_VAR("pom_c4SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c4 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_a1DDF(time, ncol) */
    DEF_VAR("so4_a1DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_a1SFWET(time, ncol) */
    DEF_VAR("so4_a1SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_a1_CLXF(time, ncol) */
    DEF_VAR("so4_a1_CLXF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for so4_a1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_a1_SRF(time, ncol) */
    DEF_VAR("so4_a1_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "so4_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_a1_sfgaex1(time, ncol) */
    DEF_VAR("so4_a1_sfgaex1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_a2DDF(time, ncol) */
    DEF_VAR("so4_a2DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_a2SFWET(time, ncol) */
    DEF_VAR("so4_a2SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_a2_CLXF(time, ncol) */
    DEF_VAR("so4_a2_CLXF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for so4_a2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_a2_SRF(time, ncol) */
    DEF_VAR("so4_a2_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "so4_a2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_a2_sfgaex1(time, ncol) */
    DEF_VAR("so4_a2_sfgaex1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a2 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_a3DDF(time, ncol) */
    DEF_VAR("so4_a3DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_a3SFWET(time, ncol) */
    DEF_VAR("so4_a3SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_a3_SRF(time, ncol) */
    DEF_VAR("so4_a3_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "so4_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_a3_sfgaex1(time, ncol) */
    DEF_VAR("so4_a3_sfgaex1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a3 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_c1DDF(time, ncol) */
    DEF_VAR("so4_c1DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_c1SFWET(time, ncol) */
    DEF_VAR("so4_c1SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_c2DDF(time, ncol) */
    DEF_VAR("so4_c2DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_c2SFWET(time, ncol) */
    DEF_VAR("so4_c2SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c2 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_c3DDF(time, ncol) */
    DEF_VAR("so4_c3DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float so4_c3SFWET(time, ncol) */
    DEF_VAR("so4_c3SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_a1DDF(time, ncol) */
    DEF_VAR("soa_a1DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_a1SFWET(time, ncol) */
    DEF_VAR("soa_a1SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_a1_SRF(time, ncol) */
    DEF_VAR("soa_a1_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "soa_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_a1_sfgaex1(time, ncol) */
    DEF_VAR("soa_a1_sfgaex1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_a2DDF(time, ncol) */
    DEF_VAR("soa_a2DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_a2SFWET(time, ncol) */
    DEF_VAR("soa_a2SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_a2_SRF(time, ncol) */
    DEF_VAR("soa_a2_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "soa_a2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_a2_sfgaex1(time, ncol) */
    DEF_VAR("soa_a2_sfgaex1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a2 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_a3DDF(time, ncol) */
    DEF_VAR("soa_a3DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_a3SFWET(time, ncol) */
    DEF_VAR("soa_a3SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_a3_SRF(time, ncol) */
    DEF_VAR("soa_a3_SRF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "soa_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_a3_sfgaex1(time, ncol) */
    DEF_VAR("soa_a3_sfgaex1", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a3 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_c1DDF(time, ncol) */
    DEF_VAR("soa_c1DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_c1SFWET(time, ncol) */
    DEF_VAR("soa_c1SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_c2DDF(time, ncol) */
    DEF_VAR("soa_c2DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_c2SFWET(time, ncol) */
    DEF_VAR("soa_c2SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c2 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_c3DDF(time, ncol) */
    DEF_VAR("soa_c3DDF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float soa_c3SFWET(time, ncol) */
    DEF_VAR("soa_c3SFWET", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    assert (varid - varids + 1 == cfg.nvars + nvars_decomp);

err_out:
    return err;
}

/*----< def_F_case_h1() >----------------------------------------------------*/
int def_F_case_h1(e3sm_io_config &cfg,
                  e3sm_io_decom  &decom,
                  e3sm_io_driver &driver,
                  int             ncid,     /* file ID */
                  int            *varids)   /* variable IDs */
{
    int one=1, two=2, three=3, *varid;
    int err, ndims, dimids[3], mdims=1, nvars_decomp=0;
    int dim_ncol, dim_time, dim_nbnd, dim_chars, dim_lev, dim_ilev;
    int nblobs_ID, nelems_D[3], max_nreqs_D[3], dimids_D[3][3];
    int fix_D[3][3], rec_D[3][3];
    float fillv = 1.e+20f, missv = 1.e+20f;

    /* add global attributes: */
    err = add_gattrs(cfg, decom, driver, ncid);
    CHECK_ERR

    PUT_GATTR_TXT("time_period_freq", "hour_1")

    /* define dimensions */
    err = def_dims(cfg, decom, driver, ncid, &dim_ncol, &dim_time, &dim_nbnd,
                   &dim_chars, &dim_lev, &dim_ilev, &nblobs_ID, nelems_D,
                   max_nreqs_D, dimids_D, fix_D, rec_D);
    CHECK_ERR

    /* define variables related to decompositions */
    if (cfg.strategy == blob && (cfg.api == pnetcdf || cfg.api == hdf5)) {
        err = def_var_decomp(cfg, decom, driver, ncid, nblobs_ID, max_nreqs_D,
                             dimids_D, varids, &nvars_decomp);
        CHECK_ERR

        varid = varids + nvars_decomp - 1;
    }
    else
        varid = varids - 1;

    /* Total 51 climate variables (15 fixed-size and 36 record variables) */

    /* below 10 are fixed-size climate variables ---------------------------*/

    /* double lat(ncol) */
    DEF_VAR("lat", NC_DOUBLE, 1, fix_D[1])
    PUT_ATTR_TXT("long_name", "latitude")
    PUT_ATTR_TXT("units", "degrees_north")
    PUT_ATTR_DECOMP(two, 1, dimids_D[1]+1)

    /* double lon(ncol) */
    DEF_VAR("lon", NC_DOUBLE, 1, fix_D[1])
    PUT_ATTR_TXT("long_name", "longitude")
    PUT_ATTR_TXT("units", "degrees_east")
    PUT_ATTR_DECOMP(two, 1, dimids_D[1]+1)

    /* double area(ncol) */
    DEF_VAR("area", NC_DOUBLE, 1, fix_D[0])
    PUT_ATTR_TXT("long_name", "gll grid areas")
    PUT_ATTR_DECOMP(one, 1, dimids_D[0]+1)

    /* double lev(lev) */
    DEF_VAR("lev", NC_DOUBLE, 1, &dim_lev)
    PUT_ATTR_TXT("long_name", "hybrid level at midpoints (1000*(A+B))")
    PUT_ATTR_TXT("units", "hPa")
    PUT_ATTR_TXT("positive", "down")
    PUT_ATTR_TXT("standard_name", "atmosphere_hybrid_sigma_pressure_coordinate")
    PUT_ATTR_TXT("formula_terms", "a: hyam b: hybm p0: P0 ps: PS")

    /* double hyam(lev) */
    DEF_VAR("hyam", NC_DOUBLE, 1, &dim_lev)
    PUT_ATTR_TXT("long_name", "hybrid A coefficient at layer midpoints")

    /* double hybm(lev) */
    DEF_VAR("hybm", NC_DOUBLE, 1, &dim_lev)
    PUT_ATTR_TXT("long_name", "hybrid B coefficient at layer midpoints")

    /* double P0 */
    DEF_VAR("P0", NC_DOUBLE, 0, NULL)
    PUT_ATTR_TXT("long_name", "reference pressure")
    PUT_ATTR_TXT("units", "Pa")

    /* double ilev(ilev) */
    DEF_VAR("ilev", NC_DOUBLE, 1, &dim_ilev)
    PUT_ATTR_TXT("long_name", "hybrid level at interfaces (1000*(A+B))")
    PUT_ATTR_TXT("units", "hPa")
    PUT_ATTR_TXT("positive", "down")
    PUT_ATTR_TXT("standard_name", "atmosphere_hybrid_sigma_pressure_coordinate")
    PUT_ATTR_TXT("formula_terms", "a: hyai b: hybi p0: P0 ps: PS")

    /* double hyai(ilev) */
    DEF_VAR("hyai", NC_DOUBLE, 1, &dim_ilev)
    PUT_ATTR_TXT("long_name", "hybrid A coefficient at layer interfaces")

    /* double hybi(ilev) */
    DEF_VAR("hybi", NC_DOUBLE, 1, &dim_ilev)
    PUT_ATTR_TXT("long_name", "hybrid B coefficient at layer interfaces")

    /* below 6 are record climate variables --------------------------------*/

    /* double time(time) */
    DEF_VAR("time", NC_DOUBLE, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "time")
    PUT_ATTR_TXT("units", "days since 0001-01-01 00:00:00")
    PUT_ATTR_TXT("calendar", "noleap")
    PUT_ATTR_TXT("bounds", "time_bnds")

    /* int date(time) */
    DEF_VAR("date", NC_INT, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "current date (YYYYMMDD)")

    /* int datesec(time) */
    DEF_VAR("datesec", NC_INT, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "current seconds of current date")

    /* double time_bnds(time, nbnd) */
    dimids[0] = dim_time;
    dimids[1] = dim_nbnd;
    DEF_VAR("time_bnds", NC_DOUBLE, 2, dimids)
    PUT_ATTR_TXT("long_name", "time interval endpoints")

    /* char date_written(time, chars) */
    dimids[0] = dim_time;
    dimids[1] = dim_chars;
    DEF_VAR("date_written", NC_CHAR, 2, dimids)

    /* char time_written(time, chars) */
    dimids[0] = dim_time;
    dimids[1] = dim_chars;
    DEF_VAR("time_written", NC_CHAR, 2, dimids)

    /* below 5 are fixed-size climate variables ----------------------------*/

    /* int ndbase */
    DEF_VAR("ndbase", NC_INT, 0, NULL)
    PUT_ATTR_TXT("long_name", "base day")

    /* int nsbase */
    DEF_VAR("nsbase", NC_INT, 0, NULL)
    PUT_ATTR_TXT("long_name", "seconds of base day")

    /* int nbdate */
    DEF_VAR("nbdate", NC_INT, 0, NULL)
    PUT_ATTR_TXT("long_name", "base date (YYYYMMDD)")

    /* int nbsec */
    DEF_VAR("nbsec", NC_INT, 0, NULL)
    PUT_ATTR_TXT("long_name", "seconds of base date")

    /* int mdt */
    DEF_VAR("mdt", NC_INT, 0, NULL)
    PUT_ATTR_TXT("long_name", "timestep")
    PUT_ATTR_TXT("units", "s")

    /* below 30 are record climate variables --------------------------------*/

    /* int ndcur(time) */
    DEF_VAR("ndcur", NC_INT, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "current day (from base day)")

    /* int nscur(time) */
    DEF_VAR("nscur", NC_INT, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "current seconds of current day")

    /* double co2vmr(time) */
    DEF_VAR("co2vmr", NC_DOUBLE, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "co2 volume mixing ratio")

    /* double ch4vmr(time) */
    DEF_VAR("ch4vmr", NC_DOUBLE, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "ch4 volume mixing ratio")

    /* double n2ovmr(time) */
    DEF_VAR("n2ovmr", NC_DOUBLE, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "n2o volume mixing ratio")

    /* double f11vmr(time) */
    DEF_VAR("f11vmr", NC_DOUBLE, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "f11 volume mixing ratio")

    /* double f12vmr(time) */
    DEF_VAR("f12vmr", NC_DOUBLE, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "f12 volume mixing ratio")

    /* double sol_tsi(time) */
    DEF_VAR("sol_tsi", NC_DOUBLE, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "total solar irradiance")
    PUT_ATTR_TXT("units", "W/m2")

    /* int nsteph(time) */
    DEF_VAR("nsteph", NC_INT, 1, &dim_time)
    PUT_ATTR_TXT("long_name", "current timestep")

    /* float CLDHGH(time, ncol) */
    DEF_VAR("CLDHGH", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Vertically-integrated high cloud")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float CLDLOW(time, ncol) */
    DEF_VAR("CLDLOW", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Vertically-integrated low cloud")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float CLDMED(time, ncol) */
    DEF_VAR("CLDMED", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Vertically-integrated mid-level cloud")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float FLNT(time, ncol) */
    DEF_VAR("FLNT", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Net longwave flux at top of model")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float LWCF(time, ncol) */
    DEF_VAR("LWCF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Longwave cloud forcing")
    PUT_ATTR_TXT("standard_name", "toa_longwave_cloud_radiative_effect")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float OMEGA500(time, ncol) */
    DEF_VAR("OMEGA500", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "Pa/s")
    PUT_ATTR_TXT("long_name", "Vertical velocity at 500 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float OMEGA850(time, ncol) */
    DEF_VAR("OMEGA850", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "Pa/s")
    PUT_ATTR_TXT("long_name", "Vertical velocity at 850 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float PRECT(time, ncol) */
    DEF_VAR("PRECT", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Total (convective and large-scale) precipitation rate (liq + ice)")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float PS(time, ncol) */
    DEF_VAR("PS", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("long_name", "Surface pressure")
    PUT_ATTR_TXT("standard_name", "surface_air_pressure")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float SWCF(time, ncol) */
    DEF_VAR("SWCF", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLOAT(_FillValue, 1, &fillv)
    PUT_ATTR_FLOAT("missing_value", 1, &missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Shortwave cloud forcing")
    PUT_ATTR_TXT("standard_name", "toa_shortwave_cloud_radiative_effect")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float T850(time, ncol) */
    DEF_VAR("T850", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Temperature at 850 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TMQ(time, ncol) */
    DEF_VAR("TMQ", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Total (vertically integrated) precipitable water")
    PUT_ATTR_TXT("standard_name", "atmosphere_mass_content_of_water_vapor")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float TS(time, ncol) */
    DEF_VAR("TS", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Surface temperature (radiative)")
    PUT_ATTR_TXT("standard_name", "surface_temperature")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float U(time, lev, ncol) */
    ndims = (cfg.strategy == blob) ? 2 : 3;
    DEF_VAR("U", NC_FLOAT, ndims, rec_D[2])
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Zonal wind")
    PUT_ATTR_TXT("standard_name", "eastward_wind")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(three, 3, dimids_D[2])

    /* float U250(time, ncol) */
    DEF_VAR("U250", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Zonal wind at 250 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float U850(time, ncol) */
    DEF_VAR("U850", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Zonal wind at 850 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float UBOT(time, ncol) */
    DEF_VAR("UBOT", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Lowest model level zonal wind")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float V250(time, ncol) */
    DEF_VAR("V250", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Meridional wind at 250 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float V850(time, ncol) */
    DEF_VAR("V850", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Meridional wind at 850 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float VBOT(time, ncol) */
    DEF_VAR("VBOT", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Lowest model level meridional wind")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    /* float Z500(time, ncol) */
    DEF_VAR("Z500", NC_FLOAT, 2, rec_D[1])
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "Geopotential Z at 500 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")
    PUT_ATTR_DECOMP(two, 2, dimids_D[1])

    assert (varid - varids + 1 == cfg.nvars + nvars_decomp);

err_out:
    return err;
}

