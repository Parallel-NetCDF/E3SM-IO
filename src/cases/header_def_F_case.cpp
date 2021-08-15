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
#include <e3sm_io_err.h>
#include <e3sm_io_case.hpp>
#include <e3sm_io_driver.hpp>

#if 0
/*-----< define_var() >------------------------------------------------------*/
static
int define_var(int              ncid,
               std::string      name,
               nc_type          xtype,
               int              can_ndims,  /* canonical num dims */
               int             *can_dimids, /* canonical dim IDs */
               int             *varid,      /* variable ID */
               e3sm_io_config  &cfg,
               e3sm_io_decom   &decom,
               e3sm_io_driver  &driver,
               int              decomid,    /* decomposition map ID */
               int              dim_time,
               int              orig_dimids[3][3],
               MPI_Datatype     itype,      /* buffer data type in memory */
               var_meta        *vars,
               io_buffers      *wr_buf)
{
    int k, err, isRecVar, ndims, *dimids;
    MPI_Offset varlen, dimlen;

    /* check if this is a record variable */
    isRecVar = (can_ndims > 0 && can_dimids[0] == dim_time) ? 1 : 0;

    if (decomid >= 0 && cfg.strategy == blob && cfg.api != adios) {
        /* blob is 1D for fixed-size and 2D for record variables */
        ndims = (can_dimids[0] == dim_time) ? 2 : 1;
        dimids = (isRecVar) ? orig_dimids[decomid] : orig_dimids[decomid]+1;
    }
    else { /* use the canonical dimensions */
        ndims  = can_ndims;
        dimids = can_dimids;
    }

    /* define the variable */
    err = driver.def_var(ncid, name.c_str(), xtype, ndims, dimids, varid);
    if (err != 0) {
        printf("Error in %s line %d: def_var %s\n", __FILE__, __LINE__, name.c_str());
        return -1;
    }

    /* add blob attributes */
    if (cfg.api == adios) {
    }
    else if (cfg.strategy == blob && decomid >= 0) {
        /* add two blob I/O attributes */
        int id = decomid+1;
        err = driver.put_att(ncid, *varid, "decomposition_ID", NC_INT, 1, &id);
        CHECK_VAR_ERR(*varid)
        err = driver.put_att(ncid, *varid, "global_dimids", NC_INT, can_ndims,
                             can_dimids);
        CHECK_VAR_ERR(*varid)
    }

    if (decomid == -1) { /* calculate the entire variable size */
        for (varlen=1, k=0; k<ndims; k++) {
            err = driver.inq_dimlen(ncid, can_dimids[k], &dimlen);
            CHECK_ERR
            varlen *= dimlen;
        }
    }
    else /* write request size has been precalculated */
        varlen = decom.count[decomid];

    /* add metadata of this variable */
    vars[*varid].vid       = *varid;
    vars[*varid].itype     = itype;
    vars[*varid].decomp_id = decomid;
    vars[*varid].isRecVar  = isRecVar;
    vars[*varid].vlen      = varlen;

    /* accumulate the write buffer size */
    if (itype == MPI_DOUBLE) {
#ifdef _DOUBLE_TYPE_
        if (isRecVar)
            wr_buf->rec_buflen += varlen + wr_buf->gap;
        else
            wr_buf->fix_buflen += varlen + wr_buf->gap;
#else
        if (isRecVar)
            wr_buf->rec_dbl_buflen += varlen + wr_buf->gap;
        else
            wr_buf->fix_dbl_buflen += varlen + wr_buf->gap;
#endif
    }
    else if (itype == MPI_INT) {
        if (isRecVar)
            wr_buf->rec_int_buflen += varlen + wr_buf->gap;
        else
            wr_buf->fix_int_buflen += varlen + wr_buf->gap;
    }
    else if (itype == MPI_CHAR) {
        if (isRecVar)
            wr_buf->rec_txt_buflen += varlen + wr_buf->gap;
        else
            wr_buf->fix_txt_buflen += varlen + wr_buf->gap;
    }
err_out:
    return err;
}
#endif

/*----< add_gattrs() >-------------------------------------------------------*/
static
int add_gattrs(e3sm_io_config &cfg,
               e3sm_io_decom  &decom,
               e3sm_io_driver &driver,
               int             ncid)
{
    int err, nprocs;

#define prefix

    if (cfg.strategy == blob && (cfg.api == pnetcdf || cfg.api == hdf5)) {
        MPI_Comm_size(cfg.io_comm, &nprocs);
        PUT_GATTR_INT("global_nprocs", nprocs)
        PUT_GATTR_INT("num_decompositions", decom.num_decomp)
        PUT_GATTR_INT("num_subfiles", cfg.num_subfiles)
    }

    /* save number of processes as global attributes */
    PUT_GATTR_INT("ne", 120)
    PUT_GATTR_INT("np", 21600)
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

    if (cfg.hist == h0)
        PUT_GATTR_TXT("time_period_freq", "day_5")
    else
        PUT_GATTR_TXT("time_period_freq", "hour_1")

err_out:
    return err;
}

/*----< def_var_decomp() >---------------------------------------------------*/
int e3sm_io_case::def_var_decomp(e3sm_io_config &cfg,
                                 e3sm_io_decom  &decom,
                                 e3sm_io_driver &driver,
                                 int             ncid,
                                 int             dim_time,
                                 int             dim_nblobs,
                                 int             dim_max_nreqs[MAX_NUM_DECOMP],
                                 int             orig_dimids[MAX_NUM_DECOMP][4])
{
    char name[128];
    int i, j, err=0, dimids[2];

    if (cfg.api == adios ) {
        for (j=0; j<decom.num_decomp; j++) {
            int ival, piodims[MAX_NUM_DECOMP];
            sprintf (name, "/__pio__/decomp/%d", (j + 512));
            err = driver.def_local_var(ncid, name, NC_INT64, 1,
                                       decom.raw_nreqs+j, &vars[j].vid);
            CHECK_ERR

            for (i=0; i<decom.ndims[j]; i++)
                 piodims[i] = (int)decom.dims[j][i];

            err = driver.put_att(ncid, vars[j].vid, "dimlen", NC_INT,
                                 decom.ndims[j], piodims);
            CHECK_ERR
            err = driver.put_att(ncid, vars[j].vid, "ndims", NC_INT, 1,
                                 decom.ndims+j);
            CHECK_ERR

            ival = 6;
            err = driver.put_att(ncid, vars[j].vid, "piotype", NC_INT, 1,
                                 &ival);
            CHECK_ERR
        }
        err = driver.def_local_var(ncid, "/__pio__/info/nproc", NC_INT, 0,
                                   NULL, &vars[j].vid);
        CHECK_ERR
    }
    else {
        std::map<int, std::string> dnames;
        var_meta *varp = vars - 1;

        for (i=0; i<decom.num_decomp; i++) {
            dimids[0] = dim_nblobs;
            dimids[1] = dim_max_nreqs[i];
            sprintf(name, "D%d.nreqs", i+1);
            DEF_VAR(name, NC_INT, 1, dimids, -1)
            PUT_ATTR_TXT("description", "Number of noncontiguous requests per blob")
            PUT_ATTR_INT("global_dimids", decom.ndims[i], orig_dimids[i]+1)

            sprintf(name, "D%d.blob_start", i+1);
            DEF_VAR(name, NC_INT64, 1, dimids, -1)
            PUT_ATTR_TXT("description", "Starting variable array index per blob")

            sprintf(name, "D%d.blob_count", i+1);
            DEF_VAR(name, NC_INT64, 1, dimids, -1)
            PUT_ATTR_TXT("description", "Number of variable array elements per blob")

            sprintf(name, "D%d.offsets", i+1);
            DEF_VAR(name, NC_INT, 2, dimids, -1)
            PUT_ATTR_TXT("description", "Flattened starting indices of noncontiguous requests")

            sprintf(name, "D%d.lengths", i+1);
            DEF_VAR(name, NC_INT, 2, dimids, -1)
            PUT_ATTR_TXT("description", "Lengths of noncontiguous requests")
        }
    }

err_out:
    return err;
}

/*----< def_F_case() >-------------------------------------------------------*/
int e3sm_io_case::def_F_case(e3sm_io_config &cfg,
                             e3sm_io_decom  &decom,
                             e3sm_io_driver &driver,
                             int             ncid)    /* file ID */
{
    int i, err, nprocs, dimids[3], mdims=1, nvars_decomp;
    int dim_ncol, dim_time, dim_nbnd, dim_chars, dim_lev, dim_ilev, dim_nblobs;
    int dim_nelems[MAX_NUM_DECOMP], dim_max_nreqs[MAX_NUM_DECOMP];
    int orig_dimids[MAX_NUM_DECOMP][4];
    int fix_dimids[MAX_NUM_DECOMP][3], rec_dimids[MAX_NUM_DECOMP][3];
    float fillv = 1.e+20f, missv = 1.e+20f;
    std::map<int, std::string> dnames;
    MPI_Offset nbnd, nchars, lev, ilev;
    var_meta *varp;

    /* add global attributes */
    err = add_gattrs(cfg, decom, driver, ncid);
    CHECK_ERR

    /* define dimensions */

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
    DEF_DIM("time",  NC_UNLIMITED,         &dim_time)
    DEF_DIM("nbnd",  2,                    &dim_nbnd)
    DEF_DIM("chars", 8,                    &dim_chars)
    DEF_DIM("lev",   decom.dims[2][0],     &dim_lev)
    DEF_DIM("ilev",  decom.dims[2][0] + 1, &dim_ilev)
    DEF_DIM("ncol",  decom.dims[2][1],     &dim_ncol)

    if (cfg.strategy == blob && cfg.api != adios) {
        for (i=0; i<decom.num_decomp; i++)
            orig_dimids[i][0] = dim_time;

        orig_dimids[0][1] = dim_ncol;
        orig_dimids[1][1] = dim_ncol;
        orig_dimids[2][1] = dim_lev;  orig_dimids[2][2] = dim_ncol;

        MPI_Comm_size(cfg.sub_comm, &nprocs);
        DEF_DIM("nblobs", (MPI_Offset)nprocs, &dim_nblobs)

        DEF_DIM("D1.nelems", decom.nelems[0], &dim_nelems[0])
        DEF_DIM("D2.nelems", decom.nelems[1], &dim_nelems[1])
        DEF_DIM("D3.nelems", decom.nelems[2], &dim_nelems[2])

        DEF_DIM("D1.max_nreqs", decom.max_nreqs[0], &dim_max_nreqs[0])
        DEF_DIM("D2.max_nreqs", decom.max_nreqs[1], &dim_max_nreqs[1])
        DEF_DIM("D3.max_nreqs", decom.max_nreqs[2], &dim_max_nreqs[2])

        for (i=0; i<decom.num_decomp; i++) {
            fix_dimids[i][0] = dim_nelems[i];
            rec_dimids[i][1] = dim_nelems[i];
        }
    }
    else { /* canonical or ADIOS blob */
        fix_dimids[0][0] = dim_ncol;
        fix_dimids[1][0] = dim_ncol;
        fix_dimids[2][0] = dim_lev;  fix_dimids[2][1] = dim_ncol;
        rec_dimids[0][1] = dim_ncol;
        rec_dimids[1][1] = dim_ncol;
        rec_dimids[2][1] = dim_lev;  rec_dimids[2][2] = dim_ncol;
    }

    for (i=0; i<decom.num_decomp; i++)
        rec_dimids[i][0] = dim_time;

    /* define variables related to decompositions */
    nvars_decomp = 0;
    if (cfg.strategy == blob) {
        if (cfg.api == adios)
            nvars_decomp = decom.num_decomp + 1;
        else
            nvars_decomp = NVARS_DECOMP * decom.num_decomp;

        err = def_var_decomp(cfg, decom, driver, ncid, dim_time, dim_nblobs,
                             dim_max_nreqs, orig_dimids);
        CHECK_ERR
    }

    /* Total 414 climate variables (15 fixed-size and 399 record variables) */
    varp = vars + nvars_decomp - 1;

    /* below 10 are fixed-size climate variables ---------------------------*/

    nbnd   = 2;
    nchars = 8;
    lev    = decom.dims[2][0];
    ilev   = decom.dims[2][0] + 1;

    /* combine DEF_VAR(), PUT_ATTR_DECOMP(), and SET_VAR_META() into a function
     * call.
    err = define_var(ncid, "lat", NC_DOUBLE, 1, fix_dimids[1], &varid, cfg, decom,
          driver, 1, dim_time, orig_dimids, MPI_DOUBLE, vars, &wr_buf);
     */

    /* double lat(ncol) */
    DEF_VAR("lat", NC_DOUBLE, 1, fix_dimids[1], 1)
    PUT_ATTR_TXT("long_name", "latitude")
    PUT_ATTR_TXT("units", "degrees_north")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, fix_buflen,     decom.count[1])
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, decom.count[1])
#endif

    /* double lon(ncol) */
    DEF_VAR("lon", NC_DOUBLE, 1, fix_dimids[1], 1)
    PUT_ATTR_TXT("long_name", "longitude")
    PUT_ATTR_TXT("units", "degrees_east")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, fix_buflen, decom.count[1])
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, decom.count[1])
#endif

    /* double area(ncol) */
    DEF_VAR("area", NC_DOUBLE, 1, fix_dimids[0], 0)
    PUT_ATTR_TXT("long_name", "gll grid areas")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, fix_buflen, decom.count[0])
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, decom.count[0])
#endif

    /* double lev(lev) */
    DEF_VAR("lev", NC_DOUBLE, 1, &dim_lev, -1)
    PUT_ATTR_TXT("long_name", "hybrid level at midpoints (1000*(A+B))")
    PUT_ATTR_TXT("units", "hPa")
    PUT_ATTR_TXT("positive", "down")
    PUT_ATTR_TXT("standard_name", "atmosphere_hybrid_sigma_pressure_coordinate")
    PUT_ATTR_TXT("formula_terms", "a: hyam b: hybm p0: P0 ps: PS")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, fix_buflen, lev)
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, lev)
#endif

    /* double hyam(lev) */
    DEF_VAR("hyam", NC_DOUBLE, 1, &dim_lev, -1)
    PUT_ATTR_TXT("long_name", "hybrid A coefficient at layer midpoints")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, fix_buflen, lev)
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, lev)
#endif

    /* double hybm(lev) */
    DEF_VAR("hybm", NC_DOUBLE, 1, &dim_lev, -1)
    PUT_ATTR_TXT("long_name", "hybrid B coefficient at layer midpoints")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, fix_buflen, lev)
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, lev)
#endif

    /* double P0 */
    DEF_VAR("P0", NC_DOUBLE, 0, NULL, -1)
    PUT_ATTR_TXT("long_name", "reference pressure")
    PUT_ATTR_TXT("units", "Pa")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, fix_buflen, 1)
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, 1)
#endif

    /* double ilev(ilev) */
    DEF_VAR("ilev", NC_DOUBLE, 1, &dim_ilev, -1)
    PUT_ATTR_TXT("long_name", "hybrid level at interfaces (1000*(A+B))")
    PUT_ATTR_TXT("units", "hPa")
    PUT_ATTR_TXT("positive", "down")
    PUT_ATTR_TXT("standard_name", "atmosphere_hybrid_sigma_pressure_coordinate")
    PUT_ATTR_TXT("formula_terms", "a: hyai b: hybi p0: P0 ps: PS")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, fix_buflen, ilev)
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, ilev)
#endif

    /* double hyai(ilev) */
    DEF_VAR("hyai", NC_DOUBLE, 1, &dim_ilev, -1)
    PUT_ATTR_TXT("long_name", "hybrid A coefficient at layer interfaces")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, fix_buflen, ilev)
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, ilev)
#endif

    /* double hybi(ilev) */
    DEF_VAR("hybi", NC_DOUBLE, 1, &dim_ilev, -1)
    PUT_ATTR_TXT("long_name", "hybrid B coefficient at layer interfaces")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, fix_buflen, ilev)
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, ilev)
#endif

    /* below 6 are record climate variables --------------------------------*/

    /* double time(time) */
    DEF_VAR("time", NC_DOUBLE, 1, &dim_time, -1)
    PUT_ATTR_TXT("long_name", "time")
    PUT_ATTR_TXT("units", "days since 0001-01-01 00:00:00")
    PUT_ATTR_TXT("calendar", "noleap")
    PUT_ATTR_TXT("bounds", "time_bnds")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, rec_buflen, 1)
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, 1)
#endif

    /* int date(time) */
    DEF_VAR("date", NC_INT, 1, &dim_time, -1)
    PUT_ATTR_TXT("long_name", "current date (YYYYMMDD)")
    SET_VAR_META(MPI_INT, rec_int_buflen, 1)

    /* int datesec(time) */
    DEF_VAR("datesec", NC_INT, 1, &dim_time, -1)
    PUT_ATTR_TXT("long_name", "current seconds of current date")
    SET_VAR_META(MPI_INT, rec_int_buflen, 1)

    /* double time_bnds(time, nbnd) */
    dimids[0] = dim_time;
    dimids[1] = dim_nbnd;
    DEF_VAR("time_bnds", NC_DOUBLE, 2, dimids, -1)
    PUT_ATTR_TXT("long_name", "time interval endpoints")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, rec_buflen, nbnd)
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, nbnd)
#endif

    /* char date_written(time, chars) */
    dimids[0] = dim_time;
    dimids[1] = dim_chars;
    DEF_VAR("date_written", NC_CHAR, 2, dimids, -1)
    SET_VAR_META(MPI_CHAR, rec_txt_buflen, nchars)

    /* char time_written(time, chars) */
    dimids[0] = dim_time;
    dimids[1] = dim_chars;
    DEF_VAR("time_written", NC_CHAR, 2, dimids, -1)
    SET_VAR_META(MPI_CHAR, rec_txt_buflen, nchars)

    /* below 5 are fixed-size climate variables ----------------------------*/

    /* int ndbase */
    DEF_VAR("ndbase", NC_INT, 0, NULL, -1)
    PUT_ATTR_TXT("long_name", "base day")
    SET_VAR_META(MPI_INT, fix_int_buflen, 1)

    /* int nsbase */
    DEF_VAR("nsbase", NC_INT, 0, NULL, -1)
    PUT_ATTR_TXT("long_name", "seconds of base day")
    SET_VAR_META(MPI_INT, fix_int_buflen, 1)

    /* int nbdate */
    DEF_VAR("nbdate", NC_INT, 0, NULL, -1)
    PUT_ATTR_TXT("long_name", "base date (YYYYMMDD)")
    SET_VAR_META(MPI_INT, fix_int_buflen, 1)

    /* int nbsec */
    DEF_VAR("nbsec", NC_INT, 0, NULL, -1)
    PUT_ATTR_TXT("long_name", "seconds of base date")
    SET_VAR_META(MPI_INT, fix_int_buflen, 1)

    /* int mdt */
    DEF_VAR("mdt", NC_INT, 0, NULL, -1)
    PUT_ATTR_TXT("long_name", "timestep")
    PUT_ATTR_TXT("units", "s")
    SET_VAR_META(MPI_INT, fix_int_buflen, 1)

    /* below 393 are record climate variables ------------------------------*/

    /* int ndcur(time) */
    DEF_VAR("ndcur", NC_INT, 1, &dim_time, -1)
    PUT_ATTR_TXT("long_name", "current day (from base day)")
    SET_VAR_META(MPI_INT, rec_int_buflen, 1)

    /* int nscur(time) */
    DEF_VAR("nscur", NC_INT, 1, &dim_time, -1)
    PUT_ATTR_TXT("long_name", "current seconds of current day")
    SET_VAR_META(MPI_INT, rec_int_buflen, 1)

    /* double co2vmr(time) */
    DEF_VAR("co2vmr", NC_DOUBLE, 1, &dim_time, -1)
    PUT_ATTR_TXT("long_name", "co2 volume mixing ratio")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, rec_buflen, 1)
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, 1)
#endif

    /* double ch4vmr(time) */
    DEF_VAR("ch4vmr", NC_DOUBLE, 1, &dim_time, -1)
    PUT_ATTR_TXT("long_name", "ch4 volume mixing ratio")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, rec_buflen, 1)
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, 1)
#endif

    /* double n2ovmr(time) */
    DEF_VAR("n2ovmr", NC_DOUBLE, 1, &dim_time, -1)
    PUT_ATTR_TXT("long_name", "n2o volume mixing ratio")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, rec_buflen, 1)
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, 1)
#endif

    /* double f11vmr(time) */
    DEF_VAR("f11vmr", NC_DOUBLE, 1, &dim_time, -1)
    PUT_ATTR_TXT("long_name", "f11 volume mixing ratio")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, rec_buflen, 1)
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, 1)
#endif

    /* double f12vmr(time) */
    DEF_VAR("f12vmr", NC_DOUBLE, 1, &dim_time, -1)
    PUT_ATTR_TXT("long_name", "f12 volume mixing ratio")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, rec_buflen, 1)
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, 1)
#endif

    /* double sol_tsi(time) */
    DEF_VAR("sol_tsi", NC_DOUBLE, 1, &dim_time, -1)
    PUT_ATTR_TXT("long_name", "total solar irradiance")
    PUT_ATTR_TXT("units", "W/m2")
#ifdef _DOUBLE_TYPE_
    SET_VAR_META(MPI_DOUBLE, rec_buflen, 1)
#else
    SET_VAR_META(MPI_DOUBLE, fix_dbl_buflen, 1)
#endif

    /* int nsteph(time) */
    DEF_VAR("nsteph", NC_INT, 1, &dim_time, -1)
    PUT_ATTR_TXT("long_name", "current timestep")
    SET_VAR_META(MPI_INT, rec_int_buflen, 1)

if (cfg.hist == h0) {
    /* float AEROD_v(time, ncol) */
    DEF_VAR("AEROD_v", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Total Aerosol Optical Depth in visible band")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ANRAIN(time, lev, ncol) */
    DEF_VAR("ANRAIN", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m-3")
    PUT_ATTR_TXT("long_name", "Average rain number conc")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float ANSNOW(time, lev, ncol) */
    DEF_VAR("ANSNOW", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m-3")
    PUT_ATTR_TXT("long_name", "Average snow number conc")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float AODABS(time, ncol) */
    DEF_VAR("AODABS", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol absorption optical depth 550 nm")
    PUT_ATTR_TXT("standard_name", "atmosphere_absorption_optical_thickness_due_to_ambient_aerosol_particles")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODABSBC(time, ncol) */
    DEF_VAR("AODABSBC", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol absorption optical depth 550 nm from BC")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODALL(time, ncol) */
    DEF_VAR("AODALL", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "AOD 550 nm for all time and species")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODBC(time, ncol) */
    DEF_VAR("AODBC", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from BC")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODDUST(time, ncol) */
    DEF_VAR("AODDUST", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from dust")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODDUST1(time, ncol) */
    DEF_VAR("AODDUST1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm model 1 from dust")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODDUST3(time, ncol) */
    DEF_VAR("AODDUST3", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm model 3 from dust")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODDUST4(time, ncol) */
    DEF_VAR("AODDUST4", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm model 4 from dust")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODMODE1(time, ncol) */
    DEF_VAR("AODMODE1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm mode 1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODMODE2(time, ncol) */
    DEF_VAR("AODMODE2", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm mode 2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODMODE3(time, ncol) */
    DEF_VAR("AODMODE3", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm mode 3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODMODE4(time, ncol) */
    DEF_VAR("AODMODE4", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm mode 4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODNIR(time, ncol) */
    DEF_VAR("AODNIR", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 850 nm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODPOM(time, ncol) */
    DEF_VAR("AODPOM", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from POM")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODSO4(time, ncol) */
    DEF_VAR("AODSO4", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from SO4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODSOA(time, ncol) */
    DEF_VAR("AODSOA", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from SOA")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODSS(time, ncol) */
    DEF_VAR("AODSS", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from seasalt")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODUV(time, ncol) */
    DEF_VAR("AODUV", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 350 nm")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AODVIS(time, ncol) */
    DEF_VAR("AODVIS", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm")
    PUT_ATTR_TXT("standard_name", "atmosphere_optical_thickness_due_to_ambient_aerosol_particles")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AQRAIN(time, lev, ncol) */
    DEF_VAR("AQRAIN", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "Average rain mixing ratio")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float AQSNOW(time, lev, ncol) */
    DEF_VAR("AQSNOW", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "Average snow mixing ratio")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float AQ_DMS(time, ncol) */
    DEF_VAR("AQ_DMS", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "DMS aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AQ_H2O2(time, ncol) */
    DEF_VAR("AQ_H2O2", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "H2O2 aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AQ_H2SO4(time, ncol) */
    DEF_VAR("AQ_H2SO4", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "H2SO4 aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AQ_O3(time, ncol) */
    DEF_VAR("AQ_O3", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "O3 aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AQ_SO2(time, ncol) */
    DEF_VAR("AQ_SO2", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "SO2 aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AQ_SOAG(time, ncol) */
    DEF_VAR("AQ_SOAG", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "SOAG aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float AREI(time, lev, ncol) */
    DEF_VAR("AREI", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "Micron")
    PUT_ATTR_TXT("long_name", "Average ice effective radius")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float AREL(time, lev, ncol) */
    DEF_VAR("AREL", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "Micron")
    PUT_ATTR_TXT("long_name", "Average droplet effective radius")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float AWNC(time, lev, ncol) */
    DEF_VAR("AWNC", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m-3")
    PUT_ATTR_TXT("long_name", "Average cloud water number conc")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float AWNI(time, lev, ncol) */
    DEF_VAR("AWNI", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m-3")
    PUT_ATTR_TXT("long_name", "Average cloud ice number conc")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float BURDEN1(time, ncol) */
    DEF_VAR("BURDEN1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Aerosol burden mode 1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float BURDEN2(time, ncol) */
    DEF_VAR("BURDEN2", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Aerosol burden mode 2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float BURDEN3(time, ncol) */
    DEF_VAR("BURDEN3", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Aerosol burden mode 3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float BURDEN4(time, ncol) */
    DEF_VAR("BURDEN4", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Aerosol burden mode 4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float CCN3(time, lev, ncol) */
    DEF_VAR("CCN3", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/cm3")
    PUT_ATTR_TXT("long_name", "CCN concentration at S=0.1%")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float CDNUMC(time, ncol) */
    DEF_VAR("CDNUMC", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "1/m2")
    PUT_ATTR_TXT("long_name", "Vertically-integrated droplet concentration")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])
}
    /* float CLDHGH(time, ncol) */
    DEF_VAR("CLDHGH", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Vertically-integrated high cloud")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

if (cfg.hist == h0) {
    /* float CLDICE(time, lev, ncol) */
    DEF_VAR("CLDICE", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged cloud ice amount")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float CLDLIQ(time, lev, ncol) */
    DEF_VAR("CLDLIQ", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged cloud liquid amount")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])
}
    /* float CLDLOW(time, ncol) */
    DEF_VAR("CLDLOW", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Vertically-integrated low cloud")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float CLDMED(time, ncol) */
    DEF_VAR("CLDMED", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Vertically-integrated mid-level cloud")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

if (cfg.hist == h0) {
    /* float CLDTOT(time, ncol) */
    DEF_VAR("CLDTOT", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Vertically-integrated total cloud")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float CLOUD(time, lev, ncol) */
    DEF_VAR("CLOUD", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Cloud fraction")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float CLOUDFRAC_CLUBB(time, lev, ncol) */
    DEF_VAR("CLOUDFRAC_CLUBB", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Cloud Fraction")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float CONCLD(time, lev, ncol) */
    DEF_VAR("CONCLD", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "fraction")
    PUT_ATTR_TXT("long_name", "Convective cloud cover")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float DCQ(time, lev, ncol) */
    DEF_VAR("DCQ", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg/s")
    PUT_ATTR_TXT("long_name", "Q tendency due to moist processes")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float DF_DMS(time, ncol) */
    DEF_VAR("DF_DMS", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float DF_H2O2(time, ncol) */
    DEF_VAR("DF_H2O2", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float DF_H2SO4(time, ncol) */
    DEF_VAR("DF_H2SO4", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float DF_O3(time, ncol) */
    DEF_VAR("DF_O3", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float DF_SO2(time, ncol) */
    DEF_VAR("DF_SO2", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float DF_SOAG(time, ncol) */
    DEF_VAR("DF_SOAG", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float DMS_SRF(time, ncol) */
    DEF_VAR("DMS_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "DMS in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float DP_KCLDBASE(time, ncol) */
    DEF_VAR("DP_KCLDBASE", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Deep conv. cloudbase level index")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float DP_MFUP_MAX(time, ncol) */
    DEF_VAR("DP_MFUP_MAX", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Deep conv. column-max updraft mass flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float DP_WCLDBASE(time, ncol) */
    DEF_VAR("DP_WCLDBASE", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Deep conv. cloudbase vertical velocity")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float DSTSFMBL(time, ncol) */
    DEF_VAR("DSTSFMBL", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Mobilization flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float DTCOND(time, lev, ncol) */
    DEF_VAR("DTCOND", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "K/s")
    PUT_ATTR_TXT("long_name", "T tendency - moist processes")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float DTENDTH(time, ncol) */
    DEF_VAR("DTENDTH", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Dynamic Tendency of Total (vertically integrated) moist static energy")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float DTENDTQ(time, ncol) */
    DEF_VAR("DTENDTQ", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Dynamic Tendency of Total (vertically integrated) specific humidity")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float EXTINCT(time, lev, ncol) */
    DEF_VAR("EXTINCT", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "/m")
    PUT_ATTR_TXT("long_name", "Aerosol extinction")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float FICE(time, lev, ncol) */
    DEF_VAR("FICE", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fractional ice content within cloud")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float FLDS(time, ncol) */
    DEF_VAR("FLDS", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Downwelling longwave flux at surface")
    PUT_ATTR_TXT("standard_name", "surface_downwelling_longwave_flux_in_air")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float FLNS(time, ncol) */
    DEF_VAR("FLNS", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Net longwave flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float FLNSC(time, ncol) */
    DEF_VAR("FLNSC", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky net longwave flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])
}
    /* float FLNT(time, ncol) */
    DEF_VAR("FLNT", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Net longwave flux at top of model")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

if (cfg.hist == h0) {
    /* float FLNTC(time, ncol) */
    DEF_VAR("FLNTC", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky net longwave flux at top of model")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float FLUT(time, ncol) */
    DEF_VAR("FLUT", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Upwelling longwave flux at top of model")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float FLUTC(time, ncol) */
    DEF_VAR("FLUTC", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky upwelling longwave flux at top of model")
    PUT_ATTR_TXT("standard_name", "toa_outgoing_longwave_flux_assuming_clear_sky")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float FREQI(time, lev, ncol) */
    DEF_VAR("FREQI", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fractional occurrence of ice")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float FREQL(time, lev, ncol) */
    DEF_VAR("FREQL", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fractional occurrence of liquid")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float FREQR(time, lev, ncol) */
    DEF_VAR("FREQR", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fractional occurrence of rain")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float FREQS(time, lev, ncol) */
    DEF_VAR("FREQS", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fractional occurrence of snow")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float FSDS(time, ncol) */
    DEF_VAR("FSDS", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Downwelling solar flux at surface")
    PUT_ATTR_TXT("standard_name", "surface_downwelling_shortwave_flux_in_air")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float FSDSC(time, ncol) */
    DEF_VAR("FSDSC", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky downwelling solar flux at surface")
    PUT_ATTR_TXT("standard_name", "surface_downwelling_shortwave_flux_in_air_assuming_clear_sky")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float FSNS(time, ncol) */
    DEF_VAR("FSNS", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Net solar flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float FSNSC(time, ncol) */
    DEF_VAR("FSNSC", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky net solar flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float FSNT(time, ncol) */
    DEF_VAR("FSNT", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Net solar flux at top of model")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float FSNTC(time, ncol) */
    DEF_VAR("FSNTC", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky net solar flux at top of model")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float FSNTOA(time, ncol) */
    DEF_VAR("FSNTOA", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Net solar flux at top of atmosphere")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float FSNTOAC(time, ncol) */
    DEF_VAR("FSNTOAC", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky net solar flux at top of atmosphere")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float FSUTOA(time, ncol) */
    DEF_VAR("FSUTOA", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Upwelling solar flux at top of atmosphere")
    PUT_ATTR_TXT("standard_name", "toa_outgoing_shortwave_flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float FSUTOAC(time, ncol) */
    DEF_VAR("FSUTOAC", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky upwelling solar flux at top of atmosphere")
    PUT_ATTR_TXT("standard_name", "toa_outgoing_shortwave_flux_assuming_clear_sky")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float F_eff(time, ncol) */
    DEF_VAR("F_eff", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Effective enrichment factor of marine organic matter")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float H2O2_SRF(time, ncol) */
    DEF_VAR("H2O2_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "H2O2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float H2SO4_SRF(time, ncol) */
    DEF_VAR("H2SO4_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "H2SO4 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float H2SO4_sfgaex1(time, ncol) */
    DEF_VAR("H2SO4_sfgaex1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "H2SO4 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ICEFRAC(time, ncol) */
    DEF_VAR("ICEFRAC", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fraction of sfc area covered by sea-ice")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ICIMR(time, lev, ncol) */
    DEF_VAR("ICIMR", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "Prognostic in-cloud ice mixing ratio")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float ICWMR(time, lev, ncol) */
    DEF_VAR("ICWMR", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "Prognostic in-cloud water mixing ratio")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float IWC(time, lev, ncol) */
    DEF_VAR("IWC", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/m3")
    PUT_ATTR_TXT("long_name", "Grid box average ice water content")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float LANDFRAC(time, ncol) */
    DEF_VAR("LANDFRAC", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fraction of sfc area covered by land")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float LHFLX(time, ncol) */
    DEF_VAR("LHFLX", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Surface latent heat flux")
    PUT_ATTR_TXT("standard_name", "surface_upward_latent_heat_flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float LINOZ_DO3(time, lev, ncol) */
    DEF_VAR("LINOZ_DO3", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/s")
    PUT_ATTR_TXT("long_name", "ozone vmr tendency by linearized ozone chemistry")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float LINOZ_DO3_PSC(time, lev, ncol) */
    DEF_VAR("LINOZ_DO3_PSC", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/s")
    PUT_ATTR_TXT("long_name", "ozone vmr loss by PSCs using Carille et al. (1990)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float LINOZ_O3CLIM(time, lev, ncol) */
    DEF_VAR("LINOZ_O3CLIM", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "climatology of ozone in LINOZ")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float LINOZ_O3COL(time, lev, ncol) */
    DEF_VAR("LINOZ_O3COL", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "DU")
    PUT_ATTR_TXT("long_name", "ozone column above")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float LINOZ_SFCSINK(time, ncol) */
    DEF_VAR("LINOZ_SFCSINK", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "Tg/yr/m2")
    PUT_ATTR_TXT("long_name", "surface o3 sink in LINOZ with an e-fold to a fixed concentration")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float LINOZ_SSO3(time, lev, ncol) */
    DEF_VAR("LINOZ_SSO3", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg")
    PUT_ATTR_TXT("long_name", "steady state ozone in LINOZ")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float LINOZ_SZA(time, ncol) */
    DEF_VAR("LINOZ_SZA", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "degrees")
    PUT_ATTR_TXT("long_name", "solar zenith angle in LINOZ")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float LND_MBL(time, ncol) */
    DEF_VAR("LND_MBL", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Soil erodibility factor")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])
}
    /* float LWCF(time, ncol) */
    DEF_VAR("LWCF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Longwave cloud forcing")
    PUT_ATTR_TXT("standard_name", "toa_longwave_cloud_radiative_effect")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

if (cfg.hist == h0) {
    /* float Mass_bc(time, lev, ncol) */
    DEF_VAR("Mass_bc", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of bc mass concentration bc_a1+bc_c1+bc_a3+bc_c3+bc_a4+bc_c4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float Mass_dst(time, lev, ncol) */
    DEF_VAR("Mass_dst", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of dst mass concentration dst_a1+dst_c1+dst_a3+dst_c3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float Mass_mom(time, lev, ncol) */
    DEF_VAR("Mass_mom", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of mom mass concentration mom_a1+mom_c1+mom_a2+mom_c2+mom_a3+mom_c3+mom_a4+mom_c4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float Mass_ncl(time, lev, ncol) */
    DEF_VAR("Mass_ncl", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of ncl mass concentration ncl_a1+ncl_c1+ncl_a2+ncl_c2+ncl_a3+ncl_c3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float Mass_pom(time, lev, ncol) */
    DEF_VAR("Mass_pom", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of pom mass concentration pom_a1+pom_c1+pom_a3+pom_c3+pom_a4+pom_c4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float Mass_so4(time, lev, ncol) */
    DEF_VAR("Mass_so4", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of so4 mass concentration so4_a1+so4_c1+so4_a2+so4_c2+so4_a3+so4_c3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float Mass_soa(time, lev, ncol) */
    DEF_VAR("Mass_soa", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of soa mass concentration soa_a1+soa_c1+soa_a2+soa_c2+soa_a3+soa_c3")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float NUMICE(time, lev, ncol) */
    DEF_VAR("NUMICE", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged cloud ice number")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float NUMLIQ(time, lev, ncol) */
    DEF_VAR("NUMLIQ", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged cloud liquid number")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float NUMRAI(time, lev, ncol) */
    DEF_VAR("NUMRAI", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged rain number")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float NUMSNO(time, lev, ncol) */
    DEF_VAR("NUMSNO", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged snow number")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float O3(time, lev, ncol) */
    DEF_VAR("O3", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("mixing_ratio", "dry")
    PUT_ATTR_TXT("long_name", "O3 concentration")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float O3_SRF(time, ncol) */
    DEF_VAR("O3_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "O3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float OCNFRAC(time, ncol) */
    DEF_VAR("OCNFRAC", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fraction of sfc area covered by ocean")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float OMEGA(time, lev, ncol) */
    DEF_VAR("OMEGA", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "Pa/s")
    PUT_ATTR_TXT("long_name", "Vertical velocity (pressure)")
    PUT_ATTR_TXT("standard_name", "lagrangian_tendency_of_air_pressure")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])
}

    /* float OMEGA500(time, ncol) */
    DEF_VAR("OMEGA500", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "Pa/s")
    PUT_ATTR_TXT("long_name", "Vertical velocity at 500 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

if (cfg.hist == h0) {
    /* float OMEGAT(time, lev, ncol) */
    DEF_VAR("OMEGAT", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "K Pa/s")
    PUT_ATTR_TXT("long_name", "Vertical heat flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float PBLH(time, ncol) */
    DEF_VAR("PBLH", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "PBL height")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float PHIS(time, ncol) */
    DEF_VAR("PHIS", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m2/s2")
    PUT_ATTR_TXT("long_name", "Surface geopotential")
    PUT_ATTR_TXT("cell_methods", "time: point")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float PRECC(time, ncol) */
    DEF_VAR("PRECC", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Convective precipitation rate (liq + ice)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float PRECL(time, ncol) */
    DEF_VAR("PRECL", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Large-scale (stable) precipitation rate (liq + ice)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float PRECSC(time, ncol) */
    DEF_VAR("PRECSC", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Convective snow rate (water equivalent)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float PRECSL(time, ncol) */
    DEF_VAR("PRECSL", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Large-scale (stable) snow rate (water equivalent)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])
}
else {
    /* float OMEGA850(time, ncol) */
    DEF_VAR("OMEGA850", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "Pa/s")
    PUT_ATTR_TXT("long_name", "Vertical velocity at 850 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float PRECT(time, ncol) */
    DEF_VAR("PRECT", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Total (convective and large-scale) precipitation rate (liq + ice)")
    PUT_ATTR_TXT("cell_methods", "time: point")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])
}

    /* float PS(time, ncol) */
    DEF_VAR("PS", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("long_name", "Surface pressure")
    PUT_ATTR_TXT("standard_name", "surface_air_pressure")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

if (cfg.hist == h0) {
    /* float PSL(time, ncol) */
    DEF_VAR("PSL", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("long_name", "Sea level pressure")
    PUT_ATTR_TXT("standard_name", "air_pressure_at_mean_sea_level")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float Q(time, lev, ncol) */
    DEF_VAR("Q", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Specific humidity")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float QFLX(time, ncol) */
    DEF_VAR("QFLX", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Surface water flux")
    PUT_ATTR_TXT("standard_name", "water_evapotranspiration_flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float QREFHT(time, ncol) */
    DEF_VAR("QREFHT", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "Reference height humidity")
    PUT_ATTR_TXT("standard_name", "specific_humidity")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float QRL(time, lev, ncol) */
    DEF_VAR("QRL", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "K/s")
    PUT_ATTR_TXT("long_name", "Longwave heating rate")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float QRS(time, lev, ncol) */
    DEF_VAR("QRS", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "K/s")
    PUT_ATTR_TXT("long_name", "Solar heating rate")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float RAINQM(time, lev, ncol) */
    DEF_VAR("RAINQM", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged rain amount")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float RAM1(time, ncol) */
    DEF_VAR("RAM1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "frac")
    PUT_ATTR_TXT("long_name", "RAM1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float RELHUM(time, lev, ncol) */
    DEF_VAR("RELHUM", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "percent")
    PUT_ATTR_TXT("long_name", "Relative humidity")
    PUT_ATTR_TXT("standard_name", "relative_humidity")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float SFDMS(time, ncol) */
    DEF_VAR("SFDMS", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "DMS surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFH2O2(time, ncol) */
    DEF_VAR("SFH2O2", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "H2O2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFH2SO4(time, ncol) */
    DEF_VAR("SFH2SO4", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "H2SO4 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFO3(time, ncol) */
    DEF_VAR("SFO3", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "O3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFSO2(time, ncol) */
    DEF_VAR("SFSO2", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "SO2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFSOAG(time, ncol) */
    DEF_VAR("SFSOAG", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "SOAG surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFbc_a1(time, ncol) */
    DEF_VAR("SFbc_a1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFbc_a3(time, ncol) */
    DEF_VAR("SFbc_a3", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFbc_a4(time, ncol) */
    DEF_VAR("SFbc_a4", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a4 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFdst_a1(time, ncol) */
    DEF_VAR("SFdst_a1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFdst_a3(time, ncol) */
    DEF_VAR("SFdst_a3", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFmom_a1(time, ncol) */
    DEF_VAR("SFmom_a1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFmom_a2(time, ncol) */
    DEF_VAR("SFmom_a2", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFmom_a3(time, ncol) */
    DEF_VAR("SFmom_a3", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFmom_a4(time, ncol) */
    DEF_VAR("SFmom_a4", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a4 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFncl_a1(time, ncol) */
    DEF_VAR("SFncl_a1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFncl_a2(time, ncol) */
    DEF_VAR("SFncl_a2", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFncl_a3(time, ncol) */
    DEF_VAR("SFncl_a3", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFnum_a1(time, ncol) */
    DEF_VAR("SFnum_a1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFnum_a2(time, ncol) */
    DEF_VAR("SFnum_a2", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFnum_a3(time, ncol) */
    DEF_VAR("SFnum_a3", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFnum_a4(time, ncol) */
    DEF_VAR("SFnum_a4", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a4 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFpom_a1(time, ncol) */
    DEF_VAR("SFpom_a1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFpom_a3(time, ncol) */
    DEF_VAR("SFpom_a3", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFpom_a4(time, ncol) */
    DEF_VAR("SFpom_a4", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a4 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFso4_a1(time, ncol) */
    DEF_VAR("SFso4_a1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFso4_a2(time, ncol) */
    DEF_VAR("SFso4_a2", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFso4_a3(time, ncol) */
    DEF_VAR("SFso4_a3", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFsoa_a1(time, ncol) */
    DEF_VAR("SFsoa_a1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFsoa_a2(time, ncol) */
    DEF_VAR("SFsoa_a2", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SFsoa_a3(time, ncol) */
    DEF_VAR("SFsoa_a3", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SHFLX(time, ncol) */
    DEF_VAR("SHFLX", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Surface sensible heat flux")
    PUT_ATTR_TXT("standard_name", "surface_upward_sensible_heat_flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SH_KCLDBASE(time, ncol) */
    DEF_VAR("SH_KCLDBASE", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Shallow conv. cloudbase level index")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SH_MFUP_MAX(time, ncol) */
    DEF_VAR("SH_MFUP_MAX", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Shallow conv. column-max updraft mass flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SH_WCLDBASE(time, ncol) */
    DEF_VAR("SH_WCLDBASE", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Shallow conv. cloudbase vertical velocity")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SNOWHICE(time, ncol) */
    DEF_VAR("SNOWHICE", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "Snow depth over ice")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SNOWHLND(time, ncol) */
    DEF_VAR("SNOWHLND", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "Water equivalent snow depth")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SNOWQM(time, lev, ncol) */
    DEF_VAR("SNOWQM", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged snow amount")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float SO2(time, lev, ncol) */
    DEF_VAR("SO2", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("mixing_ratio", "dry")
    PUT_ATTR_TXT("long_name", "SO2 concentration")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float SO2_CLXF(time, ncol) */
    DEF_VAR("SO2_CLXF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for SO2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SO2_SRF(time, ncol) */
    DEF_VAR("SO2_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "SO2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SOAG_CLXF(time, ncol) */
    DEF_VAR("SOAG_CLXF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for SOAG")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SOAG_SRF(time, ncol) */
    DEF_VAR("SOAG_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "SOAG in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SOAG_sfgaex1(time, ncol) */
    DEF_VAR("SOAG_sfgaex1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "SOAG gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SOLIN(time, ncol) */
    DEF_VAR("SOLIN", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Solar insolation")
    PUT_ATTR_TXT("standard_name", "toa_incoming_shortwave_flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SSAVIS(time, ncol) */
    DEF_VAR("SSAVIS", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol singel-scatter albedo")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SSTSFMBL(time, ncol) */
    DEF_VAR("SSTSFMBL", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Mobilization flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float SSTSFMBL_OM(time, ncol) */
    DEF_VAR("SSTSFMBL_OM", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Mobilization flux of marine organic matter at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])
}
    /* float SWCF(time, ncol) */
    DEF_VAR("SWCF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Shortwave cloud forcing")
    PUT_ATTR_TXT("standard_name", "toa_shortwave_cloud_radiative_effect")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

if (cfg.hist == h0) {
    /* float T(time, lev, ncol) */
    DEF_VAR("T", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Temperature")
    PUT_ATTR_TXT("standard_name", "air_temperature")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float TAUGWX(time, ncol) */
    DEF_VAR("TAUGWX", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "N/m2")
    PUT_ATTR_TXT("long_name", "Zonal gravity wave surface stress")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float TAUGWY(time, ncol) */
    DEF_VAR("TAUGWY", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "N/m2")
    PUT_ATTR_TXT("long_name", "Meridional gravity wave surface stress")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float TAUX(time, ncol) */
    DEF_VAR("TAUX", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "N/m2")
    PUT_ATTR_TXT("long_name", "Zonal surface stress")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float TAUY(time, ncol) */
    DEF_VAR("TAUY", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "N/m2")
    PUT_ATTR_TXT("long_name", "Meridional surface stress")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float TGCLDCWP(time, ncol) */
    DEF_VAR("TGCLDCWP", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Total grid-box cloud water path (liquid and ice)")
    PUT_ATTR_TXT("standard_name", "atmosphere_mass_content_of_cloud_condensed_water")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float TGCLDIWP(time, ncol) */
    DEF_VAR("TGCLDIWP", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Total grid-box cloud ice water path")
    PUT_ATTR_TXT("standard_name", "atmosphere_mass_content_of_cloud_ice")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float TGCLDLWP(time, ncol) */
    DEF_VAR("TGCLDLWP", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Total grid-box cloud liquid water path")
    PUT_ATTR_TXT("standard_name", "atmosphere_mass_content_of_cloud_liquid_water")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float TH7001000(time, ncol) */
    DEF_VAR("TH7001000", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Theta difference 700 mb - 1000 mb")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])
}
else {
    /* float T850(time, ncol) */
    DEF_VAR("T850", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Temperature at 850 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])
}

    /* float TMQ(time, ncol) */
    DEF_VAR("TMQ", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Total (vertically integrated) precipitable water")
    PUT_ATTR_TXT("standard_name", "atmosphere_mass_content_of_water_vapor")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

if (cfg.hist == h0) {
    /* float TREFHT(time, ncol) */
    DEF_VAR("TREFHT", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Reference height temperature")
    PUT_ATTR_TXT("standard_name", "air_temperature")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float TROP_P(time, ncol) */
    DEF_VAR("TROP_P", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("long_name", "Tropopause Pressure")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float TROP_T(time, ncol) */
    DEF_VAR("TROP_T", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_FLT(_FillValue, fillv)
    PUT_ATTR_FLT("missing_value", missv)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Tropopause Temperature")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])
}

    /* float TS(time, ncol) */
    DEF_VAR("TS", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Surface temperature (radiative)")
    PUT_ATTR_TXT("standard_name", "surface_temperature")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

if (cfg.hist == h0) {
    /* float TSMN(time, ncol) */
    DEF_VAR("TSMN", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Minimum surface temperature over output period")
    PUT_ATTR_TXT("cell_methods", "time: minimum")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float TSMX(time, ncol) */
    DEF_VAR("TSMX", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Maximum surface temperature over output period")
    PUT_ATTR_TXT("cell_methods", "time: maximum")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float TUH(time, ncol) */
    DEF_VAR("TUH", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "W/m")
    PUT_ATTR_TXT("long_name", "Total (vertically integrated) zonal MSE flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float TUQ(time, ncol) */
    DEF_VAR("TUQ", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m/s")
    PUT_ATTR_TXT("long_name", "Total (vertically integrated) zonal water flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float TVH(time, ncol) */
    DEF_VAR("TVH", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "W/m")
    PUT_ATTR_TXT("long_name", "Total (vertically integrated) meridional MSE flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float TVQ(time, ncol) */
    DEF_VAR("TVQ", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m/s")
    PUT_ATTR_TXT("long_name", "Total (vertically integrated) meridional water flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])
}

    /* float U(time, lev, ncol) */
    DEF_VAR("U", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Zonal wind")
    PUT_ATTR_TXT("standard_name", "eastward_wind")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

if (cfg.hist == h0) {
    /* float U10(time, ncol) */
    DEF_VAR("U10", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "10m wind speed")
    PUT_ATTR_TXT("standard_name", "wind_speed")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float UU(time, lev, ncol) */
    DEF_VAR("UU", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m2/s2")
    PUT_ATTR_TXT("long_name", "Zonal velocity squared")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float V(time, lev, ncol) */
    DEF_VAR("V", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Meridional wind")
    PUT_ATTR_TXT("standard_name", "northward_wind")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float VQ(time, lev, ncol) */
    DEF_VAR("VQ", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m/s kg/kg")
    PUT_ATTR_TXT("long_name", "Meridional water transport")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float VT(time, lev, ncol) */
    DEF_VAR("VT", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "K m/s")
    PUT_ATTR_TXT("long_name", "Meridional heat transport")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float VU(time, lev, ncol) */
    DEF_VAR("VU", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m2/s2")
    PUT_ATTR_TXT("long_name", "Meridional flux of zonal momentum")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float VV(time, lev, ncol) */
    DEF_VAR("VV", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m2/s2")
    PUT_ATTR_TXT("long_name", "Meridional velocity squared")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float WD_H2O2(time, ncol) */
    DEF_VAR("WD_H2O2", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/s")
    PUT_ATTR_TXT("long_name", "H2O2             wet deposition")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float WD_H2SO4(time, ncol) */
    DEF_VAR("WD_H2SO4", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/s")
    PUT_ATTR_TXT("long_name", "H2SO4            wet deposition")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float WD_SO2(time, ncol) */
    DEF_VAR("WD_SO2", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/s")
    PUT_ATTR_TXT("long_name", "SO2              wet deposition")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float WSUB(time, lev, ncol) */
    DEF_VAR("WSUB", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Diagnostic sub-grid vertical velocity")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float Z3(time, lev, ncol) */
    DEF_VAR("Z3", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "Geopotential Height (above sea level)")
    PUT_ATTR_TXT("standard_name", "geopotential_height")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float aero_water(time, lev, ncol) */
    DEF_VAR("aero_water", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "sum of aerosol water of interstitial modes wat_a1+wat_a2+wat_a3+wat_a4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float airFV(time, ncol) */
    DEF_VAR("airFV", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "frac")
    PUT_ATTR_TXT("long_name", "FV")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_a1DDF(time, ncol) */
    DEF_VAR("bc_a1DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_a1SFWET(time, ncol) */
    DEF_VAR("bc_a1SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_a1_SRF(time, ncol) */
    DEF_VAR("bc_a1_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "bc_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_a1_sfgaex1(time, ncol) */
    DEF_VAR("bc_a1_sfgaex1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_a3DDF(time, ncol) */
    DEF_VAR("bc_a3DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_a3SFWET(time, ncol) */
    DEF_VAR("bc_a3SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_a3_SRF(time, ncol) */
    DEF_VAR("bc_a3_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "bc_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_a4DDF(time, ncol) */
    DEF_VAR("bc_a4DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_a4SFWET(time, ncol) */
    DEF_VAR("bc_a4SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_a4_CLXF(time, ncol) */
    DEF_VAR("bc_a4_CLXF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for bc_a4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_a4_SRF(time, ncol) */
    DEF_VAR("bc_a4_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "bc_a4 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_a4_sfgaex1(time, ncol) */
    DEF_VAR("bc_a4_sfgaex1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a4 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_c1DDF(time, ncol) */
    DEF_VAR("bc_c1DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_c1SFWET(time, ncol) */
    DEF_VAR("bc_c1SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_c3DDF(time, ncol) */
    DEF_VAR("bc_c3DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_c3SFWET(time, ncol) */
    DEF_VAR("bc_c3SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_c4DDF(time, ncol) */
    DEF_VAR("bc_c4DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float bc_c4SFWET(time, ncol) */
    DEF_VAR("bc_c4SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c4 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float chla(time, ncol) */
    DEF_VAR("chla", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "mg L-1")
    PUT_ATTR_TXT("long_name", "ocean input data: chla")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float dst_a1DDF(time, ncol) */
    DEF_VAR("dst_a1DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float dst_a1SF(time, ncol) */
    DEF_VAR("dst_a1SF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a1 dust surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float dst_a1SFWET(time, ncol) */
    DEF_VAR("dst_a1SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float dst_a1_SRF(time, ncol) */
    DEF_VAR("dst_a1_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "dst_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float dst_a3DDF(time, ncol) */
    DEF_VAR("dst_a3DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float dst_a3SF(time, ncol) */
    DEF_VAR("dst_a3SF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a3 dust surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float dst_a3SFWET(time, ncol) */
    DEF_VAR("dst_a3SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float dst_a3_SRF(time, ncol) */
    DEF_VAR("dst_a3_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "dst_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float dst_c1DDF(time, ncol) */
    DEF_VAR("dst_c1DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float dst_c1SFWET(time, ncol) */
    DEF_VAR("dst_c1SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float dst_c3DDF(time, ncol) */
    DEF_VAR("dst_c3DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float dst_c3SFWET(time, ncol) */
    DEF_VAR("dst_c3SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float hstobie_linoz(time, lev, ncol) */
    DEF_VAR("hstobie_linoz", NC_FLOAT, 3, rec_dimids[2], 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "fraction of model time")
    PUT_ATTR_TXT("long_name", "Lowest possible Linoz level")
    PUT_ATTR_TXT("cell_methods", "time: point")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[2])

    /* float mlip(time, ncol) */
    DEF_VAR("mlip", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "uM C")
    PUT_ATTR_TXT("long_name", "ocean input data: mlip")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_a1DDF(time, ncol) */
    DEF_VAR("mom_a1DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_a1SF(time, ncol) */
    DEF_VAR("mom_a1SF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a1 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_a1SFWET(time, ncol) */
    DEF_VAR("mom_a1SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_a1_SRF(time, ncol) */
    DEF_VAR("mom_a1_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "mom_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_a1_sfgaex1(time, ncol) */
    DEF_VAR("mom_a1_sfgaex1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_a2DDF(time, ncol) */
    DEF_VAR("mom_a2DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_a2SF(time, ncol) */
    DEF_VAR("mom_a2SF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a2 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_a2SFWET(time, ncol) */
    DEF_VAR("mom_a2SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_a2_SRF(time, ncol) */
    DEF_VAR("mom_a2_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "mom_a2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_a3DDF(time, ncol) */
    DEF_VAR("mom_a3DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_a3SFWET(time, ncol) */
    DEF_VAR("mom_a3SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_a3_SRF(time, ncol) */
    DEF_VAR("mom_a3_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "mom_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_a4DDF(time, ncol) */
    DEF_VAR("mom_a4DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_a4SF(time, ncol) */
    DEF_VAR("mom_a4SF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a4 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_a4SFWET(time, ncol) */
    DEF_VAR("mom_a4SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_a4_SRF(time, ncol) */
    DEF_VAR("mom_a4_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "mom_a4 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_a4_sfgaex1(time, ncol) */
    DEF_VAR("mom_a4_sfgaex1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a4 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_c1DDF(time, ncol) */
    DEF_VAR("mom_c1DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_c1SFWET(time, ncol) */
    DEF_VAR("mom_c1SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_c2DDF(time, ncol) */
    DEF_VAR("mom_c2DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_c2SFWET(time, ncol) */
    DEF_VAR("mom_c2SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c2 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_c3DDF(time, ncol) */
    DEF_VAR("mom_c3DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_c3SFWET(time, ncol) */
    DEF_VAR("mom_c3SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_c4DDF(time, ncol) */
    DEF_VAR("mom_c4DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mom_c4SFWET(time, ncol) */
    DEF_VAR("mom_c4SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c4 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mpoly(time, ncol) */
    DEF_VAR("mpoly", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "uM C")
    PUT_ATTR_TXT("long_name", "ocean input data: mpoly")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float mprot(time, ncol) */
    DEF_VAR("mprot", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "uM C")
    PUT_ATTR_TXT("long_name", "ocean input data: mprot")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_a1DDF(time, ncol) */
    DEF_VAR("ncl_a1DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_a1SF(time, ncol) */
    DEF_VAR("ncl_a1SF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a1 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_a1SFWET(time, ncol) */
    DEF_VAR("ncl_a1SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_a1_SRF(time, ncol) */
    DEF_VAR("ncl_a1_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "ncl_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_a2DDF(time, ncol) */
    DEF_VAR("ncl_a2DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_a2SF(time, ncol) */
    DEF_VAR("ncl_a2SF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a2 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_a2SFWET(time, ncol) */
    DEF_VAR("ncl_a2SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_a2_SRF(time, ncol) */
    DEF_VAR("ncl_a2_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "ncl_a2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_a3DDF(time, ncol) */
    DEF_VAR("ncl_a3DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_a3SF(time, ncol) */
    DEF_VAR("ncl_a3SF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a3 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_a3SFWET(time, ncol) */
    DEF_VAR("ncl_a3SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_a3_SRF(time, ncol) */
    DEF_VAR("ncl_a3_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "ncl_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_c1DDF(time, ncol) */
    DEF_VAR("ncl_c1DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_c1SFWET(time, ncol) */
    DEF_VAR("ncl_c1SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_c2DDF(time, ncol) */
    DEF_VAR("ncl_c2DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_c2SFWET(time, ncol) */
    DEF_VAR("ncl_c2SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c2 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_c3DDF(time, ncol) */
    DEF_VAR("ncl_c3DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float ncl_c3SFWET(time, ncol) */
    DEF_VAR("ncl_c3SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a1DDF(time, ncol) */
    DEF_VAR("num_a1DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a1SF(time, ncol) */
    DEF_VAR("num_a1SF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "num_a1 dust surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a1SFWET(time, ncol) */
    DEF_VAR("num_a1SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a1_CLXF(time, ncol) */
    DEF_VAR("num_a1_CLXF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for num_a1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a1_SRF(time, ncol) */
    DEF_VAR("num_a1_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/kg")
    PUT_ATTR_TXT("long_name", "num_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a1_sfgaex1(time, ncol) */
    DEF_VAR("num_a1_sfgaex1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "num_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a2DDF(time, ncol) */
    DEF_VAR("num_a2DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a2SFWET(time, ncol) */
    DEF_VAR("num_a2SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a2_CLXF(time, ncol) */
    DEF_VAR("num_a2_CLXF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for num_a2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a2_SRF(time, ncol) */
    DEF_VAR("num_a2_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/kg")
    PUT_ATTR_TXT("long_name", "num_a2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a3DDF(time, ncol) */
    DEF_VAR("num_a3DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a3SF(time, ncol) */
    DEF_VAR("num_a3SF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "num_a3 dust surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a3SFWET(time, ncol) */
    DEF_VAR("num_a3SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a3_SRF(time, ncol) */
    DEF_VAR("num_a3_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/kg")
    PUT_ATTR_TXT("long_name", "num_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a4DDF(time, ncol) */
    DEF_VAR("num_a4DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a4SFWET(time, ncol) */
    DEF_VAR("num_a4SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a4_CLXF(time, ncol) */
    DEF_VAR("num_a4_CLXF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for num_a4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a4_SRF(time, ncol) */
    DEF_VAR("num_a4_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/kg")
    PUT_ATTR_TXT("long_name", "num_a4 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_a4_sfgaex1(time, ncol) */
    DEF_VAR("num_a4_sfgaex1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "num_a4 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_c1DDF(time, ncol) */
    DEF_VAR("num_c1DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_c1SFWET(time, ncol) */
    DEF_VAR("num_c1SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_c2DDF(time, ncol) */
    DEF_VAR("num_c2DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_c2SFWET(time, ncol) */
    DEF_VAR("num_c2SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c2 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_c3DDF(time, ncol) */
    DEF_VAR("num_c3DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_c3SFWET(time, ncol) */
    DEF_VAR("num_c3SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_c4DDF(time, ncol) */
    DEF_VAR("num_c4DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float num_c4SFWET(time, ncol) */
    DEF_VAR("num_c4SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c4 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_a1DDF(time, ncol) */
    DEF_VAR("pom_a1DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_a1SFWET(time, ncol) */
    DEF_VAR("pom_a1SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_a1_SRF(time, ncol) */
    DEF_VAR("pom_a1_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "pom_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_a1_sfgaex1(time, ncol) */
    DEF_VAR("pom_a1_sfgaex1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_a3DDF(time, ncol) */
    DEF_VAR("pom_a3DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_a3SFWET(time, ncol) */
    DEF_VAR("pom_a3SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_a3_SRF(time, ncol) */
    DEF_VAR("pom_a3_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "pom_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_a4DDF(time, ncol) */
    DEF_VAR("pom_a4DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_a4SFWET(time, ncol) */
    DEF_VAR("pom_a4SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_a4_CLXF(time, ncol) */
    DEF_VAR("pom_a4_CLXF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for pom_a4")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_a4_SRF(time, ncol) */
    DEF_VAR("pom_a4_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "pom_a4 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_a4_sfgaex1(time, ncol) */
    DEF_VAR("pom_a4_sfgaex1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a4 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_c1DDF(time, ncol) */
    DEF_VAR("pom_c1DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_c1SFWET(time, ncol) */
    DEF_VAR("pom_c1SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_c3DDF(time, ncol) */
    DEF_VAR("pom_c3DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_c3SFWET(time, ncol) */
    DEF_VAR("pom_c3SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_c4DDF(time, ncol) */
    DEF_VAR("pom_c4DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float pom_c4SFWET(time, ncol) */
    DEF_VAR("pom_c4SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c4 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_a1DDF(time, ncol) */
    DEF_VAR("so4_a1DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_a1SFWET(time, ncol) */
    DEF_VAR("so4_a1SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_a1_CLXF(time, ncol) */
    DEF_VAR("so4_a1_CLXF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for so4_a1")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_a1_SRF(time, ncol) */
    DEF_VAR("so4_a1_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "so4_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_a1_sfgaex1(time, ncol) */
    DEF_VAR("so4_a1_sfgaex1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_a2DDF(time, ncol) */
    DEF_VAR("so4_a2DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_a2SFWET(time, ncol) */
    DEF_VAR("so4_a2SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_a2_CLXF(time, ncol) */
    DEF_VAR("so4_a2_CLXF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for so4_a2")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_a2_SRF(time, ncol) */
    DEF_VAR("so4_a2_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "so4_a2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_a2_sfgaex1(time, ncol) */
    DEF_VAR("so4_a2_sfgaex1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a2 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_a3DDF(time, ncol) */
    DEF_VAR("so4_a3DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_a3SFWET(time, ncol) */
    DEF_VAR("so4_a3SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_a3_SRF(time, ncol) */
    DEF_VAR("so4_a3_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "so4_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_a3_sfgaex1(time, ncol) */
    DEF_VAR("so4_a3_sfgaex1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a3 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_c1DDF(time, ncol) */
    DEF_VAR("so4_c1DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_c1SFWET(time, ncol) */
    DEF_VAR("so4_c1SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_c2DDF(time, ncol) */
    DEF_VAR("so4_c2DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_c2SFWET(time, ncol) */
    DEF_VAR("so4_c2SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c2 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_c3DDF(time, ncol) */
    DEF_VAR("so4_c3DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float so4_c3SFWET(time, ncol) */
    DEF_VAR("so4_c3SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_a1DDF(time, ncol) */
    DEF_VAR("soa_a1DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_a1SFWET(time, ncol) */
    DEF_VAR("soa_a1SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_a1_SRF(time, ncol) */
    DEF_VAR("soa_a1_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "soa_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_a1_sfgaex1(time, ncol) */
    DEF_VAR("soa_a1_sfgaex1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_a2DDF(time, ncol) */
    DEF_VAR("soa_a2DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_a2SFWET(time, ncol) */
    DEF_VAR("soa_a2SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_a2_SRF(time, ncol) */
    DEF_VAR("soa_a2_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "soa_a2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_a2_sfgaex1(time, ncol) */
    DEF_VAR("soa_a2_sfgaex1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a2 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_a3DDF(time, ncol) */
    DEF_VAR("soa_a3DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_a3SFWET(time, ncol) */
    DEF_VAR("soa_a3SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_a3_SRF(time, ncol) */
    DEF_VAR("soa_a3_SRF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "soa_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_a3_sfgaex1(time, ncol) */
    DEF_VAR("soa_a3_sfgaex1", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a3 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_c1DDF(time, ncol) */
    DEF_VAR("soa_c1DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_c1SFWET(time, ncol) */
    DEF_VAR("soa_c1SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_c2DDF(time, ncol) */
    DEF_VAR("soa_c2DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_c2SFWET(time, ncol) */
    DEF_VAR("soa_c2SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c2 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_c3DDF(time, ncol) */
    DEF_VAR("soa_c3DDF", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float soa_c3SFWET(time, ncol) */
    DEF_VAR("soa_c3SFWET", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])
}
else {
    /* float U250(time, ncol) */
    DEF_VAR("U250", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Zonal wind at 250 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float U850(time, ncol) */
    DEF_VAR("U850", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Zonal wind at 850 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float UBOT(time, ncol) */
    DEF_VAR("UBOT", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Lowest model level zonal wind")
    PUT_ATTR_TXT("cell_methods", "time: point")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float V250(time, ncol) */
    DEF_VAR("V250", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Meridional wind at 250 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float V850(time, ncol) */
    DEF_VAR("V850", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Meridional wind at 850 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float VBOT(time, ncol) */
    DEF_VAR("VBOT", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Lowest model level meridional wind")
    PUT_ATTR_TXT("cell_methods", "time: point")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])

    /* float Z500(time, ncol) */
    DEF_VAR("Z500", NC_FLOAT, 2, rec_dimids[1], 1)
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "Geopotential Z at 500 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")
    SET_VAR_META(REC_ITYPE, rec_buflen, decom.count[1])
}

    assert(varp - vars + 1 == cfg.nvars + nvars_decomp);

err_out:
    return err;
}

