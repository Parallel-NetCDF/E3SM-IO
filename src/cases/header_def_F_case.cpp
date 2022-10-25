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
    int err, nprocs;

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
/* define decomposition variables. Blob I/O only */
int e3sm_io_case::def_var_decomp(e3sm_io_config &cfg,
                                 e3sm_io_decom  &decom,
                                 e3sm_io_driver &driver,
                                 case_meta      *cmeta,
                                 int ncid,
                                 int dim_time,
                                 int dim_nblobs,
                                 int dim_max_nreqs[MAX_NUM_DECOMP],
                                 int g_dimids[MAX_NUM_DECOMP][MAX_NDIMS])
{
    char name[512];
    int i, j, err=0, dimids[2];
    var_meta *varp;

    if (cfg.api == adios ) {
        MPI_Offset one = 1;
        for (j=0; j<decom.num_decomp; j++) {
            int ival, piodims[MAX_NUM_DECOMP];

            varp = vars + j;
            sprintf (name, "/__pio__/decomp/%d", (j + 512));
            err = driver.def_local_var(ncid, name, NC_INT64, 1,
                                       decom.raw_nreqs+j, &varp->vid);
            CHECK_ERR

            for (i=0; i<decom.ndims[j]; i++)
                piodims[i] = (int)decom.dims[j][i];

            PUT_ATTR_INT("dimlen", decom.ndims[j], piodims)
            PUT_ATTR_INT("ndims", 1, decom.ndims+j)
            ival = 6;
            PUT_ATTR_INT("piotype", 1, &ival)

            varp = vars + decom.num_decomp + j;
            sprintf (name, "/__pio__/track/num_decomp_block_writers/%d", (j + 512));
            err = driver.def_local_var(ncid, name, NC_INT, 1,
                                       &one, &varp->vid);
            CHECK_ERR
        }
        err = driver.def_local_var(ncid, "/__pio__/info/nproc", NC_INT, 0,
                                   NULL, &vars[2*decom.num_decomp].vid);
        CHECK_ERR
        err = driver.def_local_var(ncid, "/__pio__/info/block_nprocs", NC_INT, 0,
                                   NULL, &vars[2*decom.num_decomp+1].vid);
        CHECK_ERR
        err = driver.def_local_var(ncid, "/__pio__/info/block_list", NC_INT, 1,
                                   &one, &vars[2*decom.num_decomp+2].vid);
        CHECK_ERR
    }
    else {
        std::map<int, std::string> dnames;
        varp = vars - 1;

        for (i=0; i<decom.num_decomp; i++) {
            dimids[0] = dim_nblobs;
            dimids[1] = dim_max_nreqs[i];
            sprintf(name, "D%d.nreqs", i+1);
            DEF_VAR(name, NC_INT, 1, dimids, MPI_INT, -1)
            PUT_ATTR_TXT("description", "Number of noncontiguous requests in individual blobs")
            PUT_ATTR_INT("global_dimids", decom.ndims[i], g_dimids[i])
#ifdef GLOBAL_DIM_USE_STRING
            name[0] = '\0';
            for (j=0; j<decom.ndims[i]; j++) {
                char dim_name[128];
                ncmpi_inq_dimname(ncid, g_dimids[i][j], dim_name);
                strcat(name, dim_name);
                if (j < decom.ndims[i]-1) strcat(name,",");
            }
            PUT_ATTR_TXT("global_dims", name)
#endif

            /* For PnetCDF blob I/O, it uses variable-centric blob layout
             * strategy where there are are nprocs (or nblobs) blobs in the
             * file space occupied by each variable. Thus, given N variables
             * defined and each is partitioned among nprocs processes, there
             * will be N*nprocs blobs in the NetCDF file. A blob contains all
             * of a process's write requests to a variable.
             *
             * For HDF5 blob I/O, it uses process-centric blob layout strategy
             * where there are nprocs blobs in the entire file. A blob contains
             * all of a process's write requests to all variables.
             */

            /* For PnetCDF blob I/O, the starting file offset of a variable's
             * blob is calculated by
             *     variable's begin + D*.blob_start * variable's type size.
             * A variable's begin and type are stored in the NetCDF file
             * header.
             *
             * For HDF5 blob I/O, the calculation is more complicated, because
             * the write amount of a process to a variable can be different
             * from another process. The starting offset of a variable in a
             * blob is calculated by
             * 1. calculating the file starting offset of the blob. This
             *    requires to sum up the sizes of all blobs before it.
             * 2. calculating the starting offset of a variable inside the
             *    blob. This requires to to sum up the write amounts to all
             *    variables defined before it.
             * A less complicated approach is to additionally save the offsets
             * of all variables inside all blobs. This requires additional
             * N*nprocs integers. This is not implemented yet.
             *
             * Note NC_INT64 type is used here, but it can be NC_INT because
             * none of E3SM variables is larger than 4B elements.
             */
            sprintf(name, "D%d.blob_start", i+1);
            DEF_VAR(name, NC_INT64, 1, dimids, MPI_INT, -1)
            PUT_ATTR_TXT("description", "Starting array indices of individual blobs stored in a variable")

            sprintf(name, "D%d.blob_count", i+1);
            DEF_VAR(name, NC_INT64, 1, dimids, MPI_LONG_LONG, -1)
            PUT_ATTR_TXT("description", "Number of contiguous array elements in individual blobs stored in a variable")

            /* For PnetCDF blob I/O, the file offset of a write request to a
             * variable inside a blob is calculated by
             *     the variable blob's file offset calculated above +
             *     D*.offsets[blob_id][request_number] * variable's type size
             * A variable's begin and type are stored in NetCDF file header.
             *
             * For HDF5 blob I/O, the calculation is similar, i.e.
             *     the variable's file offset calculated above +
             *     D*.offsets[blob_id][request_number] * variable's type size
             *
             * Note NC_INT type is used here because none of E3SM variables is
             * larger than 4B elements. In general, it can be NC_INT64.
             */
            sprintf(name, "D%d.offsets", i+1);
            DEF_VAR(name, NC_INT, 2, dimids, MPI_INT, -1)
            PUT_ATTR_TXT("description", "Starting indices of flattened canonical noncontiguous requests of individual blobs")

            sprintf(name, "D%d.lengths", i+1);
            DEF_VAR(name, NC_INT, 2, dimids, MPI_INT, -1)
            PUT_ATTR_TXT("description", "Number of elements of flattened canonical noncontiguous requests of individual blobs")
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
    int g_dimids[MAX_NUM_DECOMP][MAX_NDIMS];
    float fillv = 1.e+20f, missv = 1.e+20f;
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

    DEF_DIM("time",  NC_UNLIMITED,         &dim_time)
    DEF_DIM("nbnd",  2,                    &dim_nbnd)
    DEF_DIM("chars", 8,                    &dim_chars)
    DEF_DIM("lev",   decom.dims[2][0],     &dim_lev)
    DEF_DIM("ilev",  decom.dims[2][0] + 1, &dim_ilev)
    DEF_DIM("ncol",  decom.dims[2][1],     &dim_ncol)

    if (cfg.strategy == blob && cfg.api != adios) {
        char name[64];

        g_dimids[0][0] = dim_ncol;
        g_dimids[1][0] = dim_ncol;
        g_dimids[2][0] = dim_lev;  g_dimids[2][1] = dim_ncol;

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

    /* For h0 file, there are 414 climate variables:
     *    6 scalar variables     + 408 array variables
     *   15 fixed-size variables + 399 record variables
     *   27 not partitioned      + 387 partitioned
     *
     * For h1 file, there are 51 climate variables:
     *    6 scalar variables     + 45 array variables
     *   15 fixed-size variables + 36 record variables
     *   27 not partitioned      + 24 partitioned
     */
    varp = vars + nvars_decomp - 1;

    /* below 10 are fixed-size climate variables ---------------------------*/

    /* double lat(ncol) */
    DEF_VAR("lat", NC_DOUBLE, 1, &dim_ncol, MPI_DOUBLE, 1)
    PUT_ATTR_TXT("long_name", "latitude")
    PUT_ATTR_TXT("units", "degrees_north")

    /* double lon(ncol) */
    DEF_VAR("lon", NC_DOUBLE, 1, &dim_ncol, MPI_DOUBLE, 1)
    PUT_ATTR_TXT("long_name", "longitude")
    PUT_ATTR_TXT("units", "degrees_east")

    /* double area(ncol) */
    DEF_VAR("area", NC_DOUBLE, 1, &dim_ncol, MPI_DOUBLE, 0)
    PUT_ATTR_TXT("long_name", "gll grid areas")

    /* double lev(lev) */
    DEF_VAR("lev", NC_DOUBLE, 1, &dim_lev, MPI_DOUBLE, -1)
    PUT_ATTR_TXT("long_name", "hybrid level at midpoints (1000*(A+B))")
    PUT_ATTR_TXT("units", "hPa")
    PUT_ATTR_TXT("positive", "down")
    PUT_ATTR_TXT("standard_name", "atmosphere_hybrid_sigma_pressure_coordinate")
    PUT_ATTR_TXT("formula_terms", "a: hyam b: hybm p0: P0 ps: PS")

    /* double hyam(lev) */
    DEF_VAR("hyam", NC_DOUBLE, 1, &dim_lev, MPI_DOUBLE, -1)
    PUT_ATTR_TXT("long_name", "hybrid A coefficient at layer midpoints")

    /* double hybm(lev) */
    DEF_VAR("hybm", NC_DOUBLE, 1, &dim_lev, MPI_DOUBLE, -1)
    PUT_ATTR_TXT("long_name", "hybrid B coefficient at layer midpoints")

    /* double P0 */
    DEF_VAR("P0", NC_DOUBLE, 0, NULL, MPI_DOUBLE, -1)
    PUT_ATTR_TXT("long_name", "reference pressure")
    PUT_ATTR_TXT("units", "Pa")

    /* double ilev(ilev) */
    DEF_VAR("ilev", NC_DOUBLE, 1, &dim_ilev, MPI_DOUBLE, -1)
    PUT_ATTR_TXT("long_name", "hybrid level at interfaces (1000*(A+B))")
    PUT_ATTR_TXT("units", "hPa")
    PUT_ATTR_TXT("positive", "down")
    PUT_ATTR_TXT("standard_name", "atmosphere_hybrid_sigma_pressure_coordinate")
    PUT_ATTR_TXT("formula_terms", "a: hyai b: hybi p0: P0 ps: PS")

    /* double hyai(ilev) */
    DEF_VAR("hyai", NC_DOUBLE, 1, &dim_ilev, MPI_DOUBLE, -1)
    PUT_ATTR_TXT("long_name", "hybrid A coefficient at layer interfaces")

    /* double hybi(ilev) */
    DEF_VAR("hybi", NC_DOUBLE, 1, &dim_ilev, MPI_DOUBLE, -1)
    PUT_ATTR_TXT("long_name", "hybrid B coefficient at layer interfaces")

    /* below 6 are record climate variables --------------------------------*/

    /* double time(time) */
    DEF_VAR("time", NC_DOUBLE, 1, &dim_time, MPI_DOUBLE, -1)
    PUT_ATTR_TXT("long_name", "time")
    PUT_ATTR_TXT("units", "days since 0001-01-01 00:00:00")
    PUT_ATTR_TXT("calendar", "noleap")
    PUT_ATTR_TXT("bounds", "time_bnds")

    /* int date(time) */
    DEF_VAR("date", NC_INT, 1, &dim_time, MPI_INT, -1)
    PUT_ATTR_TXT("long_name", "current date (YYYYMMDD)")

    /* int datesec(time) */
    DEF_VAR("datesec", NC_INT, 1, &dim_time, MPI_INT, -1)
    PUT_ATTR_TXT("long_name", "current seconds of current date")

    /* double time_bnds(time, nbnd) */
    dimids[0] = dim_time;
    dimids[1] = dim_nbnd;
    DEF_VAR("time_bnds", NC_DOUBLE, 2, dimids, MPI_DOUBLE, -1)
    PUT_ATTR_TXT("long_name", "time interval endpoints")

    /* char date_written(time, chars) */
    dimids[0] = dim_time;
    dimids[1] = dim_chars;
    DEF_VAR("date_written", NC_CHAR, 2, dimids, MPI_CHAR, -1)

    /* char time_written(time, chars) */
    dimids[0] = dim_time;
    dimids[1] = dim_chars;
    DEF_VAR("time_written", NC_CHAR, 2, dimids, MPI_CHAR, -1)

    /* below 5 are fixed-size climate variables ----------------------------*/

    /* int ndbase */
    DEF_VAR("ndbase", NC_INT, 0, NULL, MPI_INT, -1)
    PUT_ATTR_TXT("long_name", "base day")

    /* int nsbase */
    DEF_VAR("nsbase", NC_INT, 0, NULL, MPI_INT, -1)
    PUT_ATTR_TXT("long_name", "seconds of base day")

    /* int nbdate */
    DEF_VAR("nbdate", NC_INT, 0, NULL, MPI_INT, -1)
    PUT_ATTR_TXT("long_name", "base date (YYYYMMDD)")

    /* int nbsec */
    DEF_VAR("nbsec", NC_INT, 0, NULL, MPI_INT, -1)
    PUT_ATTR_TXT("long_name", "seconds of base date")

    /* int mdt */
    DEF_VAR("mdt", NC_INT, 0, NULL, MPI_INT, -1)
    PUT_ATTR_TXT("long_name", "timestep")
    PUT_ATTR_TXT("units", "s")

    /* below 393 are record climate variables ------------------------------*/

    /* int ndcur(time) */
    DEF_VAR("ndcur", NC_INT, 1, &dim_time, MPI_INT, -1)
    PUT_ATTR_TXT("long_name", "current day (from base day)")

    /* int nscur(time) */
    DEF_VAR("nscur", NC_INT, 1, &dim_time, MPI_INT, -1)
    PUT_ATTR_TXT("long_name", "current seconds of current day")

    /* double co2vmr(time) */
    DEF_VAR("co2vmr", NC_DOUBLE, 1, &dim_time, MPI_DOUBLE, -1)
    PUT_ATTR_TXT("long_name", "co2 volume mixing ratio")

    /* double ch4vmr(time) */
    DEF_VAR("ch4vmr", NC_DOUBLE, 1, &dim_time, MPI_DOUBLE, -1)
    PUT_ATTR_TXT("long_name", "ch4 volume mixing ratio")

    /* double n2ovmr(time) */
    DEF_VAR("n2ovmr", NC_DOUBLE, 1, &dim_time, MPI_DOUBLE, -1)
    PUT_ATTR_TXT("long_name", "n2o volume mixing ratio")

    /* double f11vmr(time) */
    DEF_VAR("f11vmr", NC_DOUBLE, 1, &dim_time, MPI_DOUBLE, -1)
    PUT_ATTR_TXT("long_name", "f11 volume mixing ratio")

    /* double f12vmr(time) */
    DEF_VAR("f12vmr", NC_DOUBLE, 1, &dim_time, MPI_DOUBLE, -1)
    PUT_ATTR_TXT("long_name", "f12 volume mixing ratio")

    /* double sol_tsi(time) */
    DEF_VAR("sol_tsi", NC_DOUBLE, 1, &dim_time, MPI_DOUBLE, -1)
    PUT_ATTR_TXT("long_name", "total solar irradiance")
    PUT_ATTR_TXT("units", "W/m2")

    /* int nsteph(time) */
    DEF_VAR("nsteph", NC_INT, 1, &dim_time, MPI_INT, -1)
    PUT_ATTR_TXT("long_name", "current timestep")

if (cfg.hist == h0) {
    /* float AEROD_v(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AEROD_v", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Total Aerosol Optical Depth in visible band")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ANRAIN(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("ANRAIN", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m-3")
    PUT_ATTR_TXT("long_name", "Average rain number conc")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ANSNOW(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("ANSNOW", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m-3")
    PUT_ATTR_TXT("long_name", "Average snow number conc")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODABS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODABS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol absorption optical depth 550 nm")
    PUT_ATTR_TXT("standard_name", "atmosphere_absorption_optical_thickness_due_to_ambient_aerosol_particles")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODABSBC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODABSBC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol absorption optical depth 550 nm from BC")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODALL(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODALL", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "AOD 550 nm for all time and species")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODBC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODBC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from BC")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODDUST(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODDUST", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from dust")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODDUST1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODDUST1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm model 1 from dust")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODDUST3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODDUST3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm model 3 from dust")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODDUST4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODDUST4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm model 4 from dust")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODMODE1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODMODE1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm mode 1")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODMODE2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODMODE2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm mode 2")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODMODE3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODMODE3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm mode 3")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODMODE4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODMODE4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm mode 4")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODNIR(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODNIR", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 850 nm")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODPOM(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODPOM", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from POM")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODSO4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODSO4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from SO4")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODSOA(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODSOA", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from SOA")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODSS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODSS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm from seasalt")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODUV(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODUV", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 350 nm")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AODVIS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AODVIS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol optical depth 550 nm")
    PUT_ATTR_TXT("standard_name", "atmosphere_optical_thickness_due_to_ambient_aerosol_particles")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AQRAIN(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("AQRAIN", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "Average rain mixing ratio")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AQSNOW(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("AQSNOW", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "Average snow mixing ratio")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AQ_DMS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AQ_DMS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "DMS aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AQ_H2O2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AQ_H2O2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "H2O2 aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AQ_H2SO4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AQ_H2SO4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "H2SO4 aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AQ_O3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AQ_O3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "O3 aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AQ_SO2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AQ_SO2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "SO2 aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AQ_SOAG(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("AQ_SOAG", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "SOAG aqueous chemistry (for gas species)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AREI(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("AREI", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "Micron")
    PUT_ATTR_TXT("long_name", "Average ice effective radius")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AREL(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("AREL", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "Micron")
    PUT_ATTR_TXT("long_name", "Average droplet effective radius")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AWNC(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("AWNC", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m-3")
    PUT_ATTR_TXT("long_name", "Average cloud water number conc")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float AWNI(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("AWNI", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m-3")
    PUT_ATTR_TXT("long_name", "Average cloud ice number conc")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float BURDEN1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("BURDEN1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Aerosol burden mode 1")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float BURDEN2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("BURDEN2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Aerosol burden mode 2")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float BURDEN3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("BURDEN3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Aerosol burden mode 3")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float BURDEN4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("BURDEN4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Aerosol burden mode 4")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float CCN3(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("CCN3", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/cm3")
    PUT_ATTR_TXT("long_name", "CCN concentration at S=0.1%")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float CDNUMC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("CDNUMC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "1/m2")
    PUT_ATTR_TXT("long_name", "Vertically-integrated droplet concentration")
    PUT_ATTR_TXT("cell_methods", "time: mean")
}
    /* float CLDHGH(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("CLDHGH", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Vertically-integrated high cloud")
    PUT_ATTR_TXT("cell_methods", "time: mean")

if (cfg.hist == h0) {
    /* float CLDICE(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("CLDICE", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged cloud ice amount")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float CLDLIQ(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("CLDLIQ", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged cloud liquid amount")
    PUT_ATTR_TXT("cell_methods", "time: mean")
}
    /* float CLDLOW(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("CLDLOW", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Vertically-integrated low cloud")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float CLDMED(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("CLDMED", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Vertically-integrated mid-level cloud")
    PUT_ATTR_TXT("cell_methods", "time: mean")

if (cfg.hist == h0) {
    /* float CLDTOT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("CLDTOT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Vertically-integrated total cloud")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float CLOUD(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("CLOUD", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Cloud fraction")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float CLOUDFRAC_CLUBB(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("CLOUDFRAC_CLUBB", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Cloud Fraction")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float CONCLD(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("CONCLD", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "fraction")
    PUT_ATTR_TXT("long_name", "Convective cloud cover")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float DCQ(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("DCQ", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg/s")
    PUT_ATTR_TXT("long_name", "Q tendency due to moist processes")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float DF_DMS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("DF_DMS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float DF_H2O2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("DF_H2O2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float DF_H2SO4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("DF_H2SO4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float DF_O3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("DF_O3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float DF_SO2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("DF_SO2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float DF_SOAG(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("DF_SOAG", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dry deposition flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float DMS_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("DMS_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "DMS in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float DP_KCLDBASE(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("DP_KCLDBASE", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Deep conv. cloudbase level index")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float DP_MFUP_MAX(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("DP_MFUP_MAX", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Deep conv. column-max updraft mass flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float DP_WCLDBASE(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("DP_WCLDBASE", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Deep conv. cloudbase vertical velocity")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float DSTSFMBL(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("DSTSFMBL", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Mobilization flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float DTCOND(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("DTCOND", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "K/s")
    PUT_ATTR_TXT("long_name", "T tendency - moist processes")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float DTENDTH(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("DTENDTH", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Dynamic Tendency of Total (vertically integrated) moist static energy")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float DTENDTQ(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("DTENDTQ", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Dynamic Tendency of Total (vertically integrated) specific humidity")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float EXTINCT(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("EXTINCT", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "/m")
    PUT_ATTR_TXT("long_name", "Aerosol extinction")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FICE(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("FICE", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fractional ice content within cloud")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FLDS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("FLDS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Downwelling longwave flux at surface")
    PUT_ATTR_TXT("standard_name", "surface_downwelling_longwave_flux_in_air")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FLNS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("FLNS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Net longwave flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FLNSC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("FLNSC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky net longwave flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
}
    /* float FLNT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("FLNT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Net longwave flux at top of model")
    PUT_ATTR_TXT("cell_methods", "time: mean")

if (cfg.hist == h0) {
    /* float FLNTC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("FLNTC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky net longwave flux at top of model")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FLUT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("FLUT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Upwelling longwave flux at top of model")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FLUTC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("FLUTC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky upwelling longwave flux at top of model")
    PUT_ATTR_TXT("standard_name", "toa_outgoing_longwave_flux_assuming_clear_sky")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FREQI(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("FREQI", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fractional occurrence of ice")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FREQL(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("FREQL", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fractional occurrence of liquid")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FREQR(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("FREQR", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fractional occurrence of rain")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FREQS(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("FREQS", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fractional occurrence of snow")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FSDS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("FSDS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Downwelling solar flux at surface")
    PUT_ATTR_TXT("standard_name", "surface_downwelling_shortwave_flux_in_air")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FSDSC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("FSDSC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky downwelling solar flux at surface")
    PUT_ATTR_TXT("standard_name", "surface_downwelling_shortwave_flux_in_air_assuming_clear_sky")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FSNS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("FSNS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Net solar flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FSNSC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("FSNSC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky net solar flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FSNT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("FSNT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Net solar flux at top of model")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FSNTC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("FSNTC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky net solar flux at top of model")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FSNTOA(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("FSNTOA", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Net solar flux at top of atmosphere")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FSNTOAC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("FSNTOAC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky net solar flux at top of atmosphere")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FSUTOA(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("FSUTOA", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Upwelling solar flux at top of atmosphere")
    PUT_ATTR_TXT("standard_name", "toa_outgoing_shortwave_flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float FSUTOAC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("FSUTOAC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Clearsky upwelling solar flux at top of atmosphere")
    PUT_ATTR_TXT("standard_name", "toa_outgoing_shortwave_flux_assuming_clear_sky")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float F_eff(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("F_eff", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Effective enrichment factor of marine organic matter")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float H2O2_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("H2O2_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "H2O2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float H2SO4_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("H2SO4_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "H2SO4 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float H2SO4_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("H2SO4_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "H2SO4 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ICEFRAC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ICEFRAC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fraction of sfc area covered by sea-ice")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ICIMR(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("ICIMR", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "Prognostic in-cloud ice mixing ratio")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ICWMR(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("ICWMR", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "Prognostic in-cloud water mixing ratio")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float IWC(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("IWC", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/m3")
    PUT_ATTR_TXT("long_name", "Grid box average ice water content")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float LANDFRAC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("LANDFRAC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fraction of sfc area covered by land")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float LHFLX(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("LHFLX", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Surface latent heat flux")
    PUT_ATTR_TXT("standard_name", "surface_upward_latent_heat_flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float LINOZ_DO3(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("LINOZ_DO3", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/s")
    PUT_ATTR_TXT("long_name", "ozone vmr tendency by linearized ozone chemistry")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float LINOZ_DO3_PSC(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("LINOZ_DO3_PSC", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/s")
    PUT_ATTR_TXT("long_name", "ozone vmr loss by PSCs using Carille et al. (1990)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float LINOZ_O3CLIM(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("LINOZ_O3CLIM", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "climatology of ozone in LINOZ")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float LINOZ_O3COL(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("LINOZ_O3COL", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "DU")
    PUT_ATTR_TXT("long_name", "ozone column above")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float LINOZ_SFCSINK(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("LINOZ_SFCSINK", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "Tg/yr/m2")
    PUT_ATTR_TXT("long_name", "surface o3 sink in LINOZ with an e-fold to a fixed concentration")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float LINOZ_SSO3(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("LINOZ_SSO3", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg")
    PUT_ATTR_TXT("long_name", "steady state ozone in LINOZ")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float LINOZ_SZA(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("LINOZ_SZA", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "degrees")
    PUT_ATTR_TXT("long_name", "solar zenith angle in LINOZ")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float LND_MBL(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("LND_MBL", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Soil erodibility factor")
    PUT_ATTR_TXT("cell_methods", "time: mean")
}
    /* float LWCF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("LWCF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Longwave cloud forcing")
    PUT_ATTR_TXT("standard_name", "toa_longwave_cloud_radiative_effect")
    PUT_ATTR_TXT("cell_methods", "time: mean")

if (cfg.hist == h0) {
    /* float Mass_bc(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("Mass_bc", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of bc mass concentration bc_a1+bc_c1+bc_a3+bc_c3+bc_a4+bc_c4")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float Mass_dst(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("Mass_dst", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of dst mass concentration dst_a1+dst_c1+dst_a3+dst_c3")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float Mass_mom(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("Mass_mom", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of mom mass concentration mom_a1+mom_c1+mom_a2+mom_c2+mom_a3+mom_c3+mom_a4+mom_c4")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float Mass_ncl(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("Mass_ncl", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of ncl mass concentration ncl_a1+ncl_c1+ncl_a2+ncl_c2+ncl_a3+ncl_c3")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float Mass_pom(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("Mass_pom", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of pom mass concentration pom_a1+pom_c1+pom_a3+pom_c3+pom_a4+pom_c4")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float Mass_so4(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("Mass_so4", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of so4 mass concentration so4_a1+so4_c1+so4_a2+so4_c2+so4_a3+so4_c3")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float Mass_soa(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("Mass_soa", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "sum of soa mass concentration soa_a1+soa_c1+soa_a2+soa_c2+soa_a3+soa_c3")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float NUMICE(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("NUMICE", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged cloud ice number")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float NUMLIQ(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("NUMLIQ", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged cloud liquid number")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float NUMRAI(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("NUMRAI", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged rain number")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float NUMSNO(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("NUMSNO", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "1/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged snow number")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float O3(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("O3", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("mixing_ratio", "dry")
    PUT_ATTR_TXT("long_name", "O3 concentration")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float O3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("O3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "O3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float OCNFRAC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("OCNFRAC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Fraction of sfc area covered by ocean")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float OMEGA(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("OMEGA", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "Pa/s")
    PUT_ATTR_TXT("long_name", "Vertical velocity (pressure)")
    PUT_ATTR_TXT("standard_name", "lagrangian_tendency_of_air_pressure")
    PUT_ATTR_TXT("cell_methods", "time: mean")
}

    /* float OMEGA500(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("OMEGA500", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "Pa/s")
    PUT_ATTR_TXT("long_name", "Vertical velocity at 500 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

if (cfg.hist == h0) {
    /* float OMEGAT(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("OMEGAT", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "K Pa/s")
    PUT_ATTR_TXT("long_name", "Vertical heat flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float PBLH(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("PBLH", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "PBL height")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float PHIS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("PHIS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m2/s2")
    PUT_ATTR_TXT("long_name", "Surface geopotential")
    PUT_ATTR_TXT("cell_methods", "time: point")

    /* float PRECC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("PRECC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Convective precipitation rate (liq + ice)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float PRECL(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("PRECL", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Large-scale (stable) precipitation rate (liq + ice)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float PRECSC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("PRECSC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Convective snow rate (water equivalent)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float PRECSL(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("PRECSL", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Large-scale (stable) snow rate (water equivalent)")
    PUT_ATTR_TXT("cell_methods", "time: mean")
}
else {
    /* float OMEGA850(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("OMEGA850", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "Pa/s")
    PUT_ATTR_TXT("long_name", "Vertical velocity at 850 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")

    /* float PRECT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("PRECT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Total (convective and large-scale) precipitation rate (liq + ice)")
    PUT_ATTR_TXT("cell_methods", "time: point")
}

    /* float PS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("PS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("long_name", "Surface pressure")
    PUT_ATTR_TXT("standard_name", "surface_air_pressure")
    PUT_ATTR_TXT("cell_methods", "time: mean")

if (cfg.hist == h0) {
    /* float PSL(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("PSL", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("long_name", "Sea level pressure")
    PUT_ATTR_TXT("standard_name", "air_pressure_at_mean_sea_level")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float Q(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("Q", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Specific humidity")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float QFLX(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("QFLX", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Surface water flux")
    PUT_ATTR_TXT("standard_name", "water_evapotranspiration_flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float QREFHT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("QREFHT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "Reference height humidity")
    PUT_ATTR_TXT("standard_name", "specific_humidity")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float QRL(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("QRL", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "K/s")
    PUT_ATTR_TXT("long_name", "Longwave heating rate")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float QRS(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("QRS", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "K/s")
    PUT_ATTR_TXT("long_name", "Solar heating rate")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float RAINQM(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("RAINQM", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged rain amount")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float RAM1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("RAM1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "frac")
    PUT_ATTR_TXT("long_name", "RAM1")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float RELHUM(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("RELHUM", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "percent")
    PUT_ATTR_TXT("long_name", "Relative humidity")
    PUT_ATTR_TXT("standard_name", "relative_humidity")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFDMS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFDMS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "DMS surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFH2O2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFH2O2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "H2O2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFH2SO4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFH2SO4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "H2SO4 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFO3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFO3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "O3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFSO2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFSO2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "SO2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFSOAG(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFSOAG", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "SOAG surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFbc_a1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFbc_a1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFbc_a3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFbc_a3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFbc_a4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFbc_a4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a4 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFdst_a1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFdst_a1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFdst_a3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFdst_a3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFmom_a1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFmom_a1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFmom_a2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFmom_a2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFmom_a3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFmom_a3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFmom_a4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFmom_a4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a4 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFncl_a1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFncl_a1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFncl_a2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFncl_a2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFncl_a3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFncl_a3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFnum_a1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFnum_a1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFnum_a2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFnum_a2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFnum_a3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFnum_a3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFnum_a4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFnum_a4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a4 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFpom_a1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFpom_a1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFpom_a3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFpom_a3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFpom_a4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFpom_a4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a4 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFso4_a1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFso4_a1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFso4_a2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFso4_a2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFso4_a3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFso4_a3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFsoa_a1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFsoa_a1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a1 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFsoa_a2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFsoa_a2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a2 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SFsoa_a3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SFsoa_a3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a3 surface flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SHFLX(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SHFLX", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Surface sensible heat flux")
    PUT_ATTR_TXT("standard_name", "surface_upward_sensible_heat_flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SH_KCLDBASE(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SH_KCLDBASE", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "1")
    PUT_ATTR_TXT("long_name", "Shallow conv. cloudbase level index")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SH_MFUP_MAX(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SH_MFUP_MAX", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Shallow conv. column-max updraft mass flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SH_WCLDBASE(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SH_WCLDBASE", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Shallow conv. cloudbase vertical velocity")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SNOWHICE(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SNOWHICE", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "Snow depth over ice")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SNOWHLND(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SNOWHLND", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "Water equivalent snow depth")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SNOWQM(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("SNOWQM", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("mixing_ratio", "wet")
    PUT_ATTR_TXT("long_name", "Grid box averaged snow amount")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SO2(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("SO2", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("mixing_ratio", "dry")
    PUT_ATTR_TXT("long_name", "SO2 concentration")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SO2_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SO2_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for SO2")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SO2_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SO2_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "SO2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SOAG_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SOAG_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for SOAG")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SOAG_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SOAG_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "mol/mol")
    PUT_ATTR_TXT("long_name", "SOAG in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SOAG_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SOAG_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "SOAG gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SOLIN(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SOLIN", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Solar insolation")
    PUT_ATTR_TXT("standard_name", "toa_incoming_shortwave_flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SSAVIS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SSAVIS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("long_name", "Aerosol singel-scatter albedo")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SSTSFMBL(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SSTSFMBL", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Mobilization flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float SSTSFMBL_OM(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SSTSFMBL_OM", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Mobilization flux of marine organic matter at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
}
    /* float SWCF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("SWCF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("Sampling_Sequence", "rad_lwsw")
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "W/m2")
    PUT_ATTR_TXT("long_name", "Shortwave cloud forcing")
    PUT_ATTR_TXT("standard_name", "toa_shortwave_cloud_radiative_effect")
    PUT_ATTR_TXT("cell_methods", "time: mean")

if (cfg.hist == h0) {
    /* float T(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("T", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Temperature")
    PUT_ATTR_TXT("standard_name", "air_temperature")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float TAUGWX(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TAUGWX", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "N/m2")
    PUT_ATTR_TXT("long_name", "Zonal gravity wave surface stress")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float TAUGWY(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TAUGWY", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "N/m2")
    PUT_ATTR_TXT("long_name", "Meridional gravity wave surface stress")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float TAUX(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TAUX", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "N/m2")
    PUT_ATTR_TXT("long_name", "Zonal surface stress")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float TAUY(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TAUY", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "N/m2")
    PUT_ATTR_TXT("long_name", "Meridional surface stress")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float TGCLDCWP(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TGCLDCWP", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Total grid-box cloud water path (liquid and ice)")
    PUT_ATTR_TXT("standard_name", "atmosphere_mass_content_of_cloud_condensed_water")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float TGCLDIWP(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TGCLDIWP", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Total grid-box cloud ice water path")
    PUT_ATTR_TXT("standard_name", "atmosphere_mass_content_of_cloud_ice")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float TGCLDLWP(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TGCLDLWP", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Total grid-box cloud liquid water path")
    PUT_ATTR_TXT("standard_name", "atmosphere_mass_content_of_cloud_liquid_water")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float TH7001000(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TH7001000", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Theta difference 700 mb - 1000 mb")
    PUT_ATTR_TXT("cell_methods", "time: mean")
}
else {
    /* float T850(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("T850", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Temperature at 850 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")
}

    /* float TMQ(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TMQ", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2")
    PUT_ATTR_TXT("long_name", "Total (vertically integrated) precipitable water")
    PUT_ATTR_TXT("standard_name", "atmosphere_mass_content_of_water_vapor")
    PUT_ATTR_TXT("cell_methods", "time: mean")

if (cfg.hist == h0) {
    /* float TREFHT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TREFHT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Reference height temperature")
    PUT_ATTR_TXT("standard_name", "air_temperature")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float TROP_P(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TROP_P", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "Pa")
    PUT_ATTR_TXT("long_name", "Tropopause Pressure")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float TROP_T(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TROP_T", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_FILL(fillv)
    PUT_ATTR_FLT1("missing_value", missv)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Tropopause Temperature")
    PUT_ATTR_TXT("cell_methods", "time: mean")
}

    /* float TS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Surface temperature (radiative)")
    PUT_ATTR_TXT("standard_name", "surface_temperature")
    PUT_ATTR_TXT("cell_methods", "time: mean")

if (cfg.hist == h0) {
    /* float TSMN(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TSMN", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Minimum surface temperature over output period")
    PUT_ATTR_TXT("cell_methods", "time: minimum")

    /* float TSMX(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TSMX", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "K")
    PUT_ATTR_TXT("long_name", "Maximum surface temperature over output period")
    PUT_ATTR_TXT("cell_methods", "time: maximum")

    /* float TUH(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TUH", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "W/m")
    PUT_ATTR_TXT("long_name", "Total (vertically integrated) zonal MSE flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float TUQ(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TUQ", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m/s")
    PUT_ATTR_TXT("long_name", "Total (vertically integrated) zonal water flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float TVH(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TVH", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "W/m")
    PUT_ATTR_TXT("long_name", "Total (vertically integrated) meridional MSE flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float TVQ(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("TVQ", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m/s")
    PUT_ATTR_TXT("long_name", "Total (vertically integrated) meridional water flux")
    PUT_ATTR_TXT("cell_methods", "time: mean")
}

    /* float U(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("U", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Zonal wind")
    PUT_ATTR_TXT("standard_name", "eastward_wind")
    PUT_ATTR_TXT("cell_methods", "time: mean")

if (cfg.hist == h0) {
    /* float U10(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("U10", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "10m wind speed")
    PUT_ATTR_TXT("standard_name", "wind_speed")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float UU(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("UU", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m2/s2")
    PUT_ATTR_TXT("long_name", "Zonal velocity squared")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float V(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("V", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Meridional wind")
    PUT_ATTR_TXT("standard_name", "northward_wind")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float VQ(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("VQ", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m/s kg/kg")
    PUT_ATTR_TXT("long_name", "Meridional water transport")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float VT(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("VT", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "K m/s")
    PUT_ATTR_TXT("long_name", "Meridional heat transport")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float VU(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("VU", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m2/s2")
    PUT_ATTR_TXT("long_name", "Meridional flux of zonal momentum")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float VV(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("VV", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m2/s2")
    PUT_ATTR_TXT("long_name", "Meridional velocity squared")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float WD_H2O2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("WD_H2O2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/s")
    PUT_ATTR_TXT("long_name", "H2O2             wet deposition")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float WD_H2SO4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("WD_H2SO4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/s")
    PUT_ATTR_TXT("long_name", "H2SO4            wet deposition")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float WD_SO2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("WD_SO2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/s")
    PUT_ATTR_TXT("long_name", "SO2              wet deposition")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float WSUB(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("WSUB", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Diagnostic sub-grid vertical velocity")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float Z3(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("Z3", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "Geopotential Height (above sea level)")
    PUT_ATTR_TXT("standard_name", "geopotential_height")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float aero_water(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("aero_water", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "sum of aerosol water of interstitial modes wat_a1+wat_a2+wat_a3+wat_a4")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float airFV(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("airFV", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "frac")
    PUT_ATTR_TXT("long_name", "FV")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_a1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_a1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_a1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_a1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_a1_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_a1_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "bc_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_a1_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_a1_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_a3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_a3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_a3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_a3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_a3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_a3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "bc_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_a4DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_a4DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_a4SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_a4SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_a4_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_a4_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for bc_a4")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_a4_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_a4_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "bc_a4 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_a4_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_a4_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_a4 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_c1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_c1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_c1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_c1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_c3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_c3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_c3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_c3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_c4DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_c4DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float bc_c4SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("bc_c4SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "bc_c4 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float chla(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("chla", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "mg L-1")
    PUT_ATTR_TXT("long_name", "ocean input data: chla")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float dst_a1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("dst_a1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float dst_a1SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("dst_a1SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a1 dust surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float dst_a1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("dst_a1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float dst_a1_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("dst_a1_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "dst_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float dst_a3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("dst_a3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float dst_a3SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("dst_a3SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_a3 dust surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float dst_a3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("dst_a3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float dst_a3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("dst_a3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "dst_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float dst_c1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("dst_c1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float dst_c1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("dst_c1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float dst_c3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("dst_c3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float dst_c3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("dst_c3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "dst_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float hstobie_linoz(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    DEF_VAR("hstobie_linoz", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    PUT_ATTR_INT("mdims", 1, &mdims)
    PUT_ATTR_TXT("units", "fraction of model time")
    PUT_ATTR_TXT("long_name", "Lowest possible Linoz level")
    PUT_ATTR_TXT("cell_methods", "time: point")

    /* float mlip(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mlip", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "uM C")
    PUT_ATTR_TXT("long_name", "ocean input data: mlip")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_a1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_a1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_a1SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_a1SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a1 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_a1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_a1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_a1_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_a1_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "mom_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_a1_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_a1_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_a2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_a2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_a2SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_a2SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a2 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_a2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_a2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_a2_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_a2_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "mom_a2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_a3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_a3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_a3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_a3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_a3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_a3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "mom_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_a4DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_a4DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_a4SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_a4SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a4 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_a4SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_a4SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_a4_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_a4_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "mom_a4 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_a4_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_a4_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_a4 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_c1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_c1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_c1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_c1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_c2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_c2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_c2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_c2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c2 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_c3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_c3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_c3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_c3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_c4DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_c4DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mom_c4SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mom_c4SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "mom_c4 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mpoly(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mpoly", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "uM C")
    PUT_ATTR_TXT("long_name", "ocean input data: mpoly")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float mprot(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("mprot", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "uM C")
    PUT_ATTR_TXT("long_name", "ocean input data: mprot")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_a1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_a1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_a1SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_a1SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a1 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_a1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_a1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_a1_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_a1_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "ncl_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_a2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_a2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_a2SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_a2SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a2 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_a2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_a2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_a2_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_a2_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "ncl_a2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_a3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_a3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_a3SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_a3SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_a3 seasalt surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_a3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_a3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_a3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_a3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "ncl_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_c1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_c1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_c1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_c1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_c2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_c2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_c2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_c2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c2 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_c3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_c3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float ncl_c3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("ncl_c3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "ncl_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a1SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a1SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "num_a1 dust surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a1_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a1_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for num_a1")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a1_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a1_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/kg")
    PUT_ATTR_TXT("long_name", "num_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a1_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a1_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "num_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a2_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a2_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for num_a2")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a2_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a2_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/kg")
    PUT_ATTR_TXT("long_name", "num_a2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a3SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a3SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "num_a3 dust surface emission")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/kg")
    PUT_ATTR_TXT("long_name", "num_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a4DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a4DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_a4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a4SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a4SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a4_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a4_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for num_a4")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a4_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a4_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/kg")
    PUT_ATTR_TXT("long_name", "num_a4 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_a4_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_a4_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "num_a4 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_c1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_c1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_c1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_c1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_c2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_c2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_c2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_c2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c2 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_c3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_c3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_c3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_c3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_c4DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_c4DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float num_c4SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("num_c4SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", " 1/m2/s")
    PUT_ATTR_TXT("long_name", "num_c4 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_a1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_a1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_a1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_a1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_a1_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_a1_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "pom_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_a1_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_a1_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_a3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_a3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_a3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_a3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_a3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_a3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "pom_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_a4DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_a4DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_a4SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_a4SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_a4_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_a4_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for pom_a4")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_a4_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_a4_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "pom_a4 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_a4_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_a4_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_a4 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_c1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_c1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_c1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_c1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_c3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_c3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_c3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_c3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_c4DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_c4DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c4 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float pom_c4SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("pom_c4SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "pom_c4 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_a1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_a1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_a1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_a1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_a1_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_a1_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for so4_a1")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_a1_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_a1_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "so4_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_a1_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_a1_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_a2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_a2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_a2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_a2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_a2_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_a2_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "molec/cm2/s")
    PUT_ATTR_TXT("long_name", "vertically intergrated external forcing for so4_a2")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_a2_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_a2_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "so4_a2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_a2_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_a2_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a2 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_a3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_a3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_a3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_a3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_a3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_a3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "so4_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_a3_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_a3_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_a3 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_c1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_c1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_c1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_c1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_c2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_c2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_c2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_c2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c2 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_c3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_c3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float so4_c3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("so4_c3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "so4_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_a1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_a1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_a1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_a1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_a1_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_a1_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "soa_a1 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_a1_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_a1_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a1 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_a2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_a2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_a2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_a2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_a2_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_a2_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "soa_a2 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_a2_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_a2_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a2 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_a3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_a3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_a3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_a3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "Wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_a3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_a3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/kg")
    PUT_ATTR_TXT("long_name", "soa_a3 in bottom layer")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_a3_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_a3_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_a3 gas-aerosol-exchange primary column tendency")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_c1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_c1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c1 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_c1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_c1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c1 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_c2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_c2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c2 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_c2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_c2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c2 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_c3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_c3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c3 dry deposition flux at bottom (grav + turb)")
    PUT_ATTR_TXT("cell_methods", "time: mean")

    /* float soa_c3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("soa_c3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "kg/m2/s")
    PUT_ATTR_TXT("long_name", "soa_c3 wet deposition flux at surface")
    PUT_ATTR_TXT("cell_methods", "time: mean")
}
else {
    /* float U250(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("U250", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Zonal wind at 250 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")

    /* float U850(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("U850", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Zonal wind at 850 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")

    /* float UBOT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("UBOT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Lowest model level zonal wind")
    PUT_ATTR_TXT("cell_methods", "time: point")

    /* float V250(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("V250", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Meridional wind at 250 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")

    /* float V850(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("V850", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Meridional wind at 850 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")

    /* float VBOT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("VBOT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m/s")
    PUT_ATTR_TXT("long_name", "Lowest model level meridional wind")
    PUT_ATTR_TXT("cell_methods", "time: point")

    /* float Z500(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    DEF_VAR("Z500", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    PUT_ATTR_TXT("units", "m")
    PUT_ATTR_TXT("long_name", "Geopotential Z at 500 mbar pressure surface")
    PUT_ATTR_TXT("cell_methods", "time: point")
}

    assert(varp - vars + 1 == cfg.nvars + nvars_decomp);

err_out:
    return err;
}

