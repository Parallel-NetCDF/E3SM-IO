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
    int err, nprocs, intBuf;

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

    /* save number of processes as global attributes */
    GET_GATTR_INT("ne", &intBuf)
    GET_GATTR_INT("np", &intBuf)

    txtBuf[0] = '\0';
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
    GET_GATTR_TXT("Conventions", txtBuf)
    GET_GATTR_TXT("institution_id", txtBuf)
    GET_GATTR_TXT("institution", txtBuf)
    GET_GATTR_TXT("contact", txtBuf)
    GET_GATTR_TXT("initial_file", txtBuf)
    GET_GATTR_TXT("topography_file", txtBuf)

    if (cfg.hist == h0)
        GET_GATTR_TXT("time_period_freq", txtBuf)
    else
        GET_GATTR_TXT("time_period_freq", txtBuf)

    assert(txtBuf[0] != '\0');
err_out:
    return err;
}

/*----< inq_var_decomp() >---------------------------------------------------*/
/* define decomposition variables. Blob I/O only */
int e3sm_io_case::inq_var_decomp(e3sm_io_config &cfg,
                                 e3sm_io_decom  &decom,
                                 e3sm_io_driver &driver,
                                 case_meta      *cmeta,
                                 int ncid,
                                 int dim_time,
                                 int dim_nblobs,
                                 int dim_max_nreqs[MAX_NUM_DECOMP],
                                 int g_dimids[MAX_NUM_DECOMP][MAX_NDIMS])
{
    char name[512], txtBuf[1024];
    int i, j, err=0, dimids[2];
    var_meta *varp;

    if (cfg.api == adios ) {
        for (j=0; j<decom.num_decomp; j++) {
            int ival, piodims[MAX_NUM_DECOMP];

            varp = vars + j;
            sprintf (name, "/__pio__/decomp/%d", (j + 512));
            err = driver.inq_varid(ncid, name, &varp->vid);
            CHECK_ERR

            for (i=0; i<decom.ndims[j]; i++)
                piodims[i] = (int)decom.dims[j][i];

            GET_ATTR_INT("dimlen", decom.ndims[j], piodims)
            GET_ATTR_INT("ndims", 1, decom.ndims+j)
            ival = 6;
            GET_ATTR_INT("piotype", 1, &ival)
        }
        err = driver.inq_varid(ncid, (char*)"/__pio__/info/nproc", &vars[j].vid);
        CHECK_ERR
    }
    else {
        std::map<int, std::string> dnames;
        varp = vars - 1;

        for (i=0; i<decom.num_decomp; i++) {
            dimids[0] = dim_nblobs;
            dimids[1] = dim_max_nreqs[i];
            sprintf(name, "D%d.nreqs", i+1);
            INQ_VAR(name, NC_INT, 1, dimids, MPI_INT, -1)
            GET_ATTR_TXT("description", txtBuf)
            GET_ATTR_INT("global_dimids", decom.ndims[i], g_dimids[i])
#ifdef GLOBAL_DIM_USE_STRING
            name[0] = '\0';
            for (j=0; j<decom.ndims[i]; j++) {
                char dim_name[128];
                ncmpi_inq_dimname(ncid, g_dimids[i][j], dim_name);
                strcat(name, dim_name);
                if (j < decom.ndims[i]-1) strcat(name,",");
            }
            GET_ATTR_TXT("global_dims", name)
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
            INQ_VAR(name, NC_INT64, 1, dimids, MPI_INT, -1)
            GET_ATTR_TXT("description", txtBuf)

            sprintf(name, "D%d.blob_count", i+1);
            INQ_VAR(name, NC_INT64, 1, dimids, MPI_LONG_LONG, -1)
            GET_ATTR_TXT("description", txtBuf)

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
            INQ_VAR(name, NC_INT, 2, dimids, MPI_INT, -1)
            GET_ATTR_TXT("description", txtBuf)

            sprintf(name, "D%d.lengths", i+1);
            INQ_VAR(name, NC_INT, 2, dimids, MPI_INT, -1)
            GET_ATTR_TXT("description", txtBuf)
        }
    }

err_out:
    return err;
}

/*----< inq_F_case() >-------------------------------------------------------*/
int e3sm_io_case::inq_F_case(e3sm_io_config &cfg,
                             e3sm_io_decom  &decom,
                             e3sm_io_driver &driver,
                             int             ncid)    /* file ID */
{
    char txtBuf[1024];
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
    err = get_gattrs(cfg, decom, driver, cmeta, ncid);
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

    INQ_DIM("time",  NC_UNLIMITED,         &dim_time)
    INQ_DIM("nbnd",  2,                    &dim_nbnd)
    INQ_DIM("chars", 8,                    &dim_chars)
    INQ_DIM("lev",   decom.dims[2][0],     &dim_lev)
    INQ_DIM("ilev",  decom.dims[2][0] + 1, &dim_ilev)
    INQ_DIM("ncol",  decom.dims[2][1],     &dim_ncol)

    if (cfg.strategy == blob && cfg.api != adios) {
        char name[64];

        g_dimids[0][0] = dim_ncol;
        g_dimids[1][0] = dim_ncol;
        g_dimids[2][0] = dim_lev;  g_dimids[2][1] = dim_ncol;

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
    INQ_VAR("lat", NC_DOUBLE, 1, &dim_ncol, MPI_DOUBLE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)

    /* double lon(ncol) */
    INQ_VAR("lon", NC_DOUBLE, 1, &dim_ncol, MPI_DOUBLE, 1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)

    /* double area(ncol) */
    INQ_VAR("area", NC_DOUBLE, 1, &dim_ncol, MPI_DOUBLE, 0)
    GET_ATTR_TXT("long_name", txtBuf)

    /* double lev(lev) */
    INQ_VAR("lev", NC_DOUBLE, 1, &dim_lev, MPI_DOUBLE, -1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("positive", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("formula_terms", txtBuf)

    /* double hyam(lev) */
    INQ_VAR("hyam", NC_DOUBLE, 1, &dim_lev, MPI_DOUBLE, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* double hybm(lev) */
    INQ_VAR("hybm", NC_DOUBLE, 1, &dim_lev, MPI_DOUBLE, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* double P0 */
    INQ_VAR("P0", NC_DOUBLE, 0, NULL, MPI_DOUBLE, -1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)

    /* double ilev(ilev) */
    INQ_VAR("ilev", NC_DOUBLE, 1, &dim_ilev, MPI_DOUBLE, -1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("positive", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("formula_terms", txtBuf)

    /* double hyai(ilev) */
    INQ_VAR("hyai", NC_DOUBLE, 1, &dim_ilev, MPI_DOUBLE, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* double hybi(ilev) */
    INQ_VAR("hybi", NC_DOUBLE, 1, &dim_ilev, MPI_DOUBLE, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* below 6 are record climate variables --------------------------------*/

    /* double time(time) */
    INQ_VAR("time", NC_DOUBLE, 1, &dim_time, MPI_DOUBLE, -1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("calendar", txtBuf)
    GET_ATTR_TXT("bounds", txtBuf)

    /* int date(time) */
    INQ_VAR("date", NC_INT, 1, &dim_time, MPI_INT, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* int datesec(time) */
    INQ_VAR("datesec", NC_INT, 1, &dim_time, MPI_INT, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* double time_bnds(time, nbnd) */
    dimids[0] = dim_time;
    dimids[1] = dim_nbnd;
    INQ_VAR("time_bnds", NC_DOUBLE, 2, dimids, MPI_DOUBLE, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* char date_written(time, chars) */
    dimids[0] = dim_time;
    dimids[1] = dim_chars;
    INQ_VAR("date_written", NC_CHAR, 2, dimids, MPI_CHAR, -1)

    /* char time_written(time, chars) */
    dimids[0] = dim_time;
    dimids[1] = dim_chars;
    INQ_VAR("time_written", NC_CHAR, 2, dimids, MPI_CHAR, -1)

    /* below 5 are fixed-size climate variables ----------------------------*/

    /* int ndbase */
    INQ_VAR("ndbase", NC_INT, 0, NULL, MPI_INT, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* int nsbase */
    INQ_VAR("nsbase", NC_INT, 0, NULL, MPI_INT, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* int nbdate */
    INQ_VAR("nbdate", NC_INT, 0, NULL, MPI_INT, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* int nbsec */
    INQ_VAR("nbsec", NC_INT, 0, NULL, MPI_INT, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* int mdt */
    INQ_VAR("mdt", NC_INT, 0, NULL, MPI_INT, -1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)

    /* below 393 are record climate variables ------------------------------*/

    /* int ndcur(time) */
    INQ_VAR("ndcur", NC_INT, 1, &dim_time, MPI_INT, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* int nscur(time) */
    INQ_VAR("nscur", NC_INT, 1, &dim_time, MPI_INT, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* double co2vmr(time) */
    INQ_VAR("co2vmr", NC_DOUBLE, 1, &dim_time, MPI_DOUBLE, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* double ch4vmr(time) */
    INQ_VAR("ch4vmr", NC_DOUBLE, 1, &dim_time, MPI_DOUBLE, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* double n2ovmr(time) */
    INQ_VAR("n2ovmr", NC_DOUBLE, 1, &dim_time, MPI_DOUBLE, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* double f11vmr(time) */
    INQ_VAR("f11vmr", NC_DOUBLE, 1, &dim_time, MPI_DOUBLE, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* double f12vmr(time) */
    INQ_VAR("f12vmr", NC_DOUBLE, 1, &dim_time, MPI_DOUBLE, -1)
    GET_ATTR_TXT("long_name", txtBuf)

    /* double sol_tsi(time) */
    INQ_VAR("sol_tsi", NC_DOUBLE, 1, &dim_time, MPI_DOUBLE, -1)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("units", txtBuf)

    /* int nsteph(time) */
    INQ_VAR("nsteph", NC_INT, 1, &dim_time, MPI_INT, -1)
    GET_ATTR_TXT("long_name", txtBuf)

if (cfg.hist == h0) {
    /* float AEROD_v(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AEROD_v", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ANRAIN(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("ANRAIN", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ANSNOW(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("ANSNOW", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODABS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODABS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODABSBC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODABSBC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODALL(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODALL", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODBC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODBC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODDUST(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODDUST", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODDUST1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODDUST1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODDUST3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODDUST3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODDUST4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODDUST4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODMODE1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODMODE1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODMODE2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODMODE2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODMODE3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODMODE3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODMODE4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODMODE4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODNIR(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODNIR", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODPOM(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODPOM", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODSO4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODSO4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODSOA(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODSOA", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODSS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODSS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODUV(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODUV", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AODVIS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AODVIS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AQRAIN(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("AQRAIN", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AQSNOW(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("AQSNOW", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AQ_DMS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AQ_DMS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AQ_H2O2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AQ_H2O2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AQ_H2SO4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AQ_H2SO4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AQ_O3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AQ_O3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AQ_SO2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AQ_SO2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AQ_SOAG(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("AQ_SOAG", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AREI(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("AREI", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AREL(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("AREL", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AWNC(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("AWNC", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float AWNI(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("AWNI", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float BURDEN1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("BURDEN1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float BURDEN2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("BURDEN2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float BURDEN3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("BURDEN3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float BURDEN4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("BURDEN4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float CCN3(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("CCN3", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float CDNUMC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("CDNUMC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
}
    /* float CLDHGH(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("CLDHGH", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

if (cfg.hist == h0) {
    /* float CLDICE(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("CLDICE", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("mixing_ratio", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float CLDLIQ(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("CLDLIQ", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("mixing_ratio", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
}
    /* float CLDLOW(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("CLDLOW", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float CLDMED(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("CLDMED", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

if (cfg.hist == h0) {
    /* float CLDTOT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("CLDTOT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float CLOUD(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("CLOUD", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float CLOUDFRAC_CLUBB(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("CLOUDFRAC_CLUBB", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float CONCLD(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("CONCLD", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float DCQ(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("DCQ", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float DF_DMS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("DF_DMS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float DF_H2O2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("DF_H2O2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float DF_H2SO4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("DF_H2SO4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float DF_O3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("DF_O3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float DF_SO2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("DF_SO2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float DF_SOAG(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("DF_SOAG", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float DMS_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("DMS_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float DP_KCLDBASE(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("DP_KCLDBASE", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float DP_MFUP_MAX(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("DP_MFUP_MAX", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float DP_WCLDBASE(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("DP_WCLDBASE", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float DSTSFMBL(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("DSTSFMBL", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float DTCOND(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("DTCOND", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float DTENDTH(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("DTENDTH", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float DTENDTQ(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("DTENDTQ", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float EXTINCT(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("EXTINCT", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FICE(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("FICE", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FLDS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("FLDS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FLNS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("FLNS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FLNSC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("FLNSC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
}
    /* float FLNT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("FLNT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

if (cfg.hist == h0) {
    /* float FLNTC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("FLNTC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FLUT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("FLUT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FLUTC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("FLUTC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FREQI(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("FREQI", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FREQL(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("FREQL", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FREQR(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("FREQR", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FREQS(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("FREQS", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FSDS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("FSDS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FSDSC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("FSDSC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FSNS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("FSNS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FSNSC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("FSNSC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FSNT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("FSNT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FSNTC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("FSNTC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FSNTOA(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("FSNTOA", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FSNTOAC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("FSNTOAC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FSUTOA(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("FSUTOA", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float FSUTOAC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("FSUTOAC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float F_eff(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("F_eff", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float H2O2_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("H2O2_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float H2SO4_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("H2SO4_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float H2SO4_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("H2SO4_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ICEFRAC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ICEFRAC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ICIMR(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("ICIMR", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ICWMR(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("ICWMR", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float IWC(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("IWC", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float LANDFRAC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("LANDFRAC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float LHFLX(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("LHFLX", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float LINOZ_DO3(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("LINOZ_DO3", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float LINOZ_DO3_PSC(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("LINOZ_DO3_PSC", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float LINOZ_O3CLIM(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("LINOZ_O3CLIM", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float LINOZ_O3COL(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("LINOZ_O3COL", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float LINOZ_SFCSINK(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("LINOZ_SFCSINK", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float LINOZ_SSO3(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("LINOZ_SSO3", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float LINOZ_SZA(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("LINOZ_SZA", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float LND_MBL(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("LND_MBL", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
}
    /* float LWCF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("LWCF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

if (cfg.hist == h0) {
    /* float Mass_bc(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("Mass_bc", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float Mass_dst(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("Mass_dst", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float Mass_mom(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("Mass_mom", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float Mass_ncl(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("Mass_ncl", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float Mass_pom(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("Mass_pom", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float Mass_so4(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("Mass_so4", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float Mass_soa(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("Mass_soa", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float NUMICE(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("NUMICE", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("mixing_ratio", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float NUMLIQ(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("NUMLIQ", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("mixing_ratio", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float NUMRAI(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("NUMRAI", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("mixing_ratio", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float NUMSNO(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("NUMSNO", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("mixing_ratio", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float O3(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("O3", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("mixing_ratio", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float O3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("O3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float OCNFRAC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("OCNFRAC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float OMEGA(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("OMEGA", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
}

    /* float OMEGA500(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("OMEGA500", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

if (cfg.hist == h0) {
    /* float OMEGAT(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("OMEGAT", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float PBLH(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("PBLH", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float PHIS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("PHIS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float PRECC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("PRECC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float PRECL(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("PRECL", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float PRECSC(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("PRECSC", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float PRECSL(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("PRECSL", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
}
else {
    /* float OMEGA850(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("OMEGA850", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float PRECT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("PRECT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
}

    /* float PS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("PS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

if (cfg.hist == h0) {
    /* float PSL(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("PSL", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float Q(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("Q", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("mixing_ratio", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float QFLX(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("QFLX", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float QREFHT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("QREFHT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float QRL(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("QRL", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float QRS(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("QRS", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float RAINQM(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("RAINQM", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("mixing_ratio", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float RAM1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("RAM1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float RELHUM(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("RELHUM", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFDMS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFDMS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFH2O2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFH2O2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFH2SO4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFH2SO4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFO3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFO3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFSO2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFSO2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFSOAG(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFSOAG", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFbc_a1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFbc_a1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFbc_a3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFbc_a3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFbc_a4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFbc_a4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFdst_a1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFdst_a1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFdst_a3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFdst_a3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFmom_a1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFmom_a1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFmom_a2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFmom_a2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFmom_a3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFmom_a3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFmom_a4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFmom_a4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFncl_a1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFncl_a1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFncl_a2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFncl_a2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFncl_a3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFncl_a3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFnum_a1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFnum_a1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFnum_a2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFnum_a2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFnum_a3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFnum_a3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFnum_a4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFnum_a4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFpom_a1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFpom_a1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFpom_a3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFpom_a3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFpom_a4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFpom_a4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFso4_a1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFso4_a1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFso4_a2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFso4_a2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFso4_a3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFso4_a3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFsoa_a1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFsoa_a1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFsoa_a2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFsoa_a2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SFsoa_a3(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SFsoa_a3", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SHFLX(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SHFLX", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SH_KCLDBASE(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SH_KCLDBASE", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SH_MFUP_MAX(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SH_MFUP_MAX", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SH_WCLDBASE(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SH_WCLDBASE", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SNOWHICE(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SNOWHICE", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SNOWHLND(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SNOWHLND", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SNOWQM(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("SNOWQM", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("mixing_ratio", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SO2(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("SO2", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("mixing_ratio", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SO2_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SO2_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SO2_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SO2_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SOAG_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SOAG_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SOAG_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SOAG_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SOAG_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SOAG_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SOLIN(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SOLIN", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SSAVIS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SSAVIS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SSTSFMBL(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SSTSFMBL", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float SSTSFMBL_OM(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SSTSFMBL_OM", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
}
    /* float SWCF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("SWCF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("Sampling_Sequence", txtBuf)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

if (cfg.hist == h0) {
    /* float T(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("T", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float TAUGWX(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TAUGWX", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float TAUGWY(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TAUGWY", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float TAUX(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TAUX", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float TAUY(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TAUY", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float TGCLDCWP(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TGCLDCWP", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float TGCLDIWP(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TGCLDIWP", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float TGCLDLWP(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TGCLDLWP", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float TH7001000(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TH7001000", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
}
else {
    /* float T850(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("T850", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
}

    /* float TMQ(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TMQ", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

if (cfg.hist == h0) {
    /* float TREFHT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TREFHT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float TROP_P(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TROP_P", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float TROP_T(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TROP_T", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_FILL(fillv)
    GET_ATTR_FLT1("missing_value", &missv)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
}

    /* float TS(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TS", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

if (cfg.hist == h0) {
    /* float TSMN(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TSMN", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float TSMX(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TSMX", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float TUH(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TUH", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float TUQ(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TUQ", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float TVH(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TVH", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float TVQ(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("TVQ", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
}

    /* float U(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("U", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

if (cfg.hist == h0) {
    /* float U10(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("U10", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float UU(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("UU", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float V(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("V", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float VQ(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("VQ", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float VT(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("VT", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float VU(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("VU", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float VV(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("VV", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float WD_H2O2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("WD_H2O2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float WD_H2SO4(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("WD_H2SO4", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float WD_SO2(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("WD_SO2", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float WSUB(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("WSUB", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float Z3(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("Z3", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("standard_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float aero_water(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("aero_water", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float airFV(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("airFV", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_a1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_a1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_a1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_a1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_a1_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_a1_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_a1_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_a1_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_a3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_a3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_a3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_a3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_a3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_a3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_a4DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_a4DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_a4SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_a4SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_a4_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_a4_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_a4_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_a4_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_a4_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_a4_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_c1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_c1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_c1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_c1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_c3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_c3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_c3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_c3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_c4DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_c4DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float bc_c4SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("bc_c4SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float chla(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("chla", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float dst_a1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("dst_a1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float dst_a1SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("dst_a1SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float dst_a1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("dst_a1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float dst_a1_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("dst_a1_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float dst_a3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("dst_a3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float dst_a3SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("dst_a3SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float dst_a3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("dst_a3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float dst_a3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("dst_a3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float dst_c1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("dst_c1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float dst_c1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("dst_c1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float dst_c3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("dst_c3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float dst_c3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("dst_c3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float hstobie_linoz(time, lev, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    INQ_VAR("hstobie_linoz", NC_FLOAT, 3, dimids, REC_ITYPE, 2)
    GET_ATTR_INT("mdims", 1, &mdims)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mlip(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mlip", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_a1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_a1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_a1SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_a1SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_a1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_a1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_a1_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_a1_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_a1_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_a1_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_a2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_a2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_a2SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_a2SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_a2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_a2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_a2_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_a2_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_a3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_a3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_a3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_a3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_a3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_a3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_a4DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_a4DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_a4SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_a4SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_a4SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_a4SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_a4_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_a4_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_a4_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_a4_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_c1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_c1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_c1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_c1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_c2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_c2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_c2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_c2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_c3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_c3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_c3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_c3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_c4DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_c4DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mom_c4SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mom_c4SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mpoly(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mpoly", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float mprot(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("mprot", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_a1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_a1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_a1SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_a1SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_a1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_a1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_a1_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_a1_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_a2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_a2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_a2SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_a2SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_a2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_a2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_a2_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_a2_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_a3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_a3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_a3SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_a3SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_a3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_a3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_a3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_a3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_c1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_c1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_c1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_c1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_c2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_c2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_c2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_c2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_c3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_c3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float ncl_c3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("ncl_c3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a1SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a1SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a1_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a1_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a1_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a1_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a1_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a1_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a2_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a2_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a2_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a2_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a3SF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a3SF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a4DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a4DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a4SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a4SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a4_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a4_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a4_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a4_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_a4_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_a4_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_c1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_c1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_c1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_c1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_c2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_c2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_c2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_c2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_c3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_c3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_c3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_c3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_c4DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_c4DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float num_c4SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("num_c4SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_a1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_a1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_a1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_a1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_a1_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_a1_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_a1_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_a1_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_a3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_a3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_a3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_a3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_a3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_a3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_a4DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_a4DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_a4SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_a4SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_a4_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_a4_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_a4_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_a4_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_a4_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_a4_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_c1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_c1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_c1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_c1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_c3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_c3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_c3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_c3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_c4DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_c4DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float pom_c4SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("pom_c4SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_a1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_a1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_a1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_a1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_a1_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_a1_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_a1_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_a1_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_a1_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_a1_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_a2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_a2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_a2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_a2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_a2_CLXF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_a2_CLXF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_a2_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_a2_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_a2_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_a2_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_a3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_a3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_a3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_a3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_a3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_a3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_a3_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_a3_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_c1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_c1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_c1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_c1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_c2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_c2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_c2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_c2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_c3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_c3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float so4_c3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("so4_c3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_a1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_a1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_a1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_a1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_a1_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_a1_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_a1_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_a1_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_a2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_a2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_a2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_a2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_a2_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_a2_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_a2_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_a2_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_a3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_a3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_a3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_a3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_a3_SRF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_a3_SRF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_a3_sfgaex1(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_a3_sfgaex1", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_c1DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_c1DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_c1SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_c1SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_c2DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_c2DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_c2SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_c2SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_c3DDF(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_c3DDF", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float soa_c3SFWET(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("soa_c3SFWET", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
}
else {
    /* float U250(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("U250", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float U850(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("U850", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float UBOT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("UBOT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float V250(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("V250", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float V850(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("V850", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float VBOT(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("VBOT", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)

    /* float Z500(time, ncol) */
    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    INQ_VAR("Z500", NC_FLOAT, 2, dimids, REC_ITYPE, 1)
    GET_ATTR_TXT("units", txtBuf)
    GET_ATTR_TXT("long_name", txtBuf)
    GET_ATTR_TXT("cell_methods", txtBuf)
}

    assert(varp - vars + 1 == cfg.nvars + nvars_decomp);

err_out:
    return err;
}

