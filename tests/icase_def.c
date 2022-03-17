/*********************************************************************
 *
 * Copyright (C) 2022, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 *********************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>

#define ERR { \
    if (err != NC_NOERR) { \
        printf("Error at %s:%d : %s\n", __FILE__,__LINE__, nc_strerror(err)); \
        exit(1); \
    } \
}

static void def(int ncid);

int main(int argc, char *argv[])
{
   int err=NC_NOERR, ncid;

   MPI_Init(&argc, &argv);
   err = nc_create_par("testfile.nc", NC_NETCDF4|NC_CLOBBER, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid); ERR
   err = nc_set_fill(ncid, NC_NOFILL, NULL); ERR

   /* define global attributes, dimensions, variables */
   def(ncid);

   err = nc_enddef(ncid); ERR
   err = nc_close(ncid); ERR

   MPI_Finalize();
   return err;
}

static void def(int ncid) {
   int err=NC_NOERR, varid, dimids[23];
   char buf[16];

err = nc_put_att(ncid, NC_GLOBAL, "title", NC_CHAR, 28, "ELM History file information"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "source", NC_CHAR, 15, "E3SM Land Model"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "source_id", NC_CHAR, 10, "6eb829238a"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "product", NC_CHAR, 12, "model-output"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "realm", NC_CHAR, 4, "land"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "case", NC_CHAR, 28, "I1850GSWCNPRDCTCBC_hcru_hcru"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "username", NC_CHAR, 4, "E3SM"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "hostname", NC_CHAR, 8, "cori-knl"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "git_version", NC_CHAR, 10, "6eb829238a"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "history", NC_CHAR, 28, "created on 07/13/21 20:37:15"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "institution_id", NC_CHAR, 12, "E3SM-Project"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "institution", NC_CHAR, 660, "LLNL (Lawrence Livermore National Laboratory, Livermore, CA 94550, USA); ANL (Argonne National Laboratory, Argonne, IL 60439, USA); BNL (Brookhaven National Laboratory, Upton, NY 11973, USA); LANL (Los Alamos National Laboratory, Los Alamos, NM 87545, USA); LBNL (Lawrence Berkeley National Laboratory, Berkeley, CA 94720, USA); ORNL (Oak Ridge National Laboratory, Oak Ridge, TN 37831, USA); PNNL (Pacific Northwest National Laboratory, Richland, WA 99352, USA); SNL (Sandia National Laboratories, Albuquerque, NM 87185, USA). Mailing address: LLNL Climate Program, c/o David C. Bader, Principal Investigator, L-103, 7000 East Avenue, Livermore, CA 94550, USA"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "contact", NC_CHAR, 35, "e3sm-data-support@listserv.llnl.gov"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "Conventions", NC_CHAR, 6, "CF-1.7"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "comment", NC_CHAR, 58, "NOTE: None of the variables are weighted by land fraction!"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "Surface_dataset", NC_CHAR, 40, "surfdata_360x720cru_simyr1850_c180216.nc"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "Initial_conditions_dataset", NC_CHAR, 24, "arbitrary initialization"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "PFT_physiological_constants_dataset", NC_CHAR, 21, "clm_params_c180524.nc"); ERR
err = nc_put_att(ncid, NC_GLOBAL, "ltype_vegetated_or_bare_soil", NC_INT, 1, buf); ERR
err = nc_put_att(ncid, NC_GLOBAL, "ltype_crop", NC_INT, 1, buf); ERR
err = nc_put_att(ncid, NC_GLOBAL, "ltype_landice", NC_INT, 1, buf); ERR
err = nc_put_att(ncid, NC_GLOBAL, "ltype_landice_multiple_elevation_classes", NC_INT, 1, buf); ERR
err = nc_put_att(ncid, NC_GLOBAL, "ltype_deep_lake", NC_INT, 1, buf); ERR
err = nc_put_att(ncid, NC_GLOBAL, "ltype_wetland", NC_INT, 1, buf); ERR
err = nc_put_att(ncid, NC_GLOBAL, "ltype_urban_tbd", NC_INT, 1, buf); ERR
err = nc_put_att(ncid, NC_GLOBAL, "ltype_urban_hd", NC_INT, 1, buf); ERR
err = nc_put_att(ncid, NC_GLOBAL, "ltype_urban_md", NC_INT, 1, buf); ERR

err = nc_def_dim(ncid, "time", 0, &dimids[0]); ERR
err = nc_def_dim(ncid, "lon", 144, &dimids[1]); ERR
err = nc_def_dim(ncid, "lat", 96, &dimids[2]); ERR
err = nc_def_dim(ncid, "gridcell", 62482, &dimids[3]); ERR
err = nc_def_dim(ncid, "topounit", 62482, &dimids[4]); ERR
err = nc_def_dim(ncid, "landunit", 271014, &dimids[5]); ERR
err = nc_def_dim(ncid, "column", 1020642, &dimids[6]); ERR
err = nc_def_dim(ncid, "pft", 2020354, &dimids[7]); ERR
err = nc_def_dim(ncid, "levgrnd", 15, &dimids[8]); ERR
err = nc_def_dim(ncid, "levurb", 5, &dimids[9]); ERR
err = nc_def_dim(ncid, "levlak", 10, &dimids[10]); ERR
err = nc_def_dim(ncid, "numrad", 2, &dimids[11]); ERR
err = nc_def_dim(ncid, "month", 12, &dimids[12]); ERR
err = nc_def_dim(ncid, "levsno", 5, &dimids[13]); ERR
err = nc_def_dim(ncid, "ltype", 9, &dimids[14]); ERR
err = nc_def_dim(ncid, "nvegwcs", 4, &dimids[15]); ERR
err = nc_def_dim(ncid, "natpft", 17, &dimids[16]); ERR
err = nc_def_dim(ncid, "string_length", 16, &dimids[17]); ERR
err = nc_def_dim(ncid, "levdcmp", 15, &dimids[18]); ERR
err = nc_def_dim(ncid, "levtrc", 10, &dimids[19]); ERR
err = nc_def_dim(ncid, "hist_interval", 2, &dimids[20]); ERR

err = nc_def_var(ncid, "levgrnd", NC_FLOAT, 1, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "coordinate soil levels"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR

err = nc_def_var(ncid, "levlak", NC_FLOAT, 1, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "coordinate lake levels"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR

err = nc_def_var(ncid, "levdcmp", NC_FLOAT, 1, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "coordinate soil levels"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR

err = nc_def_var(ncid, "time", NC_FLOAT, 1, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 4, "time"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 30, "days since 0001-01-01 00:00:00"); ERR
err = nc_put_att(ncid, varid, "calendar", NC_CHAR, 6, "noleap"); ERR
err = nc_put_att(ncid, varid, "bounds", NC_CHAR, 11, "time_bounds"); ERR

err = nc_def_var(ncid, "mcdate", NC_INT, 1, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 23, "current date (YYYYMMDD)"); ERR

err = nc_def_var(ncid, "mcsec", NC_INT, 1, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "current seconds of current date"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "s"); ERR

err = nc_def_var(ncid, "mdcur", NC_INT, 1, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "current day (from base day)"); ERR

err = nc_def_var(ncid, "mscur", NC_INT, 1, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "current seconds of current day"); ERR

err = nc_def_var(ncid, "nstep", NC_INT, 1, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 9, "time step"); ERR

err = nc_def_var(ncid, "time_bounds", NC_DOUBLE, 2, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "history time interval endpoints"); ERR

err = nc_def_var(ncid, "date_written", NC_CHAR, 2, dimids, &varid); ERR

err = nc_def_var(ncid, "time_written", NC_CHAR, 2, dimids, &varid); ERR

err = nc_def_var(ncid, "lon", NC_FLOAT, 1, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 20, "coordinate longitude"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 12, "degrees_east"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "lat", NC_FLOAT, 1, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 19, "coordinate latitude"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 13, "degrees_north"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "area", NC_FLOAT, 2, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 15, "grid cell areas"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "km^2"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "topo", NC_FLOAT, 2, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 20, "grid cell topography"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "landfrac", NC_FLOAT, 2, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 13, "land fraction"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "landmask", NC_INT, 2, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 38, "land/ocean mask (0.=ocean and 1.=land)"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_INT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_INT, 1, buf); ERR

err = nc_def_var(ncid, "pftmask", NC_INT, 2, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 40, "pft real/fake mask (0.=fake and 1.=real)"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_INT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_INT, 1, buf); ERR

err = nc_def_var(ncid, "ZSOI", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 10, "soil depth"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DZSOI", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 14, "soil thickness"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "WATSAT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 39, "saturated soil water content (porosity)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 7, "mm3/mm3"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SUCSAT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "saturated soil matric potential"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "mm"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "BSW", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 35, "slope of soil water retention curve"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "1"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "HKSAT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 32, "saturated hydraulic conductivity"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "1"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ZLAKE", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 21, "lake layer node depth"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DZLAKE", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 20, "lake layer thickness"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ACTUAL_IMMOB", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 23, "actual N immobilization"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ACTUAL_IMMOB_P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 23, "actual P immobilization"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ADSORBTION_P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 13, "adsorb P flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "AGNPP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 15, "aboveground NPP"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "AGWDNPP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 20, "aboveground wood NPP"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ALT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "current active layer thickness"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ALTMAX", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 37, "maximum annual active layer thickness"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ALTMAX_LASTYEAR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 41, "maximum prior year active layer thickness"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "AR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "autotrophic respiration (MR + GR)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "AVAILC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "C flux available for allocation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "AVAIL_RETRANSP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 42, "P flux available from retranslocation pool"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "BAF_CROP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "fractional area burned for crop"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 14, "proportion/sec"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "BAF_PEATF", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 34, "fractional area burned in peatland"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 14, "proportion/sec"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "BCDEP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 45, "total BC deposition (dry+wet) from atmosphere"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "kg/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "BGNPP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 15, "belowground NPP"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "BIOCHEM_PMIN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 36, "biochemical rate of P mineralization"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "BIOCHEM_PMIN_TO_PLANT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 44, "plant uptake of biochemical P mineralization"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "BTRAN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 25, "transpiration beta factor"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "1"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "BUILDHEAT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 56, "heat flux from urban building interior to walls and roof"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CH4PROD", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 32, "Gridcell total production of CH4"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 7, "gC/m2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CH4_SURF_AERE_SAT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 58, "aerenchyma surface CH4 flux for inundated area; (+ to atm)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "mol/m2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CH4_SURF_AERE_UNSAT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 62, "aerenchyma surface CH4 flux for non-inundated area; (+ to atm)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "mol/m2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CH4_SURF_DIFF_SAT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 64, "diffusive surface CH4 flux for inundated / lake area; (+ to atm)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "mol/m2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CH4_SURF_DIFF_UNSAT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 61, "diffusive surface CH4 flux for non-inundated area; (+ to atm)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "mol/m2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CH4_SURF_EBUL_SAT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 65, "ebullition surface CH4 flux for inundated / lake area; (+ to atm)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "mol/m2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CH4_SURF_EBUL_UNSAT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 62, "ebullition surface CH4 flux for non-inundated area; (+ to atm)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "mol/m2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "COL_PTRUNC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 34, "column-level sink for P truncation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CONC_CH4_SAT", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 48, "CH4 soil Concentration for inundated / lake area"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "mol/m3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CONC_CH4_UNSAT", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 45, "CH4 soil Concentration for non-inundated area"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "mol/m3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CONC_O2_SAT", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 47, "O2 soil Concentration for inundated / lake area"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "mol/m3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CONC_O2_UNSAT", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 44, "O2 soil Concentration for non-inundated area"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "mol/m3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CPOOL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "temporary photosynthate C pool"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CWDC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 5, "CWD C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CWDC_HR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 47, "coarse woody debris C heterotrophic respiration"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CWDC_LOSS", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 26, "coarse woody debris C loss"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CWDC_TO_LITR2C", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 46, "decomp. of coarse woody debris C to litter 2 C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CWDC_TO_LITR3C", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 46, "decomp. of coarse woody debris C to litter 3 C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CWDC_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "CWD C (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CWDN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 5, "CWD N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CWDN_TO_LITR2N", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 46, "decomp. of coarse woody debris N to litter 2 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CWDN_TO_LITR3N", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 46, "decomp. of coarse woody debris N to litter 3 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CWDN_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "CWD N (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CWDP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 5, "CWD P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CWDP_TO_LITR2P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 46, "decomp. of coarse woody debris P to litter 2 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CWDP_TO_LITR3P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 46, "decomp. of coarse woody debris P to litter 3 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "CWDP_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "CWD P (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DEADCROOTC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 18, "dead coarse root C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DEADCROOTN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 18, "dead coarse root N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DEADCROOTP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 18, "dead coarse root P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DEADSTEMC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 11, "dead stem C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DEADSTEMN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 11, "dead stem N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DEADSTEMP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 11, "dead stem P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DEFICIT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 21, "runoff supply deficit"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DENIT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "total rate of denitrification"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DESORPTION_P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 13, "desorp P flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DISPVEGC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 49, "displayed veg carbon, excluding storage and cpool"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DISPVEGN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "displayed vegetation nitrogen"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DISPVEGP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "displayed vegetation phosphorus"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DSTDEP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 47, "total dust deposition (dry+wet) from atmosphere"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "kg/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DSTFLXT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "total surface dust emission"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 7, "kg/m2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DWB", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "net change in total water mass"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DWT_CONV_CFLUX_DRIBBLED", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 71, "conversion C flux (immediate loss to atm), dribbled throughout the year"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DWT_CONV_CFLUX_GRC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 88, "conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DWT_CONV_NFLUX_GRC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 88, "conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DWT_CONV_PFLUX_GRC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 88, "conversion C flux (immediate loss to atm) (0 at all times except first timestep of year)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DWT_SLASH_CFLUX", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 46, "slash C flux to litter and CWD due to land use"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DWT_SLASH_NFLUX", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 46, "slash N flux to litter and CWD due to land use"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "DWT_SLASH_PFLUX", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 46, "slash P flux to litter and CWD due to land use"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "EFLX_DYNBAL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 48, "dynamic land cover change conversion energy flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "EFLX_GRND_LAKE", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 66, "net heat flux into lake/snow surface, excluding light transmission"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "EFLX_LH_TOT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "total latent heat flux [+ to atm]"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "EFLX_LH_TOT_R", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 23, "Rural total evaporation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "EFLX_LH_TOT_U", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 23, "Urban total evaporation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ELAI", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "exposed one-sided leaf area index"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 7, "m^2/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ER", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 56, "total ecosystem respiration, autotrophic + heterotrophic"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ERRH2O", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "total water conservation error"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "mm"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ERRH2OSNO", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 38, "imbalance in snow depth (liquid water)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "mm"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ERRSEB", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "surface energy conservation error"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ERRSOI", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 35, "soil/lake energy conservation error"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ERRSOL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 34, "solar radiation conservation error"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ESAI", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "exposed one-sided stem area index"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 7, "m^2/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FAREA_BURNED", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "timestep fractional area burned"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 10, "proportion"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FCEV", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 18, "canopy evaporation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FCH4", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 50, "Gridcell surface CH4 flux to atmosphere (+ to atm)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "kgC/m2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FCH4TOCO2", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 32, "Gridcell oxidation of CH4 to CO2"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 7, "gC/m2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FCH4_DFSAT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 66, "CH4 additional flux due to changing fsat, vegetated landunits only"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "kgC/m2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FCOV", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "fractional impermeable area"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "unitless"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FCTR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 20, "canopy transpiration"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FGEV", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 18, "ground evaporation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FGR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 79, "heat flux into soil/snow including snow melt and lake / snow light transmission"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FGR12", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 37, "heat flux between soil layers 1 and 2"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FGR_R", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 78, "Rural heat flux into soil/snow including snow melt and snow light transmission"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FGR_U", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 50, "Urban heat flux into soil/snow including snow melt"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FH2OSFC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "fraction of ground covered by surface water"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "1"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FINUNDATED", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 46, "fractional inundated area of vegetated columns"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "1"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FINUNDATED_LAG", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 51, "time-lagged inundated fraction of vegetated columns"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "1"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FIRA", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "net infrared (longwave) radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FIRA_R", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 39, "Rural net infrared (longwave) radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FIRA_U", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 39, "Urban net infrared (longwave) radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FIRE", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 37, "emitted infrared (longwave) radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FIRE_R", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "Rural emitted infrared (longwave) radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FIRE_U", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "Urban emitted infrared (longwave) radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FLDS", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "atmospheric longwave radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FPG", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 45, "fraction of potential gpp due to N limitation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 10, "proportion"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FPG_P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 45, "fraction of potential gpp due to P limitation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 10, "proportion"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FPI", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 48, "fraction of potential immobilization of nitrogen"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 10, "proportion"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FPI_P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 50, "fraction of potential immobilization of phosphorus"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 10, "proportion"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FPI_P_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 50, "fraction of potential immobilization of phosphorus"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 10, "proportion"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FPI_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 48, "fraction of potential immobilization of nitrogen"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 10, "proportion"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FPSN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 14, "photosynthesis"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "umol/m2s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FPSN_WC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "Rubisco-limited photosynthesis"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "umol/m2s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FPSN_WJ", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "RuBP-limited photosynthesis"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "umol/m2s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FPSN_WP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "Product-limited photosynthesis"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "umol/m2s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FROOTC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 11, "fine root C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FROOTC_ALLOC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "fine root C allocation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FROOTC_LOSS", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 16, "fine root C loss"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FROOTN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 11, "fine root N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FROOTP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 11, "fine root P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FROST_TABLE", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 44, "frost table depth (vegetated landunits only)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSA", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 24, "absorbed solar radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSAT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "fractional area with water table at surface"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "unitless"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSA_R", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "Rural absorbed solar radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSA_U", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "Urban absorbed solar radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSDS", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 36, "atmospheric incident solar radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSDSND", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 35, "direct nir incident solar radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSDSNDLN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 49, "direct nir incident solar radiation at local noon"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSDSNI", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 36, "diffuse nir incident solar radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSDSVD", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 35, "direct vis incident solar radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSDSVDLN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 49, "direct vis incident solar radiation at local noon"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSDSVI", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 36, "diffuse vis incident solar radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSDSVILN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 50, "diffuse vis incident solar radiation at local noon"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSH", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 13, "sensible heat"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSH_G", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 25, "sensible heat from ground"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSH_NODYNLNDUSE", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 58, "sensible heat not including correction for land use change"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSH_R", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 19, "Rural sensible heat"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSH_U", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 19, "Urban sensible heat"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSH_V", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "sensible heat from veg"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSM", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 19, "snow melt heat flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSM_R", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 25, "Rural snow melt heat flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSM_U", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 25, "Urban snow melt heat flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSNO", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 34, "fraction of ground covered by snow"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "1"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSNO_EFF", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 44, "effective fraction of ground covered by snow"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "1"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 25, "reflected solar radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSRND", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 36, "direct nir reflected solar radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSRNDLN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 50, "direct nir reflected solar radiation at local noon"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSRNI", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 37, "diffuse nir reflected solar radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSRVD", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 36, "direct vis reflected solar radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSRVDLN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 50, "direct vis reflected solar radiation at local noon"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "FSRVI", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 37, "diffuse vis reflected solar radiation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "F_CO2_SOIL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 28, "total soil-atm. CO2 exchange"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "F_CO2_SOIL_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 48, "total vertically resolved soil-atm. CO2 exchange"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "F_DENIT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 20, "denitrification flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "F_DENIT_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 20, "denitrification flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "F_N2O_DENIT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 24, "denitrification N2O flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "F_N2O_NIT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "nitrification N2O flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "F_NIT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 18, "nitrification flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "F_NIT_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 18, "nitrification flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "GC_HEAT1", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 35, "initial gridcell total heat content"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "J/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "GC_ICE1", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 34, "initial gridcell total ice content"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "mm"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "GC_LIQ1", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 34, "initial gridcell total liq content"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "mm"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "GPP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 24, "gross primary production"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "GR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 24, "total growth respiration"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "GROSS_NMIN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "gross rate of N mineralization"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "GROSS_PMIN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "gross rate of P mineralization"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "H2OCAN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 17, "intercepted water"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "mm"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "H2OSFC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 19, "surface water depth"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "mm"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "H2OSNO", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 25, "snow depth (liquid water)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "mm"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "H2OSNO_TOP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "mass of snow in top snow layer"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "kg/m2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "H2OSOI", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 48, "volumetric soil water (vegetated landunits only)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 7, "mm3/mm3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "HC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "heat content of soil/snow/lake"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "MJ/m2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "HCSOI", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 17, "soil heat content"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "MJ/m2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "HEAT_FROM_AC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 76, "sensible heat flux put into canyon due to heat removed from air conditioning"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "HR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "total heterotrophic respiration"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "HR_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 51, "total vertically resolved heterotrophic respiration"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "HTOP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 10, "canopy top"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "INT_SNOW", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 42, "accumulated swe (vegetated landunits only)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "mm"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LABILEP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 13, "soil Labile P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LABILEP_TO_SECONDP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "LABILE P TO SECONDARY MINERAL P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LABILEP_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 26, "soil labile P (vert. res.)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gp/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LAISHA", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 32, "shaded projected leaf area index"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "1"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LAISUN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 32, "sunlit projected leaf area index"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "1"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LAKEICEFRAC", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 28, "lake layer ice mass fraction"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "unitless"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LAKEICETHICK", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 64, "thickness of lake ice (including physical expansion on freezing)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LAND_UPTAKE", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 44, "NEE minus LAND_USE_FLUX, negative for update"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LAND_USE_FLUX", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 65, "total C emitted from land cover conversion and wood product pools"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LEAFC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 6, "leaf C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LEAFC_ALLOC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 17, "leaf C allocation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LEAFC_LOSS", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 11, "leaf C loss"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LEAFC_TO_LITTER", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 17, "leaf C litterfall"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LEAFN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 6, "leaf N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LEAFP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 6, "leaf P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LEAF_MR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 28, "leaf maintenance respiration"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LFC2", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 51, "conversion area fraction of BET and BDT that burned"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 7, "per sec"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITFALL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 34, "litterfall (leaves and fine roots)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITHR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 32, "litter heterotrophic respiration"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR1C", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "LITR1 C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR1C_TO_SOIL1C", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "decomp. of litter 1 C to soil 1 C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR1C_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "LITR1 C (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR1N", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "LITR1 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR1N_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 45, "litter 1 N tendency due to vertical transport"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR1N_TO_SOIL1N", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "decomp. of litter 1 N to soil 1 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR1N_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "LITR1 N (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR1P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "LITR1 P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR1P_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 45, "litter 1 P tendency due to vertical transport"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR1P_TO_SOIL1P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "decomp. of litter 1 P to soil 1 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR1P_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "LITR1 P (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR1_HR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 24, "Het. Resp. from litter 1"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR2C", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "LITR2 C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR2C_TO_SOIL2C", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "decomp. of litter 2 C to soil 2 C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR2C_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "LITR2 C (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR2N", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "LITR2 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR2N_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 45, "litter 2 N tendency due to vertical transport"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR2N_TO_SOIL2N", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "decomp. of litter 2 N to soil 2 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR2N_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "LITR2 N (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR2P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "LITR2 P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR2P_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 45, "litter 2 P tendency due to vertical transport"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR2P_TO_SOIL2P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "decomp. of litter 2 P to soil 2 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR2P_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "LITR2 P (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR2_HR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 24, "Het. Resp. from litter 2"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR3C", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "LITR3 C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR3C_TO_SOIL3C", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "decomp. of litter 3 C to soil 3 C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR3C_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "LITR3 C (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR3N", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "LITR3 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR3N_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 45, "litter 3 N tendency due to vertical transport"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR3N_TO_SOIL3N", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "decomp. of litter 3 N to soil 3 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR3N_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "LITR3 N (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR3P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "LITR3 P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR3P_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 45, "litter 3 P tendency due to vertical transport"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR3P_TO_SOIL3P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "decomp. of litter 3 P to soil 3 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR3P_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "LITR3 P (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITR3_HR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 24, "Het. Resp. from litter 3"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITTERC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 8, "litter C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITTERC_HR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 34, "litter C heterotrophic respiration"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LITTERC_LOSS", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 13, "litter C loss"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LIVECROOTC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 18, "live coarse root C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LIVECROOTN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 18, "live coarse root N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LIVECROOTP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 18, "live coarse root P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LIVESTEMC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 11, "live stem C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LIVESTEMN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 11, "live stem N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "LIVESTEMP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 11, "live stem P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "MR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 23, "maintenance respiration"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "M_LITR1C_TO_LEACHING", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 24, "litter 1 C leaching loss"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "M_LITR2C_TO_LEACHING", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 24, "litter 2 C leaching loss"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "M_LITR3C_TO_LEACHING", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 24, "litter 3 C leaching loss"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "M_SOIL1C_TO_LEACHING", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "soil 1 C leaching loss"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "M_SOIL2C_TO_LEACHING", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "soil 2 C leaching loss"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "M_SOIL3C_TO_LEACHING", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "soil 3 C leaching loss"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "M_SOIL4C_TO_LEACHING", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "soil 4 C leaching loss"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "NBP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 81, "net biome production, includes fire, landuse, and harvest flux, positive for sink"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "NDEPLOY", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "total N deployed in new growth"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "NDEP_TO_SMINN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 42, "atmospheric N deposition to soil mineral N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "NEE", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 109, "net ecosystem exchange of carbon, includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "NEM", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 68, "Gridcell net adjustment to NEE passed to atm. for methane production"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 7, "gC/m2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "NEP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 85, "net ecosystem production, excludes fire, landuse, and harvest flux, positive for sink"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "NET_NMIN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 28, "net rate of N mineralization"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "NET_PMIN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 28, "net rate of P mineralization"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "NFIRE", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "fire counts valid only in Reg.C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 14, "counts/km2/sec"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "NFIX_TO_SMINN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 49, "symbiotic/asymbiotic N fixation to soil mineral N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "NPP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "net primary production"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "OCCLP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 15, "soil occluded P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "OCCLP_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 28, "soil occluded P (vert. res.)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gp/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "OCDEP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 45, "total OC deposition (dry+wet) from atmosphere"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "kg/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "O_SCALAR", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 56, "fraction by which decomposition is reduced due to anoxia"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "1"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PARVEGLN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 40, "absorbed par by vegetation at local noon"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PBOT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 20, "atmospheric pressure"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "Pa"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PCH4", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 35, "atmospheric partial pressure of CH4"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "Pa"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PCO2", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 35, "atmospheric partial pressure of CO2"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "Pa"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PCT_LANDUNIT", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "% of each landunit on grid cell"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "%"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PCT_NAT_PFT", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 61, "% of each PFT on the natural vegetation (i.e., soil) landunit"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "%"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PDEPLOY", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "total P deployed in new growth"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PDEP_TO_SMINP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 42, "atmospheric P deposition to soil mineral P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PFT_FIRE_CLOSS", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 83, "total patch-level fire C loss for non-peat fires outside land-type converted region"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PFT_FIRE_NLOSS", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "total pft-level fire N loss"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PLANT_CALLOC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "total allocated C flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PLANT_NDEMAND", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 38, "N flux required to support initial GPP"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PLANT_NDEMAND_COL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 38, "N flux required to support initial GPP"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PLANT_PALLOC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "total allocated P flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PLANT_PDEMAND", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 38, "P flux required to support initial GPP"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PLANT_PDEMAND_COL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 38, "P flux required to support initial GPP"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "POTENTIAL_IMMOB", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 26, "potential N immobilization"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "POTENTIAL_IMMOB_P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 26, "potential P immobilization"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "POT_F_DENIT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "potential denitrification flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "POT_F_NIT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 28, "potential nitrification flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PRIMP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 14, "soil primary P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PRIMP_TO_LABILEP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "PRIMARY MINERAL P TO LABILE P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PRIMP_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "soil primary P (vert. res.)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gp/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PROD1P_LOSS", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 32, "loss from 1-yr crop product pool"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PSNSHA", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 26, "shaded leaf photosynthesis"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 13, "umolCO2/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PSNSHADE_TO_CPOOL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "C fixation from shaded canopy"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PSNSUN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 26, "sunlit leaf photosynthesis"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 13, "umolCO2/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "PSNSUN_TO_CPOOL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "C fixation from sunlit canopy"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "Q2M", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 20, "2m specific humidity"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "kg/kg"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QBOT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "atmospheric specific humidity"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "kg/kg"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QCHARGE", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 48, "aquifer recharge rate (vegetated landunits only)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QDRAI", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 20, "sub-surface drainage"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QDRAI_PERCH", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 19, "perched wt drainage"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QDRAI_XS", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 26, "saturation excess drainage"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QDRIP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 11, "throughfall"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QFLOOD", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 26, "runoff from river flooding"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QFLX_ICE_DYNBAL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 52, "ice dynamic land cover change conversion runoff flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QFLX_LIQ_DYNBAL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 52, "liq dynamic land cover change conversion runoff flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QH2OSFC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 20, "surface water runoff"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QINFL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 12, "infiltration"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QINTR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 12, "interception"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QIRRIG_GRND", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "Groundwater irrigation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QIRRIG_ORIG", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 57, "Original total irrigation water demand (surface + ground)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QIRRIG_REAL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 56, "actual water added through irrigation (surface + ground)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QIRRIG_SURF", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 24, "Surface water irrigation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QIRRIG_WM", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 49, "Surface water irrigation demand sent to MOSART/WM"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QOVER", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 14, "surface runoff"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QOVER_LAG", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "time-lagged surface runoff for soil columns"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QRGWL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 57, "surface runoff at glaciers (liquid only), wetlands, lakes"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QRUNOFF", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 48, "total liquid runoff (does not include QSNWCPICE)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QRUNOFF_NODYNLNDUSE", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 93, "total liquid runoff (does not include QSNWCPICE) not including correction for land use change"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QRUNOFF_R", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 18, "Rural total runoff"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QRUNOFF_U", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 18, "Urban total runoff"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QSNOMELT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 9, "snow melt"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QSNWCPICE", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 35, "excess snowfall due to snow capping"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QSNWCPICE_NODYNLNDUSE", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 80, "excess snowfall due to snow capping not including correction for land use change"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "mm H2O/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QSOIL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 72, "Ground evaporation (soil/snow evaporation + soil/snow sublimation - dew)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QVEGE", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 18, "canopy evaporation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "QVEGT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 20, "canopy transpiration"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "RAIN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 16, "atmospheric rain"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "RETRANSN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "plant pool of retranslocated N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "RETRANSN_TO_NPOOL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "deployment of retranslocated N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "RETRANSP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "plant pool of retranslocated P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "RETRANSP_TO_PPOOL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "deployment of retranslocated P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "RH2M", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 20, "2m relative humidity"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "%"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "RH2M_R", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 26, "Rural 2m specific humidity"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "%"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "RH2M_U", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 26, "Urban 2m relative humidity"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "%"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "RR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 47, "root respiration (fine root MR + total root GR)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SABG", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 28, "solar rad absorbed by ground"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SABG_PEN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 50, "Rural solar rad penetrating top soil or snow layer"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "watt/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SABV", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 25, "solar rad absorbed by veg"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SCALARAVG_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "average of decomposition scalar"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "fraction"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SECONDP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 16, "soil secondary P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SECONDP_TO_LABILEP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "SECONDARY MINERAL P TO LABILE P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SECONDP_TO_OCCLP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "SECONDARY MINERAL P TO OCCLUDED P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SECONDP_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "soil secondary P (vert. res.)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gp/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SEEDC_GRC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 47, "pool for seeding new PFTs via dynamic landcover"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 14, "soil mineral N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINN_TO_NPOOL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 35, "deployment of soil mineral N uptake"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINN_TO_PLANT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "plant uptake of soil mineral N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINN_TO_SOIL1N_L1", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "mineral N flux for decomp. of LITR1to SOIL1"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINN_TO_SOIL2N_L2", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "mineral N flux for decomp. of LITR2to SOIL2"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINN_TO_SOIL2N_S1", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "mineral N flux for decomp. of SOIL1to SOIL2"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINN_TO_SOIL3N_L3", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "mineral N flux for decomp. of LITR3to SOIL3"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINN_TO_SOIL3N_S2", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "mineral N flux for decomp. of SOIL2to SOIL3"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINN_TO_SOIL4N_S3", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "mineral N flux for decomp. of SOIL3to SOIL4"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 14, "soil mineral P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINP_LEACHED", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 36, "soil mineral P pool loss to leaching"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINP_TO_PLANT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "plant uptake of soil mineral P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINP_TO_PPOOL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 35, "deployment of soil mineral P uptake"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINP_TO_SOIL1P_L1", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "mineral P flux for decomp. of LITR1to SOIL1"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINP_TO_SOIL2P_L2", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "mineral P flux for decomp. of LITR2to SOIL2"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINP_TO_SOIL2P_S1", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "mineral P flux for decomp. of SOIL1to SOIL2"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINP_TO_SOIL3P_L3", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "mineral P flux for decomp. of LITR3to SOIL3"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINP_TO_SOIL3P_S2", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "mineral P flux for decomp. of SOIL2to SOIL3"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINP_TO_SOIL4P_S3", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "mineral P flux for decomp. of SOIL3to SOIL4"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMINP_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "soil mineral P (vert. res.)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gp/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMIN_NH4", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 16, "soil mineral NH4"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMIN_NH4_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "soil mineral NH4 (vert. res.)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMIN_NO3", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 16, "soil mineral NO3"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMIN_NO3_LEACHED", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "soil NO3 pool loss to leaching"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMIN_NO3_RUNOFF", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 28, "soil NO3 pool loss to runoff"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SMIN_NO3_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "soil mineral NO3 (vert. res.)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SNOBCMCL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 25, "mass of BC in snow column"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "kg/m2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SNOBCMSL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 28, "mass of BC in top snow layer"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "kg/m2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SNODSTMCL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "mass of dust in snow column"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "kg/m2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SNODSTMSL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "mass of dust in top snow layer"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "kg/m2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SNOINTABS", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 55, "Percent of incoming solar absorbed by lower snow layers"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "%"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SNOOCMCL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 25, "mass of OC in snow column"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "kg/m2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SNOOCMSL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 28, "mass of OC in top snow layer"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "kg/m2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SNOW", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 16, "atmospheric snow"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SNOWDP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 25, "gridcell mean snow height"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SNOWICE", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 8, "snow ice"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "kg/m2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SNOWLIQ", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 17, "snow liquid water"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "kg/m2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SNOW_DEPTH", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 32, "snow height of snow covered area"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SNOW_SINKS", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 25, "snow sinks (liquid water)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SNOW_SOURCES", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "snow sources (liquid water)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL1C", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "SOIL1 C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL1C_TO_SOIL2C", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "decomp. of soil 1 C to soil 2 C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL1C_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "SOIL1 C (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL1N", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "SOIL1 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL1N_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "soil 1 N tendency due to vertical transport"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL1N_TO_SOIL2N", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "decomp. of soil 1 N to soil 2 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL1N_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "SOIL1 N (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL1P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "SOIL1 P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL1P_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "soil 1 P tendency due to vertical transport"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL1P_TO_SOIL2P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "decomp. of soil 1 P to soil 2 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL1P_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "SOIL1 P (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL1_HR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "Het. Resp. from soil 1"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL2C", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "SOIL2 C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL2C_TO_SOIL3C", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "decomp. of soil 2 C to soil 3 C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL2C_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "SOIL2 C (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL2N", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "SOIL2 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL2N_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "soil 2 N tendency due to vertical transport"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL2N_TO_SOIL3N", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "decomp. of soil 2 N to soil 3 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL2N_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "SOIL2 N (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL2P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "SOIL2 P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL2P_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "soil 2 P tendency due to vertical transport"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL2P_TO_SOIL3P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "decomp. of soil 2 P to soil 3 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL2P_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "SOIL2 P (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL2_HR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "Het. Resp. from soil 2"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL3C", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "SOIL3 C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL3C_TO_SOIL4C", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "decomp. of soil 3 C to soil 4 C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL3C_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "SOIL3 C (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL3N", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "SOIL3 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL3N_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "soil 3 N tendency due to vertical transport"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL3N_TO_SOIL4N", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "decomp. of soil 3 N to soil 4 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL3N_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "SOIL3 N (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL3P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "SOIL3 P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL3P_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "soil 3 P tendency due to vertical transport"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL3P_TO_SOIL4P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "decomp. of soil 3 P to soil 4 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL3P_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "SOIL3 P (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL3_HR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "Het. Resp. from soil 3"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL4C", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "SOIL4 C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL4C_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "SOIL4 C (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL4N", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "SOIL4 N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL4N_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "soil 4 N tendency due to vertical transport"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL4N_TO_SMINN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 35, "mineral N flux for decomp. of SOIL4"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL4N_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "SOIL4 N (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL4P", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 7, "SOIL4 P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL4P_TNDNCY_VERT_TRANS", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "soil 4 P tendency due to vertical transport"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^3/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL4P_TO_SMINP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 35, "mineral P flux for decomp. of SOIL4"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL4P_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "SOIL4 P (vertically resolved)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOIL4_HR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "Het. Resp. from soil 4"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOILC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 6, "soil C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOILC_HR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 32, "soil C heterotrophic respiration"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOILC_LOSS", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 11, "soil C loss"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOILICE", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 35, "soil ice (vegetated landunits only)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "kg/m2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOILICE_ICE", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 29, "soil ice (ice landunits only)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "kg/m2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOILLIQ", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 44, "soil liquid water (vegetated landunits only)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "kg/m2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOILLIQ_ICE", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 38, "soil liquid water (ice landunits only)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "kg/m2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOILPSI", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 39, "soil water potential in each soil layer"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 3, "MPa"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOILWATER_10CM", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 64, "soil liquid water + ice in top 10cm of soil (veg landunits only)"); ERR
err = nc_put_att(ncid, varid, "standard_name", NC_CHAR, 35, "mass_content_of_water_in_soil_layer"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "kg/m2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOLUTIONP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 15, "soil solution P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOLUTIONP_vr", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 28, "soil solution P (vert. res.)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gp/m^3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOMHR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 45, "soil organic matter heterotrophic respiration"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SOM_C_LEACHED", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 46, "total flux of C from SOM pools due to leaching"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 39, "total soil respiration (HR + root resp)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "STORVEGC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 41, "stored vegetation carbon, excluding cpool"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "STORVEGN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 26, "stored vegetation nitrogen"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "STORVEGP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 28, "stored vegetation phosphorus"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SUPPLEMENT_TO_SMINN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 21, "supplemental N supply"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SUPPLEMENT_TO_SMINP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 21, "supplemental P supply"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gP/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SUPPLY", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 26, "runoff supply for land use"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 4, "mm/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SoilAlpha", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "factor limiting ground evap"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "1"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "SoilAlpha_U", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "urban factor limiting ground evap"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "1"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TAUX", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 20, "zonal surface stress"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "kg/m/s^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TAUY", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 25, "meridional surface stress"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "kg/m/s^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TBOT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "atmospheric air temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TBUILD", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 35, "internal urban building temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TCS_MONTH_BEGIN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 48, "total carbon storage at the beginning of a month"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "mm"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 11, "time: point"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TCS_MONTH_END", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 42, "total carbon storage at the end of a month"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "mm"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 11, "time: point"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TG", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 18, "ground temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TG_R", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 24, "Rural ground temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TG_U", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 24, "Urban ground temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TH2OSFC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 25, "surface water temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "THBOT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 37, "atmospheric air potential temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TKE1", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 40, "top lake level eddy thermal conductivity"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "W/(mK)"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TLAI", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "total projected leaf area index"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "1"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TLAKE", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 16, "lake temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTCOLC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 62, "total column carbon, incl veg and cpool but excl product pools"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTCOLCH4", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 57, "total belowground CH4, (0 for non-lake special landunits)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "gC/m2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTCOLN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "total column-level N but excl product pools"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTCOLP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "total column-level P but excl product pools"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTECOSYSC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 70, "total ecosystem carbon, incl veg but excl cpool but excl product pools"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTECOSYSN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 40, "total ecosystem N but excl product pools"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTECOSYSP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 40, "total ecosystem P but excl product pools"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTLITC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 19, "total litter carbon"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTLITC_1m", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 36, "total litter carbon to 1 meter depth"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTLITN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 14, "total litter N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTLITP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 14, "total litter P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTLITP_1m", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 25, "total litter P to 1 meter"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTPFTC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 41, "total patch-level carbon, including cpool"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTPFTN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 24, "total PFT-level nitrogen"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTPFTP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 26, "total PFT-level phosphorus"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTSOMC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 32, "total soil organic matter carbon"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTSOMC_1m", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 49, "total soil organic matter carbon to 1 meter depth"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTSOMN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "total soil organic matter N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTSOMP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "total soil organic matter P"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTSOMP_1m", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 38, "total soil organic matter P to 1 meter"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTVEGC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 40, "total vegetation carbon, excluding cpool"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTVEGC_ABG", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 52, "total aboveground vegetation carbon, excluding cpool"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTVEGN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 25, "total vegetation nitrogen"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gN/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TOTVEGP", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "total vegetation phosphorus"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gP/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TREFMNAV", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 40, "daily minimum of average 2-m temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TREFMNAV_R", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 46, "Rural daily minimum of average 2-m temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TREFMNAV_U", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 46, "Urban daily minimum of average 2-m temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TREFMXAV", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 40, "daily maximum of average 2-m temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TREFMXAV_R", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 46, "Rural daily maximum of average 2-m temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TREFMXAV_U", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 46, "Urban daily maximum of average 2-m temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TSA", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 18, "2m air temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TSAI", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 31, "total projected stem area index"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "1"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TSA_R", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 24, "Rural 2m air temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TSA_U", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 24, "Urban 2m air temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TSOI", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 43, "soil temperature (vegetated landunits only)"); ERR
err = nc_put_att(ncid, varid, "standard_name", NC_CHAR, 16, "soil_temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TSOI_10CM", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 36, "soil temperature in top 10cm of soil"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TSOI_ICE", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 37, "soil temperature (ice landunits only)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TV", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 22, "vegetation temperature"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "K"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TWS", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 19, "total water storage"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "mm"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TWS_MONTH_BEGIN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 47, "total water storage at the beginning of a month"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "mm"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 11, "time: point"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "TWS_MONTH_END", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 41, "total water storage at the end of a month"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "mm"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 11, "time: point"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "T_SCALAR", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 39, "temperature inhibition of decomposition"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "1"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "U10", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 9, "10-m wind"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 3, "m/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "URBAN_AC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "urban air conditioning flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "URBAN_HEAT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 18, "urban heating flux"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "VOLR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "river channel total water storage"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "m3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "VOLRMCH", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 40, "river channel main channel water storage"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "m3"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "WA", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 58, "water in the unconfined aquifer (vegetated landunits only)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 2, "mm"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "WASTEHEAT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 67, "sensible heat flux from heating/cooling sources of urban waste heat"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "W/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "WF", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 41, "soil water as frac. of whc for top 0.05 m"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 10, "proportion"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "WIND", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 35, "atmospheric wind velocity magnitude"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 3, "m/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "WOODC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 6, "wood C"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "WOODC_ALLOC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 18, "wood C eallocation"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "WOODC_LOSS", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 11, "wood C loss"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "WOOD_HARVESTC", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 38, "wood harvest carbon (to product pools)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "WOOD_HARVESTN", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 33, "wood harvest N (to product pools)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gN/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "WTGQ", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 26, "surface tracer conductance"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 3, "m/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "W_SCALAR", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 46, "Moisture (dryness) inhibition of decomposition"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "1"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "XR", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 24, "total excess respiration"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "gC/m^2/s"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "XSMRPOOL", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 30, "temporary photosynthate C pool"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 6, "gC/m^2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ZBOT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 28, "atmospheric reference height"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ZWT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 44, "water table depth (vegetated landunits only)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ZWT_CH4_UNSAT", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 70, "depth of water table for methane production used in non-inundated area"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "ZWT_PERCH", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 52, "perched water table depth (vegetated landunits only)"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 1, "m"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "cn_scalar", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 19, "N limitation factor"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 0, ""); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "cp_scalar", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 19, "P limitation factor"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 0, ""); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "leaf_npimbalance", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 57, "leaf np imbalance partial C partial P/partial C partial N"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 5, "gN/gP"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "nlim_m", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "runmean N limitation factor"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 0, ""); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "o2_decomp_depth_unsat", NC_FLOAT, 4, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 21, "o2_decomp_depth_unsat"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 8, "mol/m3/2"); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR

err = nc_def_var(ncid, "plim_m", NC_FLOAT, 3, dimids, &varid); ERR
err = nc_put_att(ncid, varid, "long_name", NC_CHAR, 27, "runmean P limitation factor"); ERR
err = nc_put_att(ncid, varid, "units", NC_CHAR, 0, ""); ERR
err = nc_put_att(ncid, varid, "cell_methods", NC_CHAR, 10, "time: mean"); ERR
err = nc_put_att(ncid, varid, "_FillValue", NC_FLOAT, 1, buf); ERR
err = nc_put_att(ncid, varid, "missing_value", NC_FLOAT, 1, buf); ERR
}
