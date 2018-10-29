/*********************************************************************
 *
 * Copyright (C) 2018, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 *
 * This program is part of the E3SM I/O benchmark.
 *
 *********************************************************************/

#include <assert.h>

#include <e3sm_io.h>

static int lat, lon, area, lev, hyam, hybm, P0, ilev, hyai, hybi, time, date,
        datesec, time_bnds, date_written, time_written, ndbase, nsbase, nbdate,
        nbsec, mdt, ndcur, nscur, co2vmr, ch4vmr, n2ovmr, f11vmr, f12vmr,
        sol_tsi, nsteph, AEROD_v, ANRAIN, ANSNOW, AODABS, AODABSBC, AODBC,
        AODDUST, AODDUST1, AODDUST3, AODDUST4, AODMODE1, AODMODE2, AODMODE3,
        AODMODE4, AODNIR, AODPOM, AODSO4, AODSOA, AODSS, AODUV, AODVIS, AQRAIN,
        AQSNOW, AQ_DMS, AQ_H2O2, AQ_H2SO4, AQ_O3, AQ_SO2, AQ_SOAG, AREI, AREL,
        AWNC, AWNI, CCN3, CDNUMC, CLDHGH, CLDICE, CLDLIQ, CLDLOW, CLDMED,
        CLDTOT, CLOUD, CLOUDFRAC_CLUBB, CONCLD, DCQ, DF_DMS, DF_H2O2, DF_H2SO4,
        DF_O3, DF_SO2, DF_SOAG, DMS_SRF, DP_KCLDBASE, DP_MFUP_MAX, DP_WCLDBASE,
        DSTSFMBL, DTCOND, DTENDTH, DTENDTQ, EXTINCT, FICE, FLDS, FLNS, FLNSC,
        FLNT, FLNTC, FLUT, FLUTC, FREQI, FREQL, FREQR, FREQS, FSDS, FSDSC,
        FSNS, FSNSC, FSNT, FSNTC, FSNTOA, FSNTOAC, FSUTOA, FSUTOAC, F_eff,
        H2O2_SRF, H2SO4_SRF, H2SO4_sfgaex1, ICEFRAC, ICIMR, ICWMR, IWC,
        LANDFRAC, LHFLX, LINOZ_DO3, LINOZ_DO3_PSC, LINOZ_O3CLIM, LINOZ_O3COL,
        LINOZ_SSO3, LINOZ_SZA, LND_MBL, LWCF, Mass_bc, Mass_dst, Mass_mom,
        Mass_ncl, Mass_pom, Mass_so4, Mass_soa, NUMICE, NUMLIQ, NUMRAI, NUMSNO,
        O3, O3_SRF, OCNFRAC, OMEGA, OMEGA500, OMEGAT, PBLH, PHIS, PRECC, PRECL,
        PRECSC, PRECSL, PS, PSL, Q, QFLX, QREFHT, QRL, QRS, RAINQM, RAM1,
        RELHUM, SFDMS, SFH2O2, SFH2SO4, SFO3, SFSO2, SFSOAG, SFbc_a1, SFbc_a3,
        SFbc_a4, SFdst_a1, SFdst_a3, SFmom_a1, SFmom_a2, SFmom_a3, SFmom_a4,
        SFncl_a1, SFncl_a2, SFncl_a3, SFnum_a1, SFnum_a2, SFnum_a3, SFnum_a4,
        SFpom_a1, SFpom_a3, SFpom_a4, SFso4_a1, SFso4_a2, SFso4_a3, SFsoa_a1,
        SFsoa_a2, SFsoa_a3, SHFLX, SH_KCLDBASE, SH_MFUP_MAX, SH_WCLDBASE,
        SNOWHICE, SNOWHLND, SNOWQM, SO2, SO2_CLXF, SO2_SRF, SOAG_CLXF,
        SOAG_SRF, SOAG_sfgaex1, SOLIN, SSAVIS, SSTSFMBL, SSTSFMBL_OM, SWCF, T,
        TAUGWX, TAUGWY, TAUX, TAUY, TGCLDCWP, TGCLDIWP, TGCLDLWP, TH7001000,
        TMQ, TREFHT, TROP_P, TROP_T, TS, TSMN, TSMX, TUH, TUQ, TVH, TVQ, U,
        U10, UU, V, VQ, VT, VU, VV, WD_H2O2, WD_H2SO4, WD_SO2, WSUB, Z3,
        aero_water, airFV, bc_a1DDF, bc_a1SFWET, bc_a1_SRF, bc_a1_sfgaex1,
        bc_a3DDF, bc_a3SFWET, bc_a3_SRF, bc_a4DDF, bc_a4SFWET, bc_a4_CLXF,
        bc_a4_SRF, bc_a4_sfgaex1, bc_c1DDF, bc_c1SFWET, bc_c3DDF, bc_c3SFWET,
        bc_c4DDF, bc_c4SFWET, chla, dst_a1DDF, dst_a1SF, dst_a1SFWET,
        dst_a1_SRF, dst_a3DDF, dst_a3SF, dst_a3SFWET, dst_a3_SRF, dst_c1DDF,
        dst_c1SFWET, dst_c3DDF, dst_c3SFWET, hstobie_linoz, mlip, mom_a1DDF,
        mom_a1SF, mom_a1SFWET, mom_a1_SRF, mom_a1_sfgaex1, mom_a2DDF, mom_a2SF,
        mom_a2SFWET, mom_a2_SRF, mom_a3DDF, mom_a3SFWET, mom_a3_SRF, mom_a4DDF,
        mom_a4SF, mom_a4SFWET, mom_a4_SRF, mom_a4_sfgaex1, mom_c1DDF,
        mom_c1SFWET, mom_c2DDF, mom_c2SFWET, mom_c3DDF, mom_c3SFWET, mom_c4DDF,
        mom_c4SFWET, mpoly, mprot, ncl_a1DDF, ncl_a1SF, ncl_a1SFWET,
        ncl_a1_SRF, ncl_a2DDF, ncl_a2SF, ncl_a2SFWET, ncl_a2_SRF, ncl_a3DDF,
        ncl_a3SF, ncl_a3SFWET, ncl_a3_SRF, ncl_c1DDF, ncl_c1SFWET, ncl_c2DDF,
        ncl_c2SFWET, ncl_c3DDF, ncl_c3SFWET, num_a1DDF, num_a1SF, num_a1SFWET,
        num_a1_CLXF, num_a1_SRF, num_a1_sfgaex1, num_a2DDF, num_a2SFWET,
        num_a2_CLXF, num_a2_SRF, num_a3DDF, num_a3SF, num_a3SFWET, num_a3_SRF,
        num_a4DDF, num_a4SFWET, num_a4_CLXF, num_a4_SRF, num_a4_sfgaex1,
        num_c1DDF, num_c1SFWET, num_c2DDF, num_c2SFWET, num_c3DDF, num_c3SFWET,
        num_c4DDF, num_c4SFWET, pom_a1DDF, pom_a1SFWET, pom_a1_SRF,
        pom_a1_sfgaex1, pom_a3DDF, pom_a3SFWET, pom_a3_SRF, pom_a4DDF,
        pom_a4SFWET, pom_a4_CLXF, pom_a4_SRF, pom_a4_sfgaex1, pom_c1DDF,
        pom_c1SFWET, pom_c3DDF, pom_c3SFWET, pom_c4DDF, pom_c4SFWET, so4_a1DDF,
        so4_a1SFWET, so4_a1_CLXF, so4_a1_SRF, so4_a1_sfgaex1, so4_a2DDF,
        so4_a2SFWET, so4_a2_CLXF, so4_a2_SRF, so4_a2_sfgaex1, so4_a3DDF,
        so4_a3SFWET, so4_a3_SRF, so4_a3_sfgaex1, so4_c1DDF, so4_c1SFWET,
        so4_c2DDF, so4_c2SFWET, so4_c3DDF, so4_c3SFWET, soa_a1DDF, soa_a1SFWET,
        soa_a1_SRF, soa_a1_sfgaex1, soa_a2DDF, soa_a2SFWET, soa_a2_SRF,
        soa_a2_sfgaex1, soa_a3DDF, soa_a3SFWET, soa_a3_SRF, soa_a3_sfgaex1,
        soa_c1DDF, soa_c1SFWET, soa_c2DDF, soa_c2SFWET, soa_c3DDF, soa_c3SFWET;

/*----< e3sm_io_header() >---------------------------------------------------*/
int e3sm_io_header(int         ncid,    /* file ID */
                   MPI_Offset  dims[2], /* dimension sizes */
                   int         nvars,   /* number of variables */
                   int        *varids)  /* variable IDs */
{
    int i, err, nerrs=0, dimids[3], iattr, mdims=1;
    int dim_ncol, dim_time, dim_nbnd, dim_chars, dim_lev, dim_ilev;
    float fillv=1.e+36f, missv = 1.e+36f;

    /* global attributes: */
    iattr = 4;
    err = ncmpi_put_att(ncid, NC_GLOBAL, "ne", NC_INT, 1, &iattr); ERR
    err = ncmpi_put_att(ncid, NC_GLOBAL, "np", NC_INT, 1, &iattr); ERR
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, "CF-1.0"); ERR
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "source", 3, "CAM"); ERR
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "case", 20, "FC5AV1C-H01B_ne4_ne4"); ERR
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "title", 5, "UNSET"); ERR
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "logname", 6, "wkliao"); ERR
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "host", 10, "compute001"); ERR
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "Version", 6, "$Name$"); ERR
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "revision_Id", 4, "$Id$"); ERR
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "initial_file", 86,
    "/home/climate1/acme/inputdata/atm/cam/inic/homme/cami_mam3_Linoz_ne4np4_L72_c160909.nc"); ERR
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "topography_file", 79,
    "/home/climate1/acme/inputdata/atm/cam/topo/USGS-gtopo30_ne4np4_16x.c20160612.nc"); ERR
    err = ncmpi_put_att_text(ncid, NC_GLOBAL, "time_period_freq", 5, "day_5"); ERR

    /* define dimensions */
    err = ncmpi_def_dim(ncid, "ncol", dims[1],      &dim_ncol); ERR
    err = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &dim_time); ERR
    err = ncmpi_def_dim(ncid, "nbnd",  2,           &dim_nbnd); ERR
    err = ncmpi_def_dim(ncid, "chars", 8,           &dim_chars); ERR
    err = ncmpi_def_dim(ncid, "lev",   dims[0],     &dim_lev); ERR
    err = ncmpi_def_dim(ncid, "ilev",  dims[0]+1,   &dim_ilev); ERR

    i = 0;

    /* define variables */
    dimids[0] = dim_ncol;
    err = ncmpi_def_var(ncid, "lat", NC_DOUBLE, 1, dimids, &lat); ERR
    err = ncmpi_put_att_text(ncid, lat, "long_name", 8, "latitude"); ERR
    err = ncmpi_put_att_text(ncid, lat, "units", 13, "degrees_north"); ERR
    varids[i++] = lat;

    dimids[0] = dim_ncol;
    err = ncmpi_def_var(ncid, "lon", NC_DOUBLE, 1, dimids, &lon); ERR
    err = ncmpi_put_att_text(ncid, lon, "long_name", 9, "longitude"); ERR
    err = ncmpi_put_att_text(ncid, lon, "units", 12, "degrees_east"); ERR
    varids[i++] = lon;

    dimids[0] = dim_ncol;
    err = ncmpi_def_var(ncid, "area", NC_DOUBLE, 1, dimids, &area); ERR
    err = ncmpi_put_att_text(ncid, area, "long_name", 14, "gll grid areas"); ERR
    varids[i++] = area;

    dimids[0] = dim_lev;
    err = ncmpi_def_var(ncid, "lev", NC_DOUBLE, 1, dimids, &lev); ERR
    err = ncmpi_put_att_text(ncid, lev, "long_name", 38, "hybrid level at midpoints (1000*(A+B))"); ERR
    err = ncmpi_put_att_text(ncid, lev, "units", 3, "hPa"); ERR
    err = ncmpi_put_att_text(ncid, lev, "positive", 4, "down"); ERR
    err = ncmpi_put_att_text(ncid, lev, "standard_name", 43, "atmosphere_hybrid_sigma_pressure_coordinate"); ERR
    err = ncmpi_put_att_text(ncid, lev, "formula_terms", 29, "a: hyam b: hybm p0: P0 ps: PS"); ERR
    varids[i++] = lev;

    dimids[0] = dim_lev;
    err = ncmpi_def_var(ncid, "hyam", NC_DOUBLE, 1, dimids, &hyam); ERR
    err = ncmpi_put_att_text(ncid, hyam, "long_name", 39, "hybrid A coefficient at layer midpoints"); ERR
    varids[i++] = hyam;

    dimids[0] = dim_lev;
    err = ncmpi_def_var(ncid, "hybm", NC_DOUBLE, 1, dimids, &hybm); ERR
    err = ncmpi_put_att_text(ncid, hybm, "long_name", 39, "hybrid B coefficient at layer midpoints"); ERR
    varids[i++] = hybm;

    dimids[0] = dim_lev;
    err = ncmpi_def_var(ncid, "P0", NC_DOUBLE, 0, NULL, &P0); ERR
    err = ncmpi_put_att_text(ncid, P0, "long_name", 18, "reference pressure"); ERR
    err = ncmpi_put_att_text(ncid, P0, "units", 2, "Pa"); ERR
    varids[i++] = P0;

    dimids[0] = dim_ilev;
    err = ncmpi_def_var(ncid, "ilev", NC_DOUBLE, 1, dimids, &ilev); ERR
    err = ncmpi_put_att_text(ncid, ilev, "long_name", 39, "hybrid level at interfaces (1000*(A+B))"); ERR
    err = ncmpi_put_att_text(ncid, ilev, "units", 3, "hPa"); ERR
    err = ncmpi_put_att_text(ncid, ilev, "positive", 4, "down"); ERR
    err = ncmpi_put_att_text(ncid, ilev, "standard_name", 43, "atmosphere_hybrid_sigma_pressure_coordinate"); ERR
    err = ncmpi_put_att_text(ncid, ilev, "formula_terms", 29, "a: hyai b: hybi p0: P0 ps: PS"); ERR
    varids[i++] = ilev;

    dimids[0] = dim_ilev;
    err = ncmpi_def_var(ncid, "hyai", NC_DOUBLE, 1, dimids, &hyai); ERR
    err = ncmpi_put_att_text(ncid, hyai, "long_name", 40, "hybrid A coefficient at layer interfaces"); ERR
    varids[i++] = hyai;

    dimids[0] = dim_ilev;
    err = ncmpi_def_var(ncid, "hybi", NC_DOUBLE, 1, dimids, &hybi); ERR
    err = ncmpi_put_att_text(ncid, hybi, "long_name", 40, "hybrid B coefficient at layer interfaces"); ERR
    varids[i++] = hybi;

    dimids[0] = dim_time;
    err = ncmpi_def_var(ncid, "time", NC_DOUBLE, 1, dimids, &time); ERR
    err = ncmpi_put_att_text(ncid, time, "long_name", 4, "time"); ERR
    err = ncmpi_put_att_text(ncid, time, "units", 30, "days since 0001-01-01 00:00:00"); ERR
    err = ncmpi_put_att_text(ncid, time, "calendar", 6, "noleap"); ERR
    err = ncmpi_put_att_text(ncid, time, "bounds", 9, "time_bnds"); ERR
    varids[i++] = time;

    dimids[0] = dim_time;
    err = ncmpi_def_var(ncid, "date", NC_INT, 1, dimids, &date); ERR
    err = ncmpi_put_att_text(ncid, date, "long_name", 23, "current date (YYYYMMDD)"); ERR
    varids[i++] = date;

    dimids[0] = dim_time;
    err = ncmpi_def_var(ncid, "datesec", NC_INT, 1, dimids, &datesec); ERR
    err = ncmpi_put_att_text(ncid, datesec, "long_name", 31, "current seconds of current date"); ERR
    varids[i++] = datesec;

    dimids[0] = dim_time;
    dimids[1] = dim_nbnd;
    err = ncmpi_def_var(ncid, "time_bnds", NC_DOUBLE, 2, dimids, &time_bnds); ERR
    err = ncmpi_put_att_text(ncid, time_bnds, "long_name", 23, "time interval endpoints"); ERR
    varids[i++] = time_bnds;

    dimids[0] = dim_time;
    dimids[1] = dim_chars;
    err = ncmpi_def_var(ncid, "date_written", NC_CHAR, 2, dimids, &date_written); ERR
    varids[i++] = date_written;

    dimids[0] = dim_time;
    dimids[1] = dim_chars;
    err = ncmpi_def_var(ncid, "time_written", NC_CHAR, 2, dimids, &time_written); ERR
    varids[i++] = time_written;

    err = ncmpi_def_var(ncid, "ndbase", NC_INT, 0, NULL, &ndbase); ERR
    err = ncmpi_put_att_text(ncid, ndbase, "long_name", 8, "base day"); ERR 
    varids[i++] = ndbase;
    err = ncmpi_def_var(ncid, "nsbase", NC_INT, 0, NULL, &nsbase); ERR
    err = ncmpi_put_att_text(ncid, nsbase, "long_name", 19, "seconds of base day"); ERR
    varids[i++] = nsbase;

    err = ncmpi_def_var(ncid, "nbdate", NC_INT, 0, NULL, &nbdate); ERR
    err = ncmpi_put_att_text(ncid, nbdate, "long_name", 20, "base date (YYYYMMDD)"); ERR
    varids[i++] = nbdate;

    err = ncmpi_def_var(ncid, "nbsec", NC_INT, 0, NULL, &nbsec); ERR
    err = ncmpi_put_att_text(ncid, nbsec, "long_name", 20, "seconds of base date"); ERR
    varids[i++] = nbsec;

    err = ncmpi_def_var(ncid, "mdt", NC_INT, 0, NULL, &mdt); ERR
    err = ncmpi_put_att_text(ncid, mdt, "long_name", 8, "timestep"); ERR
    err = ncmpi_put_att_text(ncid, mdt, "units", 1, "s"); ERR
    varids[i++] = mdt;

    dimids[0] = dim_time;
    err = ncmpi_def_var(ncid, "ndcur", NC_INT, 1, dimids, &ndcur); ERR
    err = ncmpi_put_att_text(ncid, ndcur, "long_name", 27, "current day (from base day)"); ERR
    varids[i++] = ndcur;

    dimids[0] = dim_time;
    err = ncmpi_def_var(ncid, "nscur", NC_INT, 1, dimids, &nscur); ERR
    err = ncmpi_put_att_text(ncid, nscur, "long_name", 30, "current seconds of current day"); ERR
    varids[i++] = nscur;

    dimids[0] = dim_time;
    err = ncmpi_def_var(ncid, "co2vmr", NC_DOUBLE, 1, dimids, &co2vmr); ERR
    err = ncmpi_put_att_text(ncid, co2vmr, "long_name", 23, "co2 volume mixing ratio"); ERR
    varids[i++] = co2vmr;

    dimids[0] = dim_time;
    err = ncmpi_def_var(ncid, "ch4vmr", NC_DOUBLE, 1, dimids, &ch4vmr); ERR
    err = ncmpi_put_att_text(ncid, ch4vmr, "long_name", 23, "ch4 volume mixing ratio"); ERR
    varids[i++] = ch4vmr;

    dimids[0] = dim_time;
    err = ncmpi_def_var(ncid, "n2ovmr", NC_DOUBLE, 1, dimids, &n2ovmr); ERR
    err = ncmpi_put_att_text(ncid, n2ovmr, "long_name", 23, "n2o volume mixing ratio"); ERR
    varids[i++] = n2ovmr;

    dimids[0] = dim_time;
    err = ncmpi_def_var(ncid, "f11vmr", NC_DOUBLE, 1, dimids, &f11vmr); ERR
    err = ncmpi_put_att_text(ncid, f11vmr, "long_name", 23, "f11 volume mixing ratio"); ERR
    varids[i++] = f11vmr;

    dimids[0] = dim_time;
    err = ncmpi_def_var(ncid, "f12vmr", NC_DOUBLE, 1, dimids, &f12vmr); ERR
    err = ncmpi_put_att_text(ncid, f12vmr, "long_name", 23, "f12 volume mixing ratio"); ERR
    varids[i++] = f12vmr;

    dimids[0] = dim_time;
    err = ncmpi_def_var(ncid, "sol_tsi", NC_DOUBLE, 1, dimids, &sol_tsi); ERR
    err = ncmpi_put_att_text(ncid, sol_tsi, "long_name", 22, "total solar irradiance"); ERR
    err = ncmpi_put_att_text(ncid, sol_tsi, "units", 4, "W/m2"); ERR
    varids[i++] = sol_tsi;

    dimids[0] = dim_time;
    err = ncmpi_def_var(ncid, "nsteph", NC_INT, 1, dimids, &nsteph); ERR
    err = ncmpi_put_att_text(ncid, nsteph, "long_name", 16, "current timestep"); ERR
    varids[i++] = nsteph;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AEROD_v", NC_FLOAT, 2, dimids, &AEROD_v); ERR
    err = ncmpi_put_att_float(ncid, AEROD_v, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AEROD_v, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AEROD_v, "units", 1, "1"); ERR
    err = ncmpi_put_att_text(ncid, AEROD_v, "long_name", 43, "Total Aerosol Optical Depth in visible band"); ERR
    err = ncmpi_put_att_text(ncid, AEROD_v, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AEROD_v;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "ANRAIN", NC_FLOAT, 3, dimids, &ANRAIN); ERR
    err = ncmpi_put_att_int(ncid, ANRAIN, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, ANRAIN, "units", 3, "m-3"); ERR
    err = ncmpi_put_att_text(ncid, ANRAIN, "long_name", 24, "Average rain number conc"); ERR
    err = ncmpi_put_att_text(ncid, ANRAIN, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ANRAIN;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "ANSNOW", NC_FLOAT, 3, dimids, &ANSNOW); ERR
    err = ncmpi_put_att_int(ncid, ANSNOW, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, ANSNOW, "units", 3, "m-3"); ERR
    err = ncmpi_put_att_text(ncid, ANSNOW, "long_name", 24, "Average snow number conc"); ERR
    err = ncmpi_put_att_text(ncid, ANSNOW, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ANSNOW;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODABS", NC_FLOAT, 2, dimids, &AODABS); ERR
    err = ncmpi_put_att_float(ncid, AODABS, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODABS, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODABS, "long_name", 39, "Aerosol absorption optical depth 550 nm"); ERR
    err = ncmpi_put_att_text(ncid, AODABS, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODABS;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODABSBC", NC_FLOAT, 2, dimids, &AODABSBC); ERR
    err = ncmpi_put_att_float(ncid, AODABSBC, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODABSBC, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODABSBC, "long_name", 48, "Aerosol absorption optical depth 550 nm from BC"); ERR
    err = ncmpi_put_att_text(ncid, AODABSBC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODABSBC;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODBC", NC_FLOAT, 2, dimids, &AODBC); ERR
    err = ncmpi_put_att_float(ncid, AODBC, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODBC, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODBC, "long_name", 37, "Aerosol optical depth 550 nm from BC"); ERR
    err = ncmpi_put_att_text(ncid, AODBC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODBC;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODDUST", NC_FLOAT, 2, dimids, &AODDUST); ERR
    err = ncmpi_put_att_float(ncid, AODDUST, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODDUST, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODDUST, "long_name", 38, "Aerosol optical depth 550 nm from dust"); ERR
    err = ncmpi_put_att_text(ncid, AODDUST, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODDUST;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODDUST1", NC_FLOAT, 2, dimids, &AODDUST1); ERR
    err = ncmpi_put_att_float(ncid, AODDUST1, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODDUST1, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODDUST1, "long_name", 46, "Aerosol optical depth 550 nm model 1 from dust"); ERR
    err = ncmpi_put_att_text(ncid, AODDUST1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODDUST1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODDUST3", NC_FLOAT, 2, dimids, &AODDUST3); ERR
    err = ncmpi_put_att_float(ncid, AODDUST3, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODDUST3, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODDUST3, "long_name", 46, "Aerosol optical depth 550 nm model 3 from dust"); ERR
    err = ncmpi_put_att_text(ncid, AODDUST3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODDUST3;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODDUST4", NC_FLOAT, 2, dimids, &AODDUST4); ERR
    err = ncmpi_put_att_float(ncid, AODDUST4, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODDUST4, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODDUST4, "long_name", 46, "Aerosol optical depth 550 nm model 4 from dust"); ERR
    err = ncmpi_put_att_text(ncid, AODDUST4, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODDUST4;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODMODE1", NC_FLOAT, 2, dimids, &AODMODE1); ERR
    err = ncmpi_put_att_float(ncid, AODMODE1, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODMODE1, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODMODE1, "long_name", 35, "Aerosol optical depth 550 nm mode 1"); ERR
    err = ncmpi_put_att_text(ncid, AODMODE1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODMODE1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODMODE2", NC_FLOAT, 2, dimids, &AODMODE2); ERR
    err = ncmpi_put_att_float(ncid, AODMODE2, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODMODE2, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODMODE2, "long_name", 35, "Aerosol optical depth 550 nm mode 2"); ERR
    err = ncmpi_put_att_text(ncid, AODMODE2, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODMODE2;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODMODE3", NC_FLOAT, 2, dimids, &AODMODE3); ERR
    err = ncmpi_put_att_float(ncid, AODMODE3, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODMODE3, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODMODE3, "long_name", 35, "Aerosol optical depth 550 nm mode 3"); ERR
    err = ncmpi_put_att_text(ncid, AODMODE3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODMODE3;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODMODE4", NC_FLOAT, 2, dimids, &AODMODE4); ERR
    err = ncmpi_put_att_float(ncid, AODMODE4, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODMODE4, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODMODE4, "long_name", 35, "Aerosol optical depth 550 nm mode 4"); ERR
    err = ncmpi_put_att_text(ncid, AODMODE4, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODMODE4;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODNIR", NC_FLOAT, 2, dimids, &AODNIR); ERR
    err = ncmpi_put_att_float(ncid, AODNIR, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODNIR, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODNIR, "long_name", 28, "Aerosol optical depth 850 nm"); ERR
    err = ncmpi_put_att_text(ncid, AODNIR, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODNIR;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODPOM", NC_FLOAT, 2, dimids, &AODPOM); ERR
    err = ncmpi_put_att_float(ncid, AODPOM, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODPOM, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODPOM, "long_name", 37, "Aerosol optical depth 550 nm from POM"); ERR
    err = ncmpi_put_att_text(ncid, AODPOM, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODPOM;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODSO4", NC_FLOAT, 2, dimids, &AODSO4); ERR
    err = ncmpi_put_att_float(ncid, AODSO4, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODSO4, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODSO4, "long_name", 37, "Aerosol optical depth 550 nm from SO4"); ERR
    err = ncmpi_put_att_text(ncid, AODSO4, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODSO4;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODSOA", NC_FLOAT, 2, dimids, &AODSOA); ERR
    err = ncmpi_put_att_float(ncid, AODSOA, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODSOA, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODSOA, "long_name", 37, "Aerosol optical depth 550 nm from SOA"); ERR
    err = ncmpi_put_att_text(ncid, AODSOA, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODSOA;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODSS", NC_FLOAT, 2, dimids, &AODSS); ERR
    err = ncmpi_put_att_float(ncid, AODSS, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODSS, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODSS, "long_name", 41, "Aerosol optical depth 550 nm from seasalt"); ERR
    err = ncmpi_put_att_text(ncid, AODSS, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODSS;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODUV", NC_FLOAT, 2, dimids, &AODUV); ERR
    err = ncmpi_put_att_float(ncid, AODUV, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODUV, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODUV, "long_name", 28, "Aerosol optical depth 350 nm"); ERR
    err = ncmpi_put_att_text(ncid, AODUV, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODUV;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AODVIS", NC_FLOAT, 2, dimids, &AODVIS); ERR
    err = ncmpi_put_att_float(ncid, AODVIS, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, AODVIS, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, AODVIS, "long_name", 28, "Aerosol optical depth 550 nm"); ERR
    err = ncmpi_put_att_text(ncid, AODVIS, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AODVIS;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "AQRAIN", NC_FLOAT, 3, dimids, &AQRAIN); ERR
    err = ncmpi_put_att_int(ncid, AQRAIN, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, AQRAIN, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, AQRAIN, "long_name", 25, "Average rain mixing ratio"); ERR
    err = ncmpi_put_att_text(ncid, AQRAIN, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AQRAIN;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "AQSNOW", NC_FLOAT, 3, dimids, &AQSNOW); ERR
    err = ncmpi_put_att_int(ncid, AQSNOW, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, AQSNOW, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, AQSNOW, "long_name", 25, "Average snow mixing ratio"); ERR
    err = ncmpi_put_att_text(ncid, AQSNOW, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AQSNOW;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AQ_DMS", NC_FLOAT, 2, dimids, &AQ_DMS); ERR
    err = ncmpi_put_att_text(ncid, AQ_DMS, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, AQ_DMS, "long_name", 39, "DMS aqueous chemistry (for gas species)"); ERR
    err = ncmpi_put_att_text(ncid, AQ_DMS, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AQ_DMS;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AQ_H2O2", NC_FLOAT, 2, dimids, &AQ_H2O2); ERR
    err = ncmpi_put_att_text(ncid, AQ_H2O2, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, AQ_H2O2, "long_name", 40, "H2O2 aqueous chemistry (for gas species)"); ERR
    err = ncmpi_put_att_text(ncid, AQ_H2O2, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AQ_H2O2;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AQ_H2SO4", NC_FLOAT, 2, dimids, &AQ_H2SO4); ERR
    err = ncmpi_put_att_text(ncid, AQ_H2SO4, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, AQ_H2SO4, "long_name", 41, "H2SO4 aqueous chemistry (for gas species)"); ERR
    err = ncmpi_put_att_text(ncid, AQ_H2SO4, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AQ_H2SO4;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AQ_O3", NC_FLOAT, 2, dimids, &AQ_O3); ERR
    err = ncmpi_put_att_text(ncid, AQ_O3, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, AQ_O3, "long_name", 38, "O3 aqueous chemistry (for gas species)"); ERR
    err = ncmpi_put_att_text(ncid, AQ_O3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AQ_O3;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AQ_SO2", NC_FLOAT, 2, dimids, &AQ_SO2); ERR
    err = ncmpi_put_att_text(ncid, AQ_SO2, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, AQ_SO2, "long_name", 39, "SO2 aqueous chemistry (for gas species)"); ERR
    err = ncmpi_put_att_text(ncid, AQ_SO2, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AQ_SO2;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "AQ_SOAG", NC_FLOAT, 2, dimids, &AQ_SOAG); ERR
    err = ncmpi_put_att_text(ncid, AQ_SOAG, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, AQ_SOAG, "long_name", 40, "SOAG aqueous chemistry (for gas species)"); ERR
    err = ncmpi_put_att_text(ncid, AQ_SOAG, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AQ_SOAG;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "AREI", NC_FLOAT, 3, dimids, &AREI); ERR
    err = ncmpi_put_att_int(ncid, AREI, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, AREI, "units", 6, "Micron"); ERR
    err = ncmpi_put_att_text(ncid, AREI, "long_name", 28, "Average ice effective radius"); ERR
    err = ncmpi_put_att_text(ncid, AREI, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AREI;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "AREL", NC_FLOAT, 3, dimids, &AREL); ERR
    err = ncmpi_put_att_int(ncid, AREL, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, AREL, "units", 6, "Micron"); ERR
    err = ncmpi_put_att_text(ncid, AREL, "long_name", 32, "Average droplet effective radius"); ERR
    err = ncmpi_put_att_text(ncid, AREL, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AREL;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "AWNC", NC_FLOAT, 3, dimids, &AWNC); ERR
    err = ncmpi_put_att_int(ncid, AWNC, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, AWNC, "units", 3, "m-3"); ERR
    err = ncmpi_put_att_text(ncid, AWNC, "long_name", 31, "Average cloud water number conc"); ERR
    err = ncmpi_put_att_text(ncid, AWNC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AWNC;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "AWNI", NC_FLOAT, 3, dimids, &AWNI); ERR
    err = ncmpi_put_att_int(ncid, AWNI, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, AWNI, "units", 3, "m-3"); ERR
    err = ncmpi_put_att_text(ncid, AWNI, "long_name", 29, "Average cloud ice number conc"); ERR
    err = ncmpi_put_att_text(ncid, AWNI, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = AWNI;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "CCN3", NC_FLOAT, 3, dimids, &CCN3); ERR
    err = ncmpi_put_att_int(ncid, CCN3, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, CCN3, "units", 5, "#/cm3"); ERR
    err = ncmpi_put_att_text(ncid, CCN3, "long_name", 27, "CCN concentration at S=0.1%"); ERR
    err = ncmpi_put_att_text(ncid, CCN3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = CCN3;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "CDNUMC", NC_FLOAT, 2, dimids, &CDNUMC); ERR
    err = ncmpi_put_att_text(ncid, CDNUMC, "units", 4, "1/m2"); ERR
    err = ncmpi_put_att_text(ncid, CDNUMC, "long_name", 43, "Vertically-integrated droplet concentration"); ERR
    err = ncmpi_put_att_text(ncid, CDNUMC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = CDNUMC;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "CLDHGH", NC_FLOAT, 2, dimids, &CLDHGH); ERR
    err = ncmpi_put_att_text(ncid, CLDHGH, "units", 8, "fraction"); ERR
    err = ncmpi_put_att_text(ncid, CLDHGH, "long_name", 32, "Vertically-integrated high cloud"); ERR
    err = ncmpi_put_att_text(ncid, CLDHGH, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = CLDHGH;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "CLDICE", NC_FLOAT, 3, dimids, &CLDICE); ERR
    err = ncmpi_put_att_int(ncid, CLDICE, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, CLDICE, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, CLDICE, "long_name", 34, "Grid box averaged cloud ice amount"); ERR
    err = ncmpi_put_att_text(ncid, CLDICE, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = CLDICE;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "CLDLIQ", NC_FLOAT, 3, dimids, &CLDLIQ); ERR
    err = ncmpi_put_att_int(ncid, CLDLIQ, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, CLDLIQ, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, CLDLIQ, "long_name", 37, "Grid box averaged cloud liquid amount"); ERR
    err = ncmpi_put_att_text(ncid, CLDLIQ, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = CLDLIQ;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "CLDLOW", NC_FLOAT, 2, dimids, &CLDLOW); ERR
    err = ncmpi_put_att_text(ncid, CLDLOW, "units", 8, "fraction"); ERR
    err = ncmpi_put_att_text(ncid, CLDLOW, "long_name", 31, "Vertically-integrated low cloud"); ERR
    err = ncmpi_put_att_text(ncid, CLDLOW, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = CLDLOW;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "CLDMED", NC_FLOAT, 2, dimids, &CLDMED); ERR
    err = ncmpi_put_att_text(ncid, CLDMED, "units", 8, "fraction"); ERR
    err = ncmpi_put_att_text(ncid, CLDMED, "long_name", 37, "Vertically-integrated mid-level cloud"); ERR
    err = ncmpi_put_att_text(ncid, CLDMED, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = CLDMED;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "CLDTOT", NC_FLOAT, 2, dimids, &CLDTOT); ERR
    err = ncmpi_put_att_text(ncid, CLDTOT, "units", 8, "fraction"); ERR
    err = ncmpi_put_att_text(ncid, CLDTOT, "long_name", 33, "Vertically-integrated total cloud"); ERR
    err = ncmpi_put_att_text(ncid, CLDTOT, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = CLDTOT;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "CLOUD", NC_FLOAT, 3, dimids, &CLOUD); ERR
    err = ncmpi_put_att_int(ncid, CLOUD, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, CLOUD, "units", 8, "fraction"); ERR
    err = ncmpi_put_att_text(ncid, CLOUD, "long_name", 14, "Cloud fraction"); ERR
    err = ncmpi_put_att_text(ncid, CLOUD, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = CLOUD;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "CLOUDFRAC_CLUBB", NC_FLOAT, 3, dimids, &CLOUDFRAC_CLUBB); ERR
    err = ncmpi_put_att_int(ncid, CLOUDFRAC_CLUBB, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, CLOUDFRAC_CLUBB, "units", 8, "fraction"); ERR
    err = ncmpi_put_att_text(ncid, CLOUDFRAC_CLUBB, "long_name", 14, "Cloud Fraction"); ERR
    err = ncmpi_put_att_text(ncid, CLOUDFRAC_CLUBB, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = CLOUDFRAC_CLUBB;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "CONCLD", NC_FLOAT, 3, dimids, &CONCLD); ERR
    err = ncmpi_put_att_int(ncid, CONCLD, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, CONCLD, "units", 8, "fraction"); ERR
    err = ncmpi_put_att_text(ncid, CONCLD, "long_name", 22, "Convective cloud cover"); ERR
    err = ncmpi_put_att_text(ncid, CONCLD, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = CONCLD;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "DCQ", NC_FLOAT, 3, dimids, &DCQ); ERR
    err = ncmpi_put_att_int(ncid, DCQ, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, DCQ, "units", 7, "kg/kg/s"); ERR
    err = ncmpi_put_att_text(ncid, DCQ, "long_name", 33, "Q tendency due to moist processes"); ERR
    err = ncmpi_put_att_text(ncid, DCQ, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = DCQ;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "DF_DMS", NC_FLOAT, 2, dimids, &DF_DMS); ERR
    err = ncmpi_put_att_text(ncid, DF_DMS, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, DF_DMS, "long_name", 19, "dry deposition flux"); ERR
    err = ncmpi_put_att_text(ncid, DF_DMS, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = DF_DMS;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "DF_H2O2", NC_FLOAT, 2, dimids, &DF_H2O2); ERR
    err = ncmpi_put_att_text(ncid, DF_H2O2, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, DF_H2O2, "long_name", 19, "dry deposition flux"); ERR
    err = ncmpi_put_att_text(ncid, DF_H2O2, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = DF_H2O2;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "DF_H2SO4", NC_FLOAT, 2, dimids, &DF_H2SO4); ERR
    err = ncmpi_put_att_text(ncid, DF_H2SO4, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, DF_H2SO4, "long_name", 19, "dry deposition flux"); ERR
    err = ncmpi_put_att_text(ncid, DF_H2SO4, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = DF_H2SO4;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "DF_O3", NC_FLOAT, 2, dimids, &DF_O3); ERR
    err = ncmpi_put_att_text(ncid, DF_O3, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, DF_O3, "long_name", 19, "dry deposition flux"); ERR
    err = ncmpi_put_att_text(ncid, DF_O3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = DF_O3;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "DF_SO2", NC_FLOAT, 2, dimids, &DF_SO2); ERR
    err = ncmpi_put_att_text(ncid, DF_SO2, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, DF_SO2, "long_name", 19, "dry deposition flux"); ERR
    err = ncmpi_put_att_text(ncid, DF_SO2, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = DF_SO2;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "DF_SOAG", NC_FLOAT, 2, dimids, &DF_SOAG); ERR
    err = ncmpi_put_att_text(ncid, DF_SOAG, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, DF_SOAG, "long_name", 19, "dry deposition flux"); ERR
    err = ncmpi_put_att_text(ncid, DF_SOAG, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = DF_SOAG;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "DMS_SRF", NC_FLOAT, 2, dimids, &DMS_SRF); ERR
    err = ncmpi_put_att_text(ncid, DMS_SRF, "units", 7, "mol/mol"); ERR
    err = ncmpi_put_att_text(ncid, DMS_SRF, "long_name", 19, "DMS in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, DMS_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = DMS_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "DP_KCLDBASE", NC_FLOAT, 2, dimids, &DP_KCLDBASE); ERR
    err = ncmpi_put_att_text(ncid, DP_KCLDBASE, "units", 1, "1"); ERR
    err = ncmpi_put_att_text(ncid, DP_KCLDBASE, "long_name", 32, "Deep conv. cloudbase level index"); ERR
    err = ncmpi_put_att_text(ncid, DP_KCLDBASE, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = DP_KCLDBASE;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "DP_MFUP_MAX", NC_FLOAT, 2, dimids, &DP_MFUP_MAX); ERR
    err = ncmpi_put_att_text(ncid, DP_MFUP_MAX, "units", 5, "kg/m2"); ERR
    err = ncmpi_put_att_text(ncid, DP_MFUP_MAX, "long_name", 39, "Deep conv. column-max updraft mass flux"); ERR
    err = ncmpi_put_att_text(ncid, DP_MFUP_MAX, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = DP_MFUP_MAX;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "DP_WCLDBASE", NC_FLOAT, 2, dimids, &DP_WCLDBASE); ERR
    err = ncmpi_put_att_text(ncid, DP_WCLDBASE, "units", 3, "m/s"); ERR
    err = ncmpi_put_att_text(ncid, DP_WCLDBASE, "long_name", 38, "Deep conv. cloudbase vertical velocity"); ERR
    err = ncmpi_put_att_text(ncid, DP_WCLDBASE, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = DP_WCLDBASE;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "DSTSFMBL", NC_FLOAT, 2, dimids, &DSTSFMBL); ERR
    err = ncmpi_put_att_text(ncid, DSTSFMBL, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, DSTSFMBL, "long_name", 28, "Mobilization flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, DSTSFMBL, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = DSTSFMBL;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "DTCOND", NC_FLOAT, 3, dimids, &DTCOND); ERR
    err = ncmpi_put_att_int(ncid, DTCOND, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, DTCOND, "units", 3, "K/s"); ERR
    err = ncmpi_put_att_text(ncid, DTCOND, "long_name", 28, "T tendency - moist processes"); ERR
    err = ncmpi_put_att_text(ncid, DTCOND, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = DTCOND;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "DTENDTH", NC_FLOAT, 2, dimids, &DTENDTH); ERR
    err = ncmpi_put_att_text(ncid, DTENDTH, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, DTENDTH, "long_name", 69, "Dynamic Tendency of Total (vertically integrated) moist static energy"); ERR
    err = ncmpi_put_att_text(ncid, DTENDTH, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = DTENDTH;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "DTENDTQ", NC_FLOAT, 2, dimids, &DTENDTQ); ERR
    err = ncmpi_put_att_text(ncid, DTENDTQ, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, DTENDTQ, "long_name", 67, "Dynamic Tendency of Total (vertically integrated) specific humidity"); ERR
    err = ncmpi_put_att_text(ncid, DTENDTQ, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = DTENDTQ;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "EXTINCT", NC_FLOAT, 3, dimids, &EXTINCT); ERR
    err = ncmpi_put_att_int(ncid, EXTINCT, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_float(ncid, EXTINCT, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, EXTINCT, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, EXTINCT, "units", 2, "/m"); ERR
    err = ncmpi_put_att_text(ncid, EXTINCT, "long_name", 18, "Aerosol extinction"); ERR
    err = ncmpi_put_att_text(ncid, EXTINCT, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = EXTINCT;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "FICE", NC_FLOAT, 3, dimids, &FICE); ERR
    err = ncmpi_put_att_int(ncid, FICE, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, FICE, "units", 8, "fraction"); ERR
    err = ncmpi_put_att_text(ncid, FICE, "long_name", 35, "Fractional ice content within cloud"); ERR
    err = ncmpi_put_att_text(ncid, FICE, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FICE;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "FLDS", NC_FLOAT, 2, dimids, &FLDS); ERR
    err = ncmpi_put_att_text(ncid, FLDS, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, FLDS, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, FLDS, "long_name", 36, "Downwelling longwave flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, FLDS, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FLDS;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "FLNS", NC_FLOAT, 2, dimids, &FLNS); ERR
    err = ncmpi_put_att_text(ncid, FLNS, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, FLNS, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, FLNS, "long_name", 28, "Net longwave flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, FLNS, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FLNS;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "FLNSC", NC_FLOAT, 2, dimids, &FLNSC); ERR
    err = ncmpi_put_att_text(ncid, FLNSC, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, FLNSC, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, FLNSC, "long_name", 37, "Clearsky net longwave flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, FLNSC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FLNSC;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "FLNT", NC_FLOAT, 2, dimids, &FLNT); ERR
    err = ncmpi_put_att_text(ncid, FLNT, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, FLNT, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, FLNT, "long_name", 33, "Net longwave flux at top of model"); ERR
    err = ncmpi_put_att_text(ncid, FLNT, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FLNT;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "FLNTC", NC_FLOAT, 2, dimids, &FLNTC); ERR
    err = ncmpi_put_att_text(ncid, FLNTC, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, FLNTC, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, FLNTC, "long_name", 42, "Clearsky net longwave flux at top of model"); ERR
    err = ncmpi_put_att_text(ncid, FLNTC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FLNTC;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "FLUT", NC_FLOAT, 2, dimids, &FLUT); ERR
    err = ncmpi_put_att_text(ncid, FLUT, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, FLUT, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, FLUT, "long_name", 39, "Upwelling longwave flux at top of model"); ERR
    err = ncmpi_put_att_text(ncid, FLUT, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FLUT;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "FLUTC", NC_FLOAT, 2, dimids, &FLUTC); ERR
    err = ncmpi_put_att_text(ncid, FLUTC, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, FLUTC, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, FLUTC, "long_name", 48, "Clearsky upwelling longwave flux at top of model"); ERR
    err = ncmpi_put_att_text(ncid, FLUTC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FLUTC;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "FREQI", NC_FLOAT, 3, dimids, &FREQI); ERR
    err = ncmpi_put_att_int(ncid, FREQI, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, FREQI, "units", 8, "fraction"); ERR
    err = ncmpi_put_att_text(ncid, FREQI, "long_name", 28, "Fractional occurrence of ice"); ERR
    err = ncmpi_put_att_text(ncid, FREQI, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FREQI;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "FREQL", NC_FLOAT, 3, dimids, &FREQL); ERR
    err = ncmpi_put_att_int(ncid, FREQL, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, FREQL, "units", 8, "fraction"); ERR
    err = ncmpi_put_att_text(ncid, FREQL, "long_name", 31, "Fractional occurrence of liquid"); ERR
    err = ncmpi_put_att_text(ncid, FREQL, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FREQL;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "FREQR", NC_FLOAT, 3, dimids, &FREQR); ERR
    err = ncmpi_put_att_int(ncid, FREQR, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, FREQR, "units", 8, "fraction"); ERR
    err = ncmpi_put_att_text(ncid, FREQR, "long_name", 29, "Fractional occurrence of rain"); ERR
    err = ncmpi_put_att_text(ncid, FREQR, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FREQR;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "FREQS", NC_FLOAT, 3, dimids, &FREQS); ERR
    err = ncmpi_put_att_int(ncid, FREQS, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, FREQS, "units", 8, "fraction"); ERR
    err = ncmpi_put_att_text(ncid, FREQS, "long_name", 29, "Fractional occurrence of snow"); ERR
    err = ncmpi_put_att_text(ncid, FREQS, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FREQS;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "FSDS", NC_FLOAT, 2, dimids, &FSDS); ERR
    err = ncmpi_put_att_text(ncid, FSDS, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, FSDS, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, FSDS, "long_name", 33, "Downwelling solar flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, FSDS, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FSDS;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "FSDSC", NC_FLOAT, 2, dimids, &FSDSC); ERR
    err = ncmpi_put_att_text(ncid, FSDSC, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, FSDSC, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, FSDSC, "long_name", 42, "Clearsky downwelling solar flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, FSDSC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FSDSC;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "FSNS", NC_FLOAT, 2, dimids, &FSNS); ERR
    err = ncmpi_put_att_text(ncid, FSNS, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, FSNS, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, FSNS, "long_name", 25, "Net solar flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, FSNS, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FSNS;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "FSNSC", NC_FLOAT, 2, dimids, &FSNSC); ERR
    err = ncmpi_put_att_text(ncid, FSNSC, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, FSNSC, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, FSNSC, "long_name", 34, "Clearsky net solar flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, FSNSC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FSNSC;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "FSNT", NC_FLOAT, 2, dimids, &FSNT); ERR
    err = ncmpi_put_att_text(ncid, FSNT, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, FSNT, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, FSNT, "long_name", 30, "Net solar flux at top of model"); ERR
    err = ncmpi_put_att_text(ncid, FSNT, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FSNT;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "FSNTC", NC_FLOAT, 2, dimids, &FSNTC); ERR
    err = ncmpi_put_att_text(ncid, FSNTC, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, FSNTC, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, FSNTC, "long_name", 39, "Clearsky net solar flux at top of model"); ERR
    err = ncmpi_put_att_text(ncid, FSNTC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FSNTC;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "FSNTOA", NC_FLOAT, 2, dimids, &FSNTOA); ERR
    err = ncmpi_put_att_text(ncid, FSNTOA, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, FSNTOA, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, FSNTOA, "long_name", 35, "Net solar flux at top of atmosphere"); ERR
    err = ncmpi_put_att_text(ncid, FSNTOA, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FSNTOA;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "FSNTOAC", NC_FLOAT, 2, dimids, &FSNTOAC); ERR
    err = ncmpi_put_att_text(ncid, FSNTOAC, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, FSNTOAC, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, FSNTOAC, "long_name", 44, "Clearsky net solar flux at top of atmosphere"); ERR
    err = ncmpi_put_att_text(ncid, FSNTOAC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FSNTOAC;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "FSUTOA", NC_FLOAT, 2, dimids, &FSUTOA); ERR
    err = ncmpi_put_att_text(ncid, FSUTOA, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, FSUTOA, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, FSUTOA, "long_name", 41, "Upwelling solar flux at top of atmosphere"); ERR
    err = ncmpi_put_att_text(ncid, FSUTOA, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FSUTOA;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "FSUTOAC", NC_FLOAT, 2, dimids, &FSUTOAC); ERR
    err = ncmpi_put_att_text(ncid, FSUTOAC, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, FSUTOAC, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, FSUTOAC, "long_name", 50, "Clearsky upwelling solar flux at top of atmosphere"); ERR
    err = ncmpi_put_att_text(ncid, FSUTOAC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = FSUTOAC;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "F_eff", NC_FLOAT, 2, dimids, &F_eff); ERR
    err = ncmpi_put_att_text(ncid, F_eff, "units", 1, "1"); ERR
    err = ncmpi_put_att_text(ncid, F_eff, "long_name", 52, "Effective enrichment factor of marine organic matter"); ERR
    err = ncmpi_put_att_text(ncid, F_eff, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = F_eff;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "H2O2_SRF", NC_FLOAT, 2, dimids, &H2O2_SRF); ERR
    err = ncmpi_put_att_text(ncid, H2O2_SRF, "units", 7, "mol/mol"); ERR
    err = ncmpi_put_att_text(ncid, H2O2_SRF, "long_name", 20, "H2O2 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, H2O2_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = H2O2_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "H2SO4_SRF", NC_FLOAT, 2, dimids, &H2SO4_SRF); ERR
    err = ncmpi_put_att_text(ncid, H2SO4_SRF, "units", 7, "mol/mol"); ERR
    err = ncmpi_put_att_text(ncid, H2SO4_SRF, "long_name", 21, "H2SO4 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, H2SO4_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = H2SO4_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "H2SO4_sfgaex1", NC_FLOAT, 2, dimids, &H2SO4_sfgaex1); ERR
    err = ncmpi_put_att_text(ncid, H2SO4_sfgaex1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, H2SO4_sfgaex1, "long_name", 50, "H2SO4 gas-aerosol-exchange primary column tendency"); ERR
    err = ncmpi_put_att_text(ncid, H2SO4_sfgaex1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = H2SO4_sfgaex1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ICEFRAC", NC_FLOAT, 2, dimids, &ICEFRAC); ERR
    err = ncmpi_put_att_text(ncid, ICEFRAC, "units", 8, "fraction"); ERR
    err = ncmpi_put_att_text(ncid, ICEFRAC, "long_name", 39, "Fraction of sfc area covered by sea-ice"); ERR
    err = ncmpi_put_att_text(ncid, ICEFRAC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ICEFRAC;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "ICIMR", NC_FLOAT, 3, dimids, &ICIMR); ERR
    err = ncmpi_put_att_int(ncid, ICIMR, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, ICIMR, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, ICIMR, "long_name", 36, "Prognostic in-cloud ice mixing ratio"); ERR
    err = ncmpi_put_att_text(ncid, ICIMR, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ICIMR;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "ICWMR", NC_FLOAT, 3, dimids, &ICWMR); ERR
    err = ncmpi_put_att_int(ncid, ICWMR, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, ICWMR, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, ICWMR, "long_name", 38, "Prognostic in-cloud water mixing ratio"); ERR
    err = ncmpi_put_att_text(ncid, ICWMR, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ICWMR;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "IWC", NC_FLOAT, 3, dimids, &IWC); ERR
    err = ncmpi_put_att_int(ncid, IWC, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, IWC, "units", 5, "kg/m3"); ERR
    err = ncmpi_put_att_text(ncid, IWC, "long_name", 34, "Grid box average ice water content"); ERR
    err = ncmpi_put_att_text(ncid, IWC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = IWC;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "LANDFRAC", NC_FLOAT, 2, dimids, &LANDFRAC); ERR
    err = ncmpi_put_att_text(ncid, LANDFRAC, "units", 8, "fraction"); ERR
    err = ncmpi_put_att_text(ncid, LANDFRAC, "long_name", 36, "Fraction of sfc area covered by land"); ERR
    err = ncmpi_put_att_text(ncid, LANDFRAC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = LANDFRAC;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "LHFLX", NC_FLOAT, 2, dimids, &LHFLX); ERR
    err = ncmpi_put_att_text(ncid, LHFLX, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, LHFLX, "long_name", 24, "Surface latent heat flux"); ERR
    err = ncmpi_put_att_text(ncid, LHFLX, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = LHFLX;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "LINOZ_DO3", NC_FLOAT, 3, dimids, &LINOZ_DO3); ERR
    err = ncmpi_put_att_int(ncid, LINOZ_DO3, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_DO3, "units", 2, "/s"); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_DO3, "long_name", 48, "ozone vmr tendency by linearized ozone chemistry"); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_DO3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = LINOZ_DO3;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "LINOZ_DO3_PSC", NC_FLOAT, 3, dimids, &LINOZ_DO3_PSC); ERR
    err = ncmpi_put_att_int(ncid, LINOZ_DO3_PSC, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_DO3_PSC, "units", 2, "/s"); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_DO3_PSC, "long_name", 50, "ozone vmr loss by PSCs using Carille et al. (1990)"); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_DO3_PSC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = LINOZ_DO3_PSC;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "LINOZ_O3CLIM", NC_FLOAT, 3, dimids, &LINOZ_O3CLIM); ERR
    err = ncmpi_put_att_int(ncid, LINOZ_O3CLIM, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_O3CLIM, "units", 7, "mol/mol"); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_O3CLIM, "long_name", 29, "climatology of ozone in LINOZ"); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_O3CLIM, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = LINOZ_O3CLIM;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "LINOZ_O3COL", NC_FLOAT, 3, dimids, &LINOZ_O3COL); ERR
    err = ncmpi_put_att_int(ncid, LINOZ_O3COL, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_O3COL, "units", 2, "DU"); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_O3COL, "long_name", 18, "ozone column above"); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_O3COL, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = LINOZ_O3COL;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "LINOZ_SSO3", NC_FLOAT, 3, dimids, &LINOZ_SSO3); ERR
    err = ncmpi_put_att_int(ncid, LINOZ_SSO3, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_SSO3, "units", 2, "kg"); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_SSO3, "long_name", 27, "steady state ozone in LINOZ"); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_SSO3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = LINOZ_SSO3;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "LINOZ_SZA", NC_FLOAT, 2, dimids, &LINOZ_SZA); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_SZA, "units", 7, "degrees"); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_SZA, "long_name", 27, "solar zenith angle in LINOZ"); ERR
    err = ncmpi_put_att_text(ncid, LINOZ_SZA, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = LINOZ_SZA;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "LND_MBL", NC_FLOAT, 2, dimids, &LND_MBL); ERR
    err = ncmpi_put_att_text(ncid, LND_MBL, "units", 4, "frac"); ERR
    err = ncmpi_put_att_text(ncid, LND_MBL, "long_name", 23, "Soil erodibility factor"); ERR
    err = ncmpi_put_att_text(ncid, LND_MBL, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = LND_MBL;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "LWCF", NC_FLOAT, 2, dimids, &LWCF); ERR
    err = ncmpi_put_att_text(ncid, LWCF, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, LWCF, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, LWCF, "long_name", 22, "Longwave cloud forcing"); ERR
    err = ncmpi_put_att_text(ncid, LWCF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = LWCF;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "Mass_bc", NC_FLOAT, 3, dimids, &Mass_bc); ERR
    err = ncmpi_put_att_int(ncid, Mass_bc, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, Mass_bc, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, Mass_bc, "long_name", 64, "sum of bc mass concentration bc_a1+bc_c1+bc_a3+bc_c3+bc_a4+bc_c4"); ERR
    err = ncmpi_put_att_text(ncid, Mass_bc, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = Mass_bc;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "Mass_dst", NC_FLOAT, 3, dimids, &Mass_dst); ERR
    err = ncmpi_put_att_int(ncid, Mass_dst, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, Mass_dst, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, Mass_dst, "long_name", 57, "sum of dst mass concentration dst_a1+dst_c1+dst_a3+dst_c3"); ERR
    err = ncmpi_put_att_text(ncid, Mass_dst, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = Mass_dst;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "Mass_mom", NC_FLOAT, 3, dimids, &Mass_mom); ERR
    err = ncmpi_put_att_int(ncid, Mass_mom, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, Mass_mom, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, Mass_mom, "long_name", 85, "sum of mom mass concentration mom_a1+mom_c1+mom_a2+mom_c2+mom_a3+mom_c3+mom_a4+mom_c4"); ERR
    err = ncmpi_put_att_text(ncid, Mass_mom, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = Mass_mom;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "Mass_ncl", NC_FLOAT, 3, dimids, &Mass_ncl); ERR
    err = ncmpi_put_att_int(ncid, Mass_ncl, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, Mass_ncl, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, Mass_ncl, "long_name", 71, "sum of ncl mass concentration ncl_a1+ncl_c1+ncl_a2+ncl_c2+ncl_a3+ncl_c3"); ERR
    err = ncmpi_put_att_text(ncid, Mass_ncl, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = Mass_ncl;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "Mass_pom", NC_FLOAT, 3, dimids, &Mass_pom); ERR
    err = ncmpi_put_att_int(ncid, Mass_pom, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, Mass_pom, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, Mass_pom, "long_name", 71, "sum of pom mass concentration pom_a1+pom_c1+pom_a3+pom_c3+pom_a4+pom_c4"); ERR
    err = ncmpi_put_att_text(ncid, Mass_pom, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = Mass_pom;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "Mass_so4", NC_FLOAT, 3, dimids, &Mass_so4); ERR
    err = ncmpi_put_att_int(ncid, Mass_so4, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, Mass_so4, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, Mass_so4, "long_name", 71, "sum of so4 mass concentration so4_a1+so4_c1+so4_a2+so4_c2+so4_a3+so4_c3"); ERR
    err = ncmpi_put_att_text(ncid, Mass_so4, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = Mass_so4;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "Mass_soa", NC_FLOAT, 3, dimids, &Mass_soa); ERR
    err = ncmpi_put_att_int(ncid, Mass_soa, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, Mass_soa, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, Mass_soa, "long_name", 71, "sum of soa mass concentration soa_a1+soa_c1+soa_a2+soa_c2+soa_a3+soa_c3"); ERR
    err = ncmpi_put_att_text(ncid, Mass_soa, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = Mass_soa;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "NUMICE", NC_FLOAT, 3, dimids, &NUMICE); ERR
    err = ncmpi_put_att_int(ncid, NUMICE, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, NUMICE, "units", 4, "1/kg"); ERR
    err = ncmpi_put_att_text(ncid, NUMICE, "long_name", 34, "Grid box averaged cloud ice number"); ERR
    err = ncmpi_put_att_text(ncid, NUMICE, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = NUMICE;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "NUMLIQ", NC_FLOAT, 3, dimids, &NUMLIQ); ERR
    err = ncmpi_put_att_int(ncid, NUMLIQ, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, NUMLIQ, "units", 4, "1/kg"); ERR
    err = ncmpi_put_att_text(ncid, NUMLIQ, "long_name", 37, "Grid box averaged cloud liquid number"); ERR
    err = ncmpi_put_att_text(ncid, NUMLIQ, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = NUMLIQ;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "NUMRAI", NC_FLOAT, 3, dimids, &NUMRAI); ERR
    err = ncmpi_put_att_int(ncid, NUMRAI, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, NUMRAI, "units", 4, "1/kg"); ERR
    err = ncmpi_put_att_text(ncid, NUMRAI, "long_name", 29, "Grid box averaged rain number"); ERR
    err = ncmpi_put_att_text(ncid, NUMRAI, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = NUMRAI;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "NUMSNO", NC_FLOAT, 3, dimids, &NUMSNO); ERR
    err = ncmpi_put_att_int(ncid, NUMSNO, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, NUMSNO, "units", 4, "1/kg"); ERR
    err = ncmpi_put_att_text(ncid, NUMSNO, "long_name", 29, "Grid box averaged snow number"); ERR
    err = ncmpi_put_att_text(ncid, NUMSNO, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = NUMSNO;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "O3", NC_FLOAT, 3, dimids, &O3); ERR
    err = ncmpi_put_att_int(ncid, O3, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, O3, "units", 7, "mol/mol"); ERR
    err = ncmpi_put_att_text(ncid, O3, "long_name", 16, "O3 concentration"); ERR
    err = ncmpi_put_att_text(ncid, O3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = O3;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "O3_SRF", NC_FLOAT, 2, dimids, &O3_SRF); ERR
    err = ncmpi_put_att_text(ncid, O3_SRF, "units", 7, "mol/mol"); ERR
    err = ncmpi_put_att_text(ncid, O3_SRF, "long_name", 18, "O3 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, O3_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = O3_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "OCNFRAC", NC_FLOAT, 2, dimids, &OCNFRAC); ERR
    err = ncmpi_put_att_text(ncid, OCNFRAC, "units", 8, "fraction"); ERR
    err = ncmpi_put_att_text(ncid, OCNFRAC, "long_name", 37, "Fraction of sfc area covered by ocean"); ERR
    err = ncmpi_put_att_text(ncid, OCNFRAC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = OCNFRAC;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "OMEGA", NC_FLOAT, 3, dimids, &OMEGA); ERR
    err = ncmpi_put_att_int(ncid, OMEGA, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, OMEGA, "units", 4, "Pa/s"); ERR
    err = ncmpi_put_att_text(ncid, OMEGA, "long_name", 28, "Vertical velocity (pressure)"); ERR
    err = ncmpi_put_att_text(ncid, OMEGA, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = OMEGA;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "OMEGA500", NC_FLOAT, 2, dimids, &OMEGA500); ERR
    err = ncmpi_put_att_text(ncid, OMEGA500, "units", 4, "Pa/s"); ERR
    err = ncmpi_put_att_text(ncid, OMEGA500, "long_name", 46, "Vertical velocity at 500 mbar pressure surface"); ERR
    err = ncmpi_put_att_text(ncid, OMEGA500, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = OMEGA500;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "OMEGAT", NC_FLOAT, 3, dimids, &OMEGAT); ERR
    err = ncmpi_put_att_int(ncid, OMEGAT, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, OMEGAT, "units", 6, "K Pa/s"); ERR
    err = ncmpi_put_att_text(ncid, OMEGAT, "long_name", 18, "Vertical heat flux"); ERR
    err = ncmpi_put_att_text(ncid, OMEGAT, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = OMEGAT;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "PBLH", NC_FLOAT, 2, dimids, &PBLH); ERR
    err = ncmpi_put_att_text(ncid, PBLH, "units", 1, "m"); ERR
    err = ncmpi_put_att_text(ncid, PBLH, "long_name", 10, "PBL height"); ERR
    err = ncmpi_put_att_text(ncid, PBLH, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = PBLH;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "PHIS", NC_FLOAT, 2, dimids, &PHIS); ERR
    err = ncmpi_put_att_text(ncid, PHIS, "units", 5, "m2/s2"); ERR
    err = ncmpi_put_att_text(ncid, PHIS, "long_name", 20, "Surface geopotential"); ERR
    varids[i++] = PHIS;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "PRECC", NC_FLOAT, 2, dimids, &PRECC); ERR
    err = ncmpi_put_att_text(ncid, PRECC, "units", 3, "m/s"); ERR
    err = ncmpi_put_att_text(ncid, PRECC, "long_name", 41, "Convective precipitation rate (liq + ice)"); ERR
    err = ncmpi_put_att_text(ncid, PRECC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = PRECC;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "PRECL", NC_FLOAT, 2, dimids, &PRECL); ERR
    err = ncmpi_put_att_text(ncid, PRECL, "units", 3, "m/s"); ERR
    err = ncmpi_put_att_text(ncid, PRECL, "long_name", 51, "Large-scale (stable) precipitation rate (liq + ice)"); ERR
    err = ncmpi_put_att_text(ncid, PRECL, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = PRECL;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "PRECSC", NC_FLOAT, 2, dimids, &PRECSC); ERR
    err = ncmpi_put_att_text(ncid, PRECSC, "units", 3, "m/s"); ERR
    err = ncmpi_put_att_text(ncid, PRECSC, "long_name", 39, "Convective snow rate (water equivalent)"); ERR
    err = ncmpi_put_att_text(ncid, PRECSC, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = PRECSC;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "PRECSL", NC_FLOAT, 2, dimids, &PRECSL); ERR
    err = ncmpi_put_att_text(ncid, PRECSL, "units", 3, "m/s"); ERR
    err = ncmpi_put_att_text(ncid, PRECSL, "long_name", 49, "Large-scale (stable) snow rate (water equivalent)"); ERR
    err = ncmpi_put_att_text(ncid, PRECSL, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = PRECSL;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "PS", NC_FLOAT, 2, dimids, &PS); ERR
    err = ncmpi_put_att_text(ncid, PS, "units", 2, "Pa"); ERR
    err = ncmpi_put_att_text(ncid, PS, "long_name", 16, "Surface pressure"); ERR
    err = ncmpi_put_att_text(ncid, PS, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = PS;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "PSL", NC_FLOAT, 2, dimids, &PSL); ERR
    err = ncmpi_put_att_text(ncid, PSL, "units", 2, "Pa"); ERR
    err = ncmpi_put_att_text(ncid, PSL, "long_name", 18, "Sea level pressure"); ERR
    err = ncmpi_put_att_text(ncid, PSL, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = PSL;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "Q", NC_FLOAT, 3, dimids, &Q); ERR
    err = ncmpi_put_att_int(ncid, Q, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, Q, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, Q, "long_name", 17, "Specific humidity"); ERR
    err = ncmpi_put_att_text(ncid, Q, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = Q;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "QFLX", NC_FLOAT, 2, dimids, &QFLX); ERR
    err = ncmpi_put_att_text(ncid, QFLX, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, QFLX, "long_name", 18, "Surface water flux"); ERR
    err = ncmpi_put_att_text(ncid, QFLX, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = QFLX;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "QREFHT", NC_FLOAT, 2, dimids, &QREFHT); ERR
    err = ncmpi_put_att_text(ncid, QREFHT, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, QREFHT, "long_name", 25, "Reference height humidity"); ERR
    err = ncmpi_put_att_text(ncid, QREFHT, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = QREFHT;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "QRL", NC_FLOAT, 3, dimids, &QRL); ERR
    err = ncmpi_put_att_int(ncid, QRL, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, QRL, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, QRL, "units", 3, "K/s"); ERR
    err = ncmpi_put_att_text(ncid, QRL, "long_name", 21, "Longwave heating rate"); ERR
    err = ncmpi_put_att_text(ncid, QRL, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = QRL;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "QRS", NC_FLOAT, 3, dimids, &QRS); ERR
    err = ncmpi_put_att_int(ncid, QRS, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, QRS, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, QRS, "units", 3, "K/s"); ERR
    err = ncmpi_put_att_text(ncid, QRS, "long_name", 18, "Solar heating rate"); ERR
    err = ncmpi_put_att_text(ncid, QRS, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = QRS;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "RAINQM", NC_FLOAT, 3, dimids, &RAINQM); ERR
    err = ncmpi_put_att_int(ncid, RAINQM, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, RAINQM, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, RAINQM, "long_name", 29, "Grid box averaged rain amount"); ERR
    err = ncmpi_put_att_text(ncid, RAINQM, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = RAINQM;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "RAM1", NC_FLOAT, 2, dimids, &RAM1); ERR
    err = ncmpi_put_att_text(ncid, RAM1, "units", 4, "frac"); ERR
    err = ncmpi_put_att_text(ncid, RAM1, "long_name", 4, "RAM1"); ERR
    err = ncmpi_put_att_text(ncid, RAM1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = RAM1;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "RELHUM", NC_FLOAT, 3, dimids, &RELHUM); ERR
    err = ncmpi_put_att_int(ncid, RELHUM, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, RELHUM, "units", 7, "percent"); ERR
    err = ncmpi_put_att_text(ncid, RELHUM, "long_name", 17, "Relative humidity"); ERR
    err = ncmpi_put_att_text(ncid, RELHUM, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = RELHUM;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFDMS", NC_FLOAT, 2, dimids, &SFDMS); ERR
    err = ncmpi_put_att_text(ncid, SFDMS, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFDMS, "long_name", 16, "DMS surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFDMS, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFDMS;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFH2O2", NC_FLOAT, 2, dimids, &SFH2O2); ERR
    err = ncmpi_put_att_text(ncid, SFH2O2, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFH2O2, "long_name", 17, "H2O2 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFH2O2, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFH2O2;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFH2SO4", NC_FLOAT, 2, dimids, &SFH2SO4); ERR
    err = ncmpi_put_att_text(ncid, SFH2SO4, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFH2SO4, "long_name", 18, "H2SO4 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFH2SO4, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFH2SO4;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFO3", NC_FLOAT, 2, dimids, &SFO3); ERR
    err = ncmpi_put_att_text(ncid, SFO3, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFO3, "long_name", 15, "O3 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFO3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFO3;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFSO2", NC_FLOAT, 2, dimids, &SFSO2); ERR
    err = ncmpi_put_att_text(ncid, SFSO2, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFSO2, "long_name", 16, "SO2 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFSO2, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFSO2;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFSOAG", NC_FLOAT, 2, dimids, &SFSOAG); ERR
    err = ncmpi_put_att_text(ncid, SFSOAG, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFSOAG, "long_name", 17, "SOAG surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFSOAG, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFSOAG;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFbc_a1", NC_FLOAT, 2, dimids, &SFbc_a1); ERR
    err = ncmpi_put_att_text(ncid, SFbc_a1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFbc_a1, "long_name", 18, "bc_a1 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFbc_a1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFbc_a1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFbc_a3", NC_FLOAT, 2, dimids, &SFbc_a3); ERR
    err = ncmpi_put_att_text(ncid, SFbc_a3, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFbc_a3, "long_name", 18, "bc_a3 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFbc_a3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFbc_a3;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFbc_a4", NC_FLOAT, 2, dimids, &SFbc_a4); ERR
    err = ncmpi_put_att_text(ncid, SFbc_a4, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFbc_a4, "long_name", 18, "bc_a4 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFbc_a4, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFbc_a4;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFdst_a1", NC_FLOAT, 2, dimids, &SFdst_a1); ERR
    err = ncmpi_put_att_text(ncid, SFdst_a1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFdst_a1, "long_name", 19, "dst_a1 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFdst_a1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFdst_a1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFdst_a3", NC_FLOAT, 2, dimids, &SFdst_a3); ERR
    err = ncmpi_put_att_text(ncid, SFdst_a3, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFdst_a3, "long_name", 19, "dst_a3 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFdst_a3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFdst_a3;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFmom_a1", NC_FLOAT, 2, dimids, &SFmom_a1); ERR
    err = ncmpi_put_att_text(ncid, SFmom_a1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFmom_a1, "long_name", 19, "mom_a1 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFmom_a1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFmom_a1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFmom_a2", NC_FLOAT, 2, dimids, &SFmom_a2); ERR
    err = ncmpi_put_att_text(ncid, SFmom_a2, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFmom_a2, "long_name", 19, "mom_a2 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFmom_a2, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFmom_a2;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFmom_a3", NC_FLOAT, 2, dimids, &SFmom_a3); ERR
    err = ncmpi_put_att_text(ncid, SFmom_a3, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFmom_a3, "long_name", 19, "mom_a3 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFmom_a3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFmom_a3;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFmom_a4", NC_FLOAT, 2, dimids, &SFmom_a4); ERR
    err = ncmpi_put_att_text(ncid, SFmom_a4, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFmom_a4, "long_name", 19, "mom_a4 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFmom_a4, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFmom_a4;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFncl_a1", NC_FLOAT, 2, dimids, &SFncl_a1); ERR
    err = ncmpi_put_att_text(ncid, SFncl_a1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFncl_a1, "long_name", 19, "ncl_a1 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFncl_a1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFncl_a1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFncl_a2", NC_FLOAT, 2, dimids, &SFncl_a2); ERR
    err = ncmpi_put_att_text(ncid, SFncl_a2, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFncl_a2, "long_name", 19, "ncl_a2 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFncl_a2, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFncl_a2;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFncl_a3", NC_FLOAT, 2, dimids, &SFncl_a3); ERR
    err = ncmpi_put_att_text(ncid, SFncl_a3, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFncl_a3, "long_name", 19, "ncl_a3 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFncl_a3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFncl_a3;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFnum_a1", NC_FLOAT, 2, dimids, &SFnum_a1); ERR
    err = ncmpi_put_att_text(ncid, SFnum_a1, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFnum_a1, "long_name", 19, "num_a1 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFnum_a1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFnum_a1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFnum_a2", NC_FLOAT, 2, dimids, &SFnum_a2); ERR
    err = ncmpi_put_att_text(ncid, SFnum_a2, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFnum_a2, "long_name", 19, "num_a2 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFnum_a2, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFnum_a2;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFnum_a3", NC_FLOAT, 2, dimids, &SFnum_a3); ERR
    err = ncmpi_put_att_text(ncid, SFnum_a3, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFnum_a3, "long_name", 19, "num_a3 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFnum_a3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFnum_a3;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFnum_a4", NC_FLOAT, 2, dimids, &SFnum_a4); ERR
    err = ncmpi_put_att_text(ncid, SFnum_a4, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFnum_a4, "long_name", 19, "num_a4 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFnum_a4, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFnum_a4;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFpom_a1", NC_FLOAT, 2, dimids, &SFpom_a1); ERR
    err = ncmpi_put_att_text(ncid, SFpom_a1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFpom_a1, "long_name", 19, "pom_a1 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFpom_a1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFpom_a1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFpom_a3", NC_FLOAT, 2, dimids, &SFpom_a3); ERR
    err = ncmpi_put_att_text(ncid, SFpom_a3, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFpom_a3, "long_name", 19, "pom_a3 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFpom_a3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFpom_a3;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFpom_a4", NC_FLOAT, 2, dimids, &SFpom_a4); ERR
    err = ncmpi_put_att_text(ncid, SFpom_a4, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFpom_a4, "long_name", 19, "pom_a4 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFpom_a4, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFpom_a4;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFso4_a1", NC_FLOAT, 2, dimids, &SFso4_a1); ERR
    err = ncmpi_put_att_text(ncid, SFso4_a1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFso4_a1, "long_name", 19, "so4_a1 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFso4_a1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFso4_a1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFso4_a2", NC_FLOAT, 2, dimids, &SFso4_a2); ERR
    err = ncmpi_put_att_text(ncid, SFso4_a2, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFso4_a2, "long_name", 19, "so4_a2 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFso4_a2, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFso4_a2;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFso4_a3", NC_FLOAT, 2, dimids, &SFso4_a3); ERR
    err = ncmpi_put_att_text(ncid, SFso4_a3, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFso4_a3, "long_name", 19, "so4_a3 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFso4_a3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFso4_a3;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFsoa_a1", NC_FLOAT, 2, dimids, &SFsoa_a1); ERR
    err = ncmpi_put_att_text(ncid, SFsoa_a1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFsoa_a1, "long_name", 19, "soa_a1 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFsoa_a1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFsoa_a1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFsoa_a2", NC_FLOAT, 2, dimids, &SFsoa_a2); ERR
    err = ncmpi_put_att_text(ncid, SFsoa_a2, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFsoa_a2, "long_name", 19, "soa_a2 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFsoa_a2, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFsoa_a2;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SFsoa_a3", NC_FLOAT, 2, dimids, &SFsoa_a3); ERR
    err = ncmpi_put_att_text(ncid, SFsoa_a3, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SFsoa_a3, "long_name", 19, "soa_a3 surface flux"); ERR
    err = ncmpi_put_att_text(ncid, SFsoa_a3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SFsoa_a3;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SHFLX", NC_FLOAT, 2, dimids, &SHFLX); ERR
    err = ncmpi_put_att_text(ncid, SHFLX, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, SHFLX, "long_name", 26, "Surface sensible heat flux"); ERR
    err = ncmpi_put_att_text(ncid, SHFLX, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SHFLX;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SH_KCLDBASE", NC_FLOAT, 2, dimids, &SH_KCLDBASE); ERR
    err = ncmpi_put_att_text(ncid, SH_KCLDBASE, "units", 1, "1"); ERR
    err = ncmpi_put_att_text(ncid, SH_KCLDBASE, "long_name", 35, "Shallow conv. cloudbase level index"); ERR
    err = ncmpi_put_att_text(ncid, SH_KCLDBASE, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SH_KCLDBASE;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SH_MFUP_MAX", NC_FLOAT, 2, dimids, &SH_MFUP_MAX); ERR
    err = ncmpi_put_att_text(ncid, SH_MFUP_MAX, "units", 5, "kg/m2"); ERR
    err = ncmpi_put_att_text(ncid, SH_MFUP_MAX, "long_name", 42, "Shallow conv. column-max updraft mass flux"); ERR
    err = ncmpi_put_att_text(ncid, SH_MFUP_MAX, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SH_MFUP_MAX;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SH_WCLDBASE", NC_FLOAT, 2, dimids, &SH_WCLDBASE); ERR
    err = ncmpi_put_att_text(ncid, SH_WCLDBASE, "units", 3, "m/s"); ERR
    err = ncmpi_put_att_text(ncid, SH_WCLDBASE, "long_name", 41, "Shallow conv. cloudbase vertical velocity"); ERR
    err = ncmpi_put_att_text(ncid, SH_WCLDBASE, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SH_WCLDBASE;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SNOWHICE", NC_FLOAT, 2, dimids, &SNOWHICE); ERR
    err = ncmpi_put_att_text(ncid, SNOWHICE, "units", 1, "m"); ERR
    err = ncmpi_put_att_text(ncid, SNOWHICE, "long_name", 19, "Snow depth over ice"); ERR
    err = ncmpi_put_att_text(ncid, SNOWHICE, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SNOWHICE;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SNOWHLND", NC_FLOAT, 2, dimids, &SNOWHLND); ERR
    err = ncmpi_put_att_text(ncid, SNOWHLND, "units", 1, "m"); ERR
    err = ncmpi_put_att_text(ncid, SNOWHLND, "long_name", 27, "Water equivalent snow depth"); ERR
    err = ncmpi_put_att_text(ncid, SNOWHLND, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SNOWHLND;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "SNOWQM", NC_FLOAT, 3, dimids, &SNOWQM); ERR
    err = ncmpi_put_att_int(ncid, SNOWQM, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, SNOWQM, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, SNOWQM, "long_name", 29, "Grid box averaged snow amount"); ERR
    err = ncmpi_put_att_text(ncid, SNOWQM, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SNOWQM;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "SO2", NC_FLOAT, 3, dimids, &SO2); ERR
    err = ncmpi_put_att_int(ncid, SO2, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, SO2, "units", 7, "mol/mol"); ERR
    err = ncmpi_put_att_text(ncid, SO2, "long_name", 17, "SO2 concentration"); ERR
    err = ncmpi_put_att_text(ncid, SO2, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SO2;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SO2_CLXF", NC_FLOAT, 2, dimids, &SO2_CLXF); ERR
    err = ncmpi_put_att_text(ncid, SO2_CLXF, "units", 11, "molec/cm2/s"); ERR
    err = ncmpi_put_att_text(ncid, SO2_CLXF, "long_name", 47, "vertically intergrated external forcing for SO2"); ERR
    err = ncmpi_put_att_text(ncid, SO2_CLXF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SO2_CLXF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SO2_SRF", NC_FLOAT, 2, dimids, &SO2_SRF); ERR
    err = ncmpi_put_att_text(ncid, SO2_SRF, "units", 7, "mol/mol"); ERR
    err = ncmpi_put_att_text(ncid, SO2_SRF, "long_name", 19, "SO2 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, SO2_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SO2_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SOAG_CLXF", NC_FLOAT, 2, dimids, &SOAG_CLXF); ERR
    err = ncmpi_put_att_text(ncid, SOAG_CLXF, "units", 11, "molec/cm2/s"); ERR
    err = ncmpi_put_att_text(ncid, SOAG_CLXF, "long_name", 48, "vertically intergrated external forcing for SOAG"); ERR
    err = ncmpi_put_att_text(ncid, SOAG_CLXF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SOAG_CLXF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SOAG_SRF", NC_FLOAT, 2, dimids, &SOAG_SRF); ERR
    err = ncmpi_put_att_text(ncid, SOAG_SRF, "units", 7, "mol/mol"); ERR
    err = ncmpi_put_att_text(ncid, SOAG_SRF, "long_name", 20, "SOAG in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, SOAG_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SOAG_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SOAG_sfgaex1", NC_FLOAT, 2, dimids, &SOAG_sfgaex1); ERR
    err = ncmpi_put_att_text(ncid, SOAG_sfgaex1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SOAG_sfgaex1, "long_name", 49, "SOAG gas-aerosol-exchange primary column tendency"); ERR
    err = ncmpi_put_att_text(ncid, SOAG_sfgaex1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SOAG_sfgaex1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SOLIN", NC_FLOAT, 2, dimids, &SOLIN); ERR
    err = ncmpi_put_att_text(ncid, SOLIN, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, SOLIN, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, SOLIN, "long_name", 16, "Solar insolation"); ERR
    err = ncmpi_put_att_text(ncid, SOLIN, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SOLIN;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SSAVIS", NC_FLOAT, 2, dimids, &SSAVIS); ERR
    err = ncmpi_put_att_float(ncid, SSAVIS, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, SSAVIS, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, SSAVIS, "long_name", 29, "Aerosol singel-scatter albedo"); ERR
    err = ncmpi_put_att_text(ncid, SSAVIS, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SSAVIS;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SSTSFMBL", NC_FLOAT, 2, dimids, &SSTSFMBL); ERR
    err = ncmpi_put_att_text(ncid, SSTSFMBL, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SSTSFMBL, "long_name", 28, "Mobilization flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, SSTSFMBL, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SSTSFMBL;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SSTSFMBL_OM", NC_FLOAT, 2, dimids, &SSTSFMBL_OM); ERR
    err = ncmpi_put_att_text(ncid, SSTSFMBL_OM, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, SSTSFMBL_OM, "long_name", 53, "Mobilization flux of marine organic matter at surface"); ERR
    err = ncmpi_put_att_text(ncid, SSTSFMBL_OM, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SSTSFMBL_OM;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "SWCF", NC_FLOAT, 2, dimids, &SWCF); ERR
    err = ncmpi_put_att_text(ncid, SWCF, "Sampling_Sequence", 8, "rad_lwsw"); ERR
    err = ncmpi_put_att_text(ncid, SWCF, "units", 4, "W/m2"); ERR
    err = ncmpi_put_att_text(ncid, SWCF, "long_name", 23, "Shortwave cloud forcing"); ERR
    err = ncmpi_put_att_text(ncid, SWCF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = SWCF;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "T", NC_FLOAT, 3, dimids, &T); ERR
    err = ncmpi_put_att_int(ncid, T, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, T, "units", 1, "K"); ERR
    err = ncmpi_put_att_text(ncid, T, "long_name", 11, "Temperature"); ERR
    err = ncmpi_put_att_text(ncid, T, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = T;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TAUGWX", NC_FLOAT, 2, dimids, &TAUGWX); ERR
    err = ncmpi_put_att_text(ncid, TAUGWX, "units", 4, "N/m2"); ERR
    err = ncmpi_put_att_text(ncid, TAUGWX, "long_name", 33, "Zonal gravity wave surface stress"); ERR
    err = ncmpi_put_att_text(ncid, TAUGWX, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = TAUGWX;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TAUGWY", NC_FLOAT, 2, dimids, &TAUGWY); ERR
    err = ncmpi_put_att_text(ncid, TAUGWY, "units", 4, "N/m2"); ERR
    err = ncmpi_put_att_text(ncid, TAUGWY, "long_name", 38, "Meridional gravity wave surface stress"); ERR
    err = ncmpi_put_att_text(ncid, TAUGWY, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = TAUGWY;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TAUX", NC_FLOAT, 2, dimids, &TAUX); ERR
    err = ncmpi_put_att_text(ncid, TAUX, "units", 4, "N/m2"); ERR
    err = ncmpi_put_att_text(ncid, TAUX, "long_name", 20, "Zonal surface stress"); ERR
    err = ncmpi_put_att_text(ncid, TAUX, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = TAUX;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TAUY", NC_FLOAT, 2, dimids, &TAUY); ERR
    err = ncmpi_put_att_text(ncid, TAUY, "units", 4, "N/m2"); ERR
    err = ncmpi_put_att_text(ncid, TAUY, "long_name", 25, "Meridional surface stress"); ERR
    err = ncmpi_put_att_text(ncid, TAUY, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = TAUY;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TGCLDCWP", NC_FLOAT, 2, dimids, &TGCLDCWP); ERR
    err = ncmpi_put_att_text(ncid, TGCLDCWP, "units", 5, "kg/m2"); ERR
    err = ncmpi_put_att_text(ncid, TGCLDCWP, "long_name", 48, "Total grid-box cloud water path (liquid and ice)"); ERR
    err = ncmpi_put_att_text(ncid, TGCLDCWP, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = TGCLDCWP;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TGCLDIWP", NC_FLOAT, 2, dimids, &TGCLDIWP); ERR
    err = ncmpi_put_att_text(ncid, TGCLDIWP, "units", 5, "kg/m2"); ERR
    err = ncmpi_put_att_text(ncid, TGCLDIWP, "long_name", 35, "Total grid-box cloud ice water path"); ERR
    err = ncmpi_put_att_text(ncid, TGCLDIWP, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = TGCLDIWP;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TGCLDLWP", NC_FLOAT, 2, dimids, &TGCLDLWP); ERR
    err = ncmpi_put_att_text(ncid, TGCLDLWP, "units", 5, "kg/m2"); ERR
    err = ncmpi_put_att_text(ncid, TGCLDLWP, "long_name", 38, "Total grid-box cloud liquid water path"); ERR
    err = ncmpi_put_att_text(ncid, TGCLDLWP, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = TGCLDLWP;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TH7001000", NC_FLOAT, 2, dimids, &TH7001000); ERR
    err = ncmpi_put_att_text(ncid, TH7001000, "units", 1, "K"); ERR
    err = ncmpi_put_att_text(ncid, TH7001000, "long_name", 33, "Theta difference 700 mb - 1000 mb"); ERR
    err = ncmpi_put_att_text(ncid, TH7001000, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = TH7001000;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TMQ", NC_FLOAT, 2, dimids, &TMQ); ERR
    err = ncmpi_put_att_text(ncid, TMQ, "units", 5, "kg/m2"); ERR
    err = ncmpi_put_att_text(ncid, TMQ, "long_name", 48, "Total (vertically integrated) precipitable water"); ERR
    err = ncmpi_put_att_text(ncid, TMQ, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = TMQ;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TREFHT", NC_FLOAT, 2, dimids, &TREFHT); ERR
    err = ncmpi_put_att_text(ncid, TREFHT, "units", 1, "K"); ERR
    err = ncmpi_put_att_text(ncid, TREFHT, "long_name", 28, "Reference height temperature"); ERR
    err = ncmpi_put_att_text(ncid, TREFHT, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = TREFHT;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TROP_P", NC_FLOAT, 2, dimids, &TROP_P); ERR
    err = ncmpi_put_att_float(ncid, TROP_P, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, TROP_P, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, TROP_P, "units", 2, "Pa"); ERR
    err = ncmpi_put_att_text(ncid, TROP_P, "long_name", 19, "Tropopause Pressure"); ERR
    err = ncmpi_put_att_text(ncid, TROP_P, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = TROP_P;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TROP_T", NC_FLOAT, 2, dimids, &TROP_T); ERR
    err = ncmpi_put_att_float(ncid, TROP_T, _FillValue, NC_FLOAT, 1, &fillv); ERR
    err = ncmpi_put_att_float(ncid, TROP_T, "missing_value", NC_FLOAT, 1, &missv); ERR
    err = ncmpi_put_att_text(ncid, TROP_T, "units", 1, "K"); ERR
    err = ncmpi_put_att_text(ncid, TROP_T, "long_name", 22, "Tropopause Temperature"); ERR
    err = ncmpi_put_att_text(ncid, TROP_T, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = TROP_T;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TS", NC_FLOAT, 2, dimids, &TS); ERR
    err = ncmpi_put_att_text(ncid, TS, "units", 1, "K"); ERR
    err = ncmpi_put_att_text(ncid, TS, "long_name", 31, "Surface temperature (radiative)"); ERR
    err = ncmpi_put_att_text(ncid, TS, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = TS;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TSMN", NC_FLOAT, 2, dimids, &TSMN); ERR
    err = ncmpi_put_att_text(ncid, TSMN, "units", 1, "K"); ERR
    err = ncmpi_put_att_text(ncid, TSMN, "long_name", 46, "Minimum surface temperature over output period"); ERR
    err = ncmpi_put_att_text(ncid, TSMN, "cell_methods", 13, "time: minimum"); ERR
    varids[i++] = TSMN;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TSMX", NC_FLOAT, 2, dimids, &TSMX); ERR
    err = ncmpi_put_att_text(ncid, TSMX, "units", 1, "K"); ERR
    err = ncmpi_put_att_text(ncid, TSMX, "long_name", 46, "Maximum surface temperature over output period"); ERR
    err = ncmpi_put_att_text(ncid, TSMX, "cell_methods", 13, "time: maximum"); ERR
    varids[i++] = TSMX;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TUH", NC_FLOAT, 2, dimids, &TUH); ERR
    err = ncmpi_put_att_text(ncid, TUH, "units", 3, "W/m"); ERR
    err = ncmpi_put_att_text(ncid, TUH, "long_name", 44, "Total (vertically integrated) zonal MSE flux"); ERR
    err = ncmpi_put_att_text(ncid, TUH, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = TUH;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TUQ", NC_FLOAT, 2, dimids, &TUQ); ERR
    err = ncmpi_put_att_text(ncid, TUQ, "units", 6, "kg/m/s"); ERR
    err = ncmpi_put_att_text(ncid, TUQ, "long_name", 46, "Total (vertically integrated) zonal water flux"); ERR
    err = ncmpi_put_att_text(ncid, TUQ, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = TUQ;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TVH", NC_FLOAT, 2, dimids, &TVH); ERR
    err = ncmpi_put_att_text(ncid, TVH, "units", 3, "W/m"); ERR
    err = ncmpi_put_att_text(ncid, TVH, "long_name", 49, "Total (vertically integrated) meridional MSE flux"); ERR
    err = ncmpi_put_att_text(ncid, TVH, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = TVH;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "TVQ", NC_FLOAT, 2, dimids, &TVQ); ERR
    err = ncmpi_put_att_text(ncid, TVQ, "units", 6, "kg/m/s"); ERR
    err = ncmpi_put_att_text(ncid, TVQ, "long_name", 51, "Total (vertically integrated) meridional water flux"); ERR
    err = ncmpi_put_att_text(ncid, TVQ, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = TVQ;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "U", NC_FLOAT, 3, dimids, &U); ERR
    err = ncmpi_put_att_int(ncid, U, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, U, "units", 3, "m/s"); ERR
    err = ncmpi_put_att_text(ncid, U, "long_name", 10, "Zonal wind"); ERR
    err = ncmpi_put_att_text(ncid, U, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = U;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "U10", NC_FLOAT, 2, dimids, &U10); ERR
    err = ncmpi_put_att_text(ncid, U10, "units", 3, "m/s"); ERR
    err = ncmpi_put_att_text(ncid, U10, "long_name", 14, "10m wind speed"); ERR
    err = ncmpi_put_att_text(ncid, U10, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = U10;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "UU", NC_FLOAT, 3, dimids, &UU); ERR
    err = ncmpi_put_att_int(ncid, UU, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, UU, "units", 5, "m2/s2"); ERR
    err = ncmpi_put_att_text(ncid, UU, "long_name", 22, "Zonal velocity squared"); ERR
    err = ncmpi_put_att_text(ncid, UU, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = UU;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "V", NC_FLOAT, 3, dimids, &V); ERR
    err = ncmpi_put_att_int(ncid, V, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, V, "units", 3, "m/s"); ERR
    err = ncmpi_put_att_text(ncid, V, "long_name", 15, "Meridional wind"); ERR
    err = ncmpi_put_att_text(ncid, V, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = V;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "VQ", NC_FLOAT, 3, dimids, &VQ); ERR
    err = ncmpi_put_att_int(ncid, VQ, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, VQ, "units", 8, "m/skg/kg"); ERR
    err = ncmpi_put_att_text(ncid, VQ, "long_name", 26, "Meridional water transport"); ERR
    err = ncmpi_put_att_text(ncid, VQ, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = VQ;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "VT", NC_FLOAT, 3, dimids, &VT); ERR
    err = ncmpi_put_att_int(ncid, VT, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, VT, "units", 5, "K m/s"); ERR
    err = ncmpi_put_att_text(ncid, VT, "long_name", 25, "Meridional heat transport"); ERR
    err = ncmpi_put_att_text(ncid, VT, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = VT;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "VU", NC_FLOAT, 3, dimids, &VU); ERR
    err = ncmpi_put_att_int(ncid, VU, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, VU, "units", 5, "m2/s2"); ERR
    err = ncmpi_put_att_text(ncid, VU, "long_name", 33, "Meridional flux of zonal momentum"); ERR
    err = ncmpi_put_att_text(ncid, VU, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = VU;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "VV", NC_FLOAT, 3, dimids, &VV); ERR
    err = ncmpi_put_att_int(ncid, VV, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, VV, "units", 5, "m2/s2"); ERR
    err = ncmpi_put_att_text(ncid, VV, "long_name", 27, "Meridional velocity squared"); ERR
    err = ncmpi_put_att_text(ncid, VV, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = VV;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "WD_H2O2", NC_FLOAT, 2, dimids, &WD_H2O2); ERR
    err = ncmpi_put_att_text(ncid, WD_H2O2, "units", 4, "kg/s"); ERR
    err = ncmpi_put_att_text(ncid, WD_H2O2, "long_name", 31, "H2O2             wet deposition"); ERR
    err = ncmpi_put_att_text(ncid, WD_H2O2, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = WD_H2O2;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "WD_H2SO4", NC_FLOAT, 2, dimids, &WD_H2SO4); ERR
    err = ncmpi_put_att_text(ncid, WD_H2SO4, "units", 4, "kg/s"); ERR
    err = ncmpi_put_att_text(ncid, WD_H2SO4, "long_name", 31, "H2SO4            wet deposition"); ERR
    err = ncmpi_put_att_text(ncid, WD_H2SO4, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = WD_H2SO4;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "WD_SO2", NC_FLOAT, 2, dimids, &WD_SO2); ERR
    err = ncmpi_put_att_text(ncid, WD_SO2, "units", 4, "kg/s"); ERR
    err = ncmpi_put_att_text(ncid, WD_SO2, "long_name", 31, "SO2              wet deposition"); ERR
    err = ncmpi_put_att_text(ncid, WD_SO2, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = WD_SO2;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "WSUB", NC_FLOAT, 3, dimids, &WSUB); ERR
    err = ncmpi_put_att_int(ncid, WSUB, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, WSUB, "units", 3, "m/s"); ERR
    err = ncmpi_put_att_text(ncid, WSUB, "long_name", 37, "Diagnostic sub-grid vertical velocity"); ERR
    err = ncmpi_put_att_text(ncid, WSUB, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = WSUB;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "Z3", NC_FLOAT, 3, dimids, &Z3); ERR
    err = ncmpi_put_att_int(ncid, Z3, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, Z3, "units", 1, "m"); ERR
    err = ncmpi_put_att_text(ncid, Z3, "long_name", 37, "Geopotential Height (above sea level)"); ERR
    err = ncmpi_put_att_text(ncid, Z3, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = Z3;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "aero_water", NC_FLOAT, 3, dimids, &aero_water); ERR
    err = ncmpi_put_att_int(ncid, aero_water, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, aero_water, "units", 1, "m"); ERR
    err = ncmpi_put_att_text(ncid, aero_water, "long_name", 70, "sum of aerosol water of interstitial modes wat_a1+wat_a2+wat_a3+wat_a4"); ERR
    err = ncmpi_put_att_text(ncid, aero_water, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = aero_water;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "airFV", NC_FLOAT, 2, dimids, &airFV); ERR
    err = ncmpi_put_att_text(ncid, airFV, "units", 4, "frac"); ERR
    err = ncmpi_put_att_text(ncid, airFV, "long_name", 2, "FV"); ERR
    err = ncmpi_put_att_text(ncid, airFV, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = airFV;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_a1DDF", NC_FLOAT, 2, dimids, &bc_a1DDF); ERR
    err = ncmpi_put_att_text(ncid, bc_a1DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, bc_a1DDF, "long_name", 49, "bc_a1 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, bc_a1DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_a1DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_a1SFWET", NC_FLOAT, 2, dimids, &bc_a1SFWET); ERR
    err = ncmpi_put_att_text(ncid, bc_a1SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, bc_a1SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, bc_a1SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_a1SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_a1_SRF", NC_FLOAT, 2, dimids, &bc_a1_SRF); ERR
    err = ncmpi_put_att_text(ncid, bc_a1_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, bc_a1_SRF, "long_name", 21, "bc_a1 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, bc_a1_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_a1_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_a1_sfgaex1", NC_FLOAT, 2, dimids, &bc_a1_sfgaex1); ERR
    err = ncmpi_put_att_text(ncid, bc_a1_sfgaex1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, bc_a1_sfgaex1, "long_name", 51, "bc_a1 gas-aerosol-exchange primary column tendency"); ERR
    err = ncmpi_put_att_text(ncid, bc_a1_sfgaex1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_a1_sfgaex1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_a3DDF", NC_FLOAT, 2, dimids, &bc_a3DDF); ERR
    err = ncmpi_put_att_text(ncid, bc_a3DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, bc_a3DDF, "long_name", 49, "bc_a3 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, bc_a3DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_a3DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_a3SFWET", NC_FLOAT, 2, dimids, &bc_a3SFWET); ERR
    err = ncmpi_put_att_text(ncid, bc_a3SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, bc_a3SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, bc_a3SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_a3SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_a3_SRF", NC_FLOAT, 2, dimids, &bc_a3_SRF); ERR
    err = ncmpi_put_att_text(ncid, bc_a3_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, bc_a3_SRF, "long_name", 21, "bc_a3 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, bc_a3_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_a3_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_a4DDF", NC_FLOAT, 2, dimids, &bc_a4DDF); ERR
    err = ncmpi_put_att_text(ncid, bc_a4DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, bc_a4DDF, "long_name", 49, "bc_a4 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, bc_a4DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_a4DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_a4SFWET", NC_FLOAT, 2, dimids, &bc_a4SFWET); ERR
    err = ncmpi_put_att_text(ncid, bc_a4SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, bc_a4SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, bc_a4SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_a4SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_a4_CLXF", NC_FLOAT, 2, dimids, &bc_a4_CLXF); ERR
    err = ncmpi_put_att_text(ncid, bc_a4_CLXF, "units", 11, "molec/cm2/s"); ERR
    err = ncmpi_put_att_text(ncid, bc_a4_CLXF, "long_name", 49, "vertically intergrated external forcing for bc_a4"); ERR
    err = ncmpi_put_att_text(ncid, bc_a4_CLXF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_a4_CLXF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_a4_SRF", NC_FLOAT, 2, dimids, &bc_a4_SRF); ERR
    err = ncmpi_put_att_text(ncid, bc_a4_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, bc_a4_SRF, "long_name", 21, "bc_a4 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, bc_a4_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_a4_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_a4_sfgaex1", NC_FLOAT, 2, dimids, &bc_a4_sfgaex1); ERR
    err = ncmpi_put_att_text(ncid, bc_a4_sfgaex1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, bc_a4_sfgaex1, "long_name", 50, "bc_a4 gas-aerosol-exchange primary column tendency"); ERR
    err = ncmpi_put_att_text(ncid, bc_a4_sfgaex1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_a4_sfgaex1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_c1DDF", NC_FLOAT, 2, dimids, &bc_c1DDF); ERR
    err = ncmpi_put_att_text(ncid, bc_c1DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, bc_c1DDF, "long_name", 49, "bc_c1 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, bc_c1DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_c1DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_c1SFWET", NC_FLOAT, 2, dimids, &bc_c1SFWET); ERR
    err = ncmpi_put_att_text(ncid, bc_c1SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, bc_c1SFWET, "long_name", 36, "bc_c1 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, bc_c1SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_c1SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_c3DDF", NC_FLOAT, 2, dimids, &bc_c3DDF); ERR
    err = ncmpi_put_att_text(ncid, bc_c3DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, bc_c3DDF, "long_name", 49, "bc_c3 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, bc_c3DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_c3DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_c3SFWET", NC_FLOAT, 2, dimids, &bc_c3SFWET); ERR
    err = ncmpi_put_att_text(ncid, bc_c3SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, bc_c3SFWET, "long_name", 36, "bc_c3 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, bc_c3SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_c3SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_c4DDF", NC_FLOAT, 2, dimids, &bc_c4DDF); ERR
    err = ncmpi_put_att_text(ncid, bc_c4DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, bc_c4DDF, "long_name", 49, "bc_c4 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, bc_c4DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_c4DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "bc_c4SFWET", NC_FLOAT, 2, dimids, &bc_c4SFWET); ERR
    err = ncmpi_put_att_text(ncid, bc_c4SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, bc_c4SFWET, "long_name", 36, "bc_c4 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, bc_c4SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = bc_c4SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "chla", NC_FLOAT, 2, dimids, &chla); ERR
    err = ncmpi_put_att_text(ncid, chla, "units", 6, "mg L-1"); ERR
    err = ncmpi_put_att_text(ncid, chla, "long_name", 22, "ocean input data: chla"); ERR
    err = ncmpi_put_att_text(ncid, chla, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = chla;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "dst_a1DDF", NC_FLOAT, 2, dimids, &dst_a1DDF); ERR
    err = ncmpi_put_att_text(ncid, dst_a1DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, dst_a1DDF, "long_name", 50, "dst_a1 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, dst_a1DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = dst_a1DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "dst_a1SF", NC_FLOAT, 2, dimids, &dst_a1SF); ERR
    err = ncmpi_put_att_text(ncid, dst_a1SF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, dst_a1SF, "long_name", 28, "dst_a1 dust surface emission"); ERR
    err = ncmpi_put_att_text(ncid, dst_a1SF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = dst_a1SF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "dst_a1SFWET", NC_FLOAT, 2, dimids, &dst_a1SFWET); ERR
    err = ncmpi_put_att_text(ncid, dst_a1SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, dst_a1SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, dst_a1SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = dst_a1SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "dst_a1_SRF", NC_FLOAT, 2, dimids, &dst_a1_SRF); ERR
    err = ncmpi_put_att_text(ncid, dst_a1_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, dst_a1_SRF, "long_name", 22, "dst_a1 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, dst_a1_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = dst_a1_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "dst_a3DDF", NC_FLOAT, 2, dimids, &dst_a3DDF); ERR
    err = ncmpi_put_att_text(ncid, dst_a3DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, dst_a3DDF, "long_name", 50, "dst_a3 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, dst_a3DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = dst_a3DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "dst_a3SF", NC_FLOAT, 2, dimids, &dst_a3SF); ERR
    err = ncmpi_put_att_text(ncid, dst_a3SF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, dst_a3SF, "long_name", 28, "dst_a3 dust surface emission"); ERR
    err = ncmpi_put_att_text(ncid, dst_a3SF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = dst_a3SF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "dst_a3SFWET", NC_FLOAT, 2, dimids, &dst_a3SFWET); ERR
    err = ncmpi_put_att_text(ncid, dst_a3SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, dst_a3SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, dst_a3SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = dst_a3SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "dst_a3_SRF", NC_FLOAT, 2, dimids, &dst_a3_SRF); ERR
    err = ncmpi_put_att_text(ncid, dst_a3_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, dst_a3_SRF, "long_name", 22, "dst_a3 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, dst_a3_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = dst_a3_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "dst_c1DDF", NC_FLOAT, 2, dimids, &dst_c1DDF); ERR
    err = ncmpi_put_att_text(ncid, dst_c1DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, dst_c1DDF, "long_name", 50, "dst_c1 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, dst_c1DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = dst_c1DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "dst_c1SFWET", NC_FLOAT, 2, dimids, &dst_c1SFWET); ERR
    err = ncmpi_put_att_text(ncid, dst_c1SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, dst_c1SFWET, "long_name", 37, "dst_c1 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, dst_c1SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = dst_c1SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "dst_c3DDF", NC_FLOAT, 2, dimids, &dst_c3DDF); ERR
    err = ncmpi_put_att_text(ncid, dst_c3DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, dst_c3DDF, "long_name", 50, "dst_c3 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, dst_c3DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = dst_c3DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "dst_c3SFWET", NC_FLOAT, 2, dimids, &dst_c3SFWET); ERR
    err = ncmpi_put_att_text(ncid, dst_c3SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, dst_c3SFWET, "long_name", 37, "dst_c3 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, dst_c3SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = dst_c3SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_lev;
    dimids[2] = dim_ncol;
    err = ncmpi_def_var(ncid, "hstobie_linoz", NC_FLOAT, 3, dimids, &hstobie_linoz); ERR
    err = ncmpi_put_att_int(ncid, hstobie_linoz, "mdims", NC_INT, 1, &mdims); ERR
    err = ncmpi_put_att_text(ncid, hstobie_linoz, "units", 22, "fraction of model time"); ERR
    err = ncmpi_put_att_text(ncid, hstobie_linoz, "long_name", 27, "Lowest possible Linoz level"); ERR
    varids[i++] = hstobie_linoz;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mlip", NC_FLOAT, 2, dimids, &mlip); ERR
    err = ncmpi_put_att_text(ncid, mlip, "units", 4, "uM C"); ERR
    err = ncmpi_put_att_text(ncid, mlip, "long_name", 22, "ocean input data: mlip"); ERR
    err = ncmpi_put_att_text(ncid, mlip, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mlip;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_a1DDF", NC_FLOAT, 2, dimids, &mom_a1DDF); ERR
    err = ncmpi_put_att_text(ncid, mom_a1DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_a1DDF, "long_name", 50, "mom_a1 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, mom_a1DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_a1DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_a1SF", NC_FLOAT, 2, dimids, &mom_a1SF); ERR
    err = ncmpi_put_att_text(ncid, mom_a1SF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_a1SF, "long_name", 31, "mom_a1 seasalt surface emission"); ERR
    err = ncmpi_put_att_text(ncid, mom_a1SF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_a1SF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_a1SFWET", NC_FLOAT, 2, dimids, &mom_a1SFWET); ERR
    err = ncmpi_put_att_text(ncid, mom_a1SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_a1SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, mom_a1SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_a1SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_a1_SRF", NC_FLOAT, 2, dimids, &mom_a1_SRF); ERR
    err = ncmpi_put_att_text(ncid, mom_a1_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, mom_a1_SRF, "long_name", 22, "mom_a1 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, mom_a1_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_a1_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_a1_sfgaex1", NC_FLOAT, 2, dimids, &mom_a1_sfgaex1); ERR
    err = ncmpi_put_att_text(ncid, mom_a1_sfgaex1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_a1_sfgaex1, "long_name", 51, "mom_a1 gas-aerosol-exchange primary column tendency"); ERR
    err = ncmpi_put_att_text(ncid, mom_a1_sfgaex1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_a1_sfgaex1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_a2DDF", NC_FLOAT, 2, dimids, &mom_a2DDF); ERR
    err = ncmpi_put_att_text(ncid, mom_a2DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_a2DDF, "long_name", 50, "mom_a2 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, mom_a2DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_a2DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_a2SF", NC_FLOAT, 2, dimids, &mom_a2SF); ERR
    err = ncmpi_put_att_text(ncid, mom_a2SF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_a2SF, "long_name", 31, "mom_a2 seasalt surface emission"); ERR
    err = ncmpi_put_att_text(ncid, mom_a2SF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_a2SF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_a2SFWET", NC_FLOAT, 2, dimids, &mom_a2SFWET); ERR
    err = ncmpi_put_att_text(ncid, mom_a2SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_a2SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, mom_a2SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_a2SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_a2_SRF", NC_FLOAT, 2, dimids, &mom_a2_SRF); ERR
    err = ncmpi_put_att_text(ncid, mom_a2_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, mom_a2_SRF, "long_name", 22, "mom_a2 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, mom_a2_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_a2_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_a3DDF", NC_FLOAT, 2, dimids, &mom_a3DDF); ERR
    err = ncmpi_put_att_text(ncid, mom_a3DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_a3DDF, "long_name", 50, "mom_a3 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, mom_a3DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_a3DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_a3SFWET", NC_FLOAT, 2, dimids, &mom_a3SFWET); ERR
    err = ncmpi_put_att_text(ncid, mom_a3SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_a3SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, mom_a3SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_a3SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_a3_SRF", NC_FLOAT, 2, dimids, &mom_a3_SRF); ERR
    err = ncmpi_put_att_text(ncid, mom_a3_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, mom_a3_SRF, "long_name", 22, "mom_a3 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, mom_a3_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_a3_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_a4DDF", NC_FLOAT, 2, dimids, &mom_a4DDF); ERR
    err = ncmpi_put_att_text(ncid, mom_a4DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_a4DDF, "long_name", 50, "mom_a4 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, mom_a4DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_a4DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_a4SF", NC_FLOAT, 2, dimids, &mom_a4SF); ERR
    err = ncmpi_put_att_text(ncid, mom_a4SF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_a4SF, "long_name", 31, "mom_a4 seasalt surface emission"); ERR
    err = ncmpi_put_att_text(ncid, mom_a4SF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_a4SF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_a4SFWET", NC_FLOAT, 2, dimids, &mom_a4SFWET); ERR
    err = ncmpi_put_att_text(ncid, mom_a4SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_a4SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, mom_a4SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_a4SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_a4_SRF", NC_FLOAT, 2, dimids, &mom_a4_SRF); ERR
    err = ncmpi_put_att_text(ncid, mom_a4_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, mom_a4_SRF, "long_name", 22, "mom_a4 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, mom_a4_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_a4_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_a4_sfgaex1", NC_FLOAT, 2, dimids, &mom_a4_sfgaex1); ERR
    err = ncmpi_put_att_text(ncid, mom_a4_sfgaex1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_a4_sfgaex1, "long_name", 51, "mom_a4 gas-aerosol-exchange primary column tendency"); ERR
    err = ncmpi_put_att_text(ncid, mom_a4_sfgaex1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_a4_sfgaex1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_c1DDF", NC_FLOAT, 2, dimids, &mom_c1DDF); ERR
    err = ncmpi_put_att_text(ncid, mom_c1DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_c1DDF, "long_name", 50, "mom_c1 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, mom_c1DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_c1DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_c1SFWET", NC_FLOAT, 2, dimids, &mom_c1SFWET); ERR
    err = ncmpi_put_att_text(ncid, mom_c1SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_c1SFWET, "long_name", 37, "mom_c1 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, mom_c1SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_c1SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_c2DDF", NC_FLOAT, 2, dimids, &mom_c2DDF); ERR
    err = ncmpi_put_att_text(ncid, mom_c2DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_c2DDF, "long_name", 50, "mom_c2 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, mom_c2DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_c2DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_c2SFWET", NC_FLOAT, 2, dimids, &mom_c2SFWET); ERR
    err = ncmpi_put_att_text(ncid, mom_c2SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_c2SFWET, "long_name", 37, "mom_c2 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, mom_c2SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_c2SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_c3DDF", NC_FLOAT, 2, dimids, &mom_c3DDF); ERR
    err = ncmpi_put_att_text(ncid, mom_c3DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_c3DDF, "long_name", 50, "mom_c3 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, mom_c3DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_c3DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_c3SFWET", NC_FLOAT, 2, dimids, &mom_c3SFWET); ERR
    err = ncmpi_put_att_text(ncid, mom_c3SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_c3SFWET, "long_name", 37, "mom_c3 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, mom_c3SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_c3SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_c4DDF", NC_FLOAT, 2, dimids, &mom_c4DDF); ERR
    err = ncmpi_put_att_text(ncid, mom_c4DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_c4DDF, "long_name", 50, "mom_c4 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, mom_c4DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_c4DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mom_c4SFWET", NC_FLOAT, 2, dimids, &mom_c4SFWET); ERR
    err = ncmpi_put_att_text(ncid, mom_c4SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, mom_c4SFWET, "long_name", 37, "mom_c4 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, mom_c4SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mom_c4SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mpoly", NC_FLOAT, 2, dimids, &mpoly); ERR
    err = ncmpi_put_att_text(ncid, mpoly, "units", 4, "uM C"); ERR
    err = ncmpi_put_att_text(ncid, mpoly, "long_name", 23, "ocean input data: mpoly"); ERR
    err = ncmpi_put_att_text(ncid, mpoly, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mpoly;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "mprot", NC_FLOAT, 2, dimids, &mprot); ERR
    err = ncmpi_put_att_text(ncid, mprot, "units", 4, "uM C"); ERR
    err = ncmpi_put_att_text(ncid, mprot, "long_name", 23, "ocean input data: mprot"); ERR
    err = ncmpi_put_att_text(ncid, mprot, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = mprot;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_a1DDF", NC_FLOAT, 2, dimids, &ncl_a1DDF); ERR
    err = ncmpi_put_att_text(ncid, ncl_a1DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a1DDF, "long_name", 50, "ncl_a1 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a1DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_a1DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_a1SF", NC_FLOAT, 2, dimids, &ncl_a1SF); ERR
    err = ncmpi_put_att_text(ncid, ncl_a1SF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a1SF, "long_name", 31, "ncl_a1 seasalt surface emission"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a1SF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_a1SF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_a1SFWET", NC_FLOAT, 2, dimids, &ncl_a1SFWET); ERR
    err = ncmpi_put_att_text(ncid, ncl_a1SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a1SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a1SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_a1SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_a1_SRF", NC_FLOAT, 2, dimids, &ncl_a1_SRF); ERR
    err = ncmpi_put_att_text(ncid, ncl_a1_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a1_SRF, "long_name", 22, "ncl_a1 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a1_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_a1_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_a2DDF", NC_FLOAT, 2, dimids, &ncl_a2DDF); ERR
    err = ncmpi_put_att_text(ncid, ncl_a2DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a2DDF, "long_name", 50, "ncl_a2 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a2DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_a2DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_a2SF", NC_FLOAT, 2, dimids, &ncl_a2SF); ERR
    err = ncmpi_put_att_text(ncid, ncl_a2SF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a2SF, "long_name", 31, "ncl_a2 seasalt surface emission"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a2SF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_a2SF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_a2SFWET", NC_FLOAT, 2, dimids, &ncl_a2SFWET); ERR
    err = ncmpi_put_att_text(ncid, ncl_a2SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a2SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a2SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_a2SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_a2_SRF", NC_FLOAT, 2, dimids, &ncl_a2_SRF); ERR
    err = ncmpi_put_att_text(ncid, ncl_a2_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a2_SRF, "long_name", 22, "ncl_a2 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a2_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_a2_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_a3DDF", NC_FLOAT, 2, dimids, &ncl_a3DDF); ERR
    err = ncmpi_put_att_text(ncid, ncl_a3DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a3DDF, "long_name", 50, "ncl_a3 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a3DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_a3DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_a3SF", NC_FLOAT, 2, dimids, &ncl_a3SF); ERR
    err = ncmpi_put_att_text(ncid, ncl_a3SF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a3SF, "long_name", 31, "ncl_a3 seasalt surface emission"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a3SF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_a3SF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_a3SFWET", NC_FLOAT, 2, dimids, &ncl_a3SFWET); ERR
    err = ncmpi_put_att_text(ncid, ncl_a3SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a3SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a3SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_a3SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_a3_SRF", NC_FLOAT, 2, dimids, &ncl_a3_SRF); ERR
    err = ncmpi_put_att_text(ncid, ncl_a3_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a3_SRF, "long_name", 22, "ncl_a3 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, ncl_a3_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_a3_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_c1DDF", NC_FLOAT, 2, dimids, &ncl_c1DDF); ERR
    err = ncmpi_put_att_text(ncid, ncl_c1DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, ncl_c1DDF, "long_name", 50, "ncl_c1 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, ncl_c1DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_c1DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_c1SFWET", NC_FLOAT, 2, dimids, &ncl_c1SFWET); ERR
    err = ncmpi_put_att_text(ncid, ncl_c1SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, ncl_c1SFWET, "long_name", 37, "ncl_c1 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, ncl_c1SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_c1SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_c2DDF", NC_FLOAT, 2, dimids, &ncl_c2DDF); ERR
    err = ncmpi_put_att_text(ncid, ncl_c2DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, ncl_c2DDF, "long_name", 50, "ncl_c2 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, ncl_c2DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_c2DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_c2SFWET", NC_FLOAT, 2, dimids, &ncl_c2SFWET); ERR
    err = ncmpi_put_att_text(ncid, ncl_c2SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, ncl_c2SFWET, "long_name", 37, "ncl_c2 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, ncl_c2SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_c2SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_c3DDF", NC_FLOAT, 2, dimids, &ncl_c3DDF); ERR
    err = ncmpi_put_att_text(ncid, ncl_c3DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, ncl_c3DDF, "long_name", 50, "ncl_c3 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, ncl_c3DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_c3DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "ncl_c3SFWET", NC_FLOAT, 2, dimids, &ncl_c3SFWET); ERR
    err = ncmpi_put_att_text(ncid, ncl_c3SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, ncl_c3SFWET, "long_name", 37, "ncl_c3 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, ncl_c3SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = ncl_c3SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a1DDF", NC_FLOAT, 2, dimids, &num_a1DDF); ERR
    err = ncmpi_put_att_text(ncid, num_a1DDF, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_a1DDF, "long_name", 50, "num_a1 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, num_a1DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a1DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a1SF", NC_FLOAT, 2, dimids, &num_a1SF); ERR
    err = ncmpi_put_att_text(ncid, num_a1SF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_a1SF, "long_name", 28, "num_a1 dust surface emission"); ERR
    err = ncmpi_put_att_text(ncid, num_a1SF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a1SF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a1SFWET", NC_FLOAT, 2, dimids, &num_a1SFWET); ERR
    err = ncmpi_put_att_text(ncid, num_a1SFWET, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_a1SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, num_a1SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a1SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a1_CLXF", NC_FLOAT, 2, dimids, &num_a1_CLXF); ERR
    err = ncmpi_put_att_text(ncid, num_a1_CLXF, "units", 11, "molec/cm2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_a1_CLXF, "long_name", 50, "vertically intergrated external forcing for num_a1"); ERR
    err = ncmpi_put_att_text(ncid, num_a1_CLXF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a1_CLXF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a1_SRF", NC_FLOAT, 2, dimids, &num_a1_SRF); ERR
    err = ncmpi_put_att_text(ncid, num_a1_SRF, "units", 5, " 1/kg"); ERR
    err = ncmpi_put_att_text(ncid, num_a1_SRF, "long_name", 22, "num_a1 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, num_a1_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a1_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a1_sfgaex1", NC_FLOAT, 2, dimids, &num_a1_sfgaex1); ERR
    err = ncmpi_put_att_text(ncid, num_a1_sfgaex1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_a1_sfgaex1, "long_name", 51, "num_a1 gas-aerosol-exchange primary column tendency"); ERR
    err = ncmpi_put_att_text(ncid, num_a1_sfgaex1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a1_sfgaex1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a2DDF", NC_FLOAT, 2, dimids, &num_a2DDF); ERR
    err = ncmpi_put_att_text(ncid, num_a2DDF, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_a2DDF, "long_name", 50, "num_a2 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, num_a2DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a2DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a2SFWET", NC_FLOAT, 2, dimids, &num_a2SFWET); ERR
    err = ncmpi_put_att_text(ncid, num_a2SFWET, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_a2SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, num_a2SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a2SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a2_CLXF", NC_FLOAT, 2, dimids, &num_a2_CLXF); ERR
    err = ncmpi_put_att_text(ncid, num_a2_CLXF, "units", 11, "molec/cm2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_a2_CLXF, "long_name", 50, "vertically intergrated external forcing for num_a2"); ERR
    err = ncmpi_put_att_text(ncid, num_a2_CLXF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a2_CLXF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a2_SRF", NC_FLOAT, 2, dimids, &num_a2_SRF); ERR
    err = ncmpi_put_att_text(ncid, num_a2_SRF, "units", 5, " 1/kg"); ERR
    err = ncmpi_put_att_text(ncid, num_a2_SRF, "long_name", 22, "num_a2 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, num_a2_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a2_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a3DDF", NC_FLOAT, 2, dimids, &num_a3DDF); ERR
    err = ncmpi_put_att_text(ncid, num_a3DDF, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_a3DDF, "long_name", 50, "num_a3 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, num_a3DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a3DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a3SF", NC_FLOAT, 2, dimids, &num_a3SF); ERR
    err = ncmpi_put_att_text(ncid, num_a3SF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_a3SF, "long_name", 28, "num_a3 dust surface emission"); ERR
    err = ncmpi_put_att_text(ncid, num_a3SF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a3SF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a3SFWET", NC_FLOAT, 2, dimids, &num_a3SFWET); ERR
    err = ncmpi_put_att_text(ncid, num_a3SFWET, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_a3SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, num_a3SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a3SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a3_SRF", NC_FLOAT, 2, dimids, &num_a3_SRF); ERR
    err = ncmpi_put_att_text(ncid, num_a3_SRF, "units", 5, " 1/kg"); ERR
    err = ncmpi_put_att_text(ncid, num_a3_SRF, "long_name", 22, "num_a3 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, num_a3_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a3_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a4DDF", NC_FLOAT, 2, dimids, &num_a4DDF); ERR
    err = ncmpi_put_att_text(ncid, num_a4DDF, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_a4DDF, "long_name", 50, "num_a4 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, num_a4DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a4DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a4SFWET", NC_FLOAT, 2, dimids, &num_a4SFWET); ERR
    err = ncmpi_put_att_text(ncid, num_a4SFWET, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_a4SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, num_a4SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a4SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a4_CLXF", NC_FLOAT, 2, dimids, &num_a4_CLXF); ERR
    err = ncmpi_put_att_text(ncid, num_a4_CLXF, "units", 11, "molec/cm2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_a4_CLXF, "long_name", 50, "vertically intergrated external forcing for num_a4"); ERR
    err = ncmpi_put_att_text(ncid, num_a4_CLXF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a4_CLXF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a4_SRF", NC_FLOAT, 2, dimids, &num_a4_SRF); ERR
    err = ncmpi_put_att_text(ncid, num_a4_SRF, "units", 5, " 1/kg"); ERR
    err = ncmpi_put_att_text(ncid, num_a4_SRF, "long_name", 22, "num_a4 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, num_a4_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a4_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_a4_sfgaex1", NC_FLOAT, 2, dimids, &num_a4_sfgaex1); ERR
    err = ncmpi_put_att_text(ncid, num_a4_sfgaex1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_a4_sfgaex1, "long_name", 51, "num_a4 gas-aerosol-exchange primary column tendency"); ERR
    err = ncmpi_put_att_text(ncid, num_a4_sfgaex1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_a4_sfgaex1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_c1DDF", NC_FLOAT, 2, dimids, &num_c1DDF); ERR
    err = ncmpi_put_att_text(ncid, num_c1DDF, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_c1DDF, "long_name", 50, "num_c1 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, num_c1DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_c1DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_c1SFWET", NC_FLOAT, 2, dimids, &num_c1SFWET); ERR
    err = ncmpi_put_att_text(ncid, num_c1SFWET, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_c1SFWET, "long_name", 37, "num_c1 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, num_c1SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_c1SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_c2DDF", NC_FLOAT, 2, dimids, &num_c2DDF); ERR
    err = ncmpi_put_att_text(ncid, num_c2DDF, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_c2DDF, "long_name", 50, "num_c2 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, num_c2DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_c2DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_c2SFWET", NC_FLOAT, 2, dimids, &num_c2SFWET); ERR
    err = ncmpi_put_att_text(ncid, num_c2SFWET, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_c2SFWET, "long_name", 37, "num_c2 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, num_c2SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_c2SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_c3DDF", NC_FLOAT, 2, dimids, &num_c3DDF); ERR
    err = ncmpi_put_att_text(ncid, num_c3DDF, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_c3DDF, "long_name", 50, "num_c3 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, num_c3DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_c3DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_c3SFWET", NC_FLOAT, 2, dimids, &num_c3SFWET); ERR
    err = ncmpi_put_att_text(ncid, num_c3SFWET, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_c3SFWET, "long_name", 37, "num_c3 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, num_c3SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_c3SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_c4DDF", NC_FLOAT, 2, dimids, &num_c4DDF); ERR
    err = ncmpi_put_att_text(ncid, num_c4DDF, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_c4DDF, "long_name", 50, "num_c4 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, num_c4DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_c4DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "num_c4SFWET", NC_FLOAT, 2, dimids, &num_c4SFWET); ERR
    err = ncmpi_put_att_text(ncid, num_c4SFWET, "units", 7, " 1/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, num_c4SFWET, "long_name", 37, "num_c4 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, num_c4SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = num_c4SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_a1DDF", NC_FLOAT, 2, dimids, &pom_a1DDF); ERR
    err = ncmpi_put_att_text(ncid, pom_a1DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, pom_a1DDF, "long_name", 50, "pom_a1 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, pom_a1DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_a1DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_a1SFWET", NC_FLOAT, 2, dimids, &pom_a1SFWET); ERR
    err = ncmpi_put_att_text(ncid, pom_a1SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, pom_a1SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, pom_a1SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_a1SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_a1_SRF", NC_FLOAT, 2, dimids, &pom_a1_SRF); ERR
    err = ncmpi_put_att_text(ncid, pom_a1_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, pom_a1_SRF, "long_name", 22, "pom_a1 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, pom_a1_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_a1_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_a1_sfgaex1", NC_FLOAT, 2, dimids, &pom_a1_sfgaex1); ERR
    err = ncmpi_put_att_text(ncid, pom_a1_sfgaex1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, pom_a1_sfgaex1, "long_name", 51, "pom_a1 gas-aerosol-exchange primary column tendency"); ERR
    err = ncmpi_put_att_text(ncid, pom_a1_sfgaex1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_a1_sfgaex1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_a3DDF", NC_FLOAT, 2, dimids, &pom_a3DDF); ERR
    err = ncmpi_put_att_text(ncid, pom_a3DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, pom_a3DDF, "long_name", 50, "pom_a3 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, pom_a3DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_a3DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_a3SFWET", NC_FLOAT, 2, dimids, &pom_a3SFWET); ERR
    err = ncmpi_put_att_text(ncid, pom_a3SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, pom_a3SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, pom_a3SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_a3SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_a3_SRF", NC_FLOAT, 2, dimids, &pom_a3_SRF); ERR
    err = ncmpi_put_att_text(ncid, pom_a3_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, pom_a3_SRF, "long_name", 22, "pom_a3 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, pom_a3_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_a3_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_a4DDF", NC_FLOAT, 2, dimids, &pom_a4DDF); ERR
    err = ncmpi_put_att_text(ncid, pom_a4DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, pom_a4DDF, "long_name", 50, "pom_a4 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, pom_a4DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_a4DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_a4SFWET", NC_FLOAT, 2, dimids, &pom_a4SFWET); ERR
    err = ncmpi_put_att_text(ncid, pom_a4SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, pom_a4SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, pom_a4SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_a4SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_a4_CLXF", NC_FLOAT, 2, dimids, &pom_a4_CLXF); ERR
    err = ncmpi_put_att_text(ncid, pom_a4_CLXF, "units", 11, "molec/cm2/s"); ERR
    err = ncmpi_put_att_text(ncid, pom_a4_CLXF, "long_name", 50, "vertically intergrated external forcing for pom_a4"); ERR
    err = ncmpi_put_att_text(ncid, pom_a4_CLXF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_a4_CLXF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_a4_SRF", NC_FLOAT, 2, dimids, &pom_a4_SRF); ERR
    err = ncmpi_put_att_text(ncid, pom_a4_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, pom_a4_SRF, "long_name", 22, "pom_a4 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, pom_a4_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_a4_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_a4_sfgaex1", NC_FLOAT, 2, dimids, &pom_a4_sfgaex1); ERR
    err = ncmpi_put_att_text(ncid, pom_a4_sfgaex1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, pom_a4_sfgaex1, "long_name", 51, "pom_a4 gas-aerosol-exchange primary column tendency"); ERR
    err = ncmpi_put_att_text(ncid, pom_a4_sfgaex1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_a4_sfgaex1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_c1DDF", NC_FLOAT, 2, dimids, &pom_c1DDF); ERR
    err = ncmpi_put_att_text(ncid, pom_c1DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, pom_c1DDF, "long_name", 50, "pom_c1 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, pom_c1DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_c1DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_c1SFWET", NC_FLOAT, 2, dimids, &pom_c1SFWET); ERR
    err = ncmpi_put_att_text(ncid, pom_c1SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, pom_c1SFWET, "long_name", 37, "pom_c1 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, pom_c1SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_c1SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_c3DDF", NC_FLOAT, 2, dimids, &pom_c3DDF); ERR
    err = ncmpi_put_att_text(ncid, pom_c3DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, pom_c3DDF, "long_name", 50, "pom_c3 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, pom_c3DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_c3DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_c3SFWET", NC_FLOAT, 2, dimids, &pom_c3SFWET); ERR
    err = ncmpi_put_att_text(ncid, pom_c3SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, pom_c3SFWET, "long_name", 37, "pom_c3 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, pom_c3SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_c3SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_c4DDF", NC_FLOAT, 2, dimids, &pom_c4DDF); ERR
    err = ncmpi_put_att_text(ncid, pom_c4DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, pom_c4DDF, "long_name", 50, "pom_c4 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, pom_c4DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_c4DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "pom_c4SFWET", NC_FLOAT, 2, dimids, &pom_c4SFWET); ERR
    err = ncmpi_put_att_text(ncid, pom_c4SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, pom_c4SFWET, "long_name", 37, "pom_c4 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, pom_c4SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = pom_c4SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_a1DDF", NC_FLOAT, 2, dimids, &so4_a1DDF); ERR
    err = ncmpi_put_att_text(ncid, so4_a1DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, so4_a1DDF, "long_name", 50, "so4_a1 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, so4_a1DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_a1DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_a1SFWET", NC_FLOAT, 2, dimids, &so4_a1SFWET); ERR
    err = ncmpi_put_att_text(ncid, so4_a1SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, so4_a1SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, so4_a1SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_a1SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_a1_CLXF", NC_FLOAT, 2, dimids, &so4_a1_CLXF); ERR
    err = ncmpi_put_att_text(ncid, so4_a1_CLXF, "units", 11, "molec/cm2/s"); ERR
    err = ncmpi_put_att_text(ncid, so4_a1_CLXF, "long_name", 50, "vertically intergrated external forcing for so4_a1"); ERR
    err = ncmpi_put_att_text(ncid, so4_a1_CLXF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_a1_CLXF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_a1_SRF", NC_FLOAT, 2, dimids, &so4_a1_SRF); ERR
    err = ncmpi_put_att_text(ncid, so4_a1_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, so4_a1_SRF, "long_name", 22, "so4_a1 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, so4_a1_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_a1_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_a1_sfgaex1", NC_FLOAT, 2, dimids, &so4_a1_sfgaex1); ERR
    err = ncmpi_put_att_text(ncid, so4_a1_sfgaex1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, so4_a1_sfgaex1, "long_name", 51, "so4_a1 gas-aerosol-exchange primary column tendency"); ERR
    err = ncmpi_put_att_text(ncid, so4_a1_sfgaex1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_a1_sfgaex1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_a2DDF", NC_FLOAT, 2, dimids, &so4_a2DDF); ERR
    err = ncmpi_put_att_text(ncid, so4_a2DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, so4_a2DDF, "long_name", 50, "so4_a2 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, so4_a2DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_a2DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_a2SFWET", NC_FLOAT, 2, dimids, &so4_a2SFWET); ERR
    err = ncmpi_put_att_text(ncid, so4_a2SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, so4_a2SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, so4_a2SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_a2SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_a2_CLXF", NC_FLOAT, 2, dimids, &so4_a2_CLXF); ERR
    err = ncmpi_put_att_text(ncid, so4_a2_CLXF, "units", 11, "molec/cm2/s"); ERR
    err = ncmpi_put_att_text(ncid, so4_a2_CLXF, "long_name", 50, "vertically intergrated external forcing for so4_a2"); ERR
    err = ncmpi_put_att_text(ncid, so4_a2_CLXF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_a2_CLXF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_a2_SRF", NC_FLOAT, 2, dimids, &so4_a2_SRF); ERR
    err = ncmpi_put_att_text(ncid, so4_a2_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, so4_a2_SRF, "long_name", 22, "so4_a2 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, so4_a2_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_a2_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_a2_sfgaex1", NC_FLOAT, 2, dimids, &so4_a2_sfgaex1); ERR
    err = ncmpi_put_att_text(ncid, so4_a2_sfgaex1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, so4_a2_sfgaex1, "long_name", 51, "so4_a2 gas-aerosol-exchange primary column tendency"); ERR
    err = ncmpi_put_att_text(ncid, so4_a2_sfgaex1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_a2_sfgaex1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_a3DDF", NC_FLOAT, 2, dimids, &so4_a3DDF); ERR
    err = ncmpi_put_att_text(ncid, so4_a3DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, so4_a3DDF, "long_name", 50, "so4_a3 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, so4_a3DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_a3DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_a3SFWET", NC_FLOAT, 2, dimids, &so4_a3SFWET); ERR
    err = ncmpi_put_att_text(ncid, so4_a3SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, so4_a3SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, so4_a3SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_a3SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_a3_SRF", NC_FLOAT, 2, dimids, &so4_a3_SRF); ERR
    err = ncmpi_put_att_text(ncid, so4_a3_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, so4_a3_SRF, "long_name", 22, "so4_a3 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, so4_a3_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_a3_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_a3_sfgaex1", NC_FLOAT, 2, dimids, &so4_a3_sfgaex1); ERR
    err = ncmpi_put_att_text(ncid, so4_a3_sfgaex1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, so4_a3_sfgaex1, "long_name", 51, "so4_a3 gas-aerosol-exchange primary column tendency"); ERR
    err = ncmpi_put_att_text(ncid, so4_a3_sfgaex1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_a3_sfgaex1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_c1DDF", NC_FLOAT, 2, dimids, &so4_c1DDF); ERR
    err = ncmpi_put_att_text(ncid, so4_c1DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, so4_c1DDF, "long_name", 50, "so4_c1 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, so4_c1DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_c1DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_c1SFWET", NC_FLOAT, 2, dimids, &so4_c1SFWET); ERR
    err = ncmpi_put_att_text(ncid, so4_c1SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, so4_c1SFWET, "long_name", 37, "so4_c1 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, so4_c1SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_c1SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_c2DDF", NC_FLOAT, 2, dimids, &so4_c2DDF); ERR
    err = ncmpi_put_att_text(ncid, so4_c2DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, so4_c2DDF, "long_name", 50, "so4_c2 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, so4_c2DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_c2DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_c2SFWET", NC_FLOAT, 2, dimids, &so4_c2SFWET); ERR
    err = ncmpi_put_att_text(ncid, so4_c2SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, so4_c2SFWET, "long_name", 37, "so4_c2 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, so4_c2SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_c2SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_c3DDF", NC_FLOAT, 2, dimids, &so4_c3DDF); ERR
    err = ncmpi_put_att_text(ncid, so4_c3DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, so4_c3DDF, "long_name", 50, "so4_c3 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, so4_c3DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_c3DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "so4_c3SFWET", NC_FLOAT, 2, dimids, &so4_c3SFWET); ERR
    err = ncmpi_put_att_text(ncid, so4_c3SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, so4_c3SFWET, "long_name", 37, "so4_c3 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, so4_c3SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = so4_c3SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_a1DDF", NC_FLOAT, 2, dimids, &soa_a1DDF); ERR
    err = ncmpi_put_att_text(ncid, soa_a1DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, soa_a1DDF, "long_name", 50, "soa_a1 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, soa_a1DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_a1DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_a1SFWET", NC_FLOAT, 2, dimids, &soa_a1SFWET); ERR
    err = ncmpi_put_att_text(ncid, soa_a1SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, soa_a1SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, soa_a1SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_a1SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_a1_SRF", NC_FLOAT, 2, dimids, &soa_a1_SRF); ERR
    err = ncmpi_put_att_text(ncid, soa_a1_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, soa_a1_SRF, "long_name", 22, "soa_a1 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, soa_a1_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_a1_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_a1_sfgaex1", NC_FLOAT, 2, dimids, &soa_a1_sfgaex1); ERR
    err = ncmpi_put_att_text(ncid, soa_a1_sfgaex1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, soa_a1_sfgaex1, "long_name", 51, "soa_a1 gas-aerosol-exchange primary column tendency"); ERR
    err = ncmpi_put_att_text(ncid, soa_a1_sfgaex1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_a1_sfgaex1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_a2DDF", NC_FLOAT, 2, dimids, &soa_a2DDF); ERR
    err = ncmpi_put_att_text(ncid, soa_a2DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, soa_a2DDF, "long_name", 50, "soa_a2 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, soa_a2DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_a2DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_a2SFWET", NC_FLOAT, 2, dimids, &soa_a2SFWET); ERR
    err = ncmpi_put_att_text(ncid, soa_a2SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, soa_a2SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, soa_a2SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_a2SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_a2_SRF", NC_FLOAT, 2, dimids, &soa_a2_SRF); ERR
    err = ncmpi_put_att_text(ncid, soa_a2_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, soa_a2_SRF, "long_name", 22, "soa_a2 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, soa_a2_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_a2_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_a2_sfgaex1", NC_FLOAT, 2, dimids, &soa_a2_sfgaex1); ERR
    err = ncmpi_put_att_text(ncid, soa_a2_sfgaex1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, soa_a2_sfgaex1, "long_name", 51, "soa_a2 gas-aerosol-exchange primary column tendency"); ERR
    err = ncmpi_put_att_text(ncid, soa_a2_sfgaex1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_a2_sfgaex1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_a3DDF", NC_FLOAT, 2, dimids, &soa_a3DDF); ERR
    err = ncmpi_put_att_text(ncid, soa_a3DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, soa_a3DDF, "long_name", 50, "soa_a3 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, soa_a3DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_a3DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_a3SFWET", NC_FLOAT, 2, dimids, &soa_a3SFWET); ERR
    err = ncmpi_put_att_text(ncid, soa_a3SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, soa_a3SFWET, "long_name", 30, "Wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, soa_a3SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_a3SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_a3_SRF", NC_FLOAT, 2, dimids, &soa_a3_SRF); ERR
    err = ncmpi_put_att_text(ncid, soa_a3_SRF, "units", 5, "kg/kg"); ERR
    err = ncmpi_put_att_text(ncid, soa_a3_SRF, "long_name", 22, "soa_a3 in bottom layer"); ERR
    err = ncmpi_put_att_text(ncid, soa_a3_SRF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_a3_SRF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_a3_sfgaex1", NC_FLOAT, 2, dimids, &soa_a3_sfgaex1); ERR
    err = ncmpi_put_att_text(ncid, soa_a3_sfgaex1, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, soa_a3_sfgaex1, "long_name", 51, "soa_a3 gas-aerosol-exchange primary column tendency"); ERR
    err = ncmpi_put_att_text(ncid, soa_a3_sfgaex1, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_a3_sfgaex1;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_c1DDF", NC_FLOAT, 2, dimids, &soa_c1DDF); ERR
    err = ncmpi_put_att_text(ncid, soa_c1DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, soa_c1DDF, "long_name", 50, "soa_c1 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, soa_c1DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_c1DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_c1SFWET", NC_FLOAT, 2, dimids, &soa_c1SFWET); ERR
    err = ncmpi_put_att_text(ncid, soa_c1SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, soa_c1SFWET, "long_name", 37, "soa_c1 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, soa_c1SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_c1SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_c2DDF", NC_FLOAT, 2, dimids, &soa_c2DDF); ERR
    err = ncmpi_put_att_text(ncid, soa_c2DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, soa_c2DDF, "long_name", 50, "soa_c2 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, soa_c2DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_c2DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_c2SFWET", NC_FLOAT, 2, dimids, &soa_c2SFWET); ERR
    err = ncmpi_put_att_text(ncid, soa_c2SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, soa_c2SFWET, "long_name", 37, "soa_c2 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, soa_c2SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_c2SFWET;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_c3DDF", NC_FLOAT, 2, dimids, &soa_c3DDF); ERR
    err = ncmpi_put_att_text(ncid, soa_c3DDF, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, soa_c3DDF, "long_name", 50, "soa_c3 dry deposition flux at bottom (grav + turb)"); ERR
    err = ncmpi_put_att_text(ncid, soa_c3DDF, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_c3DDF;

    dimids[0] = dim_time;
    dimids[1] = dim_ncol;
    err = ncmpi_def_var(ncid, "soa_c3SFWET", NC_FLOAT, 2, dimids, &soa_c3SFWET); ERR
    err = ncmpi_put_att_text(ncid, soa_c3SFWET, "units", 7, "kg/m2/s"); ERR
    err = ncmpi_put_att_text(ncid, soa_c3SFWET, "long_name", 37, "soa_c3 wet deposition flux at surface"); ERR
    err = ncmpi_put_att_text(ncid, soa_c3SFWET, "cell_methods", 10, "time: mean"); ERR
    varids[i++] = soa_c3SFWET;

    assert(i == nvars);

fn_exit:
    return nerrs;
}
