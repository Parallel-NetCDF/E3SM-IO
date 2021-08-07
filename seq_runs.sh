#!/bin/sh
#
# Copyright (C) 2021, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

F_DECOMP=${srcdir}/datasets/f_case_866x72_16p.nc
G_DECOMP=${srcdir}/datasets/g_case_cmpaso_16p.nc
I_DECOMP=${srcdir}/datasets/i_case_f19_g16_16p.nc
EXEC=src/e3sm_io

export LD_LIBRARY_PATH=${HDF5_LIB_PATH}:${LD_LIBRARY_PATH}

# F case
${EXEC} ${F_DECOMP} -o blob_F_out.nc -x blob -a pnetcdf      -r 25
${EXEC} ${F_DECOMP} -o can_F_out.nc  -x canonical -a pnetcdf -r 25

if test "x${ENABLE_HDF5}" = x1 ; then
   ${EXEC} ${F_DECOMP} -o blob_F_out.h5 -x blob      -a hdf5       -r 25
fi
if test "x${HDF5_HAVE_DWRITE_MULTI}" = x1 ; then
   ${EXEC} ${F_DECOMP} -o can_F_out.h5  -x canonical -a hdf5_md    -r 25
fi
if test "x${ENABLE_LOGVOL}" = x1 ; then
   ${EXEC} ${F_DECOMP} -o log_F_out.h5  -x log       -a hdf5_log   -r 25
fi
if test "x${ENABLE_ADIOS2}" = x1 ; then
   ${EXEC} ${F_DECOMP} -o blob_F_out    -x blob      -a adios -g 2 -r 25
fi

# G case
${EXEC} ${G_DECOMP} -o blob_G_out.nc -x blob -a pnetcdf
${EXEC} ${G_DECOMP} -o can_G_out.nc

if test "x${ENABLE_HDF5}" = x1 ; then
   ${EXEC} ${G_DECOMP} -o blob_G_out.h5 -x blob      -a hdf5
fi
if test "x${HDF5_HAVE_DWRITE_MULTI}" = x1 ; then
   ${EXEC} ${G_DECOMP} -o can_G_out.h5  -x canonical -a hdf5_md
fi
if test "x${ENABLE_LOGVOL}" = x1 ; then
   ${EXEC} ${G_DECOMP} -o log_G_out.h5  -x log       -a hdf5_log
fi
if test "x${ENABLE_ADIOS2}" = x1 ; then
   ${EXEC} ${G_DECOMP} -o blob_G_out    -x blob      -a adios -g 2
fi

# I case
<<<<<<< HEAD
${EXEC} ${I_DECOMP} -o blob_I_out.nc -x blob -a pnetcdf -r 2
${EXEC} ${I_DECOMP} -o can_I_out.nc -r 2

if test "x${ENABLE_HDF5}" = x1 ; then
   ${EXEC} ${I_DECOMP} -o blob_I_out.h5 -x blob      -a hdf5 -r 2
fi
if test "x${HDF5_HAVE_DWRITE_MULTI}" = x1 ; then
   ${EXEC} ${I_DECOMP} -o can_I_out.h5  -x canonical -a hdf5_md -r 2
fi
if test "x${ENABLE_LOGVOL}" = x1 ; then
   ${EXEC} ${I_DECOMP} -o log_I_out.h5  -x log       -a hdf5_log -r 2
fi
if test "x${ENABLE_ADIOS2}" = x1 ; then
   ${EXEC} ${I_DECOMP} -o blob_I_out    -x blob      -a adios -g 2 -r 2
=======
${EXEC} ${I_DECOMP} -o blob_G_out.nc -x blob -a pnetcdf
${EXEC} ${I_DECOMP} -o can_G_out.nc

if test "x${ENABLE_HDF5}" = x1 ; then
   ${EXEC} ${I_DECOMP} -o blob_G_out.h5 -x blob      -a hdf5
fi
if test "x${HDF5_HAVE_DWRITE_MULTI}" = x1 ; then
   ${EXEC} ${I_DECOMP} -o can_G_out.h5  -x canonical -a hdf5_md
fi
if test "x${ENABLE_LOGVOL}" = x1 ; then
   ${EXEC} ${I_DECOMP} -o log_G_out.h5  -x log       -a hdf5_log
fi
if test "x${ENABLE_ADIOS2}" = x1 ; then
   ${EXEC} ${I_DECOMP} -o blob_G_out    -x blob      -a adios -g 2
>>>>>>> Add I case test to make check and make ptest:
fi

