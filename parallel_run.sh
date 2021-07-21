#!/bin/sh
#
# Copyright (C) 2021, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

F_DECOMP=${srcdir}/datasets/f_case_866x72_16p.nc
G_DECOMP=${srcdir}/datasets/g_case_cmpaso_16p.nc
EXEC=src/e3sm_io

export LD_LIBRARY_PATH=${HDF5_LIB_PATH}:${LD_LIBRARY_PATH}

# F case
mpiexec -n 16 ${EXEC} ${F_DECOMP} -o blob_F_out.nc -x blob -a pnetcdf      -r 25
mpiexec -n 16 ${EXEC} ${F_DECOMP} -o can_F_out.nc  -x canonical -a pnetcdf -r 25

if test "x${ENABLE_HDF5}" = x1 ; then
   mpiexec -n 16 ${EXEC} ${F_DECOMP} -o blob_F_out.h5 -x blob      -a hdf5       -r 25
fi
if test "x${HDF5_HAVE_DWRITE_MULTI}" = x1 ; then
   mpiexec -n 16 ${EXEC} ${F_DECOMP} -o can_F_out.h5  -x canonical -a hdf5_md    -r 25
fi
if test "x${ENABLE_LOGVOL}" = x1 ; then
   mpiexec -n 16 ${EXEC} ${F_DECOMP} -o log_F_out.h5  -x log       -a hdf5_log   -r 25
fi
if test "x${ENABLE_ADIOS2}" = x1 ; then
   mpiexec -n 16 ${EXEC} ${F_DECOMP} -o blob_F_out    -x blob      -a adios -g 2 -r 25
fi

# G case
mpiexec -n 16 ${EXEC} ${G_DECOMP} -o blob_G_out.nc -x blob -a pnetcdf
mpiexec -n 16 ${EXEC} ${G_DECOMP} -o can_G_out.nc

if test "x${ENABLE_HDF5}" = x1 ; then
   mpiexec -n 16 ${EXEC} ${G_DECOMP} -o blob_G_out.h5 -x blob      -a hdf5
fi
if test "x${HDF5_HAVE_DWRITE_MULTI}" = x1 ; then
   mpiexec -n 16 ${EXEC} ${G_DECOMP} -o can_G_out.h5  -x canonical -a hdf5_md
fi
if test "x${ENABLE_LOGVOL}" = x1 ; then
   mpiexec -n 16 ${EXEC} ${G_DECOMP} -o log_G_out.h5  -x log       -a hdf5_log
fi
if test "x${ENABLE_ADIOS2}" = x1 ; then
   mpiexec -n 16 ${EXEC} ${G_DECOMP} -o blob_G_out    -x blob      -a adios -g 2
fi

