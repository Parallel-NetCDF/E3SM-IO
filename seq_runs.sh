#!/bin/sh
#
# Copyright (C) 2003, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

F_DECOMP=${srcdir}/datasets/f_case_866x72_16p.nc
G_DECOMP=${srcdir}/datasets/g_case_cmpaso_16p.nc
EXEC=src/e3sm_io

export LD_LIBRARY_PATH=${HDF5_LIB_PATH}:${LD_LIBRARY_PATH}

# F case
${EXEC} ${F_DECOMP} -o blob_F_out.nc -x blob -a pnetcdf      -r 25
${EXEC} ${F_DECOMP} -o blob_F_out.h5 -x blob -a hdf5         -r 25
${EXEC} ${F_DECOMP} -o blob_F_out    -x blob -a adios -g 2   -r 25
${EXEC} ${F_DECOMP} -o can_F_out.nc  -x canonical -a pnetcdf -r 25

# G case
${EXEC} ${G_DECOMP} -o blob_G_out.nc -x blob -a pnetcdf
${EXEC} ${G_DECOMP} -o blob_G_out.h5 -x blob -a hdf5
${EXEC} ${G_DECOMP} -o blob_G_out    -x blob -a adios -g 2
${EXEC} ${G_DECOMP} -o can_G_out.nc

