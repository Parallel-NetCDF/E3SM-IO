#!/bin/bash
#
# Copyright (C) 2021, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

EXEC="src/e3sm_io"
if test "x$#" = x0 ; then
    RUN=""
else
    RUN="${TESTMPIRUN} -n $1"
fi

CONFIGS=("f_case_866x72_16p.nc" "g_case_cmpaso_16p.nc" "i_case_f19_g16_16p.nc")

APIS=()
if test "x${ENABLE_PNC}" = x1 ; then
   APIS+=("pnetcdf canonical" "pnetcdf blob")
   export LD_LIBRARY_PATH=${PNC_LIB_PATH}:${LD_LIBRARY_PATH}
fi
if test "x${ENABLE_HDF5}" = x1 ; then
   APIS+=("hdf5 canonical" "hdf5 blob")
   export LD_LIBRARY_PATH=${HDF5_LIB_PATH}:${LD_LIBRARY_PATH}
   if test "x${ENABLE_LOGVOL}" = x1 ; then
      APIS+=("hdf5_log log")
      export LD_LIBRARY_PATH=${LOGVOL_LIB_PATH}:${LD_LIBRARY_PATH}
   fi
   # hdf5_md not yet supported
   # if test "x${HDF5_HAVE_DWRITE_MULTI}" = x1 ; then
   #    APIS+=("hdf5_md canonical")
   # fi
fi
if test "x${ENABLE_ADIOS2}" = x1 ; then
   APIS+=("adios blob")
   export LD_LIBRARY_PATH=${ADIOS2_LIB_PATH}:${LD_LIBRARY_PATH}
fi
if test "x${ENABLE_NETCDF4}" = x1 ; then
   APIS+=("netcdf4 canonical")
   export LD_LIBRARY_PATH=${HDF5_LIB_PATH}:${NETCDF4_LIB_PATH}:${LD_LIBRARY_PATH}
   if test "x${ENABLE_LOGVOL}" = x1 ; then
      APIS+=("netcdf4 log")
      export LD_LIBRARY_PATH=${LOGVOL_LIB_PATH}:${LD_LIBRARY_PATH}
   fi
fi

mkdir -p "test_output"

for API in "${APIS[@]}" ; do
    ap=($API)
    for CONFIG in "${CONFIGS[@]}" ; do
        OUT_PATH="./test_output/${ap[0]}_${ap[1]}_${CONFIG}"
        echo "${RUN} ${EXEC} -k -a ${ap[0]} -f -1 -r 2 -x ${ap[1]} -y 2 -o ${OUT_PATH} ${srcdir}/datasets/${CONFIG}"
        ${RUN} ${EXEC} -k -a ${ap[0]} -f -1 -r 2 -x ${ap[1]} -y 2 -o ${OUT_PATH} ${srcdir}/datasets/${CONFIG}
    done
done

