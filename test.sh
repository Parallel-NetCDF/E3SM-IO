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

CONFIGS=("f_case_866x72_16p" "g_case_cmpaso_16p" "i_case_f19_g16_16p")

APIS=()
if test "x${ENABLE_PNC}" = x1 ; then
   APIS+=("pnetcdf canonical" "pnetcdf blob")
   export LD_LIBRARY_PATH=${PNC_LIB_PATH}:${LD_LIBRARY_PATH}
fi

if test "x${ENABLE_HDF5}" = x1 ; then
   APIS+=("hdf5 canonical" "hdf5 blob")
   export LD_LIBRARY_PATH=${HDF5_LIB_PATH}:${LD_LIBRARY_PATH}
   if test "x${ENABLE_LOGVOL}" = x1 ; then
      APIS+=("hdf5 log" "hdf5_log log")
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
   fi
fi

OUT_PATH="./test_output"
mkdir -p ${OUT_PATH}

for API in "${APIS[@]}" ; do
    ap=($API)
    if test "x${ap[0]}" = xnetcdf4 && test "x${ap[1]}" = xlog ; then
       # This option requires the two VOL environment variables to be set.
       export HDF5_PLUGIN_PATH="$LOGVOL_LIB_PATH"
       export HDF5_VOL_CONNECTOR="LOG under_vol=0;under_info={}"
    else
       # unset the two VOL environment variables, in case they have been set
       unset HDF5_PLUGIN_PATH
       unset HDF5_VOL_CONNECTOR
    fi

    for CONFIG in "${CONFIGS[@]}" ; do
        OUT_FILE="${OUT_PATH}/${ap[0]}_${ap[1]}_${CONFIG}"
        if test "x${ap[0]}" = xnetcdf4 ; then
           OUT_FILE+=".nc4"
        elif test "x${ap[0]}" = xpnetcdf ; then
           OUT_FILE+=".nc"
        elif test "x${ap[0]}" = xhdf5 || test "x${ap[0]}" = xhdf5_log ; then
           OUT_FILE+=".h5"
        elif test "x${ap[0]}" = xadios ; then
           OUT_FILE+=".bp"
        fi

        echo "${RUN} ${EXEC} -k -a ${ap[0]} -f -1 -r 2 -x ${ap[1]} -y 2 -o ${OUT_FILE} ${srcdir}/datasets/${CONFIG}.nc"
        ${RUN} ${EXEC} -k -a ${ap[0]} -f -1 -r 2 -x ${ap[1]} -y 2 -o ${OUT_FILE} ${srcdir}/datasets/${CONFIG}.nc
    done
done

