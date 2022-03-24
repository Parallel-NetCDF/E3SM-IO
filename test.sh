#!/bin/bash
#
# Copyright (C) 2021, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

DECOMREPLAY="utils/decomreplay"
DECOMS=()

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
   DECOMS+=("hdf5 h5")
fi

if test "x${ENABLE_ADIOS2}" = x1 ; then
   APIS+=("adios blob")
   export LD_LIBRARY_PATH=${ADIOS2_LIB_PATH}:${LD_LIBRARY_PATH}
   DECOMS+=("adios bp")
fi

if test "x${ENABLE_NETCDF4}" = x1 ; then
   echo "====================================================================="
   echo "Warning: skip NetCDF-4 tests due to a bug in NetCDF-C 4.8.1 and prior"
   echo "         See https://github.com/Unidata/netcdf-c/issues/2251"
   echo "====================================================================="
#   APIS+=("netcdf4 canonical")
#   export LD_LIBRARY_PATH=${HDF5_LIB_PATH}:${NETCDF4_LIB_PATH}:${LD_LIBRARY_PATH}
#   if test "x${ENABLE_LOGVOL}" = x1 ; then
#      APIS+=("netcdf4 log")
#   fi
   DECOMS+=("netcdf4 nc4")
fi

OUT_PATH="./test_output"
mkdir -p ${OUT_PATH}

# Convert decomposition files
unset HDF5_PLUGIN_PATH
unset HDF5_VOL_CONNECTOR
for DECOM in "${DECOMS[@]}" ; do
    dc=($DECOM)
    DECOM_FORMAT=${dc[0]}
    OUT_FILE_EXT=${dc[1]}
    for CONFIG in "${CONFIGS[@]}" ; do
        if test "x${DECOM_FORMAT}" != xpnetcdf ; then
           ${DECOMREPLAY} -a ${DECOM_FORMAT} -i datasets/${CONFIG}.nc -o datasets/${CONFIG}.${OUT_FILE_EXT}
           if test "x${DECOM_FORMAT}" == xadios ; then
              cp datasets/${CONFIG}.${OUT_FILE_EXT}.dir/${CONFIG}.${OUT_FILE_EXT}.0 datasets/${CONFIG}.${OUT_FILE_EXT}
           fi
        fi
    done
done

for API in "${APIS[@]}" ; do
    ap=($API)

     if test "x${ap[0]}" = xnetcdf4 && test "x${ap[1]}" = xlog ; then
        # This option requires the two VOL environment variables to be set.
        export HDF5_PLUGIN_PATH="$LOGVOL_LIB_PATH"
        export HDF5_VOL_CONNECTOR="LOG under_vol=0;under_info={}"
        # Decomposition file must be read with native VOL, use nc file
        DECOM_EXT="nc"
     else
        # unset the two VOL environment variables, in case they have been set
        unset HDF5_PLUGIN_PATH
        unset HDF5_VOL_CONNECTOR
     fi

    for CONFIG in "${CONFIGS[@]}" ; do
        if test "x${ap[0]}" = xnetcdf4 ; then
           OUT_FILE_EXT="nc4"
           DECOM_EXT="nc4"
        elif test "x${ap[0]}" = xpnetcdf ; then
           OUT_FILE_EXT="nc"
           DECOM_EXT="nc"
        elif test "x${ap[0]}" = xhdf5 || test "x${ap[0]}" = xhdf5_log ; then
           OUT_FILE_EXT="h5"
           DECOM_EXT="h5"
        elif test "x${ap[0]}" = xadios ; then
           OUT_FILE_EXT="bp"
           DECOM_EXT="bp"
        fi
        OUT_FILE="${OUT_PATH}/${ap[0]}_${ap[1]}_${CONFIG}.${OUT_FILE_EXT}"
        
        echo "${RUN} ${EXEC} -k -a ${ap[0]} -f -1 -r 2 -x ${ap[1]} -y 2 -o ${OUT_FILE} ${srcdir}/datasets/${CONFIG}.${DECOM_EXT}"
        ${RUN} ${EXEC} -k -a ${ap[0]} -f -1 -r 2 -x ${ap[1]} -y 2 -o ${OUT_FILE} ${srcdir}/datasets/${CONFIG}.${DECOM_EXT}
    done
done

