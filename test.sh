#!/bin/bash
#
# Copyright (C) 2021, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

DECOMREPLAY="utils/decomreplay"
DECOMS=()

VERBOSE=0
EXEC="src/e3sm_io"
if test "x$#" = x0 ; then
   RUN=""
   VERBOSE=1
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
   DECOMS+=("bp bp")
fi

if test "x${ENABLE_NETCDF4}" = x1 ; then
   if test "x$#" = x0 ; then :; else
      echo "==================================================================="
      echo "Warning: skip NetCDF-4 parallel tests due to a bug in NetCDF-C"
      echo "         version 4.8.1 and prior. See bug issue"
      echo "         https://github.com/Unidata/netcdf-c/issues/2251"
      echo "==================================================================="
   fi
   APIS+=("netcdf4 canonical")
   export LD_LIBRARY_PATH=${HDF5_LIB_PATH}:${NETCDF4_LIB_PATH}:${LD_LIBRARY_PATH}
   if test "x${ENABLE_LOGVOL}" = x1 ; then
      APIS+=("netcdf4 log")
   fi
   DECOMS+=("netcdf4 nc4")
fi

OUT_PATH="./test_output"
mkdir -p ${OUT_PATH}

# Convert decomposition files
unset HDF5_PLUGIN_PATH
unset HDF5_VOL_CONNECTOR
for DECOM in "${DECOMS[@]}" ; do
    dc=($DECOM)
    api=${dc[0]}
    for CONFIG in "${CONFIGS[@]}" ; do
        IN_FILE=$srcdir/datasets/${CONFIG}.nc
        OUT_FILE=datasets/${CONFIG}.${dc[1]}
        if test -e $OUT_FILE ; then
           if test $VERBOSE = 1 ; then echo "$OUT_FILE already exist"; fi
        else
           if test $VERBOSE = 1 ; then echo "Convert $IN_FILE into $OUT_FILE"; fi
           CMD="${DECOMREPLAY} -a ${api} -i ${IN_FILE} -o ${OUT_FILE}"
           if test $VERBOSE = 1 ; then echo "CMD=$CMD"; fi
           $CMD
           if test "x${api}" = xbp ; then
              CMD="mv -f ${OUT_FILE}.dir/${CONFIG}.${dc[1]}.0 ${OUT_FILE}"
              if test $VERBOSE = 1 ; then echo "CMD=$CMD"; fi
              $CMD
              CMD="rm -rf ${OUT_FILE}.dir"
              if test $VERBOSE = 1 ; then echo "CMD=$CMD"; fi
              $CMD
           fi
        fi
    done
done

for API in "${APIS[@]}" ; do
    ap=($API)

    # unset the two VOL environment variables, in case they have been set
    unset HDF5_PLUGIN_PATH
    unset HDF5_VOL_CONNECTOR

    for CONFIG in "${CONFIGS[@]}" ; do
        IN_FILE=datasets/${CONFIG}
        OUT_FILE="${OUT_PATH}/${ap[0]}_${ap[1]}_${CONFIG}"
        if test "x${ap[0]}" = xpnetcdf ; then
           IN_FILE="${srcdir}/${IN_FILE}.nc"
           OUT_FILE+=".nc"
        elif test "x${ap[0]}" = xnetcdf4 ; then
           if test "x$#" = x0 ; then :; else continue; fi
           OUT_FILE+=".nc4"
           if test "x${ap[1]}" = xlog ; then
              # This option requires the two VOL environment variables to be set.
              export HDF5_PLUGIN_PATH="$LOGVOL_LIB_PATH"
              export HDF5_VOL_CONNECTOR="LOG under_vol=0;under_info={}"
              # Decomposition file must be read with native VOL, use nc file
              IN_FILE="${srcdir}/${IN_FILE}.nc"
           else
              IN_FILE+=".nc4"
           fi
        elif test "x${ap[0]}" = xhdf5 || test "x${ap[0]}" = xhdf5_log ; then
           IN_FILE+=".h5"
           OUT_FILE+=".h5"
        elif test "x${ap[0]}" = xadios ; then
           IN_FILE+=".bp"
           OUT_FILE+=".bp"
        fi

        CMD="${RUN} ${EXEC} -k -a ${ap[0]} -f -1 -r 2 -x ${ap[1]} -y 2 -o ${OUT_FILE} ${IN_FILE}"
        if test $VERBOSE = 1 ; then echo "${CMD}"; fi
        ${CMD}
    done
done

