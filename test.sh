#!/bin/bash
#
# Copyright (C) 2021, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

DECOMP_REPLAY="utils/decomp_copy"
BPSTAT="utils/bpstat"
DECOMPS=()

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
   export LD_LIBRARY_PATH=${PNETCDF_LIB_PATH}:${LD_LIBRARY_PATH}
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
   DECOMPS+=("hdf5 h5")
fi

if test "x${ENABLE_ADIOS2}" = x1 ; then
   APIS+=("adios blob")
   export LD_LIBRARY_PATH=${ADIOS2_LIB_PATH}:${LD_LIBRARY_PATH}
   DECOMPS+=("bp bp")
fi

if test "x${ENABLE_NETCDF4}" = x1 ; then
   APIS+=("netcdf4 canonical")
   export LD_LIBRARY_PATH=${HDF5_LIB_PATH}:${NETCDF4_LIB_PATH}:${LD_LIBRARY_PATH}
   if test "x${ENABLE_LOGVOL}" = x1 ; then
      APIS+=("netcdf4 log")
   fi
   DECOMPS+=("netcdf4 nc4")
fi

mkdir -p ${TESTOUTDIR}

# Convert decomposition files from CDF5 format to others
# unset HDF5_PLUGIN_PATH
# unset HDF5_VOL_CONNECTOR
# for DECOM in "${DECOMPS[@]}" ; do
#     dc=($DECOM)
#     api=${dc[0]}
#     for CONFIG in "${CONFIGS[@]}" ; do
#         IN_FILE=$srcdir/datasets/${CONFIG}.nc
#         OUT_FILE=datasets/${CONFIG}.${dc[1]}
#         if test -e $OUT_FILE ; then
#            if test $VERBOSE = 1 ; then echo "$OUT_FILE already exist"; fi
#         else
#            if test $VERBOSE = 1 ; then echo "Convert $IN_FILE into $OUT_FILE"; fi
#            CMD="${DECOMP_REPLAY} -a ${api} -i ${IN_FILE} -o ${OUT_FILE}"
#            if test $VERBOSE = 1 ; then echo "CMD=$CMD"; fi
#            $CMD
#            if test "x${api}" = xbp ; then
#               CMD="mv -f ${OUT_FILE}.dir/${CONFIG}.${dc[1]}.0 ${OUT_FILE}"
#               if test $VERBOSE = 1 ; then echo "CMD=$CMD"; fi
#               $CMD
#               CMD="rm -rf ${OUT_FILE}.dir"
#               if test $VERBOSE = 1 ; then echo "CMD=$CMD"; fi
#               $CMD
#            fi
#         fi
#     done
# done
# if test $VERBOSE = 1 ; then echo ""; fi

for API in "${APIS[@]}" ; do
    ap=($API)

    # unset the two VOL environment variables, in case they have been set
    unset HDF5_PLUGIN_PATH
    unset HDF5_VOL_CONNECTOR

    for CONFIG in "${CONFIGS[@]}" ; do
        IN_FILE="${srcdir}/datasets/${CONFIG}"
        OUT_FILE_BASE="${TESTOUTDIR}/${ap[0]}_${ap[1]}_${CONFIG}"
        if test "x${ap[0]}" = xpnetcdf ; then
           FILE_EXT="nc"
           IN_FILE+=".${FILE_EXT}"
        elif test "x${ap[0]}" = xnetcdf4 ; then
           FILE_EXT="nc4"
           if test "x${ap[1]}" = xlog ; then
              # This option requires the two VOL environment variables to be set.
              export HDF5_PLUGIN_PATH="$LOGVOL_LIB_PATH"
              export HDF5_VOL_CONNECTOR="LOG under_vol=0;under_info={}"
              # Decomposition file must be read with native VOL, use nc file
              IN_FILE+=".nc"
           else
              IN_FILE+=".${FILE_EXT}"
           fi
        elif test "x${ap[0]}" = xhdf5 || test "x${ap[0]}" = xhdf5_log ; then
           FILE_EXT="h5"
           IN_FILE+=".${FILE_EXT}"
        elif test "x${ap[0]}" = xadios ; then
           FILE_EXT="bp"
           IN_FILE+=".${FILE_EXT}"
        fi

        # construct real output file names
        OUT_FILE="${OUT_FILE_BASE}.${FILE_EXT}"
        if test $CONFIG = f_case_866x72_16p || test $CONFIG = i_case_f19_g16_16p ; then
           REAL_OUT_FILE="${OUT_FILE_BASE}_h0.${FILE_EXT} ${OUT_FILE_BASE}_h1.${FILE_EXT}"
           if test "x${ap[0]}" = xadios ; then
              REAL_OUT_FILE="${OUT_FILE}_h0.${FILE_EXT} ${OUT_FILE}_h1.${FILE_EXT}"
           fi
        elif test $CONFIG = g_case_cmpaso_16p ; then
           REAL_OUT_FILE="${OUT_FILE_BASE}.${FILE_EXT}"
        fi

        CMD="${RUN} ${EXEC} -k -a ${ap[0]} -r 2 -x ${ap[1]} -y 2 -o ${OUT_FILE} ${IN_FILE}"
        echo "${CMD}"
        ${CMD}

        if test "x${ap[0]}" = xadios ; then
           if test $CONFIG = f_case_866x72_16p || test $CONFIG = i_case_f19_g16_16p ; then
              CMD="${BPSTAT} ${OUT_FILE}_h0.${FILE_EXT}"
              echo "${CMD}"
              ${CMD}
              CMD="${BPSTAT} ${OUT_FILE}_h1.${FILE_EXT}"
              echo "${CMD}"
              ${CMD}
           elif test $CONFIG = g_case_cmpaso_16p ; then
              CMD="${BPSTAT} ${OUT_FILE_BASE}.${FILE_EXT}"
              echo "${CMD}"
              ${CMD}
           fi
        fi

        # for log strategy, check if the output files are log-based VOL files
        if test "x${ap[1]}" = xlog ; then
           unset HDF5_VOL_CONNECTOR
           unset HDF5_PLUGIN_PATH

           for f in $REAL_OUT_FILE
           do
             echo "$H5LDUMP -k $f"
             FILE_KIND=`$H5LDUMP -k $f`
             if test "x${FILE_KIND}" != xHDF5-LogVOL ; then
                echo "Error: Output file $f is not Log VOL, but ${FILE_KIND}"
                exit 1
             else
                echo "Success: Output file $f is ${FILE_KIND}"
                echo ""
             fi
           done
        fi

        # delete the output files/folder
        rm -rf ${REAL_OUT_FILE}*
    done
done

