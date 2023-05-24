#!/bin/bash
#
# Copyright (C) 2021, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

DECOMP_REPLAY="utils/decomp_copy"
# BPSTAT="utils/bpstat"
DECOMPS=()

VERBOSE=0
EXEC="src/e3sm_io"
if test "x$#" = x0 ; then
   RUN=""
   VERBOSE=1
else
   RUN="${TESTMPIRUN} -n $1"
fi

CONFIGS=("map_f_case_16p" "map_g_case_16p" "map_i_case_16p")

APIS=()
if test "x${ENABLE_PNC}" = x1 ; then
   APIS+=("pnetcdf canonical" "pnetcdf blob")
   if test "x${LD_LIBRARY_PATH}" = x ; then
      export LD_LIBRARY_PATH=${PNETCDF_LIB_PATH}
   else
      export LD_LIBRARY_PATH=${PNETCDF_LIB_PATH}:${LD_LIBRARY_PATH}
   fi
fi

if test "x${ENABLE_HDF5}" = x1 ; then
   APIS+=("hdf5 canonical" "hdf5 blob")
   if test "x${LD_LIBRARY_PATH}" = x ; then
      export LD_LIBRARY_PATH=${HDF5_LIB_PATH}
   else
      export LD_LIBRARY_PATH=${HDF5_LIB_PATH}:${LD_LIBRARY_PATH}
   fi
   if test "x${ENABLE_LOGVOL}" = x1 ; then
      APIS+=("hdf5_log log")
      if test "x${LOGVOL_LIB_SHARED}" = x1 ; then
         APIS+=("hdf5 log")
         export LD_LIBRARY_PATH=${LOGVOL_LIB_PATH}:${LD_LIBRARY_PATH}
      fi
   fi
   if test "x${HDF5_HAVE_MULTI_DATASET_API}" = x1 ; then
      APIS+=("hdf5_md canonical")
   fi
   DECOMPS+=("hdf5 h5")
fi

if test "x${ENABLE_ADIOS2}" = x1 ; then
   APIS+=("adios blob")
   if test "x${LD_LIBRARY_PATH}" = x ; then
      export LD_LIBRARY_PATH=${ADIOS2_LIB_PATH}
   else
      export LD_LIBRARY_PATH=${ADIOS2_LIB_PATH}:${LD_LIBRARY_PATH}
   fi
   DECOMPS+=("bp bp")
fi

if test "x${ENABLE_NETCDF4}" = x1 ; then
   APIS+=("netcdf4 canonical")
   if test "x${LD_LIBRARY_PATH}" = x ; then
      export LD_LIBRARY_PATH=${HDF5_LIB_PATH}:${NETCDF4_LIB_PATH}
   else
      export LD_LIBRARY_PATH=${HDF5_LIB_PATH}:${NETCDF4_LIB_PATH}:${LD_LIBRARY_PATH}
   fi
   if test "x${ENABLE_LOGVOL}" = x1 && test "x${LOGVOL_LIB_SHARED}" = x1 ; then
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

    for CONFIG in "${CONFIGS[@]}" ; do
        IN_FILE="${srcdir}/datasets/${CONFIG}"
        OUT_FILE_BASE="${TESTOUTDIR}/${ap[0]}_${ap[1]}_${CONFIG}"
        if test "x${ap[0]}" = xpnetcdf ; then
           FILE_EXT="nc"
           IN_FILE+=".${FILE_EXT}"
        elif test "x${ap[0]}" = xnetcdf4 ; then
           FILE_EXT="nc4"
           saved_HDF5_PLUGIN_PATH=$HDF5_PLUGIN_PATH
           saved_HDF5_VOL_CONNECTOR=$HDF5_VOL_CONNECTOR
           if test "x${ap[1]}" = xlog ; then
              # This option requires the two VOL environment variables to be set.
              export HDF5_PLUGIN_PATH="$LOGVOL_LIB_PATH"
              export HDF5_VOL_CONNECTOR="LOG under_vol=0;under_info={}"
              # Decomposition file must be read with native VOL, use nc file
              IN_FILE+=".nc"
           else
              IN_FILE+=".${FILE_EXT}"
              unset HDF5_PLUGIN_PATH
              unset HDF5_VOL_CONNECTOR
           fi
        elif test "x${ap[0]:0:4}" = xhdf5 ; then # hdf5, hdf5_log, or hdf5_md
           FILE_EXT="h5"
           IN_FILE+=".${FILE_EXT}"
        elif test "x${ap[0]}" = xadios ; then
           FILE_EXT="bp"
           IN_FILE+=".${FILE_EXT}"
        fi

        # construct real output file names
        OUT_FILE="${OUT_FILE_BASE}.${FILE_EXT}"
        if test $CONFIG = map_f_case_16p || test $CONFIG = map_i_case_16p ; then
           REAL_OUT_FILE="${OUT_FILE_BASE}_h0.${FILE_EXT} ${OUT_FILE_BASE}_h1.${FILE_EXT}"
           # if test "x${ap[0]}" = xadios ; then
           #    The line below is for BP3
           #    REAL_OUT_FILE="${OUT_FILE}_h0.${FILE_EXT} ${OUT_FILE}_h1.${FILE_EXT}"
           # fi
        elif test $CONFIG = map_g_case_16p ; then
           REAL_OUT_FILE="${OUT_FILE_BASE}.${FILE_EXT}"
        fi

        CMD="${RUN} ${EXEC} -k -a ${ap[0]} -r 2 -x ${ap[1]} -y 2 -o ${OUT_FILE} ${IN_FILE}"
        echo "CMD = ${CMD}"
        ${CMD}

        # run read operations (currently support pnetcdf, netcdf4 and canonical only)
        # Disable read for netcdf4 as it is extremely slow.
        if test "x${ap[1]}" = xcanonical && test "x${ap[0]}" == xpnetcdf ; then
           CMD="${RUN} ${EXEC} -k -a ${ap[0]} -r 2 -x ${ap[1]} -y 2 -i ${OUT_FILE} ${IN_FILE}"
           echo "CMD = ${CMD}"
           ${CMD}
        fi

        if test "x${ap[0]}" = xnetcdf4 ; then
           export HDF5_PLUGIN_PATH=$saved_HDF5_PLUGIN_PATH
           export HDF5_VOL_CONNECTOR=$saved_HDF5_VOL_CONNECTOR
        fi

        # test replay on blob files
        if test "x${ap[1]}" = xblob ; then
           if test "x${ap[0]}" = xpnetcdf ; then
              PNETCDF_REPLAY="utils/pnetcdf_blob_replay"
              if test $CONFIG = map_f_case_16p || test $CONFIG = map_i_case_16p ; then
                 rm -f ${OUT_FILE_BASE}_h0.${FILE_EXT}.can
                 CMD="${RUN} ${PNETCDF_REPLAY} -i ${OUT_FILE_BASE}_h0.${FILE_EXT} -o ${OUT_FILE_BASE}_h0.${FILE_EXT}.can"
                 echo "CMD = ${CMD}"
                 ${CMD}
                 rm -f ${OUT_FILE_BASE}_h1.${FILE_EXT}.can
                 CMD="${RUN} ${PNETCDF_REPLAY} -i ${OUT_FILE_BASE}_h1.${FILE_EXT} -o ${OUT_FILE_BASE}_h1.${FILE_EXT}.can"
                 echo "CMD = ${CMD}"
                 ${CMD}
              elif test $CONFIG = map_g_case_16p ; then
                 rm -f ${OUT_FILE_BASE}.${FILE_EXT}.can
                 CMD="${RUN} ${PNETCDF_REPLAY} -i ${OUT_FILE_BASE}.${FILE_EXT} -o ${OUT_FILE_BASE}.${FILE_EXT}.can"
                 echo "CMD = ${CMD}"
                 ${CMD}
              fi
           fi
        fi

        # check ADIOS BP files
        if test "x$BPSTAT" != x && test "x${ap[0]}" = xadios ; then
           if test $CONFIG = map_f_case_16p || test $CONFIG = map_i_case_16p ; then
              CMD="${BPSTAT} ${OUT_FILE_BASE}_h0.${FILE_EXT}"
              echo "CMD = ${CMD}"
              ${CMD}
              CMD="${BPSTAT} ${OUT_FILE_BASE}_h1.${FILE_EXT}"
              echo "CMD = ${CMD}"
              ${CMD}
           elif test $CONFIG = map_g_case_16p ; then
              CMD="${BPSTAT} ${OUT_FILE_BASE}.${FILE_EXT}"
              echo "CMD = ${CMD}"
              ${CMD}
           fi
        fi

        # for log strategy, check if the output files are log-based VOL files
        if test "x${ap[1]}" = xlog ; then

           # disable all VOLs to avoid debug message messing with the output of h5ldump
           saved_HDF5_VOL_CONNECTOR=
           if test "x$HDF5_VOL_CONNECTOR" != x ; then
              saved_HDF5_VOL_CONNECTOR=$HDF5_VOL_CONNECTOR
              unset HDF5_VOL_CONNECTOR
           fi

           for f in ${REAL_OUT_FILE} ; do
             echo "CMD = $H5LDUMP -k $f"
             FILE_KIND=`$H5LDUMP -k $f`
             if test "x${FILE_KIND}" != xHDF5-LogVOL ; then
                echo "Error: Output file $f is not Log VOL, but ${FILE_KIND}"
                exit 1
             else
                echo "Success: Output file $f is ${FILE_KIND}"
                # test h5lreplay
                CMD="${H5LREPLAY} -i ${f} -o ${f}_replay.h5"
                echo "CMD = ${CMD}"
                echo "Success: test h5lreplay on file $f"
                echo ""
             fi
           done

           # restore HDF5_VOL_CONNECTOR for next run
           if test "x$saved_HDF5_VOL_CONNECTOR" != x ; then
              export HDF5_VOL_CONNECTOR=$saved_HDF5_VOL_CONNECTOR
           fi
        fi

        # delete the output files/folder
        for f in ${REAL_OUT_FILE} ; do
            CMD="rm -rf $f*"
            echo "CMD = ${CMD}"
            ${CMD}
        done
        echo ""
        echo "================================================================"
    done
done

