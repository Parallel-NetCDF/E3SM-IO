#!/bin/bash
#
# Copyright (C) 2021, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

CONFIGS=("f_case_866x72_16p.nc" "g_case_cmpaso_16p.nc" "i_case_f19_g16_16p.nc")

APIS=("pnetcdf canonical" "pnetcdf blob")
if test "x${ENABLE_HDF5}" = x1 ; then
    APIS+=("hdf5 blob")
    if test "x${ENABLE_LOGVOL}" = x1 ; then
        APIS+=("hdf5 log")
    fi
fi
if test "x${ENABLE_ADIOS2}" = x1 ; then
    APIS+=("adios2 blob")
fi

for NP in 1 4 ; do
    for API in "${APIS[@]}" ; do
        pair=($API)
        for CONFIG in "${CONFIGS[@]}" ; do
            OUT_PATH="./test_output_${pair[0]}_${pair[1]}_${CONFIG}"
            ${TESTMPIRUN} -np 4 $1 -k -a ${pair[0]} -f -1 -r 2 -x ${pair[1]} -y 2 -o ${OUT_PATH} ${srcdir}/../datasets/${CONFIG}
        done
    done
done

