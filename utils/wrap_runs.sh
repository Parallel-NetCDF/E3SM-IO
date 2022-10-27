#!/bin/sh
#
# Copyright (C) 2022, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

# RUN="valgrind --quiet --leak-check=full"
RUN=""

unset HDF5_VOL_CONNECTOR
unset HDF5_PLUGIN_PATH

LIST_FILE=${builddir}/list_f_case.txt
echo "==== cat $LIST_FILE ===="
cat $LIST_FILE
echo ""

NC_FILE=$top_srcdir/datasets/map_f_case_16p.nc

ALL_LIBRARY_PATH=`echo $ALL_LIBRARY_PATH | ${SED} 's/ //g'`
export LD_LIBRARY_PATH=${ALL_LIBRARY_PATH}:${LD_LIBRARY_PATH}

if test "$1" = ./dat2nc ; then
   rm -f $1.nc
   CMD="${RUN} $1 -i $LIST_FILE -o $1.nc"
   echo $CMD
   $CMD
   rm -f $1.nc
elif test "$1" = ./dat2decomp ; then
   if test "x${ENABLE_HDF5}" = x1 ; then
      rm -f $1.h5
      CMD="${RUN} $1 -a hdf5 -i $LIST_FILE -o $1.h5"
      echo $CMD
      $CMD
      rm -f $1.h5
   fi
   if test "x${ENABLE_NETCDF4}" = x1 ; then
      rm -f $1.nc4
      CMD="${RUN} $1 -a netcdf4 -i $LIST_FILE -o $1.nc4"
      echo $CMD
      $CMD
      rm -f $1.nc4
   fi
   if test "x${ENABLE_ADIOS2}" = x1 ; then
      rm -f $1.nc4
      CMD="${RUN} $1 -a bp -i $LIST_FILE -o $1.bp"
      echo $CMD
      $CMD
      rm -rf $1.bp
   fi
elif test "$1" = ./decomp_copy ; then
   if test "x${ENABLE_PNC}" = x1 && test "x${ENABLE_HDF5}" = x1 ; then
      rm -f $1.h5
      CMD="${RUN} $1 -a hdf5 -i $NC_FILE -o $1.h5"
      echo $CMD
      $CMD
      rm -f $1.h5
   fi
   if test "x${ENABLE_PNC}" = x1 && test "x${ENABLE_NETCDF4}" = x1 ; then
      rm -f $1.nc4
      CMD="${RUN} $1 -a netcdf4 -i $NC_FILE -o $1.nc4"
      echo $CMD
      $CMD
      rm -f $1.nc4
   fi
   if test "x${ENABLE_PNC}" = x1 && test "x${ENABLE_ADIOS2}" = x1 ; then
      rm -f $1.nc4
      CMD="${RUN} $1 -a bp -i $NC_FILE -o $1.bp"
      echo $CMD
      $CMD
      rm -rf $1.bp
   fi
elif test "$1" = ./datstat ; then
   CMD="${RUN} $1 -d $top_srcdir/datasets/piodecomp16tasks16io02dims_ioid_548.dat"
   echo $CMD
   $CMD
fi


