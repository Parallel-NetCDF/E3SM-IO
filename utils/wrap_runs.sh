#!/bin/sh
#
# Copyright (C) 2022, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

unset HDF5_VOL_CONNECTOR
unset HDF5_PLUGIN_PATH

rm -f list_f_case.txt
ls -1 $top_srcdir/datasets/piodecomp16tasks*.dat > list_f_case.txt
echo "==== cat list_f_case.txt ===="
cat list_f_case.txt
echo ""

NC_FILE=$top_srcdir/datasets/f_case_866x72_16p.nc

ALL_LIBRARY_PATH=`echo $ALL_LIBRARY_PATH | ${SED} 's/ //g'`
export LD_LIBRARY_PATH=${ALL_LIBRARY_PATH}:${LD_LIBRARY_PATH}

if test "$1" = ./dat2nc ; then
   rm -f $1.nc
   CMD="$1 -i list_f_case.txt -o $1.nc"
   echo $CMD
   $CMD
   rm -f $1.nc
elif test "$1" = ./dat2decomp ; then
   if test "x${ENABLE_HDF5}" = x1 ; then
      rm -f $1.h5
      CMD="$1 -a hdf5 -i list_f_case.txt -o $1.h5"
      echo $CMD
      $CMD
      rm -f $1.h5
   fi
   if test "x${ENABLE_NETCDF4}" = x1 ; then
      rm -f $1.nc4
      CMD="$1 -a netcdf4 -i list_f_case.txt -o $1.nc4"
      echo $CMD
      $CMD
      rm -f $1.nc4
   fi
   if test "x${ENABLE_ADIOS2}" = x1 ; then
      rm -f $1.nc4
      CMD="$1 -a bp -i list_f_case.txt -o $1.bp"
      echo $CMD
      $CMD
      rm -rf $1.bp.dir
   fi
elif test "$1" = ./decomp_copy ; then
   if test "x${ENABLE_PNC}" = x1 && test "x${ENABLE_HDF5}" = x1 ; then
      rm -f $1.h5
      CMD="$1 -a hdf5 -i $NC_FILE -o $1.h5"
      echo $CMD
      $CMD
      rm -f $1.h5
   fi
   if test "x${ENABLE_PNC}" = x1 && test "x${ENABLE_NETCDF4}" = x1 ; then
      rm -f $1.nc4
      CMD="$1 -a netcdf4 -i $NC_FILE -o $1.nc4"
      echo $CMD
      $CMD
      rm -f $1.nc4
   fi
   if test "x${ENABLE_PNC}" = x1 && test "x${ENABLE_ADIOS2}" = x1 ; then
      rm -f $1.nc4
      CMD="$1 -a bp -i $NC_FILE -o $1.bp"
      echo $CMD
      $CMD
      rm -rf $1.bp.dir
   fi
elif test "$1" = ./datstat ; then
   CMD="$1 -d $top_srcdir/datasets/piodecomp16tasks16io02dims_ioid_548.dat"
   echo $CMD
   $CMD
fi

rm -f list_f_case.txt

