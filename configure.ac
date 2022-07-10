dnl
dnl Copyright (C) 2020, Northwestern University and Argonne National Laboratory
dnl See COPYRIGHT notice in top-level directory.
dnl

dnl -*- Mode: shell-script-mode; -*-
dnl Process this file with GNU autoconf(1) to produce a configure script.
dnl

dnl autoconf v2.69 was released in 2012-04-24
AC_PREREQ([2.69])
AC_INIT([e3sm_io],[1.3.0],[],[e3sm_io],[])

AM_EXTRA_RECURSIVE_TARGETS([tests])

AC_CONFIG_HEADERS([config.h])
AH_TOP([/*
 * Copyright (C) 2020, Northwestern University and Argonne National Laboratory
 * See COPYRIGHT notice in top-level directory.
 */
])
AC_CONFIG_SRCDIR([src/e3sm_io.h])

AC_CONFIG_AUX_DIR([./scripts])

m4_ifdef([AC_CONFIG_MACRO_DIRS], [AC_CONFIG_MACRO_DIRS([m4])], [AC_CONFIG_MACRO_DIR([m4])])

AM_INIT_AUTOMAKE([foreign])
AM_SILENT_RULES([yes])
AM_MAINTAINER_MODE([enable])

dnl parse the version numbers to 3 env variables
E3SM_IO_VERSION_MAJOR=`echo ${PACKAGE_VERSION} | cut -d. -f1`
E3SM_IO_VERSION_MINOR=`echo ${PACKAGE_VERSION} | cut -d. -f2`
E3SM_IO_VERSION_SUB=`echo ${PACKAGE_VERSION} | cut -d. -f3`
E3SM_IO_VERSION=${PACKAGE_VERSION}
AC_SUBST(E3SM_IO_VERSION_MAJOR)
AC_SUBST(E3SM_IO_VERSION_MINOR)
AC_SUBST(E3SM_IO_VERSION_SUB)
AC_SUBST(E3SM_IO_VERSION)

dnl Note that command 'date' is not portable across Unix platforms.
dnl But release date matters only to e3sm_io developers who make the releases.
E3SM_IO_RELEASE_DATE="`date '+%B %-d, %Y'`"
AC_SUBST(E3SM_IO_RELEASE_DATE)
E3SM_IO_RELEASE_DATE_FULL="`date '+%Y-%m-%d'`"
AC_SUBST(E3SM_IO_RELEASE_DATE_FULL)
AC_SUBST(PACKAGE_VERSION)

AH_TOP([#ifndef _CONFIG_H
#define _CONFIG_H])
AH_BOTTOM([#endif])

AC_PROG_SED
dnl check sed command option -i and set SED_I (this requires RM defined)
UD_PROG_SED_I
AC_PROG_EGREP

MPI_INSTALL=
AC_ARG_WITH(mpi,
   [AS_HELP_STRING([--with-mpi=DIR],
                   [Provide the installation prefix path for MPI implementation.])],[ dnl this clause is run when --with-mpi or --without-mpi is used
   if test "x${withval}" = xno ; then
      AC_MSG_ERROR([
      -----------------------------------------------------------------------
        e3sm_io is built on top of MPI. Configure option --without-mpi or
        --with-mpi=no should not be used. Abort.
      -----------------------------------------------------------------------])
   elif test "x${withval}" = x ; then
      AC_MSG_ERROR(--with-mpi is set but the value is NULL)
   elif test "x${withval}" != xyes && test ! -d "${withval}" ; then
      # user may use --with-mpi without an argument, which results in withval
      # being "yes". This case is OK and we simply take no action, as e3sm_io
      # requires MPI compilers and will check them.
      AC_MSG_ERROR(Directory '${withval}' specified in --with-mpi does not exist or is not a directory)
   fi
   MPI_INSTALL=${withval}
   ]
)
if test "x$TESTMPIRUN" = x ; then
   if test "x$MPI_INSTALL" = x ; then
      TESTMPIRUN=mpiexec
   else
      TESTMPIRUN=$MPI_INSTALL/bin/mpiexec
   fi
fi

AC_ARG_VAR(MPICC,  [MPI C compiler, @<:@default: CC@:>@])
AC_ARG_VAR(MPICXX, [MPI C++ compiler, @<:@default: CXX@:>@])

ac_user_MPICC=$MPICC
ac_user_MPICXX=$MPICXX
if test "x$MPICC"  = x && test "x$CC"  != x ; then ac_user_MPICC=$CC   ; fi
if test "x$MPICXX" = x && test "x$CXX" != x ; then ac_user_MPICXX=$CXX ; fi

CANDIDATE_MPICC="${MPICC} mpicc mpicc_r"
CANDIDATE_MPICXX="${MPICXX} mpicxx mpic++ mpiCC mpcxx mpc++ mpicxx_r mpiCC_r mpcxx_r mpic++_r mpc++_r"
dnl add GNU MPI compilers
CANDIDATE_MPICC="$CANDIDATE_MPICC mpigcc mpgcc mpigcc_r mpgcc_r"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX mpig++ mpg++ mpig++_r mpg++_r"
dnl add IBM MPI compilers
CANDIDATE_MPICC="$CANDIDATE_MPICC mpcc_r mpcc mpixlc_r mpixlc"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX mpCC_r mpCC mpixlcxx_r mpixlcxx mpixlC_r mpixlC"
dnl add IBM BGL MPI compilers
CANDIDATE_MPICC="$CANDIDATE_MPICC blrts_xlc mpxlc_r mpxlc"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX blrts_xlC mpxlC_r mpxlC mpixlc++ mpxlcxx mpxlc++ mpxlCC mpixlc++_r mpxlcxx_r mpxlc++_r mpxlCC_r"
dnl add Fujitsu MPI compilers
CANDIDATE_MPICC="$CANDIDATE_MPICC mpifccpx"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX mpiFCCpx"
dnl add Cray MPI compiler wrappers
CANDIDATE_MPICC="$CANDIDATE_MPICC cc"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX CC"
dnl add Intel MPI compiler wrappers
CANDIDATE_MPICC="$CANDIDATE_MPICC mpiicc icc"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX mpiicpc mpiicxx mpiic++ mpiiCC icpc"
dnl add PGI MPI compiler wrappers
CANDIDATE_MPICC="$CANDIDATE_MPICC mpipgcc mppgcc"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX mpipgCC mppgCC"

dnl find the full path of MPICC from CANDIDATE_MPICC and MPI_INSTALL
if test "x${ac_user_MPICC}" = x ; then
   dnl if MPICC or CC has not been set by users, then search from
   dnl CANDIDATE_MPICC, and find the full path of MPICC
   UD_MPI_PATH_PROGS([MPICC], [$CANDIDATE_MPICC])
else
   dnl check whether user specified MPICC is valid
   UD_MPI_PATH_PROG([MPICC], [$ac_user_MPICC])
fi

if test "x${MPICC}" = x ; then
   if test "x$ac_user_MPICC" = x ; then
      ERR_MSG="No MPI C compiler can be found"
   else
      ERR_MSG="Specified MPI C compiler \"$ac_user_MPICC\" cannot be found"
   fi
   if test "x$MPI_INSTALL" != x ; then
      ERR_MSG="$ERR_MSG under $MPI_INSTALL"
   fi
   AC_MSG_ERROR([
   -----------------------------------------------------------------------
     $ERR_MSG
     e3sm_io requires a working MPI C compiler. Please specify the
     location of an MPI C compiler, either in the MPICC environment
     variable (not CC variable) or through --with-mpi configure flag.
     Abort.
   -----------------------------------------------------------------------])
fi
CC=${MPICC}
AC_PROG_CC

dnl find the full path of MPICXX from CANDIDATE_MPICXX and MPI_INSTALL
if test "x${ac_user_MPICXX}" = x ; then
   dnl if MPICXX or CXX has not been set by users, then search from
   dnl CANDIDATE_MPICXX, and find the full path of MPICXX
   UD_MPI_PATH_PROGS([MPICXX], [$CANDIDATE_MPICXX])
else
   dnl check whether user specified MPICXX is valid
   UD_MPI_PATH_PROG([MPICXX], [$ac_user_MPICXX])
fi
if test "x${MPICXX}" = x ; then
   if test "x$ac_user_MPICXX" = x ; then
      ERR_MSG="No MPI C++ compiler can be found"
   else
      ERR_MSG="Specified MPI C++ compiler \"$ac_user_MPICXX\" cannot be found"
   fi
   if test "x$MPI_INSTALL" != x ; then
      ERR_MSG="$ERR_MSG under $MPI_INSTALL"
   fi
   AC_MSG_ERROR([
   -----------------------------------------------------------------------
     $ERR_MSG
     e3sm_io requires a working MPI C++ compiler. Please specify the
     location of an MPI C++ compiler, either in the MPICXX environment
     variable (not CXX variable) or through --with-mpi configure flag.
     Abort.
   -----------------------------------------------------------------------])
fi
CXX=${MPICXX}
AC_PROG_CXX

AX_CXX_COMPILE_STDCXX(17, , optional)
AM_CONDITIONAL(HAVE_CXX17, [test "x${HAVE_CXX17}" = x1])

AX_CXX_COMPILE_STDCXX(11, , mandatory)

dnl Set output variable CPP to a command that runs the C preprocessor.
dnl Some C compilers require -E to be used as C preprocessor.
AC_PROG_CPP

dnl check if MPICXX works for basic MPI call: MPI_Comm_rank()
AC_LANG_PUSH(C++)
AC_CHECK_FUNC([MPI_Comm_rank], [],
   dnl maybe -lmpi is needed at link stage
   [AC_SEARCH_LIBS([MPI_Comm_rank], [mpi mpi++ mpich mpichcxx mpi_cxx], [],
                   [AC_MSG_ERROR([
   -----------------------------------------------------------------------
     Invalid MPI compiler specified or detected: "${MPICXX}"
     A working MPI C++ compiler is required. Please specify the location
     of one either in the MPICXX environment variable (not CXX variable)
     or through --with-mpi configure flag. Abort.
   -----------------------------------------------------------------------])
])])

AC_CHECK_FUNC([MPI_File_open], [],
   dnl maybe -lmpi++ is needed at link stage
   [AC_SEARCH_LIBS([MPI_File_open], [mpio], [],
                   [AC_MSG_ERROR([
   -----------------------------------------------------------------------
     The underneath MPI implementation does not support MPI-IO.
     e3sm_io requires MPI-IO support to work properly. Abort.
   -----------------------------------------------------------------------])
])])
AC_LANG_POP(C++)

AC_PROG_RANLIB

AC_HEADER_STAT
AC_HEADER_STDC
AC_CHECK_HEADERS([fts.h])
if test "x$ac_cv_header_fts_h" != xyes ; then
   # Force to use MPICXX to check C++ header
   saved_CC=$CC
   CC=$MPICXX
   AC_CHECK_HEADERS_ONCE([filesystem])
   CC=$saved_CC
fi
AM_CONDITIONAL(HAVE_FILESYSTEM_HPP, [test "x$ac_cv_header_filesystem" = xyes])

AC_CHECK_DECL([access], [], [], [[#include <unistd.h>]])
if test "x$ac_cv_have_decl_access" = xyes ; then
   AC_CHECK_FUNCS([access])
fi
AC_CHECK_DECL([unlink], [], [], [[#include <unistd.h>]])
if test "x$ac_cv_have_decl_unlink" = xyes ; then
   AC_CHECK_FUNCS([unlink])
fi

AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
                    [Enable e3sm_io internal debug mode.
                     @<:@default: disabled@:>@])],
    [debug=${enableval}], [debug=no]
)
AM_CONDITIONAL(E3SM_IO_DEBUG, [test "x$debug" = xyes])

if test "x${debug}" = xyes; then
   dnl add -g flag if not presented
   dnl remove all -O and -fast flags
   dnl add -O0 to all flags
   # check exit status of grep command is more portable than using -q
   str_found=`echo "${CXXFLAGS}" | ${EGREP} -- "-g"`
   if test "x$?" != x0 ; then
      CXXFLAGS="$CXXFLAGS -g"
   fi
   CXXFLAGS=`echo $CXXFLAGS | ${SED} 's/-O. *//g' | ${SED} 's/-fast *//g'`
   CXXFLAGS="$CXXFLAGS -O0"

   str_found=`echo "${CFLAGS}" | ${EGREP} -- "-g"`
   if test "x$?" != x0 ; then
      CFLAGS="$CFLAGS -g"
   fi
   CFLAGS=`echo $CFLAGS | ${SED} 's/-O. *//g' | ${SED} 's/-fast *//g'`
   CFLAGS="$CFLAGS -O0"

   unset str_found
   AC_DEFINE([E3SM_IO_DEBUG], [1], [Enable debug mode])
fi

UD_PROG_M4
M4FLAGS=""

# Check for pnc
have_pnc=no
AC_ARG_WITH([pnetcdf],
   [AS_HELP_STRING([--with-pnetcdf@<:@=INC,LIB | =DIR@:>@],
                   [Provide the PnetCDF installation path(s):
    --with-pnetcdf=INC,LIB for include and lib paths separated by a comma.
    --with-pnetcdf=DIR for the path containing include/ and lib/ subdirectories.
    @<:@default: enabled@:>@
   ])], [
   case $withval in
     *,*)
        pnc_inc="`echo $withval |cut -f1 -d,`"
        pnc_lib="`echo $withval |cut -f2 -d, -s`"
        ;;
     *)
        if test -n "$withval"; then
          pnc_inc="$withval/include"
          pnc_lib="$withval/lib"
        fi
        ;;
   esac
   if test "x$pnc_inc" != x ; then
      if test "x$CPPFLAGS" = x ; then
         CPPFLAGS="-I$pnc_inc"
      elif ! echo "${CPPFLAGS}" | ${EGREP} -q -w -- "-I$pnc_inc" ; then
         # append only if not already appear in CPPFLAGS
         CPPFLAGS="$CPPFLAGS -I$pnc_inc"
      fi
   fi
   if test "x$pnc_lib" != x ; then
      if test "x$LDFLAGS" = x ; then
         LDFLAGS="-L$pnc_lib"
      elif ! echo "${LDFLAGS}" | ${EGREP} -q -w -- "-L$pnc_lib" ; then
         # append only if not already appear in LDFLAGS
         LDFLAGS="$LDFLAGS -L$pnc_lib"
      fi
   fi

   AC_CHECK_HEADER([pnetcdf.h], [have_pnc=yes], [have_pnc=no])
   if test "x$have_pnc" = xno ; then
         AC_MSG_ERROR([
         -----------------------------------------------------------------------
         Missing PnetCDF-header files 'pnetcdf.h' required to build E3SM_IO. Use
         configure command-line option --with-pnetcdf=/path/to/implementation
         to specify the location of PnetCDF installation. Abort.
         -----------------------------------------------------------------------])
   else
      dnl Check if PnetCDF version is 1.10.0 or later
      AC_MSG_CHECKING([whether PnetCDF version is 1.10.0 or later])
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      #include <pnetcdf.h>
      #if (PNETCDF_VERSION_MAJOR*1000000 + PNETCDF_VERSION_MINOR*1000 + PNETCDF_VERSION_PATCH < 1010000)
      #error PnetCDF version is older than 1.10.0
      #endif
      ]])], [pnc_ge_1_10_0=yes], [pnc_ge_1_10_0=no])
      AC_MSG_RESULT([$pnc_ge_1_10_0])

      if test x$pnc_ge_1_10_0 = xno; then
         AC_MSG_ERROR([
         -----------------------------------------------------------------------
         PnetCDF version 1.10.0 and later is required. Abort.
         -----------------------------------------------------------------------])
      fi
      AC_SUBST(PNETCDF_LIB_PATH, [$pnc_lib])
   fi
])
AM_CONDITIONAL(ENABLE_PNC, [test "x$have_pnc" = xyes])

# Check for HDF5
have_hdf5=no
have_multi_dset=no
AC_ARG_WITH([hdf5],
   [AS_HELP_STRING([--with-hdf5@<:@=INC,LIB | =DIR@:>@],
                   [Enable HDF5 feature and provide the HDF5 installation path(s):
    --with-hdf5=INC,LIB for include and lib paths separated by a comma.
    --with-hdf5=DIR for the path containing include/ and lib/ subdirectories.
    @<:@default: disabled@:>@
   ])], [
   case $withval in
     *,*)
        hdf5_inc="`echo $withval |cut -f1 -d,`"
        hdf5_lib="`echo $withval |cut -f2 -d, -s`"
        ;;
     *)
        if test -n "$withval"; then
          hdf5_inc="$withval/include"
          hdf5_lib="$withval/lib"
        fi
        ;;
   esac
   if test "x$hdf5_inc" != x ; then
      if test "x$CPPFLAGS" = x ; then
         CPPFLAGS="-I$hdf5_inc"
      elif ! echo "${CPPFLAGS}" | ${EGREP} -q -w -- "-I$hdf5_inc" ; then
         # append only if not already appear in CPPFLAGS
         CPPFLAGS="$CPPFLAGS -I$hdf5_inc"
      fi
   fi
   if test "x$hdf5_lib" != x ; then
      if test "x$LDFLAGS" = x ; then
         LDFLAGS="-L$hdf5_lib"
      elif ! echo "${LDFLAGS}" | ${EGREP} -q -w -- "-L$hdf5_lib" ; then
         # append only if not already appear in LDFLAGS
         LDFLAGS="$LDFLAGS -L$hdf5_lib"
      fi
   fi

   AC_CHECK_HEADER([hdf5.h], [have_hdf5=yes], [have_hdf5=no])
   if test "x$have_hdf5" = xno ; then
      AC_MSG_ERROR([
      -----------------------------------------------------------------------
      Missing HDF5-header files 'hdf5.h' required to build E3SM_IO. Use
      configure command-line option --with-hdf5=/path/to/implementation
      to specify the location of HDF5 installation. Abort.
      -----------------------------------------------------------------------])
   else
      AC_MSG_CHECKING([HDF5 library version])
      AC_COMPUTE_INT([hdf5_major], [H5_VERS_MAJOR],   [[#include <hdf5.h>]])
      AC_COMPUTE_INT([hdf5_minor], [H5_VERS_MINOR],   [[#include <hdf5.h>]])
      AC_COMPUTE_INT([hdf5_patch], [H5_VERS_RELEASE], [[#include <hdf5.h>]])
      hdf5_version=${hdf5_major}.${hdf5_minor}.${hdf5_patch}
      AC_MSG_RESULT([$hdf5_version])

      AC_MSG_CHECKING([whether HDF5 is configured with parallel I/O enabled])
      AC_COMPUTE_INT([h5_parallel], [H5_HAVE_PARALLEL], [[#include <hdf5.h>]], [h5_parallel=0])
      if test "x$h5_parallel" = x1 ; then
         hdf5_parallel_io=yes
      else
         hdf5_parallel_io=no
      fi
      AC_MSG_RESULT([$hdf5_parallel_io])
      if test "x$hdf5_parallel_io" = xno ; then
         AC_MSG_ERROR([
         -----------------------------------------------------------------------
         E3SM_IO requires parallel I/O feature enabled in HDF5. However, This
         version of HDF5 library, $hdf5_version, is not configured to support
         parallel I/O. Abort.
         -----------------------------------------------------------------------])
      fi
   fi
   if test "x$have_hdf5" = xyes ; then
      # Check whether HDF5 supprots H5S_BLOCK
      AC_CHECK_DECLS([H5S_BLOCK], [], [], [#include <hdf5.h>])
      if test "x$ac_cv_have_decl_H5S_BLOCK" = xyes ; then
         AC_DEFINE([HDF5_HAVE_H5SBLOCK], [1], [Whether HDF5 defines H5S_BLOCK])
      fi
      AC_SEARCH_LIBS([H5Dwrite_multi], [hdf5], [have_multi_dset=yes], [have_multi_dset=no])
      if test "x${have_multi_dset}" = xno ; then
         AC_MSG_WARN([
         ---------------------------------------------------------------
         This version of HDF5 library, $hdf5_version, does not support
         'multi-dataset' APIs. The 'hdf5_md' I/O option is thus disabled.
         ---------------------------------------------------------------])
      fi
      AC_SUBST(HDF5_LIB_PATH, [$hdf5_lib])
   fi
])
AM_CONDITIONAL(ENABLE_HDF5, [test "x$have_hdf5" = xyes])
AM_CONDITIONAL(HDF5_HAVE_DWRITE_MULTI, [test "x$have_multi_dset" = xyes])

# Check for NETCDF4
have_netcdf4=no
AC_ARG_WITH([netcdf4],
   [AS_HELP_STRING([--with-netcdf4@<:@=INC,LIB | =DIR@:>@],
                   [Enable ADIOS-2 feature and provide the ADIOS-2 installation path(s):
    --with-netcdf4=INC,LIB for include and lib paths separated by a comma.
    --with-netcdf4=DIR for the path containing include/ and lib/ subdirectories.
    @<:@default: disabled@:>@
   ])], [
   if test "x${withval}" != xyes && test "x${withval}" != xno ; then
      NETCDF4_INSTALL=${withval}
   fi
   case $withval in
     *,*)
        netcdf4_inc="`echo $withval |cut -f1 -d,`"
        netcdf4_lib="`echo $withval |cut -f2 -d, -s`"
        netcdf4_bindir="$netcdf4_inc/../bin"
        ;;
     *)
        if test -n "$withval"; then
          netcdf4_inc="$withval/include"
          netcdf4_lib="$withval/lib"
          netcdf4_bindir="$withval/bin"
        fi
        ;;
   esac
   if test "x$NETCDF4_INSTALL" != x ; then
      AC_PATH_PROG([netcdf4_config],[nc-config],,[$netcdf4_bindir])
   else
      dnl Check netcdf4-config under $PATH
      AC_PATH_PROG([netcdf4_config],[nc-config])
   fi
   if test "x$netcdf4_config" != x ; then
      netcdf4_cflags=-I`$netcdf4_config --includedir`
      netcdf4_lflags=-L`$netcdf4_config --libdir`
   else
      netcdf4_cflags="-I$netcdf4_inc"
      netcdf4_lflags="-L$netcdf4_lib"
   fi
   if test "x$netcdf4_cflags" != x ; then
      if test "x$CPPFLAGS" = x ; then
         CPPFLAGS="${netcdf4_cflags}"
      elif ! echo "${CPPFLAGS}" | ${EGREP} -q -w -- "${netcdf4_cflags}" ; then
         # append only if not already appear in CPPFLAGS
         CPPFLAGS="$CPPFLAGS ${netcdf4_cflags}"
      fi
   fi
   if test "x$netcdf4_lflags" != x ; then
      if test "x$LDFLAGS" = x ; then
         LDFLAGS="$netcdf4_lflags"
      elif ! echo "${LDFLAGS}" | ${EGREP} -q -w -- "${netcdf4_lflags}" ; then
         # append only if not already appear in LDFLAGS
         LDFLAGS="$LDFLAGS $netcdf4_lflags"
      fi
   fi
   AC_CHECK_HEADERS([netcdf.h netcdf_meta.h], [have_netcdf4=yes], [have_netcdf4=no])
   if test "x$have_netcdf4" = xno ; then
      AC_MSG_ERROR([
      -----------------------------------------------------------------------
      Missing NetCDF-header files 'netcdf.h' or 'netcdf_meta.h' which are
      required to build E3SM_IO. Use configure command-line option
      --with-netcdf4=/path/to/implementation to specify the location of
      NetCDF-4 installation. Abort.
      -----------------------------------------------------------------------])
   fi

   if test "x$have_netcdf4" = xyes ; then
      if test "x$netcdf4_config" != x ; then
         netcdf4_version=`$netcdf4_config --version | cut -d' ' -f2 | cut -d'-' -f1`
      else
         AC_COMPUTE_INT([netcdf4_major], [NC_VERSION_MAJOR], [[#include <netcdf_meta.h>]])
         AC_COMPUTE_INT([netcdf4_minor], [NC_VERSION_MINOR], [[#include <netcdf_meta.h>]])
         AC_COMPUTE_INT([netcdf4_patch], [NC_VERSION_PATCH], [[#include <netcdf_meta.h>]])
         netcdf4_version=${netcdf4_major}.${netcdf4_minor}.${netcdf4_patch}
      fi
      AC_MSG_CHECKING([NetCDF-C version])
      AC_MSG_RESULT([${netcdf4_version}])

      AC_CHECK_DECLS([NC_NODIMSCALE_ATTACH], [], [], [#include <netcdf.h>])
      if test "x${ac_cv_have_decl_NC_NODIMSCALE_ATTACH}" = xno ; then
         AC_MSG_ERROR([
         ------------------------------------------------------------
         NetCDF-C library version 4.9.0 and later is required. The
         supplied NetCDF-C library is $netcdf4_version, which does not
         include the option to disable HDF5 dimension sale feature.
         Without this option, F and I cases will faile with error of
         'NetCDF: HDF error'. Abort
         ------------------------------------------------------------])
      fi

      AC_MSG_CHECKING([if NetCDF-C $netcdf4_version is configured with parallel I/O enabled])
      if test "x$netcdf4_config" != x ; then
         netcdf4_has_par=`$netcdf4_config --has-parallel4`
      else
         AC_CHECK_HEADERS([netcdf_par.h], [netcdf4_has_par=yes], [netcdf4_has_par=no], [[#include <netcdf.h>]])
      fi
      AC_MSG_RESULT($netcdf4_has_par)
      if test "x${netcdf4_has_par}" != xyes; then
         AC_MSG_ERROR([
         ------------------------------------------------------------
         This NetCDF-C library was not built with parallel NetCDF4
         support. Abort.
         ------------------------------------------------------------])
         have_netcdf4 = no
      elif test "x$netcdf4_lib" != x ; then
         AC_SUBST(NETCDF4_LIB_PATH, [$netcdf4_lib])
      fi
   fi
])
AM_CONDITIONAL(ENABLE_NETCDF4, [test "x$have_netcdf4" = xyes])

# Check for logvol
have_logvol=no
AC_ARG_WITH([logvol],
   [AS_HELP_STRING([--with-logvol@<:@=INC,LIB | =DIR@:>@],
                   [Enable HDF5 log-based VOL feature and provide the VOL installation path(s):
    --with-logvol=INC,LIB for include and lib paths separated by a comma.
    --with-logvol=DIR for the path containing include/ and lib/ subdirectories.
    @<:@default: disabled@:>@
   ])], [
   if test "x$have_hdf5" = xno ; then
      AC_MSG_ERROR([
      -----------------------------------------------------------------------
      Enabling the log-based VOL option requires HDF5. Please add command-line
      option '--with-hdf5=DIR' to indicate the HDF5 installation path. Abort.
      -----------------------------------------------------------------------])
   fi
   case $withval in
     *,*)
        logvol_inc="`echo $withval |cut -f1 -d,`"
        logvol_lib="`echo $withval |cut -f2 -d, -s`"
        logvol_bindir="$logvol_inc/../bin"
        ;;
     *)
        if test -n "$withval"; then
          logvol_inc="$withval/include"
          logvol_lib="$withval/lib"
          logvol_bindir="$withval/bin"
        fi
        ;;
   esac
   if test "x$logvol_inc" != x ; then
      if test "x$CPPFLAGS" = x ; then
         CPPFLAGS="-I$logvol_inc"
      elif ! echo "${CPPFLAGS}" | ${EGREP} -q -w -- "-I$logvol_inc" ; then
         # append only if not already appear in CPPFLAGS
         CPPFLAGS="$CPPFLAGS -I$logvol_inc"
      fi
   fi
   if test "x$logvol_lib" != x ; then
      if test "x$LDFLAGS" = x ; then
         LDFLAGS="-L$logvol_lib"
      elif ! echo "${LDFLAGS}" | ${EGREP} -q -w -- "-L$logvol_lib" ; then
         # append only if not already appear in LDFLAGS
         LDFLAGS="$LDFLAGS -L$logvol_lib"
      fi
   fi
   AC_CHECK_HEADER([H5VL_log.h], [have_logvol=yes], [have_logvol=no])
   if test "x$have_logvol" = xno ; then
      AC_MSG_ERROR([
      -----------------------------------------------------------------------
      Missing logvol-header files 'H5VL_log.h' required to build E3SM_IO. Use
      configure command-line option --with-logvol=/path/to/implementation
      to specify the location of logvol installation. Abort.
      -----------------------------------------------------------------------])
   fi

   LIBS_BAK=${LIBS}
   LIBS="$LIBS -lH5VL_log"
   AC_CHECK_FUNC([H5Pset_buffered], [have_h5pset_buffered=yes], [have_h5pset_buffered=no])
   if test "x$have_h5pset_buffered" = xyes ; then
      AC_DEFINE([LOGVOL_HAVE_H5PSET_BUFFERED], [1], [Whether log-based VOL have H5Pset_buffered])
   fi
   LIBS=${LIBS_BAK}

   dnl Check if logvol version is 1.3.0 or later
   AC_MSG_CHECKING([whether logvol version is 1.3.0 or later])
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
   #include <H5VL_log.h>
   #if (H5VL_LOG_VERSION_MAJOR*1000000 + H5VL_LOG_VERSION_MINOR*1000 + H5VL_LOG_VERSION_PATCH < 1003000)
   #error logvol version is older than 1.3.0
   #endif
   ]])], [logvol_ge_1_3_0=yes], [logvol_ge_1_3_0=no])
   AC_MSG_RESULT([$logvol_ge_1_3_0])

   if test x$logvol_ge_1_3_0 = xyes; then
         AC_DEFINE([LOGVOL_GE_1003000], [1], [whether logvol version is 1.3.0 or later])
   fi

   if test "x$have_logvol" = xyes ; then
      AC_SUBST(LOGVOL_LIB_PATH, [$logvol_lib])
   fi
   if test "x$logvol_bindir" != x ; then
      AC_PATH_PROG([h5ldump],[h5ldump],,[$logvol_bindir])
   else
      dnl Check h5ldump under $PATH
      AC_PATH_PROG([h5ldump],[h5ldump])
   fi
   if test "x$h5ldump" != x ; then
      AC_SUBST(H5LDUMP, [$h5ldump])
   fi
])
AM_CONDITIONAL(ENABLE_LOGVOL, [test "x$have_logvol" = xyes])

# Check for ADIOS2
have_adios2=no
AC_ARG_WITH([adios2],
   [AS_HELP_STRING([--with-adios2@<:@=INC,LIB | =DIR@:>@],
                   [Enable ADIOS-2 feature and provide the ADIOS-2 installation path(s):
    --with-adios2=INC,LIB for include and lib paths separated by a comma.
    --with-adios2=DIR for the path containing include/ and lib/ subdirectories.
    @<:@default: disabled@:>@
   ])], [
   if test "x${withval}" != xyes && test "x${withval}" != xno ; then
      ADIOS2_INSTALL=${withval}
   fi
   case $withval in
     *,*)
        adios2_inc="`echo $withval |cut -f1 -d,`"
        adios2_lib="`echo $withval |cut -f2 -d, -s`"
        adios2_bindir="$adios2_inc/../bin"
        ;;
     *)
        if test -n "$withval"; then
           adios2_inc="$withval/include"
           if test -d $withval/lib ; then
              adios2_lib="$withval/lib"
           elif test -d $withval/lib64 ; then
              adios2_lib="$withval/lib64"
           fi
           adios2_bindir="$ADIOS2_INSTALL/bin"
        fi
        ;;
   esac
   adios2_config=
   if test "x$ADIOS2_INSTALL" != x ; then
      AC_PATH_PROG([adios2_config],[adios2-config],,[$adios2_bindir])
   else
      dnl Check adios2-config under $PATH
      AC_PATH_PROG([adios2_config],[adios2-config])
   fi
   if test "x$adios2_config" != x ; then
      adios2_cflags=`$adios2_config --c-flags -m`
      adios2_lflags=`$adios2_config --c-libs --cxx-libs -m`
   else
      adios2_cflags="-I$adios2_inc"
      adios2_lflags="-L$adios2_lib"
   fi
   if test "x$adios2_cflags" != x ; then
      if test "x$CPPFLAGS" = x ; then
         CPPFLAGS="${adios2_cflags}"
      elif ! echo "${CPPFLAGS}" | ${EGREP} -q -w -- "${adios2_cflags}" ; then
         # append only if not already appear in CPPFLAGS
         CPPFLAGS="$CPPFLAGS ${adios2_cflags}"
      fi
   fi
   if test "x$adios2_lflags" != x ; then
      if test "x$LDFLAGS" = x ; then
         LDFLAGS="$adios2_lflags"
      elif ! echo "${LDFLAGS}" | ${EGREP} -q -w -- "${adios2_lflags}" ; then
         # append only if not already appear in LDFLAGS
         LDFLAGS="$LDFLAGS $adios2_lflags"
      fi
   fi
   AC_CHECK_HEADER([adios2_c.h], [have_adios2=yes], [have_adios2=no])
   if test "x$have_adios2" = xno ; then
      AC_MSG_ERROR([
      -----------------------------------------------------------------------
      Missing ADIOS2-header files 'adios2_c.h' required to build E3SM_IO. Use
      configure command-line option --with-adios2=/path/to/implementation
      to specify the location of ADIOS2 installation. Abort.
      -----------------------------------------------------------------------])
   fi
   if test "x$have_adios2" = xyes ; then
      AC_SUBST(ADIOS2_LIB_PATH, [$adios2_lib])
   fi
])
AM_CONDITIONAL(ENABLE_ADIOS2, [test "x$have_adios2" = xyes])

if test "x${have_pnc}" = xno && test "x${have_nc4}" = xno && test "x${have_hdf5}" = xno  && test "x${have_adios2}" = xno ; then
   AC_MSG_ERROR([
   -----------------------------------------------------------------------
    At least one API driver must be enabled.
    Use configure command-line option --with-pnetcdf=/path/to/implementation
    to specify the location of PnetCDF installation.
    Use configure command-line option --with-netcdf4=/path/to/implementation
    to specify the location of NetCDF installation.
    Use configure command-line option --with-hdf5=/path/to/implementation
    to specify the location of HDF5 installation.
    Use configure command-line option --with-adios2=/path/to/implementation
    to specify the location of ADIOS2 installation.
    Abort.
   -----------------------------------------------------------------------])
fi

AC_ARG_VAR(TESTMPIRUN, [MPI run command for "make ptest", @<:@default: mpiexec -n NP@:>@])
if test "x${TESTMPIRUN}" = x ; then
   dnl if TESTMPIRUN has not been set by users, then
   dnl set default to "mpiexec -n NP"
   UD_MPI_PATH_PROGS([TESTMPIRUN], [mpiexec mpirun srun])
fi

AC_ARG_VAR(TESTOUTDIR, [Output file directory for "make check" and "make ptest", @<:@default: ./test_output@:>@])
if test "x${TESTOUTDIR}" = x ; then
   dnl set default to current directory
   TESTOUTDIR=./test_output
fi

# Enable C++ 17 if available
if test "x${HAVE_CXX17}" = x1 ; then
   CXXFLAGS="$CXXFLAGS -std=c++17"
else
   CXXFLAGS="$CXXFLAGS -std=c++11"
fi

if test "x$have_pnc" = xyes ; then
   if test "x$LIBS" = x ; then
      LIBS="-lpnetcdf"
   else
      LIBS="$LIBS -lpnetcdf"
   fi
   AC_DEFINE([ENABLE_PNC], [1], [Enable PnetCDF I/O method])
fi

if test "x$have_netcdf4" = xyes ; then
   if test "x$LIBS" = x ; then
      LIBS="-lnetcdf"
   else
      LIBS="$LIBS -lnetcdf"
   fi
   AC_DEFINE([ENABLE_NETCDF4], [1], [Enable NetCDF 4 I/O method])
fi

if test "x$have_hdf5" = xyes ; then
   if test "x$LIBS" = x ; then
      LIBS="-lhdf5"
   else
      LIBS="$LIBS -lhdf5"
   fi
   AC_DEFINE([ENABLE_HDF5], [1], [Enable HDF5 I/O method])
   if test "x$have_multi_dset" = xyes ; then
      AC_DEFINE([HDF5_HAVE_DWRITE_MULTI], [1], [HDF5 has multi-dataset I/O support])
   fi
fi

if test "x$have_logvol" = xyes ; then
   if test "x$LIBS" = x ; then
      LIBS="-lH5VL_log"
   else
      LIBS="$LIBS -lH5VL_log"
   fi
   AC_DEFINE([ENABLE_LOGVOL], [1], [Enable HDF5 Log-layout based VOL I/O method])
fi

if test "x$have_adios2" = xyes ; then
   if test "x$LIBS" = x ; then
      LIBS="${adios2_lflags}"
   else
      LIBS="$LIBS ${adios2_lflags}"
   fi
   AC_DEFINE([ENABLE_ADIOS2], [1], [Enable ADIOS I/O method])

   # Check if type adios2_varinfo is defined, introduced in 2.8.0.
   # Use C++ compiler to check, due to ADIOSMacros.h included in adios2.h
   # includes <string>, a C++ header
   saved_CC=$CC
   CC=$MPICXX
   AC_CHECK_TYPES([adios2_varinfo], [], [], [[#include <adios2.h>
                                              #include <adios2_c.h>]])
   CC=$saved_CC
fi

AC_ARG_ENABLE([utilities],
    [AS_HELP_STRING([--disable-utilities],
                    [skip building utility programs
                     @<:@default: enabled@:>@])],
    [build_utils=${enableval}], [build_utils=yes]
)
AM_CONDITIONAL(BUILD_UTILS, [test "x$build_utils" = xyes])

AC_CONFIG_FILES(  Makefile \
                  src/Makefile \
                  src/drivers/Makefile \
                  src/cases/Makefile \
                  utils/Makefile \
                  datasets/Makefile
)

AC_OUTPUT
echo "------------------------------------------------------------------------------"
echo \
"
   ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

   Features:  Internal debug mode         - ${debug}
              PnetCDF                     - ${have_pnc}
              NetCDF-4                    - ${have_netcdf4}
              HDF5                        - ${have_hdf5}
              HDF5 log-layout based VOL   - ${have_logvol}
              HDF5 multi-dataset APIs     - ${have_multi_dset}
              ADIOS-2                     - ${have_adios2}
              Build utility programs      - ${build_utils}

   Compilers: MPICC    = ${MPICC}
              MPICXX   = ${MPICXX}"
if test "x${CPPFLAGS}" != x ; then
   echo "\
              CPPFLAGS = ${CPPFLAGS}"
fi
echo "\
              CFLAGS   = ${CFLAGS}"
if test "${CXXFLAGS}" != x ; then
   echo "\
              CXXFLAGS = ${CXXFLAGS}"
fi
if test "x${LDFLAGS}" != x ; then
   echo "\
              LDFLAGS  = ${LDFLAGS}"
fi
if test "x${LIBS}" != x ; then
   echo "\
              LIBS     = ${LIBS}"
fi
echo "\

   Now run 'make' to build the executable.

------------------------------------------------------------------------------"

