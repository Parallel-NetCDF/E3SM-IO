dnl
dnl Copyright (C) 2020, Northwestern University and Argonne National Laboratory
dnl See COPYRIGHT notice in top-level directory.
dnl

dnl -*- Mode: shell-script-mode; -*-
dnl Process this file with GNU autoconf(1) to produce a configure script.
dnl

dnl autoconf v2.69 was released in 2012-04-24
AC_PREREQ([2.69])
AC_INIT([e3sm_io],[1.0.0],[],[e3sm_io],[])

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
])
AC_CHECK_HEADER([pnetcdf.h], [have_pnc=yes], [have_pnc=no])
if test "x$have_pnc" = xno ; then
      AC_MSG_ERROR([
      -----------------------------------------------------------------------
      Missing PnetCDF-header files 'pnetcdf.h' required to build E3SM_IO. Use
      configure command-line option --with-pnetcdf=/path/to/implementation
      to specify the location of PnetCDF installation. Abort.
      -----------------------------------------------------------------------])
fi
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
      PnetCDF version too old. Abort.
      -----------------------------------------------------------------------])
fi

# Check for HDF5
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
      AC_MSG_WARN([
      -----------------------------------------------------------------------
      Missing HDF5-header files 'hdf5.h' required to build E3SM_IO. Use
      configure command-line option --with-hdf5=/path/to/implementation
      to specify the location of HDF5 installation. Disable HDF5 support.
      -----------------------------------------------------------------------])
   else
      AC_MSG_CHECKING([whether HDF5 parallel I/O is enabled])
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
      #include <hdf5.h>
      #if !defined(H5_HAVE_PARALLEL) || H5_HAVE_PARALLEL != 1
      #error Parallel I/O is not enabled in HDF5
      #endif
      ]])], [hdf5_parallel_io=yes], [hdf5_parallel_io=no])
      AC_MSG_RESULT([$hdf5_parallel_io])
      if test "x$hdf5_parallel_io" = xno ; then
          AC_MSG_WARN([
          -----------------------------------------------------------------------
              E3SM_IO requires parallel I/O feature enabled in HDF5. Disable HDF5 support.
          -----------------------------------------------------------------------])
          have_hdf5=no
      fi
   fi
   if test "x$have_hdf5" = xyes ; then
      AC_SEARCH_LIBS([H5Dwrite_multi], [hdf5], [have_multi_dset=yes], [have_multi_dset=no])
   fi
])
AM_CONDITIONAL(ENABLE_HDF5, [test "x$have_hdf5" = xyes])
AM_CONDITIONAL(HDF5_HAVE_DWRITE_MULTI, [test "x$have_multi_dset" = xyes])

# Check for logvol
have_logvol=no
AC_ARG_WITH([logvol],
   [AS_HELP_STRING([--with-logvol@<:@=INC,LIB | =DIR@:>@],
                   [Enable HDF5 log-based VOL feature and provide the VOL installation path(s):
    --with-logvol=INC,LIB for include and lib paths separated by a comma.
    --with-logvol=DIR for the path containing include/ and lib/ subdirectories.
    @<:@default: disabled@:>@
   ])], [
   case $withval in
     *,*)
        logvol_inc="`echo $withval |cut -f1 -d,`"
        logvol_lib="`echo $withval |cut -f2 -d, -s`"
        ;;
     *)
        if test -n "$withval"; then
          logvol_inc="$withval/include"
          logvol_lib="$withval/lib"
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
      AC_MSG_WARN([
      -----------------------------------------------------------------------
      Missing logvol-header files 'H5VL_log.h' required to build E3SM_IO. Use
      configure command-line option --with-logvol=/path/to/implementation
      to specify the location of logvol installation. Disable logvol support.
      -----------------------------------------------------------------------])
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
          adios2_lib="$withval/lib"
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
      adios2_lflags=`$adios2_config --c-libs -m`
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
      AC_MSG_WARN([
      -----------------------------------------------------------------------
      Missing ADIOS2-header files 'adios2_c.h' required to build E3SM_IO. Use
      configure command-line option --with-adios2=/path/to/implementation
      to specify the location of ADIOS2 installation. Disable ADIOS2 support.
      -----------------------------------------------------------------------])
   fi
])
AM_CONDITIONAL(ENABLE_ADIOS2, [test "x$have_adios2" = xyes])

if test "x$LIBS" = x ; then
LIBS="-lpnetcdf"
else
LIBS="$LIBS -lpnetcdf"
fi

if test "x$have_hdf5" = xyes ; then
    if test "x$LIBS" = x ; then
    LIBS="-lhdf5"
    else
    LIBS="$LIBS -lhdf5"
    fi
    AC_DEFINE([ENABLE_HDF5], [1], [Build HDF5 benchmark])
    if test "x$have_multi_dset" = xyes ; then
      AC_DEFINE([HDF5_HAVE_DWRITE_MULTI], [1], [Have multi-dataset I/O support])
    fi
fi

if test "x$have_logvol" = xyes ; then
    if test "x$LIBS" = x ; then
    LIBS="-lH5VL_log"
    else
    LIBS="$LIBS -lH5VL_log"
    fi
    AC_DEFINE([ENABLE_LOGVOL], [1], [Build logvol benchmark])
fi

if test "x$have_adios2" = xyes ; then
    if test "x$LIBS" = x ; then
    LIBS="${adios2_lflags}"
    else
    LIBS="$LIBS ${adios2_lflags}"
    fi
    AC_DEFINE([ENABLE_ADIOS2], [1], [Build ADIOS2 benchmark])
fi

AC_CONFIG_FILES(  Makefile \
                  src/Makefile \
                  src/drivers/Makefile \
                  src/cases/Makefile \
                  utils/Makefile
)

AC_OUTPUT
echo "------------------------------------------------------------------------------"
echo \
"
   ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

   Features:  Internal debug mode         - ${debug}
              ADIOS-2 support             - ${have_adios2}
              HDF5 log-based VOL support  - ${have_logvol}

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

