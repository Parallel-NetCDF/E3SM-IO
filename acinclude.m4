#
# Copyright (C) 2019, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
#

dnl Convert a string to all uppercase.
dnl
define([uppercase],
[translit($1, abcdefghijklmnopqrstuvwxyz, ABCDEFGHIJKLMNOPQRSTUVWXYZ)])

dnl
dnl check the availability of one MPI executable in $2
dnl
dnl $2 can be a single command, This is the case when user set the environment.
dnl The variable may contain the executable name followed by zeor or more
dnl command-line options. In the latter case, we check the first string token,
dnl the command name, and ignore the rest command-line options. For example,
dnl UD_MPI_PATH_PROG([MPICC], [mpicc -O2])
dnl
dnl In addition, the first token of $2 may contain the full path of
dnl the command. For example, UD_MPI_PATH_PROG([MPICC], [/usr/bin/mpicc -O2])
dnl
AC_DEFUN([UD_MPI_PATH_PROG], [
   if test "x$2" = x ; then
      AC_MSG_ERROR("2nd argument cannot be NULL")
   fi

   dnl 1st token in $2 must be the program name, rests are command-line options
   ac_first_token=`echo $2 | cut -d" " -f1`
   ac_rest_tokens=`echo $2 | cut -d" " -s -f2-`

   dnl First check if ac_first_token contain a full path
   dnl If yes, check, check if the file exists. Need not check MPI_INSTALL.
   ac_mpi_prog_path=`AS_DIRNAME(["$ac_first_token"])`
   if test "x$ac_mpi_prog_path" != "x." ; then
      AC_MSG_CHECKING([whether $ac_first_token exists and is executable])
      if test -x "$ac_first_token" ; then
         AC_MSG_RESULT([yes])
         $1="$2"
      else
         AC_MSG_RESULT([no])
         $1=
      fi
   else
      dnl ac_first_token does not contain a full path
      ac_mpi_prog_$1=
      if test "x$MPI_INSTALL" != x ; then
         dnl First, check if it can be found under $MPI_INSTALL, i.e.
         dnl --with-mpi is used on configure command line
         if test -d "${MPI_INSTALL}/bin" ; then
            AC_MSG_CHECKING([$ac_first_token under ${MPI_INSTALL}/bin])
            if test -x "$MPI_INSTALL/bin/$ac_first_token" ; then
               AC_MSG_RESULT([yes])
               ac_mpi_prog_$1=$MPI_INSTALL/bin/$ac_first_token
            else
               AC_MSG_RESULT([no])
            fi
         else
            dnl ${MPI_INSTALL}/bin does not exist, search $MPI_INSTALL
            AC_MSG_CHECKING([$ac_first_token under ${MPI_INSTALL}])
            if test -x "$MPI_INSTALL/$ac_first_token" ; then
               AC_MSG_RESULT([yes])
               ac_mpi_prog_$1=$MPI_INSTALL/$ac_first_token
            else
               AC_MSG_RESULT([no])
            fi
         fi
         if test "x$ac_mpi_prog_$1" != x ; then
            $1="${ac_mpi_prog_$1} $ac_rest_tokens"
         else
            $1=
         fi
      else
         dnl MPI_INSTALL is not set, i.e. --with-mpi is not used
         AC_PATH_PROG([ac_mpi_prog_$1], [$ac_first_token])
         if test "x$ac_mpi_prog_$1" != x ; then
            $1="${ac_mpi_prog_$1} $ac_rest_tokens"
         else
            $1=
         fi
      fi
   fi
])

dnl
dnl check the availability of a list of MPI executables
dnl
dnl Note $2 can be a list of executable commands to be searched with each
dnl command being the executable file name without command-line option. This is
dnl the case when user does not set the environment variable, for example
dnl MPICC, and we must search one from the candidate list. For example,
dnl UD_MPI_PATH_PROGS([MPICC], [mpicc mpixlc mpifccpx mpipgcc])
dnl
AC_DEFUN([UD_MPI_PATH_PROGS], [
   ac_mpi_prog_$1=
   if test "x$MPI_INSTALL" != x ; then
      if test -d "${MPI_INSTALL}/bin" ; then
         AC_PATH_PROGS([ac_mpi_prog_$1], [$2], [], [$MPI_INSTALL/bin])
      else
         dnl ${MPI_INSTALL}/bin does not exist, search $MPI_INSTALL
         AC_PATH_PROGS([ac_mpi_prog_$1], [$2], [], [$MPI_INSTALL])
      fi
   else
      AC_PATH_PROGS([ac_mpi_prog_$1], [$2])
   fi
   if test "x${ac_mpi_prog_$1}" = x ; then
      dnl AC_CHECK_FILES fails when $2 is not found in cross compile
      dnl AC_CHECK_FILES([$2], [ac_mpi_prog_$1=$2])
      AC_PATH_PROGS([ac_mpi_prog_$1], [$2])
      dnl AC_CHECK_PROGS([ac_mpi_prog_$1], [$2])
      dnl AC_CHECK_PROGS([ac_mpi_prog_$1], [$2], [], [/])
      dnl ac_first_token=`echo $2 | cut -d" " -f1`
      dnl if test -f $ac_first_token ; then
         dnl ac_mpi_prog_$1=$2
      dnl fi
   fi
   $1=${ac_mpi_prog_$1}
])

dnl
dnl Check how sed command handling in-place option -i and define SED_I
dnl
AC_DEFUN([UD_PROG_SED_I],
[
   AC_REQUIRE([AC_PROG_SED])
   AC_CACHE_CHECK([for sed handling option -i ], ac_cv_SED_I,[
   cat > conftest.sed_i <<EOF
   test str1
EOF
   ac_cv_err=`$SED -i '' -e 's|str1|str2|g' conftest.sed_i 2>&1`
   if test "x$ac_cv_err" = x ; then
      ac_cv_SED_I="$SED -i ''"
   else
      ac_cv_err=`sed -i'' -e 's|str1|str2|g' conftest.sed_i 2>&1`
      if test "x$ac_cv_err" = x ; then
         ac_cv_SED_I="$SED -i''"
      else
         AC_MSG_ERROR("No proper sed -i option found")
      fi
   fi
   AS_UNSET(ac_cv_err)])
   SED_I="$ac_cv_SED_I"
   AC_SUBST(SED_I)
   rm -f conftest.sed_i
])

