#
# Copyright (C) 2021, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# @configure_input@

SUFFIXES = .o .c .cpp

AM_CPPFLAGS  = -I${top_srcdir}/src
AM_CPPFLAGS += -I${top_srcdir}/src/drivers
AM_CPPFLAGS += -I${top_srcdir}/src/cases

bin_PROGRAMS = datstat dat2decomp decomp_copy
TESTPROGRAMS = datstat dat2decomp decomp_copy

if ENABLE_ADIOS2
if HAVE_CXX17
   bin_PROGRAMS += bpstat
endif
endif

if ENABLE_PNC
TESTPROGRAMS += dat2nc

bin_PROGRAMS += dat2nc \
                pnetcdf_blob_replay

dat2nc_SOURCES = dat2nc.c
dat2nc_LDADD = -lpnetcdf

pnetcdf_blob_replay_SOURCES = pnetcdf_blob_replay.c
pnetcdf_blob_replay_LDADD = -lpnetcdf
endif

dat2decomp_SOURCES = dat2decomp.cpp 
dat2decomp_LDADD = ../src/drivers/libe3sm_io_drivers.la
if E3SM_IO_PROFILING
   dat2decomp_LDADD += ../src/libe3sm_io_profiling.la
endif

decomp_copy_SOURCES = decomp_copy.cpp 
decomp_copy_LDADD = ../src/drivers/libe3sm_io_drivers.la
if E3SM_IO_PROFILING
   decomp_copy_LDADD += ../src/libe3sm_io_profiling.la
endif

datstat_SOURCES = datstat.cpp

bpstat_SOURCES = bpstat.cpp
if HAVE_FILESYSTEM_HPP
   bpstat_LDADD = -lstdc++fs
endif

BUILT_SOURCES = list_f_case.txt

list_f_case.txt:
	ls -1 $(top_srcdir)/datasets/piodecomp16tasks*.dat > $@

AM_TESTS_ENVIRONMENT  = export TESTPROGRAMS="$(TESTPROGRAMS)";
AM_TESTS_ENVIRONMENT += export top_srcdir="$(top_srcdir)";
AM_TESTS_ENVIRONMENT += export builddir="$(builddir)";
AM_TESTS_ENVIRONMENT += export SED="$(SED)";

ALL_LIBRARY_PATH = .

if ENABLE_PNC
   ALL_LIBRARY_PATH += ":@PNETCDF_LIB_PATH@"
   AM_TESTS_ENVIRONMENT += export ENABLE_PNC=1;
endif
if ENABLE_NETCDF4
   ALL_LIBRARY_PATH += ":@NETCDF4_LIB_PATH@"
   AM_TESTS_ENVIRONMENT += export ENABLE_NETCDF4=1;
endif
if ENABLE_HDF5
   ALL_LIBRARY_PATH += ":@HDF5_LIB_PATH@"
   AM_TESTS_ENVIRONMENT += export ENABLE_HDF5=1;
endif
if HDF5_HAVE_MULTI_DATASET_API
   AM_TESTS_ENVIRONMENT += export HDF5_HAVE_MULTI_DATASET_API=1;
endif
if ENABLE_LOGVOL
   ALL_LIBRARY_PATH += ":@LOGVOL_LIB_PATH@"
   AM_TESTS_ENVIRONMENT += export ENABLE_LOGVOL=1;
endif
if ENABLE_ADIOS2
   ALL_LIBRARY_PATH += ":@ADIOS2_LIB_PATH@"
   AM_TESTS_ENVIRONMENT += export ENABLE_ADIOS2=1;
endif
AM_TESTS_ENVIRONMENT += export ALL_LIBRARY_PATH="$(ALL_LIBRARY_PATH)";

TESTS = $(TESTPROGRAMS)
TEST_EXTENSIONS = .sh
LOG_COMPILER = $(srcdir)/wrap_runs.sh
SH_LOG_COMPILER =

EXTRA_DIST = wrap_runs.sh

DISTCLEANFILES = list_f_case.txt

