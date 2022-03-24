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
dat2decomp_LDADD = ${top_builddir}/src/drivers/libe3sm_io_drivers.a ${top_builddir}/src/e3sm_io_profile.o

decomp_copy_SOURCES = decomp_copy.cpp 
decomp_copy_LDADD = ${top_builddir}/src/drivers/libe3sm_io_drivers.a ${top_builddir}/src/e3sm_io_profile.o

datstat_SOURCES = datstat.cpp

bpstat_SOURCES = bpstat.cpp
bpstat_LDADD = -lstdc++fs

AM_TESTS_ENVIRONMENT  = export TESTPROGRAMS="$(TESTPROGRAMS)";
AM_TESTS_ENVIRONMENT += export top_srcdir="$(top_srcdir)";
if ENABLE_PNC
   AM_TESTS_ENVIRONMENT += export ENABLE_PNC=1;
endif
if ENABLE_NETCDF4
   AM_TESTS_ENVIRONMENT += export ENABLE_NETCDF4=1;
endif
if ENABLE_HDF5
   AM_TESTS_ENVIRONMENT += export ENABLE_HDF5=1;
endif
if ENABLE_ADIOS2
   AM_TESTS_ENVIRONMENT += export ENABLE_ADIOS2=1;
endif

TESTS = $(TESTPROGRAMS)
TEST_EXTENSIONS = .sh
LOG_COMPILER = $(srcdir)/wrap_runs.sh
SH_LOG_COMPILER =

EXTRA_DIST = wrap_runs.sh

