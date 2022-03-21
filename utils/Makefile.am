# @configure_input@

SUFFIXES = .o .c .cpp

AM_CPPFLAGS  = -I${top_srcdir}/src
AM_CPPFLAGS += -I${top_srcdir}/src/drivers
AM_CPPFLAGS += -I${top_srcdir}/src/cases

bin_PROGRAMS = datstat dat2decom

if ENABLE_ADIOS2
if HAVE_CXX17
   bin_PROGRAMS += bpstat
endif
endif

if ENABLE_PNC
bin_PROGRAMS += dat2nc \
                pnetcdf_blob_replay

dat2nc_SOURCES = dat2nc.c
dat2nc_LDADD = -lpnetcdf

pnetcdf_blob_replay_SOURCES = pnetcdf_blob_replay.c
pnetcdf_blob_replay_LDADD = -lpnetcdf
endif

dat2decom_SOURCES = dat2decom.cpp 
dat2decom_LDADD = ${top_builddir}/src/drivers/libe3sm_io_drivers.a ${top_builddir}/src/e3sm_io_profile.o

datstat_SOURCES = datstat.cpp

bpstat_SOURCES = bpstat.cpp
bpstat_LDADD = -lstdc++fs
