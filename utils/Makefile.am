# @configure_input@

SUFFIXES = .o .c .cpp

AM_CPPFLAGS  = -I${top_srcdir}/src
AM_CPPFLAGS += -I${top_srcdir}/src/drivers
AM_CPPFLAGS += -I${top_srcdir}/src/cases

bin_PROGRAMS = datstat

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

datstat_SOURCES = datstat.cpp

bpstat_SOURCES = bpstat.cpp
bpstat_LDADD = -lstdc++fs
