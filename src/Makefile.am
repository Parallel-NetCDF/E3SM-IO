AM_CPPFLAGS  = -I${top_srcdir}/src
AM_CPPFLAGS += -I${top_srcdir}/src/drivers
AM_CPPFLAGS += -I${top_srcdir}/src/cases

bin_PROGRAMS = e3sm_io

e3sm_io_SOURCES = e3sm_io.c read_decomp.c e3sm_io_core.cpp
e3sm_io_SOURCES += e3sm_io.h
e3sm_io_LDADD = cases/libe3sm_io_cases.a \
                drivers/libe3sm_io_drivers.a \
                -lpnetcdf

if ENABLE_LOGVOL
e3sm_io_LDADD += -lH5VL_log
endif

SUBDIRS = cases drivers
DIST_SUBDIRS = cases drivers