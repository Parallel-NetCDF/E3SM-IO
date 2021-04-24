AM_CPPFLAGS  = -I${top_srcdir}
AM_CPPFLAGS += -I${top_srcdir}/drivers
AM_CPPFLAGS += -I${top_srcdir}/cases

bin_PROGRAMS = e3sm_io dat2nc

e3sm_io_SOURCES = e3sm_io.cpp read_decomp.cpp 
e3sm_io_SOURCES += e3sm_io.hpp
e3sm_io_LDADD = drivers/libe3sm_io_drivers.a \
                cases/libe3sm_io_cases.a \
                -lpnetcdf

if ENABLE_LOGVOL
e3sm_io_LDADD += -lH5VL_log
endif

dat2nc_SOURCES = dat2nc.c 
dat2nc_LDADD = -lpnetcdf

SUBDIRS = drivers cases
DIST_SUBDIRS = drivers cases