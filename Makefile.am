AM_CPPFLAGS = -I$(top_srcdir)

bin_PROGRAMS = e3sm_io

e3sm_io_SOURCES = e3sm_io.cpp read_decomp.cpp 
e3sm_io_SOURCES += e3sm_io.hpp
e3sm_io_LDADD = drivers/libe3sm_io_drivers.a \
                cases/libe3sm_io_cases.a \
                -lpnetcdf

if ENABLE_LOGVOL
e3sm_io_LDADD += -lH5VL_log
endif

SUBDIRS = drivers cases
DIST_SUBDIRS = drivers cases