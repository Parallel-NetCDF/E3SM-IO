AM_CPPFLAGS  = -I${top_srcdir}/src
AM_CPPFLAGS += -I${top_srcdir}/src/drivers
AM_CPPFLAGS += -I${top_srcdir}/src/cases

M4FLAGS += -I${top_srcdir}/m4

M4_SRCS = e3sm_io_profile.m4 

M4H_SRCS = e3sm_io_profile.m4h

bin_PROGRAMS = e3sm_io

e3sm_io_SOURCES = e3sm_io.c read_decomp.c e3sm_io_core.cpp
e3sm_io_SOURCES += e3sm_io.h
e3sm_io_SOURCES += $(M4_SRCS:.m4=.cpp) $(M4H_SRCS:.m4h=.hpp)
e3sm_io_LDADD = cases/libe3sm_io_cases.a \
                drivers/libe3sm_io_drivers.a \
                -lpnetcdf

if ENABLE_LOGVOL
e3sm_io_LDADD += -lH5VL_log
endif

$(M4_SRCS:.m4=.cpp): Makefile
$(M4H_SRCS:.m4h=.hpp): Makefile

.m4.cpp:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

.m4h.hpp:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

BUILT_SOURCES = $(M4_SRCS:.m4=.cpp) $(M4H_SRCS:.m4h=.hpp)

SUBDIRS = cases drivers
DIST_SUBDIRS = cases drivers