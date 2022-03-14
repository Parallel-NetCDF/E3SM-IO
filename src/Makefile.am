AM_CPPFLAGS  = -I${top_srcdir}/src
AM_CPPFLAGS += -I${top_srcdir}/src/drivers
AM_CPPFLAGS += -I${top_srcdir}/src/cases

AM_M4FLAGS = -I${srcdir}

M4FLAGS += -I${top_srcdir}/m4

M4_SRCS = e3sm_io_profile.m4 \
          e3sm_io_profile_timers.m4

M4H_SRCS = e3sm_io_profile.m4h

bin_PROGRAMS = e3sm_io

H_SRCS          = e3sm_io.h \
                  e3sm_io_err.h

e3sm_io_SOURCES = e3sm_io.c \
                  e3sm_io_core.cpp \
                  calc_metadata.c

e3sm_io_SOURCES += $(H_SRCS) \
                   $(M4_SRCS:.m4=.cpp) \
                   $(M4H_SRCS:.m4h=.hpp)

e3sm_io_LDADD = cases/libe3sm_io_cases.a \
                drivers/libe3sm_io_drivers.a

if ENABLE_LOGVOL
e3sm_io_LDADD += -lH5VL_log
endif

if ENABLE_PNC
e3sm_io_SOURCES += read_decomp_pnc.c 
endif
if ENABLE_NETCDF4
e3sm_io_SOURCES += read_decomp_nc4.c 
endif

$(M4_SRCS:.m4=.cpp): Makefile e3sm_io_profile_timers.m4
$(M4H_SRCS:.m4h=.hpp): Makefile e3sm_io_profile_timers.m4

.m4.cpp: e3sm_io_profile_timers.m4
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

.m4h.hpp: e3sm_io_profile_timers.m4
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

BUILT_SOURCES = $(M4_SRCS:.m4=.cpp) $(M4H_SRCS:.m4h=.hpp)

SUBDIRS = cases drivers
DIST_SUBDIRS = cases drivers

EXTRA_DIST = $(M4H_SRCS) $(M4_SRCS)

CLEANFILES = $(M4_SRCS:.m4=.cpp) $(M4H_SRCS:.m4h=.hpp)
