#
# Copyright (C) 2021, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# @configure_input@

AM_CPPFLAGS  = -I${top_srcdir}/src
AM_CPPFLAGS += -I${top_srcdir}/src/drivers
AM_CPPFLAGS += -I${top_srcdir}/src/cases

AM_M4FLAGS = -I${srcdir}

M4FLAGS += -I${top_srcdir}/m4

bin_PROGRAMS = e3sm_io

H_SRCS          = e3sm_io.h \
                  e3sm_io_err.h

M4H_SRCS        = e3sm_io_profile.m4h

e3sm_io_SOURCES = e3sm_io.c \
                  e3sm_io_core.cpp \
                  read_decomp.cpp \
                  calc_metadata.c

e3sm_io_LDADD   = cases/libe3sm_io_cases.la \
                  drivers/libe3sm_io_drivers.la

if ENABLE_LOGVOL
   e3sm_io_LDADD += -lH5VL_log
endif

BUILT_SOURCES = $(M4H_SRCS:.m4h=.hpp)

M4_SRCS       = e3sm_io_profile.m4 \
                e3sm_io_profile_timers.m4

if E3SM_IO_PROFILING
   noinst_LTLIBRARIES = libe3sm_io_profiling.la
   libe3sm_io_profiling_la_SOURCES = $(M4_SRCS:.m4=.cpp)
   e3sm_io_LDADD   += libe3sm_io_profiling.la
endif

.m4.cpp:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

.m4h.hpp:
	$(M4) $(AM_M4FLAGS) $(M4FLAGS) $< >$@

SUBDIRS = drivers cases
DIST_SUBDIRS = drivers cases

EXTRA_DIST = $(H_SRCS) $(M4H_SRCS) $(M4_SRCS)

CLEANFILES = $(M4_SRCS:.m4=.cpp) $(M4H_SRCS:.m4h=.hpp)

