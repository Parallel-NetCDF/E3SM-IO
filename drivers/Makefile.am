#
# Copyright (C) 2021, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

SUFFIXES = .a .o .c .m4 .h

AM_CPPFLAGS  = -I${top_srcdir}
AM_CPPFLAGS += -I${top_srcdir}/drivers
AM_CPPFLAGS += -I${top_srcdir}/cases

noinst_LIBRARIES = libe3sm_io_drivers.a

C_SRCS = e3sm_io_pnc.cpp
H_SRCS = e3sm_io_pnc.hpp

if ENABLE_HDF5
C_SRCS += e3sm_io_hdf5.cpp
H_SRCS += e3sm_io_hdf5.hpp
endif


libe3sm_io_drivers_a_SOURCES = $(C_SRCS) $(H_SRCS)

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out

