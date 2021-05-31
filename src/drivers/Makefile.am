#
# Copyright (C) 2021, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

SUFFIXES = .a .o .c .m4 .h

AM_CPPFLAGS  = -I${top_srcdir}/src
AM_CPPFLAGS += -I${top_srcdir}/src/drivers
AM_CPPFLAGS += -I${top_srcdir}/src/cases

noinst_LIBRARIES = libe3sm_io_drivers.a

C_SRCS = e3sm_io_driver_pnc.cpp
H_SRCS = e3sm_io_driver_pnc.hpp

if ENABLE_HDF5
C_SRCS +=   e3sm_io_driver_hdf5.cpp \
            e3sm_io_driver_hdf5_agg.cpp
H_SRCS +=   e3sm_io_driver_hdf5.hpp \
            e3sm_io_driver_hdf5_int.hpp
endif

if ENABLE_ADIOS2
C_SRCS += e3sm_io_driver_adios2.cpp
H_SRCS += e3sm_io_driver_adios2.hpp
endif

libe3sm_io_drivers_a_SOURCES = $(C_SRCS) $(H_SRCS)

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out

