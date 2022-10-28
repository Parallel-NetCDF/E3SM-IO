# Copyright (C) 2021, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# @configure_input@

SUFFIXES = .o .c .h .cpp .hpp

AM_CPPFLAGS  = -I${top_srcdir}/src
AM_CPPFLAGS += -I${top_builddir}/src
AM_CPPFLAGS += -I${top_srcdir}/src/drivers
AM_CPPFLAGS += -I${top_srcdir}/src/cases

noinst_LTLIBRARIES = libe3sm_io_drivers.la

C_SRCS = e3sm_io_driver.cpp 

H_SRCS = e3sm_io_driver.hpp 

if ENABLE_HDF5
C_SRCS +=   e3sm_io_driver_hdf5.cpp \
            e3sm_io_driver_h5blob.cpp \
            e3sm_io_driver_hdf5_md.cpp \
            blob_ncmpio.c

H_SRCS +=   e3sm_io_driver_hdf5.hpp \
            e3sm_io_driver_h5blob.hpp \
            blob_ncmpio.h
endif

if ENABLE_PNC
C_SRCS +=   e3sm_io_driver_pnc.cpp
H_SRCS +=   e3sm_io_driver_pnc.hpp 
endif

if ENABLE_NETCDF4
C_SRCS +=   e3sm_io_driver_nc4.cpp
H_SRCS +=   e3sm_io_driver_nc4.hpp 
endif

if ENABLE_LOGVOL
C_SRCS +=   e3sm_io_driver_hdf5_log.cpp
H_SRCS +=   e3sm_io_driver_hdf5_log.hpp 
endif

if ENABLE_ADIOS2
C_SRCS += e3sm_io_driver_adios2.cpp
H_SRCS += e3sm_io_driver_adios2.hpp
endif

libe3sm_io_drivers_la_SOURCES = $(C_SRCS) $(H_SRCS)

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out

