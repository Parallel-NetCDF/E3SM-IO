#
# Copyright (C) 2021, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

SUFFIXES = .a .o .c .m4 .h

AM_CXXFLAGS  = -I${top_srcdir}
AM_CXXFLAGS += -I${top_srcdir}/drivers
AM_CXXFLAGS += -I${top_srcdir}/cases

noinst_LIBRARIES = libe3sm_io_cases.a

C_SRCS = var_io_F_case.cpp var_io_G_case.cpp header_io_F_case.cpp header_io_G_case.cpp
H_SRCS = e3sm_io_cases.hpp

libe3sm_io_cases_a_SOURCES = $(C_SRCS) $(H_SRCS)

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out

