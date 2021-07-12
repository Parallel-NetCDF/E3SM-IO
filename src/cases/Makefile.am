#
# Copyright (C) 2021, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#

SUFFIXES = .a .o .c .m4 .h .cpp .hpp

AM_CPPFLAGS  = -I${top_srcdir}/src
AM_CPPFLAGS += -I${top_srcdir}/src/drivers
AM_CPPFLAGS += -I${top_srcdir}/src/cases

noinst_LIBRARIES = libe3sm_io_cases.a

CXX_SRCS =  var_io_F_case.cpp \
            var_io_G_case.cpp \
            e3sm_io_case_F.cpp \
            e3sm_io_case_G.cpp \
            header_def_F_case.cpp \
            header_inq_F_case.cpp \
            header_def_G_case.cpp \
            header_inq_G_case.cpp \
            pnetcdf_blob_var_io_F_case.cpp \
            pnetcdf_blob_var_io_G_case.cpp

H_SRCS =    e3sm_io_case.hpp \
            e3sm_io_F_case.hpp \
            e3sm_io_G_case.hpp 

if ENABLE_ADIOS2
CXX_SRCS += var_io_F_case_scorpio.cpp \
            e3sm_io_case_F_scorpio.cpp \
            header_io_F_case_scorpio.cpp \
            var_io_G_case_scorpio.cpp \
            e3sm_io_case_G_scorpio.cpp \
            header_io_G_case_scorpio.cpp

H_SRCS +=   e3sm_io_case_scorpio.hpp \
            e3sm_io_F_case_scorpio.hpp \
            e3sm_io_G_case_scorpio.hpp
endif

libe3sm_io_cases_a_SOURCES = $(CXX_SRCS) $(H_SRCS)

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out

