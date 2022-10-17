#
# Copyright (C) 2021, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# @configure_input@

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src . datasets
DIST_SUBDIRS = src datasets utils

if BUILD_UTILS
   SUBDIRS += utils
endif

TESTS_ENVIRONMENT  = export srcdir="$(srcdir)";
TESTS_ENVIRONMENT += export TESTMPIRUN="$(TESTMPIRUN)";
TESTS_ENVIRONMENT += export TESTOUTDIR="$(TESTOUTDIR)";

if ENABLE_PNC
   TESTS_ENVIRONMENT += export ENABLE_PNC=1;
   TESTS_ENVIRONMENT += export PNETCDF_LIB_PATH="@PNETCDF_LIB_PATH@";
endif
if ENABLE_NETCDF4
   TESTS_ENVIRONMENT += export ENABLE_NETCDF4=1;
   TESTS_ENVIRONMENT += export NETCDF4_LIB_PATH="@NETCDF4_LIB_PATH@";
endif
if ENABLE_ADIOS2
   TESTS_ENVIRONMENT += export ENABLE_ADIOS2=1;
   TESTS_ENVIRONMENT += export ADIOS2_LIB_PATH="@ADIOS2_LIB_PATH@";
endif
if ENABLE_HDF5
   TESTS_ENVIRONMENT += export ENABLE_HDF5=1;
   TESTS_ENVIRONMENT += export HDF5_LIB_PATH="@HDF5_LIB_PATH@";
endif
if HDF5_HAVE_MULTI_DATASET_API
   TESTS_ENVIRONMENT += export HDF5_HAVE_MULTI_DATASET_API=1;
endif
if ENABLE_LOGVOL
   TESTS_ENVIRONMENT += export ENABLE_LOGVOL=1;
   TESTS_ENVIRONMENT += export LOGVOL_LIB_PATH="@LOGVOL_LIB_PATH@";
   TESTS_ENVIRONMENT += export H5LDUMP="@H5LDUMP@";
   TESTS_ENVIRONMENT += export H5LREPLAY="@H5LREPLAY@";
   TESTS_ENVIRONMENT += export LOGVOL_LIB_SHARED="@LOGVOL_LIB_SHARED@";
endif

TESTS = test.sh

EXTRA_DIST =   test.sh \
               m4/foreach_idx.m4 \
               m4/foreach.m4 \
               m4/list_len.m4 \
               m4/utils.m4 \
               README.md \
               COPYRIGHT \
               AUTHORS \
               docs

ptest: src/e3sm_io
	@echo "==========================================================="
	@echo "    Parallel testing on 16 MPI processes"
	@echo "==========================================================="
	@$(TESTS_ENVIRONMENT) \
	$(srcdir)/test.sh 16 || exit 1

clean-local:
	-rm -rf test_output/

.PHONY: ptest

