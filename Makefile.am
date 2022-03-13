SUBDIRS = src utils datasets

TESTS_ENVIRONMENT  = export srcdir="$(srcdir)";
TESTS_ENVIRONMENT += export TESTMPIRUN="@TESTMPIRUN@";

if ENABLE_PNC
   TESTS_ENVIRONMENT += export ENABLE_PNC=1;
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
if HDF5_HAVE_DWRITE_MULTI
   TESTS_ENVIRONMENT += export HDF5_HAVE_DWRITE_MULTI=1;
endif
if ENABLE_LOGVOL
   TESTS_ENVIRONMENT += export ENABLE_LOGVOL=1;
   TESTS_ENVIRONMENT += export LOGVOL_LIB_PATH="@LOGVOL_LIB_PATH@";
endif

TESTS = test.sh

EXTRA_DIST =   test.sh \
               m4/foreach_idx.m4 \
               m4/foreach.m4 \
               m4/list_len.m4 \
               m4/utils.m4 \
               README.md \
               INSTALL.md \
               COPYRIGHT \
               AUTHORS \
               slurm.knl

ptest: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    Parallel testing on 16 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/test.sh 16 || exit 1

clean-local:
	-rm -rf test_output/

.PHONY: ptest


