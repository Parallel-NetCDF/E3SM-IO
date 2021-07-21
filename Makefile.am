SUBDIRS = src utils
DIST_SUBDIRS = src utils

TESTS_ENVIRONMENT  = export srcdir="$(srcdir)";
TESTS_ENVIRONMENT += export HDF5_LIB_PATH="@HDF5_LIB_PATH@";
TESTS = seq_runs.sh

EXTRA_DIST = seq_runs.sh parallel_run.sh

CLEANFILES = can_F_out_h0.nc         can_F_out_h1.nc \
             blob_F_out_h0.nc.0000   blob_F_out_h0.nc.0001 \
             blob_F_out_h1.nc.0000   blob_F_out_h1.nc.0001 \
             blob_F_out_h0.h5.0000   blob_F_out_h0.h5.0001 \
             blob_F_out_h1.h5.0000   blob_F_out_h1.h5.0001 \
             can_G_out.nc \
             blob_G_out.nc.0000   blob_F_out.nc.0001 \
             blob_G_out.h5.0000   blob_F_out.h5.0001

ptest: $(check_PROGRAMS)
	@echo "==========================================================="
	@echo "    Parallel testing on 16 MPI processes"
	@echo "==========================================================="
	$(TESTS_ENVIRONMENT) \
	$(srcdir)/parallel_run.sh 16 || exit 1

clean-local:
	-rm -rf blob_F_out_h0.bp.dir blob_F_out_h1.bp.dir blob_G_out.bp.dir

.PHONY: ptest


