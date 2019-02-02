
MPICC		= mpicc
CFLAGS          = -O2

PnetCDF_DIR	= $(HOME)/PnetCDF/1.11.0

INCLUDES	= -I$(PnetCDF_DIR)/include -I.
LDFLAGS		= -L$(PnetCDF_DIR)/lib
LIBS		= -lpnetcdf $(shell $(PnetCDF_DIR)/bin/pnetcdf-config --libs)

.c.o:
	$(MPICC) $(CFLAGS) $(INCLUDES) -c $<

all: e3sm_io

dat2nc: dat2nc.o
	$(MPICC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

e3sm_io.o: e3sm_io.h
def_header.o: e3sm_io.h

e3sm_io: e3sm_io.o def_header.o
	$(MPICC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

# romio_patch.c contains fix in https://github.com/pmodels/mpich/pull/3089
# ROMIO_PATCH	= -Wl,--wrap=ADIOI_Type_create_hindexed_x -l:libmpich_intel.a
# ROMIO_PATCH	= -Wl,--wrap=ADIOI_Type_create_hindexed_x -l:libmpi.a
ROMIO_PATCH	= -Wl,--wrap=ADIOI_Type_create_hindexed_x

e3sm_io.romio_patch: e3sm_io.o def_header.o romio_patch.o
	$(MPICC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS) $(ROMIO_PATCH)

romio_patch.o: romio_patch.c
	$(MPICC) $(CFLAGS) -o $@ -c $^

clean:
	rm -f core.* *.o dat2nc e3sm_io e3sm_io.romio_patch
	rm -f testfile_h0_varn.nc testfile_h1_varn.nc
	rm -f testfile_h0_vard.nc testfile_h1_vard.nc

.PHONY: clean

