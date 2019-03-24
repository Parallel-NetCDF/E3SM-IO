
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

e3sm_io.o: e3sm_io.c e3sm_io.h
read_decomp.o: read_decomp.c e3sm_io.h

header_io_F_case.o: header_io_F_case.c e3sm_io.h
var_io_F_case.o: var_io_F_case.c e3sm_io.h
header_io_G_case.o: header_io_G_case.c e3sm_io.h
var_io_G_case.o: var_io_G_case.c e3sm_io.h

OBJS = read_decomp.o header_io_F_case.o var_io_F_case.o header_io_G_case.o var_io_G_case.o

e3sm_io: e3sm_io.o $(OBJS)
	$(MPICC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

# romio_patch.c contains fix in https://github.com/pmodels/mpich/pull/3089
# ROMIO_PATCH	= -Wl,--wrap=ADIOI_Type_create_hindexed_x -l:libmpich_intel.a
# ROMIO_PATCH	= -Wl,--wrap=ADIOI_Type_create_hindexed_x -l:libmpi.a
ROMIO_PATCH	= -Wl,--wrap=ADIOI_Type_create_hindexed_x

e3sm_io.romio_patch: e3sm_io.o header_io_case_F.o romio_patch.o
	$(MPICC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS) $(ROMIO_PATCH)

romio_patch.o: romio_patch.c
	$(MPICC) $(CFLAGS) -o $@ -c $^

clean:
	rm -f core.* *.o dat2nc e3sm_io e3sm_io.romio_patch
	rm -f f_case_h0_varn.nc f_case_h1_varn.nc
	rm -f f_case_h0_vard.nc f_case_h1_vard.nc
	rm -f g_case_hist_varn.nc

.PHONY: clean

