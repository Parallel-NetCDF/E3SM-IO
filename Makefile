
MPICC		= mpicc
CFLAGS          = -O2

PnetCDF_DIR	= $(HOME)/PnetCDF/1.10.0
INCLUDES	= -I$(PnetCDF_DIR)/include -I.
LDFLAGS		= -L$(PnetCDF_DIR)/lib
LIBS		= -lpnetcdf

.c.o:
	$(MPICC) $(CFLAGS) $(INCLUDES) -c $<

all:

e3sm_io: e3sm_io.o
	$(MPICC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	rm -f core.* *.o dat2nc e3sm_io
