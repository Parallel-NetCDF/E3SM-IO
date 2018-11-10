
MPICC		= mpicc

PnetCDF_DIR	= $(HOME)/PnetCDF/1.10.0

CFLAGS          = -O2
INCLUDES	= -I$(PnetCDF_DIR)/include -I.
LDFLAGS		= -L$(PnetCDF_DIR)/lib
LIBS		= -lpnetcdf

.c.o:
	$(MPICC) $(CFLAGS) $(INCLUDES) -c $<

noncontig_buf: noncontig_buf.o
	$(MPICC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	rm -f core.* *.o noncontig_buf testfile

.PHONY: clean

