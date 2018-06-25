This repository contains the benchmark program for evaluating the performance
of [PnetCDF](https://github.com/Parallel-NetCDF/PnetCDF) library using the I/O
patterns from the [E3SM](https://github.com/E3SM-Project/E3SM).

* Prepare the I/O decomposition file
  * PIO generates I/O decomposition data file in a text format. It needs to be
    converted to a NetCDF file first, as this benchmark program reads the
    decomposition NetCDF file in parallel.
  * Build the conversion utility program, dat2nc.c, by running command
    `make dat2nc`.
  * Then run command `./dat2nc inputfile.dat -o outputfile.nc`.

* Compile
  * Run command `make e3sm_io` to generate the benchmark program e3sm_io.

* Run
  * example run command:
    `mpiexec -n 8 ./e3sm_io -q datasets/48602x72_64p.nc -n 10`
  * Below shows the command-line options:
```
     Usage: e3sm_io [-h] [-q] [-n nvars] [-o output_file] input_file
            [-h] Print help
            [-q] Quiet mode
            [-n nvars]: number of variables (default 1)
            [-o output_file]: output file name (default ./testfile.nc)
            input_file: intput netCDF file name
```
* Example outputs on screen
```
     ---- benchmarking vard API -----------------------
     -----------------------------------------------------------
     MAX heap memory allocated by PnetCDF internally is 0.27 MiB
     Total number of variables          = 10
     Total write amount                 = 133.49 MiB = 0.13 GiB
     Max number of requests             = 22144
     Max Time of open + metadata define = 0.0005 sec
     Max Time of I/O preparing          = 0.0354 sec
     Max Time of ncmpi_put_vard         = 0.4521 sec
     Max Time of close                  = 0.0037 sec
     Time of TOTAL                      = 0.5385 sec
     I/O bandwidth                      = 247.9046 MiB/sec

     ---- benchmarking varn API -----------------------
     -----------------------------------------------------------
     MAX heap memory allocated by PnetCDF internally is 20.47 MiB
     Total number of variables          = 10
     Total write amount                 = 133.49 MiB = 0.13 GiB
     Max number of requests             = 22144
     Max Time of open + metadata define = 0.0004 sec
     Max Time of I/O preparing          = 0.0211 sec
     Max Time of ncmpi_iput_varn        = 0.0497 sec
     Max Time of ncmpi_wait_all         = 0.4593 sec
     Max Time of close                  = 0.0079 sec
     Max Time of TOTAL                  = 0.5385 sec
     I/O bandwidth                      = 247.9046 MiB/sec

     ---- benchmarking vara API -----------------------
     -----------------------------------------------------------
     MAX heap memory allocated by PnetCDF internally is 45.78 MiB
     Total number of variables          = 10
     Total write amount                 = 133.49 MiB = 0.13 GiB
     Max number of requests             = 22144
     Max Time of open + metadata define = 0.0004 sec
     Max Time of I/O preparing          = 0.0210 sec
     Max Time of ncmpi_iput_vara        = 0.9343 sec
     Max Time of ncmpi_wait_all         = 0.6547 sec
     Max Time of close                  = 0.0141 sec
     Max Time of TOTAL                  = 1.6247 sec
     I/O bandwidth                      = 82.1640 MiB/sec
     Max Time of TOTAL                  = 0.4916 sec
     I/O bandwidth                      = 271.5230 MiB/sec
```
