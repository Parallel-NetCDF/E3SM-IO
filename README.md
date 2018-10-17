## Parallel I/O kernel case study for E3SM

This repository contains a case study of parallel I/O kernel from the
[E3SM](https://github.com/E3SM-Project/E3SM) climate simulation model. The
benchmark program, e3sm_io.c, can be used to evaluate the performance of
[PnetCDF](https://github.com/Parallel-NetCDF/PnetCDF) library for handling the
I/O patterns used in E3SM. This study focuses on the most challenging I/O
pattern in E3SM that writes the cubed sphere variables in a long list of
short and noncontiguous requests on every MPI process.

The benchmark program evaluates three implementations of different PnetCDF
APIs: blocking vard, nonblocking varn, and nonblocking vara. All these three
APIs can aggregate multiple requests across more than one variable.

The I/O patterns (data decompositions among MPI processes) used in this case
study were captured by the [PIO](https://github.com/NCAR/ParallelIO) library.
A data decomposition file records the data access patterns at the array element
level for each of the MPI processes. The patterns are stored in a text file,
referred by PIO as the `decomposition file`.

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
    `mpiexec -n 8 ./e3sm_io -q -n 10 datasets/48602x72_64p.nc`
  * Below shows the command-line options:
```
     Usage: ./e3sm_io [-h] [-q] [-k] [-n nvars] [-o output_dir] input_file
            [-h] Print help
            [-q] Quiet mode
            [-k] Keep the output files when program exits
            [-n nvars]: number of variables (default 1)
            [-o output_dir]: output directory name (default ./)
            input_file: name of input netCDF file describing the decomposition
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
     Max Time of TOTAL                  = 0.4916 sec
     I/O bandwidth                      = 271.5230 MiB/sec

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
```
