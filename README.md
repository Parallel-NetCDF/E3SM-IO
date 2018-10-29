## Parallel I/O kernel case study for E3SM

This repository contains a case study of parallel I/O kernel from the
[E3SM](https://github.com/E3SM-Project/E3SM) climate simulation model. The
benchmark program, e3sm_io.c, can be used to evaluate the performance of
[PnetCDF](https://github.com/Parallel-NetCDF/PnetCDF) library for handling the
I/O patterns used in E3SM. This study focuses on the most challenging I/O
pattern in E3SM that writes the cubed sphere variables in a long list of
short and noncontiguous requests on every MPI process.

The benchmark program evaluates two implementations of different PnetCDF
APIs: blocking vard and nonblocking varn. Both APIs can aggregate multiple
requests across variables.

The I/O patterns (data decompositions among MPI processes) used in this case
study were captured by the [PIO](https://github.com/NCAR/ParallelIO) library.
A data decomposition file records the data access patterns at the array element
level for each of the MPI processes. The patterns are stored in a text file,
referred by PIO as the `decomposition file`.

* Prepare the I/O decomposition file
  * PIO generates I/O decomposition data files in a text format. It needs to be
    converted to a NetCDF file first, as this benchmark program reads the
    decomposition NetCDF file in parallel.
  * Build the conversion utility program, dat2nc.c, by running command
    `make dat2nc`.
  * Then run command
```
     ./dat2nc -q -o outputfile.nc -1 decomp_1D.dat -2 decomp_2D.dat -3 decomp_3D.dat`.
```
  * Command options of `./dat2nc`:
```
     Usage: ./dat2nc [OPTION]... [FILE]...
            -h               Print help
            -q               Quiet mode (reports when fail)
            -l num           max number of characters per line in input file
            -o out_file      name of output netCDF file
            -1 input_file    name of 1st 1D decomposition file
            -2 input_file    name of 2nd 1D decomposition file
            -3 input_file    name of     2D decomposition file
```
  * Three example input decomposition files are available in directory datasets.
    * piodecomp16tasks16io01dims_ioid_514.dat
    * piodecomp16tasks16io02dims_ioid_548.dat
    * piodecomp16tasks16io01dims_ioid_516.dat
  * Example of an output file generated from the 3 input files is provided in
    datasets/866x72_16p.nc.gz. The metadata of the file is shown below.
```
    % cd ./datasets
    % ncmpidump -h 866x72_16p.nc
    netcdf 866x72_16p {
    // file format: CDF-1
    dimensions:
        num_procs = 16 ;
        D3.max_nreqs = 4032 ;
        D2.max_nreqs = 56 ;
        D1.max_nreqs = 70 ;
    variables:
        int D3.nreqs(num_procs) ;
                D3.nreqs:description = "Number of noncontiguous subarray requests by each MPI process" ;
        int D3.offsets(num_procs, D3.max_nreqs) ;
                D3.offsets:description = "Flattened starting indices of noncontiguous requests. Each row corresponds to requests by an MPI process." ;
        int D2.nreqs(num_procs) ;
                D2.nreqs:description = "Number of noncontiguous subarray requests by each MPI process" ;
        int D2.offsets(num_procs, D2.max_nreqs) ;
                D2.offsets:description = "Flattened starting indices of noncontiguous requests. Each row corresponds to requests by an MPI process." ;
        int D1.nreqs(num_procs) ;
                D1.nreqs:description = "Number of noncontiguous subarray requests by each MPI process" ;
        int D1.offsets(num_procs, D1.max_nreqs) ;
                D1.offsets:description = "Flattened starting indices of noncontiguous requests. Each row corresponds to requests by an MPI process." ;

    // global attributes:
        :dim_len_Y = 72 ;
        :dim_len_X = 866 ;
        :D3.max_nreqs = 4032 ;
        :D3.min_nreqs = 3744 ;
        :D2.max_nreqs = 56 ;
        :D2.min_nreqs = 52 ;
        :D1.max_nreqs = 70 ;
        :D1.min_nreqs = 40 ;
    }
```

* Compile
  * Run command `make e3sm_io` to generate the benchmark program e3sm_io.

* Run
  * example run command:
    `mpiexec -n 8 ./e3sm_io -q datasets/866x72_16p.nc`
  * Below shows the command-line options:
```
    Usage: e3sm_io [OPTION]... [FILE]...
       [-h] Print help
       [-q] Quiet mode
       [-k] Keep the output files when program exits
       [-o output_dir]: output directory name (default ./)
       input_file: name of input netCDF file describing data decompositions
```
* Example outputs on screen
```
    % mpiexec -n 8 e3sm_io -q -k datasets/866x72_16p.nc 

    ---- benchmarking vard API -----------------------
    -----------------------------------------------------------
    MAX heap memory allocated by PnetCDF internally is 2.16 MiB
    Total number of variables          = 408
    Total write amount                 = 16.13 MiB = 0.02 GiB
    Max number of requests             = 325153
    Max Time of open + metadata define = 0.0258 sec
    Max Time of I/O preparing          = 0.0489 sec
    Max Time of ncmpi_put_vard         = 0.6843 sec
    Max Time of close                  = 0.0095 sec
    Max Time of TOTAL                  = 0.7686 sec
    I/O bandwidth                      = 20.9924 MiB/sec

    ---- benchmarking varn API -----------------------
    -----------------------------------------------------------
    MAX heap memory allocated by PnetCDF internally is 36.59 MiB
    Total number of variables          = 408
    Total write amount                 = 16.13 MiB = 0.02 GiB
    Max number of requests             = 325153
    Max Time of open + metadata define = 0.0246 sec
    Max Time of I/O preparing          = 0.0019 sec
    Max Time of ncmpi_iput_varn        = 0.0988 sec
    Max Time of ncmpi_wait_all         = 0.7720 sec
    Max Time of close                  = 0.0159 sec
    Max Time of TOTAL                  = 0.9131 sec
    I/O bandwidth                      = 17.6687 MiB/sec
```
## Questions/Comments:
email: wkliao@eecs.northwestern.edu

Copyright (C) 2018, Northwestern University.

See COPYRIGHT notice in top-level directory.

