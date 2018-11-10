## Parallel I/O Kernel Case Study For E3SM

This repository contains a case study of parallel I/O kernel from the
[E3SM](https://github.com/E3SM-Project/E3SM) climate simulation model. The
benchmark program, e3sm_io.c, can be used to evaluate
[PnetCDF](https://github.com/Parallel-NetCDF/PnetCDF) library for its
performance on the I/O patterns used by E3SM. One of the most challenging I/O
patterns used in E3SM is the cubed sphere variables whose I/O pattern consists
of a long list of small and noncontiguous requests on every MPI process.

The benchmark program evaluates two implementations of using PnetCDF blocking
vard and nonblocking varn APIs separately. Both APIs can aggregate multiple
requests across variables. In addition to timings and I/O bandwidths, the
benchmark reports the PnetCDF internal memory footprint (high water mark).

The I/O patterns (data decompositions among MPI processes) used in this case
study were captured by the [PIO](https://github.com/NCAR/ParallelIO) library.
A data decomposition file records the data access patterns at the array element
level for each MPI process. The patterns are stored in a text file, referred by
PIO as the `decomposition file`. Note E3SM uses three unique data decomposition
patterns shared by its 381 variables.

* Prepare an combined data decomposition file
  * PIO generates I/O decomposition data files in a text format with file
    extension `.dat`. The three unique decomposition files need first be
    combined and converted into a NetCDF file, as this benchmark program reads
    the combined decomposition NetCDF file in parallel. Using NetCDF file makes
    read faster.
  * Build the conversion utility program, `dat2nc.c`, by running command
    `make dat2nc`.
  * The command to combine and convert the three dat files to a NetCDF file:
    ```
      % ./dat2nc -q -o outputfile.nc -1 decomp_1.dat -2 decomp_2.dat -3 decomp_3.dat
    ```
  * Command-line options of `./dat2nc`:
    ```
      % ./dat2nc -h
      Usage: ./dat2nc [OPTION]... [FILE]...
            -h               Print help
            -q               Quiet mode (reports when fail)
            -l num           max number of characters per line in input file
            -o out_file      name of output netCDF file
            -1 input_file    name of 1st 1D decomposition file
            -2 input_file    name of 2nd 1D decomposition file
            -3 input_file    name of     2D decomposition file
    ```
  * Three small input decomposition files are provide in directory `datasets/`.
    * `piodecomp16tasks16io01dims_ioid_514.dat`  (decomposition along the fastest dimensions)
    * `piodecomp16tasks16io01dims_ioid_516.dat`  (decomposition along the fastest dimensions)
    * `piodecomp16tasks16io02dims_ioid_548.dat`  (decomposition along the fastest two dimensions)
  * The combined NetCDF output file from these 3 decomposition files is
    provided in `datasets/866x72_16p.nc`. Its metadata is shown below.
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

* Compile command:
  * Edit `Makefile` to customize the compiler, compile options, location of
    PnetCDF library, etc.
  * The minimum required PnetCDF version is 1.10.0.
  * Run command `make e3sm_io` to generate the executable program named
    `e3sm_io`.

* Run command:
  * Example run command using `mpiexec` and 8 MPI processes:
    `mpiexec -n 8 ./e3sm_io -q datasets/866x72_16p.nc`
  * The number of MPI processes can be different from the value set in the
    variable `num_procs` in the decomposition NetCDF file. For example, in the
    case of file `866x72_16p.nc`, `num_procs` is 16, which is the number of MPI
    processes originally used to produce the decomposition dat files. When
    using less number of MPI processes to run this benchmark, the workload will
    be divided among all the processes. When using more number of processes,
    those processes with MPI ranks greater than or equal to 16 will have no
    data to write but simply participate the collective I/O subroutines.
  * Command-line options:
    ```
      % ./e3sm_io -h
      Usage: ./e3sm_io [OPTION]... [FILE]...
         [-h] Print help
         [-q] Quiet mode
         [-k] Keep the output files when program exits
         [-o output_dir]: output directory name (default ./)
         input_file: name of input netCDF file describing data decompositions
    ```
* Example outputs on screen
  ```
    % mpiexec -n 8 ./e3sm_io -q -k datasets/866x72_16p.nc

    Total number of MPI processes      = 8
    Input decomposition file           = datasets/866x72_16p.nc
    Output file directory              = .
    Variable dimensions (C order)      = 72 x 866

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
* Output files
  * The above example command uses command-line option `-k` to keep the output
    files (otherwise the default is to delete them when program exits.) Each
    run of `e3sm_io` produces two output netCDF files named `testfile_vard.nc`
    and `testfile_varn.nc`. The contents of two files should be the same. Their
    file header (metadata) obtainable by command `ncdump -h testfile_vard.nc`
    from running the provided decomposition file `866x72_16p.nc` is available
    in [datasets/outputfile_header.txt](datasets/outputfile_header.txt).

## Questions/Comments:
email: wkliao@eecs.northwestern.edu

Copyright (C) 2018, Northwestern University.

See [COPYRIGHT](COPYRIGHT) notice in top-level directory.

