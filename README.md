## Parallel I/O Kernel Case Study -- E3SM

This repository contains a case study of parallel I/O kernel from the
[E3SM](https://github.com/E3SM-Project/E3SM) climate simulation model. The
E3SM I/O module makes use of [PIO](https://github.com/NCAR/ParallelIO)
library which is built on top of
[PnetCDF](https://github.com/Parallel-NetCDF/PnetCDF) and
[NetCDF-4](http://www.unidata.ucar.edu/software/netcdf). PnetCDF is a
parallel I/O library for accessing the classic NetCDF files, i.e., CDF-1,
CDF-2, and CDF-5 formats. NetCDF-4 provides parallel I/O capability for
[HDF5](https://www.hdfgroup.org/solutions/hdf5) based NetCDF file format.
The benchmark program in this repository, e3sm_io.c, is designed to evaluate
the E3SM I/O kernel when using the PnetCDF library to perform the I/O task.
In particular, it studies one of the E3SM's most challenging I/O patterns
when the problem domain is represented by cubed sphere grids which produce
long lists of small and noncontiguous requests in each of all MPI processes.

The I/O kernel is implemented by using PnetCDF nonblocking varn APIs, which
can aggregate multiple requests to a variable or across variables. In addition
to timings and I/O bandwidths, the benchmark reports the PnetCDF internal
memory footprint (high water mark).

The I/O patterns (data decompositions among MPI processes) used in this case
study were captured by the [PIO](https://github.com/NCAR/ParallelIO) library.
A data decomposition file records the data access patterns at the array element
level for each MPI process. The access offsets are stored in a text file,
referred to by PIO as the `decomposition file. This benchmark currently studies
two cases from E3SM, namely F and G cases. The F case uses three unique data
decomposition patterns shared by 388 2D and 3D variables (2 sharing
Decomposition 1, 323 sharing Decomposition 2, and 63 sharing Decomposition 3).
The G case uses 6 data decompositions shared by 52 variables (6 sharing
Decomposition 1, 2 sharing Decomposition 2, 25 sharing Decomposition 3, 2
sharing Decomposition 4, 2 sharing Decomposition 5, and 4 sharing Decomposition
6).

### Software Requirements
* [PnetCDF 1.12.2](https://parallel-netcdf.github.io/Release/pnetcdf-1.12.2.tar.gz)
* MPI C and C++ compilers
  + The plugin uses the constant initializer; a C++ compiler supporting std 11 is required
* Autotools utility
  + autoconf 2.69
  + automake 1.16.1
  + libtoolize 2.4.6
  + m4 1.4.18
* (Optional) [HDF5 1.12.0](https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_1_12_0/source/hdf5-1.12.0.tar.gz)
  + Parallel I/O support (--enable-parallel) is required
* (Optional) [Log-based VOL](https://github.com/DataLib-ECP/vol-log-based.git)
  + Experimental
  + Must have HDF5
* (Optional) [ADIOS2 2.7.1](https://github.com/ornladios/ADIOS2.git)
  + Parallel I/O support (-DADIOS2_USE_MPI=ON) is required

### Building Steps
* Build PnetCDF
  + Download and extract the PnetCDF source code
  + Configure PnetCDF with mpi C compiler
  + Run make install
  + Example build commands are given below. This example will install
    the PnetCDF library under the folder `${HOME}/PnetCDF/1.12.2`.
    ```
    % wget https://parallel-netcdf.github.io/Release/pnetcdf-1.12.2.tar.gz
    % tar -zxf pnetcdf-1.12.2.tar.gz
    % cd pnetcdf-1.12.2
    % ./configure --prefix=${HOME}/PnetCDF/1.12.2 CC=mpicc
    % make -j 16 install
    ```
    The PnetCDF library is now installed under the folder `${HOME}/PnetCDF/1.12.2.`
* (Optional) Build HDF5 with parallel I/O support
  + Download and extract the HDF5 source code
  + Configure HDF5 with parallel I/O enabled
  + Run make install
  + Example build commands are given below. This example will install
    the HD5 library under the folder `${HOME}/HDF5/1.12.0`.
    ```
    % wget https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_1_12_0/source/hdf5-1.12.0.tar.gz
    % tar -zxf hdf5-1.12.0.tar.gz 
    % cd hdf5-1.12.0
    % ./configure --prefix=${HOME}/HDF5/1.12.0 --enable-parallel CC=mpicc
    % make -j 16 install
    ```
    The HDF5 library is now installed under the folder `${HOME}/HDF5/1.12.0.`
* (Optional) Build log-based VOL plugin.
  + Clone the source code from the log-based VOL repository
  + Run command autoreconf -i
  + Configure log-based VOL 
    + Shared library is required to enable log-based VOL by environment variables
    + Compile with zlib library to enable metadata compression
  + Example commands are given below.
    ```
    % git clone https://github.com/DataLib-ECP/vol-log-based.git
    % cd log_io_vol
    % autoreconf -i
    % ./configure --prefix=${HOME}/Log_IO_VOL --with-hdf5=${HOME}/HDF5/1.12.0 --enable-shared --enable-zlib
    % make -j 16 install
    ```
    The VOL plugin library is now installed under the folder `${HOME}/Log_IO_VOL.`
* Build the E3SM-I/O benchmark
  + Clone this E3SM-I/O benchmark repository
  + Run command autoreconf -i
  + Configure the E3SM-I/O benchmark with mpi compilers
    + Add HDF5 installation path (--with-hdf5=/path/to/implementation) to enable HDF5 API support
      + E3SM-I/O benchmark only supports sequential run using HDF5 native VOL due to some HDF5 API limitations.
        + Parallel I/O is supported when using the log-based VOL.
    + Add log-based VOL installation path (--with-logvol=/path/to/implementation) to enable HDF5 API support
    + Add ADIOS2 installation path (--with-adios2=/path/to/implementation) to enable HDF5 API support
  + Run make install
  + Example commands are given below.
    ```
    % git clone https://github.com/Parallel-NetCDF/E3SM-IO.git
    % cd E3SM-IO
    % autoreconf -i
    % ./configure --with-pnetcdf=${HOME}/PnetCDF/1.12.2 --with-hdf5=${HOME}/HDF5/1.12.0 --with-logvol=${HOME}/Log_IO_VOL --with-adios2=${HOME}/ADIOS2 CC=mpicc CXX=mpicxx
    % make -j 16 install
    ```

### Prepare the data decomposition file in NetCDF file format
* For the F case, the three data decomposition files generated by the PIO library
  are in text format with file extension name `.dat`. The decomposition files
  must first be combined and converted into a NetCDF file to be read in
  parallel as the input file to this benchmark program. Similarly, for the G
  case, the six decomposition files need to be converted first.
* A utility program, `dat2nc.c`, is included to convert the text files. To
  build this utility program, run command
  ```
    % make dat2nc
  ```
* The command to combine the three .dat files to a NetCDF file for the F case
  as an example, is:
  ```
    % ./dat2nc -o outputfile.nc -1 decomp_1.dat -2 decomp_2.dat -3 decomp_3.dat
  ```
* Command-line options of `./dat2nc`:
  ```
    % ./dat2nc -h
    Usage: ./dat2nc [OPTION]... [FILE]...
          -h               Print help
          -v               Verbose mode
          -l num           Max number of characters per line in input file
          -o out_file      Name of output netCDF file
          -1 input_file    name of 1st decomposition file
          -2 input_file    name of 2nd decomposition file
          -3 input_file    name of 3rd decomposition file
          -4 input_file    name of 4th decomposition file
          -5 input_file    name of 5th decomposition file
          -6 input_file    name of 6th decomposition file
  ```
* Three small input decomposition files for an F case are provide in
  directory `datasets/`.
  * `piodecomp16tasks16io01dims_ioid_514.dat`  (decomposition along the fastest dimensions)
  * `piodecomp16tasks16io01dims_ioid_516.dat`  (decomposition along the fastest dimensions)
  * `piodecomp16tasks16io02dims_ioid_548.dat`  (decomposition along the fastest two dimensions)
* The NetCDF file converted from these 3 decomposition .dat files is provided
  in folder `datasets` named `f_case_866x72_16p.nc`. Its metadata is shown
  below.
  ```
    % cd ./datasets
    % ncmpidump -h f_case_866x72_16p.nc
    netcdf f_case_866x72_16p {
    // file format: CDF-1
    dimensions:
        num_decomp = 3 ;
        decomp_nprocs = 16 ;
        D1.total_nreqs = 47 ;
        D2.total_nreqs = 407 ;
        D3.total_nreqs = 29304 ;
    variables:
        int D1.nreqs(decomp_nprocs) ;
            D1.nreqs:description = "Number of noncontiguous requests per process" ;
        int D1.offsets(D1.total_nreqs) ;
            D1.offsets:description = "Flattened starting indices of noncontiguous requests" ;
        int D1.lengths(D1.total_nreqs) ;
            D1.lengths:description = "Lengths of noncontiguous requests" ;
            D1.lengths:max = 36 ;
            D1.lengths:min = 9 ;
        int D2.nreqs(decomp_nprocs) ;
            D2.nreqs:description = "Number of noncontiguous requests per process" ;
        int D2.offsets(D2.total_nreqs) ;
            D2.offsets:description = "Flattened starting indices of noncontiguous requests" ;
        int D2.lengths(D2.total_nreqs) ;
            D2.lengths:description = "Lengths of noncontiguous requests" ;
            D2.lengths:max = 4 ;
            D2.lengths:min = 1 ;
        int D3.nreqs(decomp_nprocs) ;
            D3.nreqs:description = "Number of noncontiguous requests per process" ;
        int D3.offsets(D3.total_nreqs) ;
            D3.offsets:description = "Flattened starting indices of noncontiguous requests" ;
        int D3.lengths(D3.total_nreqs) ;
            D3.lengths:description = "Lengths of noncontiguous requests" ;
            D3.lengths:max = 4 ;
            D3.lengths:min = 1 ;

    // global attributes:
        :command_line = "./dat2nc -o f_case_866x72_16p.nc -1 datasets/piodecomp16tasks16io01dims_ioid_514.dat -2 datasets/piodecomp16tasks16io01dims_ioid_516.dat -3 datasets/piodecomp16tasks16io02dims_ioid_548.dat " ;
        :D1.ndims = 1 ;
        :D1.dims = 866 ;
        :D1.max_nreqs = 4 ;
        :D1.min_nreqs = 2 ;
        :D2.ndims = 1 ;
        :D2.dims = 866 ;
        :D2.max_nreqs = 39 ;
        :D2.min_nreqs = 13 ;
        :D3.ndims = 2 ;
        :D3.dims = 72, 866 ;
        :D3.max_nreqs = 2808 ;
        :D3.min_nreqs = 936 ;
    }
  ```
* A NetCDF file containing 6 decompositions from a small G case is also
  included in folder `datasets` named `g_case_cmpaso_16p.nc`.

### Run command:
* An example run command using `mpiexec` and 16 MPI processes is:
  ```
    % mpiexec -n 16 ./e3sm_io -c datasets/f_case_866x72_16p.nc
  ```
* The number of MPI processes used to run this benchmark can be different
  from the value of the variable `decomp_nprocs` stored in the decomposition
  NetCDF file. For example, in file `f_case_866x72_16p.nc`, `decomp_nprocs`
  is 16, the number of MPI processes originally used to produce the decomposition .dat files. When running this benchmark using less number of
  MPI processes, the I/O workload will be divided among all the allocated MPI
  processes. When using more processes than `decomp_nprocs`, those processes
  with MPI ranks greater than or equal to `decomp_nprocs` will have no data
  to write but still participate the collective I/O in the benchmark.
* Command-line options:
  ```
    % ./e3sm_io -h
    Usage: src/e3sm_io [OPTION]... FILE
       [-h] Print help
       [-v] Verbose mode
       [-k] Keep the output files when program exits
       [-d] Run test that uses PnetCDF vard API
       [-n] Run test that uses PnetCDF varn API
       [-m] Run test using noncontiguous write buffer
       [-t] Write 2D variables followed by 3D variables
       [-R] Test reading performance
       [-W] Test writing performance
       [-f num] File number to run in F case (-1 (both) (default), 0, 1)
       [-r num] Number of records (default 1)
       [-s num] Stride between IO tasks (default 1)
       [-o output_dir] Output directory name (default ./)
       [-i target_dir] Path to directory containing the input files
       [-a api] Underlying API to test (pnc (default), hdf5, adios2)
       [-c chunk_size] Use chunked storage layout with chunk_size (0 (no chunking) (default))
       [-z filter] Apply the filter if supported by the underlying API (none (default), deflate)
       FILE: Name of input netCDF file describing data decompositions
  ```
* An example batch script file for running a job on Cori @NERSC with 8 KNL
  nodes, 64 MPI processes per node, is provided in `./slurm.knl`.
* A median-size decomposition file `datasets/f_case_48602x72_512p.nc` contains
  the I/O pattern from a bigger problem size used in an E3SM experiment ran on
  512 MPI processes.
* Two large decomposition files `f_case_72x777602_21632p.nc.gz` (60 MB) and
  `g_case_11135652x80_9600p.nc.gz` (126 MB) from high-resolution simulations of
  F and G cases running on 21632 and 9600 MPI processes respectively are
  available upon request.

### Environment variables
* E3SM_IO_HDF5_ENABLE_LOGVOL 
  + 1: Use Log I/O VOL if available
  + 0: Use the native VOL (default)
* E3SM_IO_HDF5_USE_LOGVOL_WRITEN
  + 1: Use the H5Dwrite_N API in Log I/O VOL
  + 0: Use the HDF5 driver varn implementation (default)
  + Only effective when E3SM_IO_HDF5_ENABLE_LOGVOL is 1
* E3SM_IO_HDF5_MERGE_VARN
  + 1: Merge varn hyper-slabs into one dataspace selection
  + 0: Call H5Dwrite per hyper-slab (default)

### example outputs shown on screen
```
  % mpiexec -n 512 ./e3sm_io -k -r 3 -o $SCRATCH/FS_1M_64 -c datasets/f_case_48602x72_512p.nc

  Total number of MPI processes      = 512
  Input decomposition file           = datasets/f_case_48602x72_512p.nc
  Output file directory              = $SCRATCH/FS_1M_64
  Variable dimensions (C order)      = 72 x 48602
  Write number of records (time dim) = 3
  Using noncontiguous write buffer   = no

  ==== benchmarking vard API ================================
  Variable written order: same as variables are defined

  History output file                = f_case_h0_vard.nc
  No. variables use no decomposition =  26
  No. variables use decomposition D1 =   2
  No. variables use decomposition D2 = 323
  No. variables use decomposition D3 =  63
  Total number of variables          = 414
  MAX heap memory allocated by PnetCDF internally is 2.22 MiB
  Total write amount                 = 2699.36 MiB = 2.64 GiB
  Max number of requests             = 310598
  Max Time of open + metadata define = 0.0533 sec
  Max Time of I/O preparing          = 0.1156 sec
  Max Time of ncmpi_put_vard         = 5.4311 sec
  Max Time of close                  = 0.0306 sec
  Max Time of TOTAL                  = 5.6385 sec
  I/O bandwidth (open-to-close)      = 478.7341 MiB/sec
  I/O bandwidth (write-only)         = 496.9981 MiB/sec
  -----------------------------------------------------------
  History output file                = f_case_h1_vard.nc
  No. variables use no decomposition =  26
  No. variables use decomposition D1 =   2
  No. variables use decomposition D2 =  22
  No. variables use decomposition D3 =   1
  Total number of variables          =  51
  MAX heap memory allocated by PnetCDF internally is 2.22 MiB
  Total write amount                 = 52.30 MiB = 0.05 GiB
  Max number of requests             = 6022
  Max Time of open + metadata define = 0.0338 sec
  Max Time of I/O preparing          = 0.0014 sec
  Max Time of ncmpi_put_vard         = 0.2489 sec
  Max Time of close                  = 0.0055 sec
  Max Time of TOTAL                  = 0.2907 sec
  I/O bandwidth (open-to-close)      = 179.8902 MiB/sec
  I/O bandwidth (write-only)         = 210.1002 MiB/sec
  -----------------------------------------------------------

  ==== benchmarking varn API ================================
  Variable written order: same as variables are defined

  History output file                = f_case_h0_varn.nc
  No. variables use no decomposition =  26
  No. variables use decomposition D1 =   2
  No. variables use decomposition D2 = 323
  No. variables use decomposition D3 =  63
  Total number of variables          = 414
  MAX heap memory allocated by PnetCDF internally is 35.07 MiB
  Total write amount                 = 2699.36 MiB = 2.64 GiB
  Max number of requests             = 310464
  Max Time of open + metadata define = 0.0635 sec
  Max Time of I/O preparing          = 0.0018 sec
  Max Time of ncmpi_iput_varn        = 0.2468 sec
  Max Time of ncmpi_wait_all         = 5.8602 sec
  Max Time of close                  = 0.0190 sec
  Max Time of TOTAL                  = 6.2001 sec
  I/O bandwidth (open-to-close)      = 435.3753 MiB/sec
  I/O bandwidth (write-only)         = 460.6144 MiB/sec
  -----------------------------------------------------------
  History output file                = f_case_h1_varn.nc
  No. variables use no decomposition =  26
  No. variables use decomposition D1 =   2
  No. variables use decomposition D2 =  22
  No. variables use decomposition D3 =   1
  Total number of variables          =  51
  MAX heap memory allocated by PnetCDF internally is 35.07 MiB
  Total write amount                 = 52.30 MiB = 0.05 GiB
  Max number of requests             = 5888
  Max Time of open + metadata define = 0.0370 sec
  Max Time of I/O preparing          = 0.0005 sec
  Max Time of ncmpi_iput_varn        = 0.0048 sec
  Max Time of ncmpi_wait_all         = 0.2423 sec
  Max Time of close                  = 0.0058 sec
  Max Time of TOTAL                  = 0.2925 sec
  I/O bandwidth (open-to-close)      = 178.7747 MiB/sec
  I/O bandwidth (write-only)         = 215.7512 MiB/sec
  -----------------------------------------------------------
```
### Output files
* The above example command uses command-line option `-k` to keep the output
  files (otherwise the default is to delete them when program exits.) For the
  F case, each run of `e3sm_io` produces four output netCDF files named
  `f_case_h0_vard.nc`, `f_case_h1_vard.nc`, `f_case_h0_varn.nc`, and
  `f_case_h1_varn.nc`. Their file header (metadata) obtainable by command
  `ncdump -h` from running the provided decomposition file
  `f_case_866x72_16p.nc` is available in
  [datasets/f_case_h0.txt](datasets/f_case_h0.txt), and
  [datasets/f_case_h1.txt](datasets/f_case_h1.txt).
* For the G case, there is one output file, namely `g_case_hist_varn.nc`, and
  its file header can be found in
  [datasets/g_case_hist.txt](datasets/g_case_hist.txt).

### Current build status
* [Travis CI ![Build Status](https://travis-ci.org/Parallel-NetCDF/E3SM-IO.svg?branch=master)](https://travis-ci.org/Parallel-NetCDF/E3SM-IO)

### Developers
* Wei-keng Liao <wkliao@northwestern.edu>
* Kai-yuan Hou <kai-yuanhou2020@u.northwestern.edu>

Copyright (C) 2021, Northwestern University.

See [COPYRIGHT](COPYRIGHT) notice in top-level directory.
