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
  + The plugin uses the constant initializer; a C++ compiler supporting std 11
    is required
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
* (Optional) [ADIOS2 2.7.1](https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.7.1.tar.gz)
  + Parallel I/O support (-DADIOS2_USE_MPI=ON) is required

### Building Steps
* Build PnetCDF
  + Download and extract the PnetCDF source code
  + Configure PnetCDF with mpi C compiler
  + Run make install
  + Example build commands are given below. This example will install
    the PnetCDF library under folder `${HOME}/PnetCDF/1.12.2`.
    ```
    % wget https://parallel-netcdf.github.io/Release/pnetcdf-1.12.2.tar.gz
    % tar -zxf pnetcdf-1.12.2.tar.gz
    % cd pnetcdf-1.12.2
    % ./configure --prefix=${HOME}/PnetCDF/1.12.2 CC=mpicc
    % make -j 16 install
    ```
    The PnetCDF library is now installed under folder `${HOME}/PnetCDF/1.12.2.`
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
    The HDF5 library is now installed under folder `${HOME}/HDF5/1.12.0.`
* (Optional) Build log-based VOL plugin.
  + Clone the source code from the log-based VOL repository
  + Run command autoreconf -i
  + Configure log-based VOL 
    + Shared library is required to enable log-based VOL by environment
      variables
    + Compile with zlib library to enable metadata compression
  + Example commands are given below.
    ```
    % git clone https://github.com/DataLib-ECP/vol-log-based.git
    % cd log_io_vol
    % autoreconf -i
    % ./configure --prefix=${HOME}/Log_IO_VOL --with-hdf5=${HOME}/HDF5/1.12.0 --enable-shared --enable-zlib
    % make -j 16 install
    ```
    The VOL plugin library is now installed under folder `${HOME}/Log_IO_VOL.`
* (Optional) Build ADIOS2 with parallel I/O support
  + Download and extract the ADIOS2 source code
  + Configure ADIOS2 with MPI support enabled
  + Run make install
  + Example build commands are given below. This example will install
    the ADIOS2 library under the folder `${HOME}/ADIOS2/2.7.1`.
    ```
    % wget https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.7.1.tar.gz
    % tar -zxf v2.7.1.tar.gz
    % mkdir ADIOS2_BUILD
    % cd ADIOS2_BUILD
    % cmake -DCMAKE_INSTALL_PREFIX=${HOME}/ADIOS2/2.7.1 -DADIOS2_USE_MPI=ON ../ADIOS2-2.7.1
    % make -j 16
    % make -j 16 install
    ```
    The ADIOS2 library is now installed under folder `${HOME}/ADIOS2/2.7.1.`
* Build the E3SM-I/O benchmark
  + Clone this E3SM-I/O benchmark repository
  + Run command autoreconf -i
  + Configure the E3SM-I/O benchmark with mpi compilers
    + Add HDF5 installation path (--with-hdf5=/path/to/implementation) to
      enable HDF5 API support. Using native HDF5 is supported only when running
      the benchmark with command-line option '-a hdf5_ra'.
    + Add HDF5 installation path (--with-hdf5=/path/to/implementation) that
      contains the implementation of multi-dataset APIs. This is required when
      running the benchmark with command-line option '-a hdf5_md'.
    + Add log-based VOL installation path
      (--with-logvol=/path/to/implementation). This is required when running
      the benchmark with command-line option '-a hdf5_log'.
    + Add ADIOS2 installation path (--with-adios2=/path/to/implementation) to
      enable ADIOS API support. This is required when running the benchmark
      with command-line option '-a adios'.
  + Run 'make install'
  + Example commands are given below.
    ```
    % git clone https://github.com/Parallel-NetCDF/E3SM-IO.git
    % cd E3SM-IO
    % autoreconf -i
    % ./configure --with-pnetcdf=${HOME}/PnetCDF/1.12.2 \
                  --with-hdf5=${HOME}/HDF5/1.12.0 \
                  --with-logvol=${HOME}/Log_IO_VOL \
                  --with-adios2=${HOME}/ADIOS2 \
                  CC=mpicc CXX=mpicxx
    % make -j 16 install
    ```

### Prepare the data decomposition file in NetCDF file format
* For the F case, there are three data decomposition files generated by the PIO
  library in text format with file extension name `.dat`. The decomposition
  files must first be combined and converted into a NetCDF file to be read in
  parallel as the input file to this benchmark program. Similarly, for the G
  case, there are six decomposition files that need to be converted first.
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
          -o out_file      Name of output NetCDF file
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
* Example run commands using `mpiexec` and 16 MPI processes:
  + Run write test with the default settings.
    ```
      % mpiexec -n 16 ./e3sm_io -o ${Output_File_Path} datasets/f_case_866x72_16p.nc
    ```
  + Using ADIOS2 APIs (if enabled) to run write test.
    ```
      % mpiexec -n 16 ./e3sm_io -a adios -o ${Output_File_Path} datasets/f_case_866x72_16p.nc
    ```
  + Run read the test using HDF5 API with rearranger in the E3SM benchmark.
    ```
      % ln -s ${Path_to_E3SM_F_Case_H0_File} ${Input_File_Path}_h0.h5
      % ln -s ${Path_to_E3SM_F_Case_H1_File} ${Input_File_Path}_h1.h5
      % mpiexec -n 16 ./e3sm_io -a hdf5_ra -i ${Input_File_Path}.h5 datasets/f_case_866x72_16p.nc
  + Run read the test with PnetCDF and use the data read to run write test with ADIOS2.
    ```
      % ln -s ${Path_to_E3SM_F_Case_H0_File} ${Input_File_Path}_h0.nc
      % ln -s ${Path_to_E3SM_F_Case_H1_File} ${Input_File_Path}_h1.nc
      % mpiexec -n 16 ./e3sm_io -a adios -i ${Input_File_Path}.nc -o ${Output_File_Path} datasets/f_case_866x72_16p.nc
    ```
* The number of MPI processes used to run this benchmark can be different from
  the value of variable `decomp_nprocs` stored in the decomposition NetCDF
  file. For example, in file `f_case_866x72_16p.nc`, `decomp_nprocs` is 16, the
  number of MPI processes originally used to generate the decomposition .dat
  files. When running this benchmark using less number of MPI processes, the
  I/O workload will be divided among all the allocated MPI processes. When
  using more processes than `decomp_nprocs`, the processes with MPI ranks
  greater than or equal to `decomp_nprocs` will have no data to write but still
  participate the collective I/O in the benchmark.
* Command-line options:
  ```
    % ./e3sm_io -h
    Usage: ./src/e3sm_io [OPTION]... FILE
       [-h] Print this help message
       [-v] Verbose mode
       [-k] Keep the output files when program exits
       [-d] Run test that uses PnetCDF vard API
       [-n] Run test that uses PnetCDF varn API
       [-m] Run test using noncontiguous write buffer
       [-t] Write 2D variables followed by 3D variables
       [-f num] Set history output files h0 and/or h1: 0 for h0 only, 1 for h1
                only, -1 for both (default: -1)
       [-r num] Number of records/time steps (default: 1)
       [-s num] MPI rank stride for selecting processes to perform I/O tasks
                (default: 1)
       [-g num] Number of subfiles, used in blob I/O only (default: 1)
       [-i path] Enable read performance evaluation and set the input file
                 (folder) path
       [-o path] Enable write performance evaluation and set the output file
                 (folder) path
       [-a api]  I/O library name to perform write operation
           pnetcdf:   PnetCDF library (default)
           hdf5_ra:   HDF5 library with request rearranger on top of it
           hdf5_log:  HDF5 library with Log-based VOL
           hdf5_md:   HDF5 library with multi-dataset APIs
           adios:     ADIOS2 library using BP3 format
       [-x strategy] I/O strategy to write
           canonical: Store E3SM variables in the canonical layout (default)
           log:       Store E3SM variables as is in log-based storage layout
           blob:      Write data is stored in a contiguous block (blob),
                      ignoring variable's canonical order
       [-c size] Data chunk size to be used when compression is enabled.
                 (default 0, i.e. no chunking)
       [-z filter] Enable data compression in write and use the supplied the
                 filter name (default: none)
       FILE: Name of input NetCDF file describing data decompositions
  ```
* For option '-i', the API to be used in the read test is determined by the
  input file format, which is detected internally. The API selected in option
  '-a' is for write only, except for when the input file is an HDF5 and option
  '-a' is either hdf5_ra or hdf5_md, in which case the HDF5 rearranger API or
  multi-dataset API is used to read the input file.
* The F write case will create two history files. The supplied file name in
  option '-o' will be used to create two new files name by inserting/appending
  strings '_h0' and '_h1' to indicate the two history files.  If the input path
  contains file extension '.nc', '_h0' and '_h1' will be inserted before the
  file extension. Otherwise, they will be appended at the end.  See examples in
  "Output files" section below.
* If both read and write options are enabled, i.e. both '-i' and '-o' are set.
  The benchmark will read the input file first and use the data read to write
  to the output file. If read option is not set, the benchmark will write
  random data.
* Current supported APIs (option '-a') and I/O strategies (option '-x')
  + **pnetcdf + canonical**
    * A single NetCDF file in the classic CDF5 format will be created. All data
      objects stored in the file are in the canonical order and understandable
      by NetCDF and its third-party software.
    * If the output file systems allow users to customize the file striping
      configuration, such as Lustre, users are recommended to write to a folder
      with a high file striping count to obtain a good I/O performance.
  + **pnetcdf + blob**
    * Multiple subfiles in NetCDF format will be created. The files conform
      with NetCDF file format specification.
    * There will be one subfile per compute node used.
    * Input file name provided in option '-i' will be used as a base to create
      the subfile name. The subfile names will have the numerical IDs as the
      suffix.
    * Because all variables are stored in a blob fashion in the files, the
      subfiles altogether can only be understood by the conversion utility
      tool, [utils/pnetcdf_blob_replay.c](utils/pnetcdf_blob_replay.c), which
      is to be run off-line to convert the subfiles into a single regular
      NetCDF file in CDF5 format.
  + **hdf5_ra + canonical**
    * hdf5_ra re-arranges the write requests among all MPI processes into less
      but large contiguous requests before calling HDF5 to write the data. This
      is essentially the same as the BOX data rearrangment implemented in
      Scorpio.
    * The output file is a regular HDF5, which is understandable by regular
      HDF5 and its third-party software.
  + **hdf5_md + canonical**
    * hdf5_md reads/writes data using the multi-dataset APIs, which is a new
      HDF5 feature and currently under development. The APIs allow users to
      read/write multiple requests in a single API call and thus achieve a
      better I/O performance.
    * This option requires to configure this benchmark with the
      [develop branch](https://bitbucket.hdfgroup.org/projects/HDFFV/repos/hdf5/browse)
      of HDF5 that implements the multi-dataset APIs.
    * The output file is a regular HDF5, which is understandable by regular
      HDF5 and its third-party software.
  + **hdf5_log + log**
    * hdf5_log reads/writes data using the log-based VOL, which stores data in
      a log layout, rather than a canonical layout. The output file is a valid
      HDF5 file, but requires the log-based VOL to read and understand the data
      structures.
    * The output file is a single HDF5 file.
  + **adios + blob**
    * Multiple subfiles in BP format will be created.
    * The number of subfile is determined by command-line option '-g'.
    * The input file name provided in option '-i' will be used as a base to
      create folder names that store the subfiles. The folder names will have
      suffix ".bp.dir" appended. Each subfile name in its folder will have
      ".bp" and a numerical ID appended.
    * Because all variables are stored in a blob fashion in the files, the
      subfiles can only be understood by the Scorpio's conversion utility tool,
      [adios2pio-nm](https://github.com/E3SM-Project/scorpio/tree/master/tools/adios2pio-nm),
      which is to be run off-line to convert the subfiles into a single regular
      NetCDF file.

### Example files
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
* E3SM_IO_HDF5_USE_LOGVOL_WRITEN
  + 1: Use the H5Dwrite_N API in Log I/O VOL
  + 0: Use the HDF5 driver varn implementation (default)
  + Only effective when E3SM_IO_HDF5_ENABLE_LOGVOL is 1
* E3SM_IO_HDF5_MERGE_VARN
  + 1: Merge varn hyper-slabs into one dataspace selection
  + 0: Call H5Dwrite per hyper-slab (default)

### example outputs shown on screen
```
  % mpiexec -n 512 ./e3sm_io -n -k -r 3 -o $SCRATCH/FS_1M_64/can_F_out.nc datasets/f_case_48602x72_512p.nc

  Total number of MPI processes      = 512
  Input decomposition file           = datasets/f_case_48602x72_512p.nc
  Output file/directory              = $SCRATCH/FS_1M_64/can_F_out.nc
  Variable dimensions (C order)      = 72 x 48602
  Write number of records (time dim) = 3
  Using noncontiguous write buffer   = no

  ==== benchmarking varn API ================================
  Variable written order: same as variables are defined

  History output file                = $SCRATCH/FS_1M_64/can_F_out_h0.nc
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
  History output file                = $SCRATCH/FS_1M_64/can_F_out_h1.nc
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
  files (otherwise the default is to delete them when the program exits.) For
  the F case, each run of `e3sm_io` produces two history output files whose
  names are created by inserting 'h0', and 'h1' to user-supplied file name.
  The header of F case files from running the provided decomposition file
  `f_case_866x72_16p.nc` using PnetCDF obtainable by command `ncdump -h` is
  available in [datasets/f_case_h0.txt](datasets/f_case_h0.txt), and
  [datasets/f_case_h1.txt](datasets/f_case_h1.txt).
* The G case only creates is one output file. The header of G case file running
  the provided decomposition file `g_case_cmpaso_16p.nc` using PnetCDF can be
  found in [datasets/g_case_hist.txt](datasets/g_case_hist.txt).
* The ADIOS2 API option automatically appends '.bp.dir' extension to the
  user-provided input path and creates two folders for F case (one for G case.)
  + When using ADIOS2, the names of output subfiles will be appended with
    file extension `.bp.dir`.

### Current build status
* [Travis CI ![Build Status](https://travis-ci.org/Parallel-NetCDF/E3SM-IO.svg?branch=master)](https://travis-ci.org/Parallel-NetCDF/E3SM-IO)

### Developers
* Wei-keng Liao <<wkliao@northwestern.edu>>
* Kai-yuan Hou <<kai-yuanhou2020@u.northwestern.edu>>

Copyright (C) 2021, Northwestern University.
See [COPYRIGHT](COPYRIGHT) notice in top-level directory.

### Project funding supports:
This research was supported by the Exascale Computing Project (17-SC-20-SC), a
joint project of the U.S. Department of Energy's Office of Science and National
Nuclear Security Administration, responsible for delivering a capable exascale
ecosystem, including software, applications, and hardware technology, to
support the nation's exascale computing imperative.

