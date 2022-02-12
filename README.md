## Parallel I/O Kernel Case Study -- E3SM

This repository contains a case study of parallel I/O kernel from the
[E3SM](https://github.com/E3SM-Project/E3SM) climate simulation model. The E3SM
I/O module, [Scorpio](https://github.com/E3SM-Project/scorpio), can be built on
top of [PnetCDF](https://github.com/Parallel-NetCDF/PnetCDF),
[NetCDF-4](http://www.unidata.ucar.edu/software/netcdf),
[HDF5](https://www.hdfgroup.org/solutions/hdf5), and
[ADIOS](https://github.com/ornladios/ADIOS2).
The benchmark program in this repository is designed to evaluate the E3SM I/O
kernel when using the above mentioned libraries to perform the I/O task. The
challenge of E3SM I/O patterns is that the problem domain is represented by
cubed sphere grids which produce long lists of small and noncontiguous requests
in each of all MPI processes.

The I/O patterns (describing the data decomposition of multi-dimensional arrays
representing the problem domain among MPI processes) used in this case study
were captured by the [PIO](https://github.com/NCAR/ParallelIO) library during
the E3SM production runs. A data decomposition file records the data access
patterns at the array element level by each MPI process. The access offsets are
stored in a text file, referred to by PIO as the "decomposition map file". This
benchmark currently studies three cases from E3SM, namely F, G and I cases.

In the F case, there are 414 climate variables stored in the 'h0' file.
Among the 414 variables, 6 are scalar variables and 408 are array variables.
Among the 414 variables, 15 are fixed-size variables and 399 are record variables.
Among the 414 variables, 387 are partitioned and 27 are not.
Among the 387 partitioned variables,
1 uses decomposition map 0,
323 use decomposition map 1, and
63 use decomposition map 2.
The 'h1' file stores 51 climate variables.
Among the 51 variables, 6 are scalar variables and 45 are array variables.
Among the 51 variables, 15 are fixed-size variables and 36 are record variables.
Among the 51 variables, 24 are partitioned and 27 are not partitioned.
Among the 24 partitioned variables,
1 uses decomposition map 0,
22 use decomposition map 1,
and 1 uses decomposition map 2.

In the G case, there are 52 climate variables stored in the output file.
All 52 variables are array variables. None is scalar.
Among the 52 variables, 11 are fixed-size variables and 41 are record variables.
Among the 52 variables, 41 are partitioned and 11 are not.
Among the 41 partitioned variables,
6 use decomposition map 0,
2 use decomposition map 1,
25 use decomposition map 2,
2 use decomposition map 3,
2 use decomposition map 4, and
4 use decomposition map 5.

In the I case, there are 560 climate variables stored in the 'h0' file.
All 560 variables are array variables. None is scalar.
Among the 560 variables, 18 are fixed-size variables and 542 are record variables.
Among the 560 variables, 546 are partitioned and 14 are not.
Among the 546 partitioned variables,
465 use decomposition map 0,
75 use decomposition map 1,
4 use decomposition map 2,
1 uses decomposition map 3, and
1 uses decomposition map 4.
The 'h1' file stores 552 climate variables.
All 552 variables are array variables. None is scalar.
Among the 552 variables, 10 are fixed-size variables and 542 are record variables.
Among the 552 variables, 538 are partitioned and 14 are not partitioned.
Among the 538 partitioned variables,
465 use decomposition map 0,
69 use decomposition map 1,
2 use decomposition map 2,
1 uses decomposition map 3, and
1 uses decomposition map 4.

## Compile and Run Instructions for E3SM-IO

### Software Requirements
* Autotools utility
  + autoconf 2.69
  + automake 1.16.1
  + libtoolize 2.4.6
  + m4 1.4.18
* MPI C and C++ compilers
  + Configured with a std 11 C++ compiler (supporting constant initializer)
* [PnetCDF 1.12.2](https://parallel-netcdf.github.io/Release/pnetcdf-1.12.2.tar.gz)
* (Optional) [HDF5 1.13.0](https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.13/hdf5-1.13.0/src/hdf5-1.13.0.tar.gz)
  + Configured with parallel I/O support (--enable-parallel is required)
* (Optional) [HDF5 Log-based VOL](https://github.com/DataLib-ECP/vol-log-based.git)
  + Experimental software developed as part of the Datalib project
* (Optional) [ADIOS 2.7.1](https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.7.1.tar.gz)
  + Configured with parallel I/O support (-DADIOS2_USE_MPI=ON is required)

### Instructions for Building Dependent I/O Libraries
* Build PnetCDF
  + Download a PnetCDF official released software
  + Configure PnetCDF with MPI C compiler
  + Run `make install`
  + Example build commands are given below. This example will install
    the PnetCDF library under folder `${HOME}/PnetCDF/1.12.2`.
    ```
    % wget https://parallel-netcdf.github.io/Release/pnetcdf-1.12.2.tar.gz
    % tar -zxf pnetcdf-1.12.2.tar.gz
    % cd pnetcdf-1.12.2
    % ./configure --prefix=${HOME}/PnetCDF/1.12.2 CC=mpicc
    % make -j 16 install
    ```
* (Optional) Build HDF5 with parallel I/O support
  + Download an HDF5 official released software.
  + Configure HDF5 with parallel I/O enabled.
  + Run `make install`
  + Example build commands are given below. This example will install
    the HD5 library under the folder `${HOME}/HDF5/1.13.0`.
    ```
    % wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.13/hdf5-1.13.0/src/hdf5-1.13.0.tar.gz
    % tar -zxf hdf5-1_13_0.tar.gz
    % cd hdf5-1.13.0
    % ./configure --prefix=${HOME}/HDF5/1.13.0 --enable-parallel CC=mpicc
    % make -j 16 install
    ```
* (Optional) Build HDF5 log-based VOL plugin.
  + Download the official released software.
  + Configure log-based VOL 
    + Enable shared library support (--enable-shared)
    + Compile with zlib library to enable metadata compression (--enable-zlib)
  + Example commands are given below.
    ```
    % wget https://github.com/DataLib-ECP/vol-log-based/archive/refs/tags/logvol.1.1.0.tar.gz
    % tar -zxf logvol.1.1.0.tar.gz
    % cd vol-log-based-logvol.1.1.0
    % ./configure --prefix=${HOME}/Log_VOL/1.1.0 --with-hdf5=${HOME}/HDF5/1.13.0 --enable-shared CC=mpicc
    % make -j 16 install
    ```
* (Optional) Build ADIOS with parallel I/O support
  + Download and extract the ADIOS source codes
  + Configure ADIOS with MPI support enabled (-DADIOS2_USE_MPI=ON)
  + Run `make install`
  + Example build commands are given below. This example will install
    the ADIOS2 library under the folder `${HOME}/ADIOS2/2.7.1`.
    ```
    % wget https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.7.1.tar.gz
    % tar -zxf v2.7.1.tar.gz
    % mkdir ADIOS2_BUILD
    % cd ADIOS2_BUILD
    % cmake -DCMAKE_INSTALL_PREFIX=${HOME}/ADIOS2/2.7.1 -DADIOS2_USE_MPI=ON ../ADIOS2-2.7.1
    % make -j 16 install
    ```

### Build E3SM-I/O benchmark
  + Clone this E3SM-I/O benchmark repository
  + Run command `autoreconf -i`
  + Configure the E3SM-I/O benchmark with MPI C and C++ compilers
    + Add HDF5 installation path (--with-hdf5=/path/to/implementation) that
      contains the HDF5 library. This is required when running the benchmark
      with HDF5 based I/O methods.
    + Add HDF5 log-based VOL installation path
      (--with-logvol=/path/to/implementation). This is required when running
      the benchmark with command-line option `-a hdf5_log -x log`.
    + Add ADIOS installation path (--with-adios2=/path/to/implementation) to
      enable ADIOS API support. This is required when running the benchmark
      with command-line option `-a adios -x log`.
  + Run `make`
  + Example commands are given below.
    ```
    % git clone https://github.com/Parallel-NetCDF/E3SM-IO.git
    % cd E3SM-IO
    % autoreconf -i
    % ./configure --with-pnetcdf=${HOME}/PnetCDF/1.12.2 \
                  --with-hdf5=${HOME}/HDF5/1.13.0 \
                  --with-logvol=${HOME}/Log_VOL/1.1.0 \
                  --with-adios2=${HOME}/ADIOS2/2.7.1 \
                  CC=mpicc CXX=mpicxx
    % make
    ```
  + The executable file, named 'e3sm_io', is created in folder 'src'.

### Prepare the Data Decomposition Map File
* Data decomposition maps generated by the PIO library are in text format (with
  file extension name ".dat". The decomposition maps must first be combined and
  converted into a NetCDF file to be read in parallel as the input file by this
  benchmark program. For the F, G, and I cases, there are 3, 6, and 5 data
  decomposition text files, respectively.
* See [utils/README](./utils) for instructions to run a utility program named
  `dat2nc.c` to convert the decomposition map files.

### Run command:
* Example run commands using `mpiexec` and 16 MPI processes are given below.
  + Run the write test with default settings, i.e. using PnetCDF library and
    producing files storing variables in a canonical data layout.
    ```
    % mpiexec -n 16 src/e3sm_io -o can_F_out.nc datasets/f_case_866x72_16p.nc
    ```
* The number of MPI processes used to run this benchmark can be smaller than
  the one used when creating the decomposition maps, i.e. the value of variable
  `decomp_nprocs` stored in the decomposition NetCDF file. For example, in file
  `datasets/f_case_866x72_16p.nc`, the value of scalar variable `decomp_nprocs`
  is 16, which is the number of MPI processes originally used to generate the
  decomposition `.dat` files. When running this benchmark using a smaller
  number of MPI processes, the I/O workload will be divided among all the
  allocated MPI processes. When using more processes than `decomp_nprocs`, the
  processes with MPI ranks greater than or equal to `decomp_nprocs` will have
  no data to write but still participate the collective I/O in the benchmark.

* Command-line Options:
  ```
    % ./e3sm_io -h
    Usage: ./e3sm_io [OPTION] FILE
       [-h] Print this help message
       [-v] Verbose mode
       [-k] Keep the output files when program exits (default: deleted)
       [-m] Run test using noncontiguous write buffer (default: contiguous)
       [-f num] Output history files h0 or h1: 0 for h0 only, 1 for h1 only,
                -1 for both. Affect only F and I cases. (default: -1)
       [-r num] Number of time records/steps written in F case h1 file and I
                case h0 file (default: 1)
       [-y num] Data flush frequency. (1: flush every time step, the default,
		and -1: flush once for all time steps. (No effect on ADIOS
		and HDF5 blob I/O options, which always flushes at file close).
       [-s num] Stride interval of ranks for selecting MPI processes to perform
                I/O tasks (default: 1, i.e. all MPI processes).\n\
       [-g num] Number of subfiles, used by ADIOS I/O only (default: 1).
       [-o path] Output file path (folder name when subfiling is used, file
                 name otherwise).
       [-a api]  I/O library name
           pnetcdf:   PnetCDF library (default)
           hdf5:      HDF5 library
           hdf5_log:  HDF5 library with Log-based VOL
           adios:     ADIOS library using BP3 format
       [-x strategy] I/O strategy
           canonical: Store variables in the canonical layout (default).
           log:       Store variables in the log-based storage layout.
	   blob:      Pack and store all data written locally in a contiguous
	              block (blob), ignoring variable's canonical order.
       FILE: Name of input file storing data decomposition maps
  ```
* Both F and I cases create two history files, referred to as 'h0' and 'h1'
  files. The supplied file name in option `-o` will be used to construct the
  output file names by inserting/appending strings "_h0" and "_h1" to indicate
  the two history files. If the input path contains file extension `.nc` or
  `.h5`, "_h0" and "_h1" will be inserted before the file extension. Otherwise,
  they will be appended at the end. See examples in "Output files" section
  below.
* When using HDF5 API (i.e. "-a hdf5" or "-a hdf5_log"), the environment
  variable `HDF5_VOL_CONNECTOR`, if used, must be set to match the I/O strategy
  used.
  + If I/O strategy is canonical, i.e. "-a hdf5 -x canonical"),
    `HDF5_VOL_CONNECTOR` must not be set to use log-based VOL.
  + If I/O strategy is log ("-x log"), `HDF5_VOL_CONNECTOR` must be set to use
    log-based VOL.  

### Current supported APIs (option `-a`) and I/O strategies (option `-x`)
  + **-a pnetcdf -x canonical**
    * A single NetCDF file in the classic CDF5 format will be created. All
      variables stored in the file are in the canonical order and
      understandable by NetCDF and its third-party software.
    * If the output file systems allow users to customize the file striping
      configuration, such as Lustre, users are recommended to write to a folder
      with a high file striping count to obtain a good I/O performance.
    * Example run command:
      ```
      mpiexec -n 16 src/e3sm_io datasets/f_case_866x72_16p.nc -k -o can_F_out.nc -a pnetcdf -x canonical -r 25
      ```
  + **-a pnetcdf -x blob**
    * Multiple subfiles in the NetCDF format will be created. The files conform
      with NetCDF file format specification.
    * There will be one subfile per compute node used.
    * File name provided in option `-o` will be used as a base to create the
      subfile names, which have the numerical IDs appended as the suffix.
    * Because all variables are stored in a blob fashion in the files, the
      subfiles altogether can only be understood by the conversion utility
      tool [utils/pnetcdf_blob_replay.c](utils/pnetcdf_blob_replay.c), which
      is developed to run off-line after the completion of an E3SM run to
      convert the subfiles into a single regular NetCDF file.
    * The blobs are per-record based, which means all write requests to the
      same variable by different MPI processes are packed and stored in a
      contiguous file space, called blob. Within that blob, data layout follows
      the process rank order.
    * Example run command:
      ```
      mpiexec -n 16 src/e3sm_io datasets/f_case_866x72_16p.nc -k -o blob_F_out.nc -a pnetcdf -x blob -r 25
      ```
  + **-a hdf5 -x blob**
    * This is the blob I/O implementation using HDF5 library. Different from
      the PnetCDF blob I/O, the implementation of uses the per-process based
      blob I/O strategy, in which each process writes only one blob at file
      close time, no matter how many data sets/variables are written. All
      write requests to all variables by a process are first cached in memory
      until file close time, in which they are packed into a contiguous buffer,
      a blob, to be flushed out in a single write call. There is an additional
      write for the header data blob written by the root process only. This
      per-process based strategy is the same one used by ADIOS.
    * Multiple subfiles in HDF5 format will be created.
    * There will be one subfile per compute node used.
    * File name provided in option `-o` will be used as a base to create the
      subfile names, which have the numerical IDs appended as the suffix.
    * The HDF5 subfiles cannot be understood by the traditional HDF5 software.
      A utility tool program will be developed in the future to convert the
      subfiles into a single regular HDF5 filet.
    * Example run command:
      ```
      mpiexec -n 16 src/e3sm_io datasets/f_case_866x72_16p.nc -k -o blob_F_out.h5 -a hdf5 -x blob -r 25
      ```
  + **-a hdf5 -x canonical**
    * This option writes data using the HDF5 library.
    * The data layout of datasets store in the output file is in a canonical
      order.
    * The VOL to be used is determined by the environment variable
      `HDF5_VOL_CONNECTOR`.
      + If `HDF5_VOL_CONNECTOR` is not set, HDF5 will use the native VOL. The
	output file of the native VOL is a regular HDF5, which is
	understandable by regular HDF5 and its third-party software.
      + This combination does not allow the use of log-based VOL, E3SM-I/O will
	fail if `HDF5_VOL_CONNECTOR` is set to use log-based VOL. Use options
	"-a hdf5 -x log" instead.
    * Example run command:
      ```
      mpiexec -n 16 src/e3sm_io datasets/f_case_866x72_16p.nc -k -o can_F_out.h5 -a hdf5 -x canonical -r 25
      ```
  + **-a hdf5 -x log**
    * This option writes data using the HDF5 log-based VOL.
    * The log-based VOL stores data in a log layout, rather than a canonical
      layout. The output file is a valid HDF5 file but requires the log-based
      VOL to read and understand the data structures.
    * The `HDF5_VOL_CONNECTOR` environment variable must be set to use log-based VOL.
      + E3SM-I/O will not run if `HDF5_VOL_CONNECTOR` is not properly set.
    * Example run command:
      ```
      mpiexec -n 16 src/e3sm_io datasets/f_case_866x72_16p.nc -k -o can_F_out.h5 -a hdf5 -x log -r 25
      ```
  + **-a hdf5_log -x log**
    * This option writes data using the HDF5 log-based VOL and specifically
      makes use of the new API, "H5Dwrite_n", to allow writing multiple
      subarrays of a dataset in a single API call.
    * Example run command:
      ```
      mpiexec -n 16 src/e3sm_io datasets/f_case_866x72_16p.nc -k -o log_F_out.h5 -a hdf5_log -x log -r 25
      ```
  + **-a adios -x blob**
    * This option writes data using the ADIOS library.
    * Multiple subfiles in BP format will be created.
    * The number of subfile is determined by command-line option `-g`.
    * File name provided in option `-o` will be used as a base to create the
      folder names which store the subfiles. The folder names will have suffix
      ".bp.dir" appended. Each subfile name in its folder will have ".bp" and a
      numerical ID appended.
    * Because all variables are stored in a blob fashion in the files, the
      subfiles can only be understood by the Scorpio's conversion utility tool,
      [adios2pio-nm](https://github.com/E3SM-Project/scorpio/tree/master/tools/adios2pio-nm),
      which is developed to run off-line to convert the subfiles into a single
      regular NetCDF file.
    * This option requires the original PIO decomposition maps in the text
      format. They can be included in the converted NetCDF decomposition file
      by adding a command-line option `-r` when running dat2nc. See README file
      in the utils folder for more instructions. If the original decomposition
      map is not in the decomposition file, the E3SM benchmark will create it
      by expanding the offset and length pairs in the converted decomposition
      map into list of offsets accessed.
    * Example run command:
      ```
      mpiexec -n 16 src/e3sm_io datasets/f_case_866x72_16p.nc -k -o blob_F_out -a adios -x blob -r 25
      ```

### Example input and run script files
* A median-size decomposition file `datasets/f_case_48602x72_512p.nc` contains
  the I/O pattern from an E3SM experiment ran on 512 MPI processes.
* Two large decomposition files `f_case_72x777602_21632p.nc.gz` (60 MB) and
  `g_case_11135652x80_9600p.nc.gz` (126 MB) from high-resolution simulations of
  F and G cases running on 21632 and 9600 MPI processes respectively are
  available upon request.
* An example batch script file for running a job on Cori @NERSC with 8 KNL
  nodes, 64 MPI processes per node, is provided in [slurm.knl](./slurm.knl).

### Example Output Shown on Screen
```
  % mpiexec -n 16 ./e3sm_io -o can_F_out.nc datasets/f_case_866x72_16p.nc
  ==== Benchmarking F case =============================
  Total number of MPI processes      = 16
  Number of IO processes             = 16
  Input decomposition file           = datasets/f_case_866x72_16p.nc
  Number of decompositions           = 3
  Output file/directory              = can_F_out.nc
  Using noncontiguous write buffer   = no
  Variable write order: same as variables are defined
  ==== PnetCDF canonical I/O using varn API ============
  History output file                = can_F_out_h0.nc
  No. variables use no decomposition =     27
  No. variables use decomposition D0 =      1
  No. variables use decomposition D1 =    323
  No. variables use decomposition D2 =     63
  Total no. climate variables        =    414
  Total no. noncontiguous requests   = 1977687
  Max   no. noncontiguous requests   = 189503
  Write no. records (time dim)       =      1
  I/O flush frequency                =      1
  No. I/O flush calls                =      1
  -----------------------------------------------------------
  Total write amount                 = 16.16 MiB = 0.02 GiB
  Max Time of I/O preparing          = 0.0012 sec
  Max Time of file open/create       = 0.0003 sec
  Max Time of define variables       = 0.0674 sec
  Max Time of posting write requests = 0.0331 sec
  Max Time of write flushing         = 0.7664 sec
  Max Time of close                  = 0.0044 sec
  Max end-to-end time                = 0.8729 sec
  I/O bandwidth (write-only)         = 21.0894 MiB/sec
  I/O bandwidth (open-to-close)      = 18.5172 MiB/sec
  -----------------------------------------------------------
  ==== Benchmarking F case =============================
  Total number of MPI processes      = 16
  Number of IO processes             = 16
  Input decomposition file           = datasets/f_case_866x72_16p.nc
  Number of decompositions           = 3
  Output file/directory              = can_F_out.nc
  Using noncontiguous write buffer   = no
  Variable write order: same as variables are defined
  ==== PnetCDF canonical I/O using varn API ============
  History output file                = can_F_out_h1.nc
  No. variables use no decomposition =     27
  No. variables use decomposition D0 =      1
  No. variables use decomposition D1 =     22
  No. variables use decomposition D2 =      1
  Total no. climate variables        =     51
  Total no. noncontiguous requests   =  38332
  Max   no. noncontiguous requests   =   3668
  Write no. records (time dim)       =      1
  I/O flush frequency                =      1
  No. I/O flush calls                =      1
  -----------------------------------------------------------
  Total write amount                 = 0.34 MiB = 0.00 GiB
  Max Time of I/O preparing          = 0.0000 sec
  Max Time of file open/create       = 0.0048 sec
  Max Time of define variables       = 0.0063 sec
  Max Time of posting write requests = 0.0007 sec
  Max Time of write flushing         = 0.0099 sec
  Max Time of close                  = 0.0006 sec
  Max end-to-end time                = 0.0223 sec
  I/O bandwidth (write-only)         = 34.0947 MiB/sec
  I/O bandwidth (open-to-close)      = 15.1156 MiB/sec
  -----------------------------------------------------------
  read_decomp=0.00 e3sm_io_core=0.90 MPI init-to-finalize=0.90
  -----------------------------------------------------------
```
### Output Files
* The above example command uses command-line option `-k` to keep the output
  files (otherwise the default is to delete them when the program exits.) For
  the F case, each run of `e3sm_io` produces two history output files whose
  names are created by inserting "_h0", and "_h1" to user-supplied file name.
  The header of F case files from running the provided decomposition file
  `f_case_866x72_16p.nc` using PnetCDF obtainable by command `ncmpidump -h` is
  available in [datasets/f_case_h0.txt](datasets/f_case_h0.txt), and
  [datasets/f_case_h1.txt](datasets/f_case_h1.txt).
* The G case only creates one output file. When using the PnetCDF I/O method
  and the provided decomposition file `g_case_cmpaso_16p.nc` to run, the header
  of output file can be found in [datasets/g_case_hist.txt](datasets/g_case_hist.txt).
* The option '-a adios' automatically appends ".bp.dir" extension to the
  user-provided input path and creates two folders for F  and I cases (one for
  G case.)
  + The names of output subfiles will be appended with file extension
    ".bp.dir".

### Current build status
* [Github Action](https://github.com/Parallel-NetCDF/E3SM-IO/actions)

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

