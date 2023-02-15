# Build Instructions for E3SM-IO Benchmark

* [Software Requirements](#software-requirements)
* [Instructions for Building Dependent I/O Libraries](#instructions-for-building-dependent-io-libraries)
* [Build E3SM-IO](#build-e3sm-io)
* [Prepare the Data Decomposition Map Files](#prepare-the-data-decomposition-map-files)
* [Run command](#run-command)
* [Example input and job script files](#example-input-and-job-script-files)
* [Example Output Shown on Screen](#example-output-shown-on-screen)
* [Output Files](#output-files)

## Software Requirements
* Autotools utility
  + autoconf 2.69
  + automake 1.16.1
  + libtool 2.4.6
  + m4 1.4.18
* MPI C and C++ compilers
  + Configured with a std 11 C++ compiler (supporting constant initializer)
* (Optional) [PnetCDF 1.12.3](https://parallel-netcdf.github.io/Release/pnetcdf-1.12.3.tar.gz)
* (Optional) [HDF5 1.14.0](https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.0/src/hdf5-1.14.0.tar.gz)
  + Configured with parallel I/O support (configured with `--enable-parallel` is required)
* (Optional) [HDF5 Log VOL connector](https://github.com/DataLib-ECP/vol-log-based.git) 1.4.0
  + Software developed as part of the Datalib project
* (Optional) [ADIOS 2.8.3](https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.8.3.tar.gz)
  + Configured with parallel I/O support (cmake with `-DADIOS2_USE_MPI=ON` is required)
* (Optional) [NetCDF-C 4.9.0](https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.9.0.tar.gz)
  + Configured with parallel HDF5 support (i.e. `--enable-netcdf4`)
  + Note using NetCDF-C versions prior to 4.9.0 will fail to run due to a
    [bug](https://github.com/Unidata/netcdf-c/issues/2251) related to dimension
    scales.

## Instructions for Building Dependent I/O Libraries
* Build PnetCDF
  + Download a PnetCDF official released software
  + Configure PnetCDF with MPI C compiler
  + Run `make install`
  + Example build commands are given below. This example will install
    the PnetCDF library under folder `${HOME}/PnetCDF/1.12.3`.
    ```console
    % wget https://parallel-netcdf.github.io/Release/pnetcdf-1.12.3.tar.gz
    % tar -zxf pnetcdf-1.12.3.tar.gz
    % cd pnetcdf-1.12.3
    % ./configure --prefix=${HOME}/PnetCDF/1.12.3 CC=mpicc
    % make -j 4 install
    ```
* Build HDF5 with parallel I/O support
  + Download an HDF5 official released software (version 1.13.0 and later is required).
  + Configure HDF5 with parallel I/O enabled.
  + Run `make install`
  + Example build commands are given below. This example will install
    the HDF5 library under the folder `${HOME}/HDF5/1.14.0`.
    ```console
    % wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.0/src/hdf5-1.14.0.tar.gz
    % tar -zxf hdf5-1.14.0.tar.gz
    % cd hdf5-1.14.0
    % ./configure --prefix=${HOME}/HDF5/1.14.0 --enable-parallel CC=mpicc
    % make -j 4 install
    ```
* Build HDF5 Log VOL connector
  + Download the official released software.
  + Configure  Log VOL connector
    + Enable shared library support (--enable-shared)
    + Compile with zlib library to enable metadata compression (--enable-zlib)
  + Example commands are given below.
    ```console
    % wget https://github.com/DataLib-ECP/vol-log-based/archive/refs/tags/logvol.1.4.0.tar.gz
    % tar -zxf logvol.1.4.0.tar.gz
    % cd vol-log-based-logvol.1.4.0
    % ./configure --prefix=${HOME}/Log_VOL/1.4.0 --with-hdf5=${HOME}/HDF5/1.14.0 --enable-shared CC=mpicc
    % make -j 4 install
    ```
* Build ADIOS with parallel I/O support
  + Download and extract the ADIOS source codes
  + Configure ADIOS with MPI support enabled (-DADIOS2_USE_MPI=ON)
  + Run `make install`
  + Example build commands are given below. This example will install
    the ADIOS2 library under the folder `${HOME}/ADIOS2/2.8.3`.
    ```console
    % wget https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.8.3.tar.gz
    % tar -zxf v2.8.3.tar.gz
    % mkdir ADIOS2_BUILD
    % cd ADIOS2_BUILD
    % cmake -DCMAKE_INSTALL_PREFIX=${HOME}/ADIOS2/2.8.3 -DADIOS2_USE_MPI=ON ../ADIOS2-2.8.3
    % make -j 4 install
    ```
* Build NetCDF-C
  + Download a NetCDF-C official released software (version 4.9.0 later is required).
  + Configure NetCDF-C with parallel HDF5 I/O enabled.
  + Run `make install`
  + Example build commands are given below. This example will install
    the NetCDF library under the folder `${HOME}/NetCDF/4.9.0`.
    ```console
    % wget https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.9.0.tar.gz
    % tar -zxf v4.9.0.tar.gz
    % cd netcdf-c-4.9.0
    % ./configure --prefix=${HOME}/NetCDF/4.9.0 \
                  CC=mpicc \
                  CPPFLAGS=-I${HOME}/HDF5/1.14.0/include \
                  LDFLAGS=-L${HOME}/HDF5/1.14.0/lib \
                  LIBS=-lhdf5
    % make -j 4 install
    ```

## Build E3SM-IO
  + Clone this E3SM-I/O benchmark repository
  + Run command `autoreconf -i`
  + Configure the E3SM-I/O benchmark with MPI C and C++ compilers
    + Add PnetCDF installation path (--with-pnetcdf=/path/to/implementation)
      that contains the PnetCDF library. This is required when running the
      benchmark with PnetCDF I/O methods.
    + Add HDF5 installation path (--with-hdf5=/path/to/implementation) that
      contains the HDF5 library. This is required when running the benchmark
      with HDF5 based I/O methods.
    + Add the installation path of the Log VOL connector
      (--with-logvol=/path/to/implementation). This is required when running
      the benchmark with command-line option `-a hdf5_log -x log`.
    + Add ADIOS installation path (--with-adios2=/path/to/implementation) to
      enable ADIOS API support. This is required when running the benchmark
      with command-line option `-a adios -x log`.
    + Add NetCDF4 installation path (--with-netcdf4=/path/to/implementation)
      that contains the NetCDF4 library. This is required when running the
      benchmark with NetCDF4 I/O methods.
  + Run `make`
  + Example commands are given below.
    ```console
    % git clone https://github.com/Parallel-NetCDF/E3SM-IO.git
    % cd E3SM-IO
    % autoreconf -i
    % ./configure --with-pnetcdf=${HOME}/PnetCDF/1.12.3 \
                  --with-hdf5=${HOME}/HDF5/1.14.0 \
                  --with-logvol=${HOME}/Log_VOL/1.4.0 \
                  --with-adios2=${HOME}/ADIOS2/2.8.3 \
                  --with-netcdf4=${HOME}/NetCDF/4.9.0 \
                  CC=mpicc CXX=mpicxx
    % make -j 8
    ```
  + The executable file, named 'e3sm_io', is created in folder 'src'.
  + Note the make command can take long to finish, as there is a total of about
    1000 climate variables across all F/G/I cases to be defined and each has
    several attributes.

## Prepare the Data Decomposition Map Files
* Data decomposition maps generated by the PIO library are in text format (with
  file extension name ".dat". The decomposition maps must first be combined and
  converted into a NetCDF file to be read in parallel as the input file by this
  benchmark program. For the F, G, and I cases, there are 3, 6, and 5 data
  decomposition text files, respectively.
* See [utils/README.md](../utils/README.md) for instructions to run utility programs
  + `dat2nc` converts the decomposition map .dat files to NetCDF CDF-5 files.
  + `dat2decomp` is more general utility program that can convert the
    decomposition map .dat files in text format to a CDF5/HDF5/NetCDF-4/BP
    file.
  + `decomp_copy` copies and converts a decomposition map file in an
    HDF5/NetCDF-4/BP format to a different format.

## Run command
* Example run commands using `mpiexec` and 16 MPI processes are given below.
  + Run the write test with default settings, i.e. using PnetCDF library and
    producing files storing variables in a canonical data layout.
    ```console
    % mpiexec -n 16 src/e3sm_io -o can_F_out.nc datasets/map_f_case_16p.nc
    ```
* The number of MPI processes used to run this benchmark can be smaller than
  the one used when creating the decomposition maps, i.e. the value of variable
  `decomp_nprocs` stored in the decomposition NetCDF file. For example, in file
  `datasets/map_f_case_16p.nc`, the value of scalar variable `decomp_nprocs`
  is 16, which is the number of MPI processes originally used to generate the
  decomposition `.dat` files. When running this benchmark using a smaller
  number of MPI processes, the I/O workload will be divided among all the
  allocated MPI processes. When using more processes than `decomp_nprocs`, the
  processes with MPI ranks greater than or equal to `decomp_nprocs` will have
  no data to write but still participate the collective I/O in the benchmark.

* Command-line Options:
  ```console
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
       [-t time] Add sleep time to emulate the computation in order to 
                 overlapping I/O when Async VOL is used.
       [-o path] Output file path (folder name when subfiling is used, file
                 name otherwise).
       [-a api]  I/O library name
           pnetcdf:   PnetCDF library (default)
           netcdf4:   NetCDF-4 library
           hdf5:      HDF5 library
           hdf5_md:   HDF5 library using multi-dataset I/O APIs
           hdf5_log:  HDF5 library with Log VOL connector
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
  variable `HDF5_VOL_CONNECTOR` can be set to use a VOL connector, e.g. Cache
  VOL, Async VOL, or Log VOL.
  + If I/O strategy is canonical, i.e. "-a hdf5 -x canonical") and
    `HDF5_VOL_CONNECTOR` is set to use Log VOL connector, then the output file
    will be in log layout.
  + If I/O strategy is log ("-x log"), `HDF5_VOL_CONNECTOR` does not have to
    set to use the Log VOL connector. Internally, E3SM-IO will explicitly call
    `H5Pset_vol()` to enable Log VOL.

### Current supported APIs (option `-a`) and I/O strategies (option `-x`)
  + Table below lists the supported combinations.
     |           | pnetcdf | hdf5 | hdf5_log | hdf5_md  | netcdf4* | adios |
     |-----------|:-------:|:----:|:--------:|:--------:|:--------:|:-----:|
     | canonical | yes     | yes  | no       | yes      | yes      | no    |
     | log       | no      | yes  | yes      | no       | yes      | no    |
     | blob      | yes     | yes  | no       | no       | no       | yes   |
     
     `*` NetCDF-C version 4.9.0 or newer is required.

  + **-a pnetcdf -x canonical**
    * A single NetCDF file in the classic CDF5 format will be created. All
      variables stored in the file are in the canonical order and
      understandable by NetCDF and its third-party software.
    * If the parallel file systems allow users to customize the file striping
      configuration, such as Lustre, users are recommended to configure the
      output folder with a high file striping count to obtain a good I/O
      performance.
    * Example run command:
      ```console
      mpiexec -n 16 src/e3sm_io datasets/map_f_case_16p.nc -k -o can_F_out.nc -a pnetcdf -x canonical -r 25
      ```
  + **-a pnetcdf -x blob**
    * Multiple subfiles in the NetCDF format will be created. The files conform
      with NetCDF file format specification.
    * There will be one subfile per compute node used.
    * File name provided in option `-o` will be used as a base to create the
      subfile names, which have the numerical IDs appended as the suffix.
    * Because all variables are stored in a blob fashion in the files, the
      subfiles altogether can only be understood by the conversion utility
      tool [utils/pnetcdf_blob_replay.c](../utils/pnetcdf_blob_replay.c), which
      is developed to run off-line after the completion of an E3SM run to
      convert the subfiles into a single regular NetCDF file.
    * The blobs are per-record based, which means all write requests to the
      same variable by different MPI processes are packed and stored in a
      contiguous file space, called blob. Within that blob, data layout follows
      the process rank order.
    * Example run command:
      ```console
      mpiexec -n 16 src/e3sm_io datasets/map_f_case_16p.nc -k -o blob_F_out.nc -a pnetcdf -x blob -r 25
      ```
  + **-a hdf5 -x canonical**
    * This option writes/reads data using HDF5 APIs `H5Dwrite`/`H5Dread`.
    * If the environment variable `HDF5_VOL_CONNECTOR` is unset or set without
      Log VOL, then the output file will be in the canonical layout and command
      `h5ldump -k` will show the file kind of `HDF5`.
    * If the environment variable `HDF5_VOL_CONNECTOR` is set to use Log VOL,
      then the output file will be in the log layout.  Running command
      `h5ldump -k` will show the file kind of `HDF5-LogVOL`.
    * Example run command:
      ```console
      mpiexec -n 16 src/e3sm_io datasets/map_f_case_16p.h5 -k -o can_F_out.h5 -a hdf5 -x canonical -r 25
      ```
  + **-a hdf5 -x blob**
    * This is the blob I/O implementation using HDF5 library. Different from
      the PnetCDF blob I/O, it uses the per-process based blob I/O strategy, in
      which each process writes only one blob in the output file at file close
      time, no matter how many data sets/variables are written. All write
      requests to all variables by a process are first cached in memory until
      file close time, at which time they are packed into a contiguous buffer,
      and flushed out by a single write call. There is an additional write for
      the header data blob written by the root process only. This per-process
      based strategy is the same one used by ADIOS.
    * Multiple subfiles will be created. Each subfile is also an HDF5 file.
    * There will be one subfile per compute node used.
    * File name provided in the command-line option `-o` will be used as a base
      to create the subfile names, which have the numerical IDs appended as the
      suffix.
    * The HDF5 subfiles cannot be understood by the traditional HDF5 software.
      A utility tool program will be developed in the future to convert the
      subfiles into a single regular HDF5 file.
    * Example run command:
      ```console
      mpiexec -n 16 src/e3sm_io datasets/map_f_case_16p.h5 -k -o blob_F_out.h5 -a hdf5 -x blob -r 25
      ```
    * If the environment variable HDF5_VOL_CONNECTOR is set to use Log VOL,
      then the subfiles will also be in the log layout. Running command
      `h5ldump -k` will show the file kind of `HDF5-LogVOL`.
  + **-a hdf5 -x log**
    * This option requires the Log VOL feature enabled at the configure time,
      i.e. "`--with-logvol=${LOGVOL_DIR}`" used at the configure command line.
    * All datasets stored in the files will be in the log layout. Running
      command `h5ldump -k` will show the file kind of `HDF5-LogVOL`.
    * E3SM-IO will write/read data using HDF5 APIs `H5Dwrite`/`H5Dread`.
    * If the environment variable `HDF5_VOL_CONNECTOR` is unset or set without
      Log VOL, then E3SM-IO will explicitly call `H5Pset_vol()` to enable the
      HDF5 Log VOL connector.
    * If the environment variable `HDF5_VOL_CONNECTOR` is set to use other VOL
      connectors, such as Cache and Async VOLS, then E3SM-IO will stack the Log
      VOL on top of those connectors.
    * The output file is a valid HDF5 file but requires the Log VOL connector
      to read and understand the data structures.
    * Example run command:
      ```console
      mpiexec -n 16 src/e3sm_io datasets/map_f_case_16p.h5 -k -o can_F_out.h5 -a hdf5 -x log -r 25
      ```
  + **-a hdf5_md -x canonical**
    * This option writes/reads data using HDF5 multi-dataset APIs
      `H5Dwrite_multi`/`H5Dread_multi`. Command `h5ldump -k` will show the file
      kind of `HDF5`.
    * If the environment variable `HDF5_VOL_CONNECTOR` containing Log VOL, the
      environment variable will be unset.
    * Warning! HDF5 versions 1.13.3 and 1.14.0 will switch collective I/O mode
      to independent internally when one of the datasets requires data type
      conversion. See https://github.com/HDFGroup/hdf5/issues/1859
    * Example run command:
      ```console
      mpiexec -n 16 src/e3sm_io datasets/map_f_case_16p.h5 -k -o can_F_out.h5 -a hdf5_md -x canonical -r 25
      ```
  + **-a hdf5_md -x log**
    * This option is not supported.
  + **-a hdf5_md -x blob**
    * This option is not supported.
  + **-a hdf5_log -x log**
    * This option writes data using the HDF5 Log VOL connector by explicitly
      calling `H5Pset_vol()` internally.
    * For dataset I/O, this option calls the APIs `H5Dwrite_n()` and
      `H5Dread_n()` created in the Log VOL, rather than the HDF5 `H5Dwrite()`
      or `H5Dread()`. The two new APIs allow to write and read multiple
      subarrays of a dataset in a single API call. They are expected to perform
      better, as their computational costs and memory footprints for metadata
      operations are smaller.
    * Datasets stored in the output will be in the log layout. Running command
      `h5ldump -k` will show the file kind of `HDF5-LogVOL`.
    * Example run command:
      ```console
      mpiexec -n 16 src/e3sm_io datasets/map_f_case_16p.h5 -k -o log_F_out.h5 -a hdf5_log -x log -r 25
      ```
  + **-a hdf5_log -x canonical**
    * This option is not supported.
  + **-a hdf5_log -x blob**
    * This option is not supported.
  + **-a netcdf4 -x canonical**
    * This option writes data using the NetCDF-4 library.
    * The output files are in the HDF5 format. Running command `h5ldump -k`
      will show the file kind of `NetCDF-4`.
    * The data layout of datasets store in the output file is in a canonical
      order.
    * Because the number of write requests are different among processes, the
      independent I/O mode is used when writing the data to files.
    * Example run command:
      ```console
      mpiexec -n 16 src/e3sm_io datasets/map_f_case_16p.nc4 -k -o can_F_out.nc4 -a netcdf4 -x canonical -r 25
      ```
    * If environment variables `HDF5_VOL_CONNECTOR` and `HDF5_PLUGIN_PATH` are
      set to use Log VOL, then the execution will abort, as this option is
      equivalent to `-a netcdf4 -x log`.
  + **-a netcdf4 -x log**
    * This option writes data using the NetCDF-4 library which calls the HDF5
      Log VOL connector underneath.
    * **Requirements** - The two environment variables `HDF5_VOL_CONNECTOR` and
      `HDF5_PLUGIN_PATH` must be set to use Log VOL connector in order to run.
      The e3sm_io program will check and error out if they are not set.
    * The Log VOL stores data in a log layout, rather than a canonical layout.
      The output file is a valid HDF5 file but requires the Log VOL to read and
      understand the data structures.
    * Running command `h5ldump -k` will show the file kind of `HDF5-LogVOL`.
    * Example run command:
      ```console
      export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/LOG_VOL/lib
      export HDF5_PLUGIN_PATH=${HOME}/LOG_VOL/lib
      export HDF5_VOL_CONNECTOR="LOG under_vol=0;under_info={}"
      mpiexec -n 16 src/e3sm_io datasets/map_f_case_16p.nc4 -k -o log_F_out.nc4 -a netcdf4 -x log -r 25
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
      in folder `utils` for more instructions. If the original decomposition
      map is not in the decomposition file, the E3SM benchmark will create it
      by expanding the offset and length pairs in the converted decomposition
      map into list of offsets accessed.
    * Example run command:
      ```console
      mpiexec -n 16 src/e3sm_io datasets/map_f_case_16p.bp -k -o blob_F_out -a adios -x blob -r 25
      ```
### Run E3SM-IO with Cache VOL and Async VOL:
*  When using HDF5 as the I/O method, [Cache VOL](https://github.com/hpc-io/vol-cache) and [Async VOL](https://github.com/hpc-io/vol-async) can be enabled. See [cache_async_vol.md](./cache_async_vol.md).

## Example input and job script files
* Three small-size decomposition map files are available for testing. They
  are generated from E3SM runs on 16 MPI processes.
  + F case uses 3 decomposition maps.
    + File `datasets/map_f_case_16p.nc` is in NetCDF classic CDF-5 format
    + File `datasets/map_f_case_16p.h5` is in HDF5 format
    + File `datasets/map_f_case_16p.nc4` is in NetCDF4 format
    + File `datasets/map_f_case_16p.bp` is in ADIOS BP format
  + G case uses 6 decomposition maps.
    + File `datasets/map_g_case_16p.nc` is in NetCDF classic CDF-5 format
    + File `datasets/map_g_case_16p.h5` is in HDF5 format
    + File `datasets/map_g_case_16p.nc4` is in NetCDF4 format
    + File `datasets/map_g_case_16p.bp` is in ADIOS BP format
  + I case uses 5 decomposition maps.
    + File `datasets/map_i_case_16p.nc` is in NetCDF classic CDF-5 format
    + File `datasets/map_i_case_16p.h5` is in HDF5 format
    + File `datasets/map_i_case_16p.nc4` is in NetCDF4 format
    + File `datasets/map_i_case_16p.bp` is in ADIOS BP format
* File `datasets/f_case_48602x72_512p.nc` contains 3 decomposition maps for a
  median-size F case produced from a 512-process run.
* Three large decomposition files are available upon request.
  + `f_case_21600p.nc` (266 MB) for F case produced from 21600 processes.
  + `g_case_9600p.nc` (303 MB) for G case produced from 9600 processes.
  + `i_case_1344p.nc` (12 MB)for I case produced from 1344 processes.
* An example batch script file for running a job on Cori @NERSC with 8 KNL
  nodes, 64 MPI processes per node, is provided in [slurm.knl](./slurm.knl).

## Example Output Shown on Screen
```console
  % mpiexec -n 16 src/e3sm_io -o can_F_out.nc datasets/map_f_case_16p.nc
  ==== Benchmarking F case =============================
  Total number of MPI processes      = 16
  Number of IO processes             = 16
  Input decomposition file           = datasets/map_f_case_16p.nc
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
  Total no. attributes               =   1421
  Total no. noncontiguous requests   = 1977687
  Max   no. noncontiguous requests   = 189503
  Min   no. noncontiguous requests   =  63170
  Write no. records (time dim)       =      1
  I/O flush frequency                =      1
  No. I/O flush calls                =      1
  -----------------------------------------------------------
  Total write amount                         = 16.16 MiB = 0.02 GiB
  Time of I/O preparing              min/max =   0.0008 /   0.0013
  Time of file open/create           min/max =   0.0005 /   0.0006
  Time of define variables           min/max =   0.0031 /   0.0033
  Time of posting write requests     min/max =   0.0124 /   0.0257
  Time of write flushing             min/max =   0.2817 /   0.2837
  Time of close                      min/max =   0.0029 /   0.0029
  end-to-end time                    min/max =   0.3175 /   0.3176
  Emulate computation time (sleep)   min/max =   0.0000 /   0.0000
  I/O bandwidth in MiB/sec (write-only)      = 56.9648
  I/O bandwidth in MiB/sec (open-to-close)   = 50.8962
  -----------------------------------------------------------
  ==== Benchmarking F case =============================
  Total number of MPI processes      = 16
  Number of IO processes             = 16
  Input decomposition file           = datasets/map_f_case_16p.nc
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
  Total no. attributes               =    142
  Total no. noncontiguous requests   =  38332
  Max   no. noncontiguous requests   =   3668
  Min   no. noncontiguous requests   =   1225
  Write no. records (time dim)       =      1
  I/O flush frequency                =      1
  No. I/O flush calls                =      1
  -----------------------------------------------------------
  Total write amount                         = 0.34 MiB = 0.00 GiB
  Time of I/O preparing              min/max =   0.0000 /   0.0000
  Time of file open/create           min/max =   0.0005 /   0.0005
  Time of define variables           min/max =   0.0002 /   0.0003
  Time of posting write requests     min/max =   0.0002 /   0.0004
  Time of write flushing             min/max =   0.0034 /   0.0034
  Time of close                      min/max =   0.0002 /   0.0003
  end-to-end time                    min/max =   0.0049 /   0.0049
  Emulate computation time (sleep)   min/max =   0.0000 /   0.0000
  I/O bandwidth in MiB/sec (write-only)      = 98.0321
  I/O bandwidth in MiB/sec (open-to-close)   = 68.3491
  -----------------------------------------------------------
  read_decomp=0.00 e3sm_io_core=0.32 MPI init-to-finalize=0.33
  -----------------------------------------------------------
```
## Output Files
* The above example command uses command-line option `-k` to keep the output
  files (otherwise the default is to delete them when the program exits.) For
  the F case, each run of `e3sm_io` produces two history output files whose
  names are created by inserting "_h0", and "_h1" to user-supplied file name.
  The header of F case files from running the provided decomposition file
  `map_f_case_16p.nc` using PnetCDF obtainable by command `ncmpidump -h` is
  available in [datasets/f_case_h0.txt](../datasets/f_case_h0.txt), and
  [datasets/f_case_h1.txt](../datasets/f_case_h1.txt).
* The G case only creates one output file. When using the PnetCDF I/O method
  and the provided decomposition file `map_g_case_16p.nc` to run, the header
  of output file can be found in [datasets/g_case_hist.txt](../datasets/g_case_hist.txt).
* The option '-a adios' automatically appends ".bp.dir" extension to the
  user-provided input path and creates two folders for F  and I cases (one for
  G case.)
  + The names of output subfiles will be appended with file extension
    ".bp.dir".

---
* Copyright (C) 2021, Northwestern University.
* See [COPYRIGHT](../COPYRIGHT) notice in top-level directory.

