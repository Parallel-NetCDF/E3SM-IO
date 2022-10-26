# Running E3SM-IO with HDF5 Cache VOL and Async VOL
* [Build Instructions](#build-instructions)
* [Run Instructions](#run-instructions)
* [Example Output](#example-output)

This case study runs E3SM-IO using HDF5 with [Cache VOL](https://github.com/hpc-io/vol-cache) and [Async VOL](https://github.com/hpc-io/vol-async) enabled.

Cache VOL is an HDF5 plugin that incorporates fast storage layers (e.g, burst buffer, node-local storage) into parallel I/O workflow for caching and staging data to improve the I/O efficiency.

Async VOL is another HDF5 plugin that takes advantage of an asynchronous interface by scheduling I/O as early as possible and overlaps computation or communication with I/O operations, which hides the cost associated with I/O and improves the overall performance.

Both Cache VOL and Async VOL can be enabled by directly setting the environment variables without modifying E3SM-IO source codes. This case study gives an instruction on how to install Cache VOL and Async VOL, and gives an demo of how to run E3SM-IO with them.

## Build Instructions
### Prerequisite
+ Set up environment:

    ```shell
    export HDF5_DIR=#the dir you want to install HDF5 to
    export ABT_DIR=#the dir you want to install argobots to
    export ASYNC_DIR=#the dir you want to install Async VOL to
    export CACHE_DIR=#the dir you want to install Cache VOL to
 
    export HDF5_ROOT=$HDF5_DIR
    ```

+ HDF5 1.13.2: `--enable-parallel`, `--enable-threadsafe`, and `--enable-unsupported` are [required by Cache VOL](https://vol-cache.readthedocs.io/en/latest/gettingstarted.html#installation) at configure time.

    ```shell
    # the following env variable will be used:
    # HDF5_DIR

    % export $HDF5_DIR=#path to the directory to install HDF5
    % wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.13/hdf5-1.13.2/src/hdf5-1.13.2.tar.gz
    % tar -xf hdf5-1.13.2.tar.gz
    % cd hdf5-1.13.2
    
    % ./configure --prefix=$HDF5_DIR --enable-parallel --enable-threadsafe --enable-unsupported CC=mpicc CXX=mpicxx
    % make
    % make install
    ```

+ Argobots, [required by Async VOL](https://hdf5-vol-async.readthedocs.io/en/latest/gettingstarted.html#build-async-i-o-vol):

    ```shell
    # the following env variable will be used:
    # ABT_DIR

    % git clone https://github.com/pmodels/argobots.git
    % cd argobots
    % ./autogen.sh

    % ./configure --prefix=$ABT_DIR CC=mpicc CXX=mpicxx
    % make
    % make install
    ```

+ Async VOL

    ```shell
    # the following env variables will be used:
    # HDF5_DIR, ABT_DIR, ASYNC_DIR, HDF5_ROOT

    % git clone https://github.com/hpc-io/vol-async.git
    % cd vol-async
    % mkdir build
    % cd build
    % CC=mpicc CXX=mpicxx cmake .. -DCMAKE_INSTALL_PREFIX=$ASYNC_DIR
    % make
    % make install
    ```

+ Cache VOL

    ```shell
    # the following env variables will be used:
    # ABT_DIR, ASYNC_DIR, CAHCE_DIR

    % git clone https://github.com/hpc-io/vol-cache.git
    % cd vol-cache
    % mkdir build
    % cd build
    % export LD_LIBRARY_PATH="$ABT_DIR/lib:$LD_LIBRARY_PATH"
    % CC=mpicc CXX=mpicxx HDF5_VOL_DIR=$ASYNC_DIR cmake .. -DCMAKE_INSTALL_PREFIX=$CAHCE_DIR
    % make
    % make install
    ```
 

### Installing E3SM-IO
`--enable-threading` option is required at configure time. Using
this option will initialize MPI with multi-thread support, i.e. initialize with `MPI_Init_thread` instead of `MPI_Init`. Using
`MPI_Init_thread` is [required by Cache VOL](https://vol-cache.readthedocs.io/en/latest/gettingstarted.html#examples).

```shell
# clone this repo
% # git clone git@github.com:Parallel-NetCDF/E3SM-IO.git
% git clone -b thread git@github.com:Parallel-NetCDF/E3SM-IO.git

% cd E3SM-IO
% autoreconf -i
% ./configure \
    --with-hdf5=$HDF5_DIR \
    --disable-profiling \
    --enable-threading \
    CC=mpicc \
    CXX=mpicxx
% make
% make install
```

## Run Instructions
1. Set up environment:

    ```shell
    # the followings are already set during installation.
    export HDF5_DIR=#path to hdf5 install dir
    export HDF5_ROOT=$HDF5_DIR
    export ABT_DIR=#path to argobots install dir
    export ASYNC_DIR=#path to Async VOL install dir
    export CACHE_DIR=#path to Cache VOL install dir

    # the followings are newly added env variables.
    export HDF5_PLUGIN_PATH=$CACHE_DIR/lib:$ASYNC_DIR/lib
    export LD_LIBRARY_PATH=$HDF5_PLUGIN_PATH:$ABT_DIR/lib:$HDF5_DIR/lib:$LD_LIBRARY_PATH
    export HDF5_VOL_CONNECTOR="cache_ext config=cache_1.cfg;under_vol=512;under_info={under_vol=0;under_info={}}"
    ```


1. Run commands
    ```shell
    % cd E3SM-IO
    % mpiexec -n 16 src/e3sm_io datasets/map_f_case_16p.h5 -k -o can_F_out.h5 -a hdf5 -x canonical -r 5 -t 2
    ```
## Example Output
```
==== Benchmarking F case =============================
Total number of MPI processes      = 16
Number of IO processes             = 16
Input decomposition file           = /files2/scratch/zhd1108/E3SM-IO/E3SM-IO/datasets/map_f_case_16p.h5
Number of decompositions           = 3
Output file/directory              = can_F_out.h5
Using noncontiguous write buffer   = no
Variable write order: same as variables are defined
==== HDF5 canonical I/O ==============================
History output file                = can_F_out_h0.h5
No. variables use no decomposition =     27
No. variables use decomposition D0 =      1
No. variables use decomposition D1 =    323
No. variables use decomposition D2 =     63
Total no. climate variables        =    414
Total no. attributes               =   1421
Total no. noncontiguous requests   = 4207968
Max   no. noncontiguous requests   = 272135
Min   no. noncontiguous requests   = 252670
Write no. records (time dim)       =      1
I/O flush frequency                =      1
No. I/O flush calls                =      1
-----------------------------------------------------------
Total write amount                         = 16.97 MiB = 0.02 GiB
Time of I/O preparing              min/max =   0.0010 /   0.0011
Time of file open/create           min/max =   0.5526 /   0.5896
Time of define variables           min/max =  65.7061 /  65.7251
Time of posting write requests     min/max =  33.5788 /  33.5797
Time of write flushing             min/max =   0.0000 /   0.0000
Time of close                      min/max =   5.3138 /   5.3703
end-to-end time                    min/max = 105.2095 / 105.2659
Emulate computation time (sleep)   min/max =   2.0000 /   2.0000
I/O bandwidth in MiB/sec (write-only)      = 0.5053
I/O bandwidth in MiB/sec (open-to-close)   = 0.1612
-----------------------------------------------------------
==== Benchmarking F case =============================
Total number of MPI processes      = 16
Number of IO processes             = 16
Input decomposition file           = /files2/scratch/zhd1108/E3SM-IO/E3SM-IO/datasets/map_f_case_16p.h5
Number of decompositions           = 3
Output file/directory              = can_F_out.h5
Using noncontiguous write buffer   = no
Variable write order: same as variables are defined
==== HDF5 canonical I/O ==============================
History output file                = can_F_out_h1.h5
No. variables use no decomposition =     27
No. variables use decomposition D0 =      1
No. variables use decomposition D1 =     22
No. variables use decomposition D2 =      1
Total no. climate variables        =     51
Total no. attributes               =    142
Total no. noncontiguous requests   = 400226
Max   no. noncontiguous requests   =  25963
Min   no. noncontiguous requests   =  24026
Write no. records (time dim)       =      5
I/O flush frequency                =      1
No. I/O flush calls                =      5
-----------------------------------------------------------
Total write amount                         = 1.88 MiB = 0.00 GiB
Time of I/O preparing              min/max =   0.0000 /   0.0000
Time of file open/create           min/max =   0.3636 /   0.3906
Time of define variables           min/max =   6.8050 /   6.8171
Time of posting write requests     min/max =  24.1418 /  24.1848
Time of write flushing             min/max =   0.0000 /   0.0000
Time of close                      min/max =   1.3208 /   1.3250
end-to-end time                    min/max =  32.7165 /  32.7208
Emulate computation time (sleep)   min/max =   2.0000 /   2.0000
I/O bandwidth in MiB/sec (write-only)      = 0.0776
I/O bandwidth in MiB/sec (open-to-close)   = 0.0574
-----------------------------------------------------------
read_decomp=1.56 e3sm_io_core=137.99 MPI init-to-finalize=139.56
-----------------------------------------------------------
```