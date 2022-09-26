# Case Study - Running with Cache Vol and Async Vol
This case study runs E3SM-IO using HDF5 with [Cache Vol](https://github.com/hpc-io/vol-cache) and [Async Vol](https://github.com/hpc-io/vol-async) enabled. 

## E3SM-IO with HDF5 using Cache Vol and Async Vol
### Installing HDF5, Async Vol and Cache Vol

1. <details> <summary>Install HDF5 1.13.2 (enable threadsafe, parallel, unsupported):</summary>

    ```shell
    % wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.13/hdf5-1.13.2/src/hdf5-1.13.2.tar.gz
    % tar -xf hdf5-1.13.2.tar.gz
    % cd hdf5-1.13.2
    
    # $HDF5_DIR points to HDF5 install dir.
    % ./configure --prefix=$HDF5_DIR --enable-parallel --enable-threadsafe --enable-unsupported CC=mpicc CXX=mpicxx
    % make; make install 
    ```

    </details>
1. <details> <summary>Install Argobots (required by Async Vol):</summary>

    ```shell
    % git clone https://github.com/pmodels/argobots.git
    % cd argobots
    % ./autogen.sh
    # $ABT_DIR points to Argobots install dir.
    % ./configure --prefix=$ABT_DIR CC=mpicc CXX=mpicxx
    % make; make install
    ```
    </details>

1. <details> <summary>Async Vol</summary>

    ```shell
    % export ABT_DIR=#path to argobots install dir
    % export HDF5_DIR=#path to hdf5 install dir
    % export HDF5_ROOT=${HDF5_DIR}
    % git clone https://github.com/hpc-io/vol-async.git
    % cd vol-async; mkdir build; cd build
    # $ASYNC_DIR points to Async Vol install dir
    % CC=mpicc CXX=mpicxx cmake .. -DCMAKE_INSTALL_PREFIX=$ASYNC_DIR
    % make; make install
    ```
    </details>

1. <details> <summary>Cache Vol.</summary>

    ```shell
    % git clone https://github.com/hpc-io/vol-cache.git
    % cd vol-cache; mkdir build; cd build
    % export LD_LIBRARY_PATH="$ABT_DIR/lib:$LD_LIBRARY_PATH"
    # $CAHCE_DIR points to Cache Vol install dir
    % CC=mpicc CXX=mpicxx HDF5_VOL_DIR=$ASYNC_DIR cmake .. -DCMAKE_INSTALL_PREFIX=$CAHCE_DIR
    % make; make install
    ```
    </details>

### Installing E3SM-IO
```shell
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
% make; make install
```

## Running E3SM-IO
1. <details> <summary>Set up environment:</summary>

    ```shell
    export HDF5_DIR=#path to hdf5 install dir
    export HDF5_ROOT=$HDF5_DIR
    export ABT_DIR=#path to argobots install dir
    export ASYNC_DIR=#path to async vol install dir
    export CACHE_DIR=#path to cache vol install dir

    export HDF5_PLUGIN_PATH=$CACHE_DIR/lib:$ASYNC_DIR/lib
    export LD_LIBRARY_PATH=$HDF5_PLUGIN_PATH:$ABT_DIR/lib:$HDF5_DIR/lib:$LD_LIBRARY_PATH
    export HDF5_VOL_CONNECTOR="cache_ext config=cache_1.cfg;under_vol=512;under_info={under_vol=0;under_info={}}"
    ```
    </details>

2. Run
    ```shell
    % cd E3SM-IO
    % mpiexec -n 16 src/e3sm_io datasets/f_case_866x72_16p.h5 -k -o can_F_out.h5 -a hdf5 -x canonical -r 5 -t 2
    ==== Benchmarking F case =============================
    Total number of MPI processes      = 16
    Number of IO processes             = 16
    Input decomposition file           = /files2/scratch/zhd1108/E3SM-IO/E3SM-IO/datasets/f_case_866x72_16p.h5
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
    Total no. noncontiguous requests   = 1977687
    Max   no. noncontiguous requests   = 189503
    Min   no. noncontiguous requests   =  63170
    Write no. records (time dim)       =      1
    I/O flush frequency                =      1
    No. I/O flush calls                =      1
    -----------------------------------------------------------
    Total write amount                         = 16.97 MiB = 0.02 GiB
    Time of I/O preparing              min/max =   0.0016 /   0.0017
    Time of file open/create           min/max =   0.3652 /   0.3887
    Time of define variables           min/max =  52.5990 /  65.7430
    Time of posting write requests     min/max =  28.7611 /  28.7702
    Time of write flushing             min/max =   0.4261 /   0.4370
    Time of close                      min/max =   4.9671 /   4.9838
    end-to-end time                    min/max = 100.3077 / 100.3243
    Simulated computation time         min/max =   2.0001 /   2.0001
    I/O bandwidth in MiB/sec (write-only)      = 0.5898
    I/O bandwidth in MiB/sec (open-to-close)   = 0.1691
    -----------------------------------------------------------
    ==== Benchmarking F case =============================
    Total number of MPI processes      = 16
    Number of IO processes             = 16
    Input decomposition file           = /files2/scratch/zhd1108/E3SM-IO/E3SM-IO/datasets/f_case_866x72_16p.h5
    Number of decompositions           = 3
    Output file/directory              = can_F_out.h5
    Using noncontiguous write buffer   = no
    Variable write order: same as variables are defined
    ==== HDF5 canonical I/O ==============================
    ```