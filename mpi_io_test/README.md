## MPI Collective Writes Using Contiguous and Noncontiguous User Buffers

This repository contains a program designed to evaluate the performance of MPI
collective write given a fileview consisting of a long list of noncontiguous,
small regions in the file. Under the same fileview, it compares two cases. One
uses a contiguous user buffer and the other a noncontiguous buffer.

The noncontiguous I/O pattern is used by the production runs of
[E3SM](https://github.com/E3SM-Project/E3SM) climate simulation model. The
offset-length pairs of individual write requests from all MPI processes are
traced and saved in the provided NetCDF files. PnetCDF library is used only to
read in the offset-length pairs. Once the inputs are read, the remaining
program makes only MPI-IO calls. Two noncontiguous MP derived datatypes are
constructed by calling `MPI_Type_create_hindexed()`. One is used to set the
fileview and the other to describe a noncontiguous user buffer when calling
`MPI_File_write_all()`. The timings of two collective write calls are measured
reported at the end. Both collective writes use the same fileview and data
amount.

When running this program on
[Cori](http://www.nersc.gov/users/computational-systems/cori) KNL nodes, we
found that the collective write performs significantly poorer when using a
noncontiguous user buffer than a contiguous one. Preliminary time profiling
shows a significant amount of time spent on posting `MPI_Irecv()`, and
`MPI_Isend()` in subroutine `ADIOI_LUSTRE_W_Exchange_data()` in ROMIO.
This occurs especially when a high number of MPI processes per KNL node
is used. We are still investigating the real cause.

* Compile command:
  * Edit file `Makefile` to customize the C compiler, compile options,
    location of PnetCDF library, etc.
  * [PnetCDF](https://github.com/Parallel-NetCDF/PnetCDF) library is used
    only to read the decomposition file.
  * Run command `make` to compile and generate the executable program named
    `noncontig_buf`.
  * To compile for KNL nodes on Cori @NERSC, run command model swap below,
    before running `make`.
    ```
      % module swap craype-haswell craype-mic-knl
    ```
  * A PnetCDF library built for KNL nodes is available on Cori in
    `/global/u2/w/wkliao/PnetCDF/1.10.0.KNL`.

* Example data decomposition file:
  * Three decomposition files in NetCDF format are provided, which store the
    data access patterns in form of offset-length pairs.
    * `48602x72_512p_D1.nc`, `48602x72_512p_D2.nc` and `48602x72_512p_D3.nc`
    * The file name `48602x72_512p_D3.nc` explains the global array of size
      48602x72 (in Fortran dimension order), decomposition done among 512
      processes, and pattern 3. The file header can be shown using command
      `ncdump -h`.
      ```
        % ncdump -h 48602x72_512p_D3.nc
        netcdf \48602x72_512p_D3 {
        dimensions:
              num_procs = 512 ;
              array_size = 3499344 ;
              nreq_00000 = 1152 ;
              nreq_00001 = 1152 ;
              nreq_00002 = 1152 ;
              ...
        variables:
              int off_00000(nreq_00000) ;
              int len_00000(nreq_00000) ;
              int off_00001(nreq_00001) ;
              int len_00001(nreq_00001) ;
              int off_00002(nreq_00002) ;
              int len_00002(nreq_00002) ;
              ...
      ```
* Run command:
  * Example run command using `mpiexec` and 512 MPI processes:
    ```
      % mpiexec -n 512 ./noncontig_buf -q 48602x72_512p_D3.nc -o output_file
    ```
  * It is recommended to use the same number of MPI processes as the value set
    in the dimension `num_procs` in the decomposition NetCDF file.
  * Command-line options:
    ```
      % ./noncontig_buf -h
      Usage: noncontig_buf [OPTION]... [FILE]...
         [-h] Print help
         [-q] Quiet mode
         [-n] number of variables
         [-o] output file name (default "./testfile")
         input_file: name of input netCDF file describing data decompositions
    ```
  * An example batch script file for job running on Cori @NERSC is given in
    `./pbs.knl`.

* Example outputs on screen
  ```
    % srun -n 512 -c 4 --cpu_bind=cores ./noncontig_buf -q -n 63 -o $SCRATCH/FS_1M_8/testfile $SCRATCH/FS_1M_8/48602x72_512p_D3.nc

    input  file name = $SCRATCH/FS_1M_8/48602x72_512p_D3.nc
    output file name = $SCRATCH/FS_1M_8/testfile
    -----------------------------------------------------------
    Total number of MPI processes        = 512
    Total number of variables            = 63
    Total write amount                   = 840.98 MiB = 0.82 GiB
    Max no. noncontig requests per var   = 4608
    Min no. noncontig requests per var   = 1080
    Max length of contig request         = 24 bytes
    Min length of contig request         = 4 bytes
    Max write amount per variable        = 0.03 MiB
    Min write amount per variable        = 0.02 MiB
    Max write time when buf is contig    = 8.9165 sec
    Max write time when buf is noncontig = 69.1925 sec
    Max time of MPI_Pack()               = 0.0061 sec
    -----------------------------------------------------------


    % srun -n 512 -c 4 --cpu_bind=cores ./noncontig_buf -q -n 1000 -o $SCRATCH/FS_1M_8/testfile $SCRATCH/FS_1M_8/48602x72_512p_D2.nc

    input  file name = $SCRATCH/FS_1M_8/48602x72_512p_D2.nc
    output file name = $SCRATCH/FS_1M_8/testfile
    -----------------------------------------------------------
    Total number of MPI processes        = 512
    Total number of variables            = 1000
    Total write amount                   = 185.40 MiB = 0.18 GiB
    Max no. noncontig requests per var   = 64
    Min no. noncontig requests per var   = 15
    Max length of contig request         = 24 bytes
    Min length of contig request         = 4 bytes
    Max write amount per variable        = 0.00 MiB
    Min write amount per variable        = 0.00 MiB
    Max write time when buf is contig    = 1.9888 sec
    Max write time when buf is noncontig = 4.8056 sec
    Max time of MPI_Pack()               = 0.0015 sec
    -----------------------------------------------------------
  ```

## Questions/Comments:
email: wkliao@eecs.northwestern.edu

Copyright (C) 2018, Northwestern University.

See [COPYRIGHT](../COPYRIGHT) notice in top-level directory.

