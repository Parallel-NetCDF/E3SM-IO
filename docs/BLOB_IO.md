## Design of Blob I/O Strategies

- [PnetCDF Blob I/O](#pnetcdf-blob-io)
- [HDF5 Blob I/O](#hdf5-blob-io)

### PnetCDF Blob I/O
* This I/O operation is triggered when the E3SM-IO command-line options "-a
  pnetcdf -x blob" are used.
* Additional global attributes, dimensions, and variables will be created to
  describe the data decomposition and layout, including the number of MPI
  processes used to create the file, decomposition maps, and flags indicating
  the map used by a variable.
  See [src/cases/header_def_F_case.cpp](../src/cases/header_def_F_case.cpp).
  These additional data objects are described below.
  + Additional global attributes (See subroutine `add_gattrs()`)
    * **global_nprocs** - a 4-byte integer storing the number of MPI processes
      used when creating the file.
    * **num_decompositions** - a 4-byte integer storing the number of
      decomposition maps.
    * **num_subfiles** - a 4-byte integer storing the number of subfiles. MPI
      processes are divided into groups and processes in each group access only
      one subfile. The group division is exclusive.
  + Additional dimensions (See subroutine `e3sm_io_case::def_F_case()`)
    * **nblobs** - number of blobs (number of processes sharing a subfile. Same
      value as global_nprocs when num_subfiles is 1).
    * There are 2 additional dimensions per decomposition map.
      + **D1.nelems** - number of array elements in decomposition map D1 in a
        subfile. See subroutine `blob_metadata()` in
        [calc_metadata.c](../src/calc_metadata.c). The sum of this dimension
        across all subfiles is equal to the size of original dimension
        decomposed by map D1.
      + **D1.max_nreqs** - max number of flattened noncontiguous requests
        (offset-length pairs) among processes sharing a subfile for
        decomposition map D1. Note the numers of noncontiguous requests
        assigned to processes can be different. See subroutine
        `blob_metadata()` in [calc_metadata.c](../src/calc_metadata.c).
      + If there are more decomposition maps, the additional map dimensions
        will be D2.nelems, D2.max_nreqs, D3.nelems, D3.max_nreqs, and so on.
  + Additional variables and their attributes (See subroutine
    `e3sm_io_case::def_var_decomp()`)
    * There are 5 additional variables per decomposition map. Below uses an
      example of map D3.
      + **int D3.nreqs(nblobs)**
        * D3.nreqs:description = "Number of noncontiguous requests in
          individual blobs"
        * D3.nreqs:global_dimids = 3, 5 ;
        * **Note**
          + `nblobs` is the number of blobs (processes) sharing this subfile.
          + Each element of D3.nreqs is the number of noncontiguous requests in
            offset-length pairs assigned to a process.
          + Attribute global_dimids contains the IDs of global dimensions.  In
            this example, map D3 decomposes along global dimensions 3 and 5
            (`lev` and `ncol` respectively).
          + Global dimensions are referred to the original dimensions
            regardless of subfiling.
      + **int64 D3.blob_start(nblobs)**
        * D3.blob_start:description = "Starting array indices of individual
          blobs stored in a variable"
      + **int64 D3.blob_count(nblobs)**
        * D3.blob_count:description = "Number of contiguous array elements in
          individual blobs stored in a variable"
      + **int D3.offsets(nblobs, D3.max_nreqs)**
        * D3.offsets:description = "Starting indices of flattened canonical
          noncontiguous requests of individual blobs"
        * **Note**
          + As the numers of noncontiguous requests assigned to processes can
            be different, the 2nd dimension of this 2D array actually has
            staggered lengths.
      + **int D3.lengths(nblobs, D3.max_nreqs)**
        * D3.lengths:description = "Number of elements of flattened canonical
          noncontiguous requests of individual blobs"
        * **Note**
          + Similar to D3.offsets, the 2nd dimension of this 2D array has
            staggered lengths.
    * Note that the contents of decomposition variables are not used when
      writing the climate variables. Thus, they can be defined together with
      climate variables in the same define mode and written in the same data
      mode as climate variables.
* Changes in variable definitions (See C macro `DEF_VAR` in
  [src/cases/e3sm_io_case.hpp](../src/cases/e3sm_io_case.hpp).)
  + The dimensions of a decomposed variable are changed to use decomposition
    map dimensions. (See
    [src/cases/header_def_F_case.cpp](../src/cases/header_def_F_case.cpp)). For
    example, given 6 dimensions defined in a NetCDF file, variable `CLOUD` is
    originally defined as a 3D array of dimension `time` x `lev` x `ncol`.
    ```c
    time = UNLIMITED ; // (1 currently)
    nbnd = 2 ;
    chars = 8 ;
    lev = 72 ;
    ilev = 73 ;
    ncol = 866 ;

    float CLOUD(time, lev, ncol) ;
    ```
    Because variable `CLOUD` is decomposed using map D3 and D3 decomposes along
    dimension `lev` and `ncol`, definition of `CLOUD`is now changed to
    ```c
    float CLOUD(time, D3.nelems) ;
          CLOUD:decomposition_ID = 3 ;
          CLOUD:global_dimids = 0, 3, 5 ;
    ```
    The decomposed dimensions, `lev` and `ncol`, are replaced with the size of
    decomposition map ID, D3. In addition, two attributes are added per
    decomposed variable. The first one is a 4-byte integer, indicating the
    decomposition ID and the second a 4-byte integer array storing the
    (original) global dimension IDs. In this case, 3 corresponds to dimension
    `lev` and 5 `ncol`. Note this example shows dimension `time` is not
    decomposed.
* Changes of arguments `start` and `count` in PnetCDF put API calls
  + See subroutine `e3sm_io_case::var_wr_case()` in
    [src/cases/var_wr_case.cpp](../src/cases/var_wr_case.cpp) and
    `blob_metadata()` in [src/calc_metadata.c](../src/calc_metadata.c).
  + A blob is a contiguous space in file, so argument `start` and `count`
    always describe the starting offsets to a variable and number of array
    elements in the blob. `count` is first calculated based on the number of
    elements written by a process.
    ```c
    for (i=0; i<decom->num_decomp; i++) {
        for (j=0; j<decom->contig_nreqs[i]; j++)
            decom->count[i] += decom->blocklens[i][j];
    ```
    `count` is then used to calculate `start`. Note `start` of a process
    depends on the amounts decomposed to processes of lower ranks.
    ```c
    err = MPI_Exscan(decom->count, decom->start, decom->num_decomp, MPI_OFFSET,
                     MPI_SUM, cfg->sub_comm);
    ```
    These calculations are done before `ncmpi_enddef(`) is called.
  + In this design, blob I/O is for each individual variable. Thus there are
    `P` blobs for each fixed-size variable and each record of a record variable
    in a file. `P` is the number of processes sharing a subfile. In other
    words, in the file space occupied by a fixed-size variable or a record of a
    record variable, there are `P` blobs, each written by a process. This
    design is referred to as **variable-centric** data layout.
* Only nonblocking `ncmpi_iput_vara` APIs are used.
  + See subroutine `e3sm_io_case::var_wr_case()` in
    [src/cases/var_wr_case.cpp](../src/cases/var_wr_case.cpp).
  + One `ncmpi_iput_vara()` is called per variable.
  + All write requests are pending until the call to `ncmpi_wait_all`.
* **Advantages of variable-centric data layout**
  + Pending nonblocking requests can be flushed at any time. This can
    effectively reduce the memory footprint.
  + New time records can be added without paying an expensive cost of moving
    any existing data in files.

### HDF5 Blob I/O
* This I/O operation is triggered when the E3SM-IO command-line options "-a
  hdf5 -x blob" are used.
* Each HDF5 (sub)file contains only two datasets. For example,
  ```console
  % h5ls -r blob_F_out_h0.h5.0000
  /                        Group
  /data_blob               Dataset {9612472}
  /header_blob             Dataset {134552}
  ```
  Dataset `/header_blob` stores the metadata and dataset `/data_blob` stores
  all the variables, including decomposition maps and climate variables.
* The CDF-5 header format is borrowed to store all the metadata. The
  implementation of metadata operations is also borrowed from PnetCDF. HDF
  utility tools, such as h5dump, are not able to understand the metadata.
  See `src/drivers/blob_ncmpio.h` and `src/drivers/blob_ncmpio.c`.
* The same additional global attributes, dimensions, decomposition variables,
  and their attributes as the PnetCDF blob I/O are created in dataset
  `header_blob`.
* Dataset `header_blob` is written by rank 0 only.
* All write requests are cached into internally allocated buffers. See
  `e3sm_io_driver_h5blob::put_vara()` in file
  `src/drivers/e3sm_io_driver_h5blob.cpp`. All cached write data is flushed out
  to the file only when closing the file. See `e3sm_io_driver_h5blob::close()`.
  The flush makes only one call to MPI collective write function.
* The I/O pattern of this design is that each process writes to a contiguous
  file space no matter how many variables are defined and written in the file.
  In other words, there is only one blob per MPI process in the file. This
  design is referred to as **process-centric** data layout.
* In contrast to PnetCDF blob I/O design, there is one blob per process per
  variable in the file in the variable-centric data layout.
* ADIOS and its BP file format also use the process-centric data layout.
* **Drawbacks of process-centric data layout**
  + Memory footprint can be large, as all write requests are cached in
    internally allocated buffers until file close time.
  + It will be very expensive to add new time records to the existing variables
    in a previous closed file. This is because the data layout is process
    centric, i.e. all data written by a process must be packed into a blob and
    appended to another blob in the file. Expanding a process's blob is
    required to move all blobs of the processes with higher ranks to higher
    file offset locations. In PnetCDF blob I/O which uses a variable-centric
    data layout, there is no such penalty when adding new time records to
    variables.


