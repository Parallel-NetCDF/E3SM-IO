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
* See [INSTALL.md](./INSTALL.md)

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

