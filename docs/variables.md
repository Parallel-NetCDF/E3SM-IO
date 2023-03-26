## Variables and Their Decomposition Maps Used in E3SM
This E3SM I/O kernel benchmark includes three common simulation cases, namely
F, G, and I cases. Both F and I cases create two history files, namely 'h0' and
'h1'.

|                 Case/History |    F/h0 |    F/h1 |     G  |   I/h0 |   I/h1 |
|-----------------------------:|--------:|--------:|-------:|-------:|-------:|
|           # vars use Map 1   |       1 |      1  |      6 |    465 |    465 |
|           # vars use Map 2   |     323 |      22 |      2 |     75 |     69 |
|           # vars use Map 3   |      63 |       1 |     25 |      4 |      2 |
|           # vars use Map 4   |         |         |      2 |      1 |      1 |
|           # vars use Map 5   |         |         |      2 |      1 |      1 |
|           # vars use Map 6   |         |         |      4 |        |        |
| Total # partitioned vars     |     387 |      24 |     41 |    546 |    538 |
| Total # non-partitioned vars |      27 |      27 |     11 |     14 |     14 |
|       Total # variables      |     414 |      51 |     52 |    560 |    552 |
|       MAX # noncontig reqs   | 184,644 | 173,206 | 21,110 | 41,400 | 38,650 |
|       MIN # noncontig reqs   | 174,926 | 164,090 | 18,821 | 33,120 | 30,920 |

`MAX # noncontig reqs` is the maximum number of noncontiguous write requests among all processes.

`MIN # noncontig reqs` is the minimum number of noncontiguous write requests among all processes.

## F case
* **F case** uses 3 decomposition maps and produces 2 files.
  + **h0** file: There are 414 climate variables stored in the 'h0' file.
    + Among the 414 variables, 6 are scalar variables and 408 are array variables.
    + Among the 414 variables, 15 are fixed-size variables and 399 are record variables.
      * Among 15 fixed-size variables,
        + 5 are scalars of type int
        + 1 is scalar of type double
        + 9 are 1D array of type double
      * Among 399 record variables (the most significant dimension is time.)
        + 5 are 1D array of type int
        + 8 are 1D array of type double, 3 of them are partitioned
        + 2 are 2D array of type char
        + 321 are 2D array of type float, all are partitioned
        + 63 are 3D array of type float, all are partitioned
    + Among the 414 variables, 387 are partitioned and 27 are not.
    + Among the 387 partitioned variables,
      * 1 uses decomposition map 0,
      * 323 use decomposition map 1, and
      * 63 use decomposition map 2.
  + **h1** file: There are 51 climate variables stored in the 'h1' file.
    * Among the 51 variables, 6 are scalar variables and 45 are array variables.
    * Among the 51 variables, 15 are fixed-size variables and 36 are record variables.
      * Among 15 fixed-size variables,
        + 5 are scalars of type int
        + 1 is scalar of type double
        + 9 are 1D arrays of type double, 3 are partitioned
      * Among 36 record variables (the most significant dimension is time.)
        + 5 are 1D arrays of type int
        + 8 are 1D arrays of type double
        + 2 are 2D arrays of type char
        + 20 are 2D arrays of type float, all are partitioned
        + 1 is 3D array of type float, partitioned
    * Among the 51 variables, 24 are partitioned and 27 are not partitioned.
    * Among the 24 partitioned variables,
      * 1 uses decomposition map 0,
      * 22 use decomposition map 1, and
      * 1 uses decomposition map 2.

## G case
* **G case** uses 6 decomposition maps and produces one file.
  + There are 52 climate variables stored in the output file.
  + All 52 variables are array variables. None is scalar.
  + Among the 52 variables, 11 are fixed-size variables and 41 are record variables.
    * Among 11 fixed-size variables,
      + 3 are 1D arrays type int, all are partitioned
      + 3 are 2D arrays type int, all are partitioned
      + 5 are 1D arrays of type double, 1 is partitioned
    * Among 41 record variables (the most significant dimension is time.)
      + 1 is 2D array of type char
      + 6 is 1D arrays of type double
      + 4 are 2D arrays of type double, all are partitioned
      + 30 are 3D arrays of type double, all are partitioned
  + Among the 52 variables, 41 are partitioned and 11 are not.
  + Among the 41 partitioned variables,
    * 6 use decomposition map 0,
    * 2 use decomposition map 1,
    * 25 use decomposition map 2,
    * 2 use decomposition map 3,
    * 2 use decomposition map 4, and
    * 4 use decomposition map 5.

## I case
* **I case** uses 5 decomposition maps and produces 2 files.
  + **h0** file: There are 560 climate variables stored in the 'h0' file.
    + All 560 variables are array variables. None is scalar.
    + Among the 560 variables, 18 are fixed-size variables and 542 are record variables.
      * Among 18 fixed-size variables,
        + 2 are 2D arrays type int, both are partitioned
        + 5 are 1D arrays of type float
        + 3 are 2D arrays of type float, all are partitioned
        + 8 are 3D arrays of type float, all are is partitioned
      * Among 542 record variables (the most significant dimension is time.)
        + 2 are 2D arrays of type char
        + 5 are 1D arrays of type int
        + 1 is 1D arrays of type float
        + 1 is 2D array of type double
        + 73 are 4D arrays of type float, all are partitioned
        + 460 are 3D arrays of type float, all are partitioned
    + Among the 560 variables, 546 are partitioned and 14 are not.
    + Among the 546 partitioned variables,
      * 465 use decomposition map 0,
      * 75 use decomposition map 1,
      * 4 use decomposition map 2,
      * 1 uses decomposition map 3, and
      * 1 uses decomposition map 4.
  + **h1** file: There are 552 climate variables stored in the 'h1' file.
    + All 552 variables are array variables. None is scalar.
    + Among the 552 variables, 10 are fixed-size variables and 542 are record variables.
      * Among 10 fixed-size variables,
        + 2 are 2D arrays type int, both are partitioned
        + 5 are 1D arrays of type float
        + 3 are 2D arrays of type float, all are partitioned
      * Among 542 record variables (the most significant dimension is time.)
        + 2 are 2D arrays of type char
        + 5 are 1D arrays of type int
        + 1 is 1D arrays of type float
        + 1 is 2D array of type double
        + 73 are 4D arrays of type float, all are partitioned
        + 460 are 3D arrays of type float, all are partitioned
    + Among the 552 variables, 538 are partitioned and 14 are not partitioned.
    + Among the 538 partitioned variables,
      * 465 use decomposition map 0,
      * 69 use decomposition map 1,
      * 2 use decomposition map 2,
      * 1 uses decomposition map 3, and
      * 1 uses decomposition map 4.

## The decomposition maps
* The decomposition map files from the E3SM production runs are available upon request.
  + `piodecomp21600tasks_F_case.nc` (266 MB) for F case produced from 21600
    processes.
  + `GMPAS-NYF_T62_oRRS18to6v3_9600p.nc` (303 MB) for G case produced from 9600
    processes.
  + `i_case_1344p.nc` (12 MB)for I case produced from 1344 processes.

* The partitioning dimension sizes of each decomposition map are given below.

| Case  |           F |              G  |               I |
|------:|:-----------:|:---------------:|:---------------:|
| Map 1 |      48,602 |  3,693,225      |       360 x 720 |
| Map 2 |      48,602 | 11,135,652      |  15 x 360 x 720 |
| Map 3 | 72 x 48,602 |  3,693,225 x 80 |  10 x 360 x 720 |
| Map 4 |             | 11,135,652 x 80 |  19 x 360 x 720 |
| Map 5 |             |  7,441,216 x 80 |  17 x 360 x 720 |
| Map 6 |             |  3,693,225 x 81 |                 |

