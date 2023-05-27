## Small-scale Decomposition Map Files for Testing

* Three small-scale cases of E3SM decomposition maps are available for testing.
  They are F, G, and I cases. The decomposition maps were generated from E3SM
  runs on 16 MPI processes.
  + F case uses 3 decomposition maps.
    + File `map_f_case_16p.nc` is in NetCDF classic CDF-5 format
    + File `map_f_case_16p.h5` is in HDF5 format
    + File `map_f_case_16p.nc4` is in NetCDF4 format
    + File `map_f_case_16p.bp` is in ADIOS BP format
  + G case uses 6 decomposition maps.
    + File `map_g_case_16p.nc` is in NetCDF classic CDF-5 format
    + File `map_g_case_16p.h5` is in HDF5 format
    + File `map_g_case_16p.nc4` is in NetCDF4 format
    + File `map_g_case_16p.bp` is in ADIOS BP format
  + I case uses 5 decomposition maps.
    + File `map_i_case_16p.nc` is in NetCDF classic CDF-5 format
    + File `map_i_case_16p.h5` is in HDF5 format
    + File `map_i_case_16p.nc4` is in NetCDF4 format
    + File `map_i_case_16p.bp` is in ADIOS BP format

* Metadata (file header) of the three cases are available in text format.
  + F case produces 2 output files, namely **h0** and **h1**.
    + File `f_case_h0.txt`
    + File `f_case_h1.txt`
  + G case produces 1 output file
    + File `g_case_hist.txt`
  + I case produces 2 output files, namely **h0** and **h1**.
    + File `i_case_h0.txt`
    + File `i_case_h1.txt`

### Dimensions in Decomposition Maps
  |    | F case     | G case                  | I case                           |
  |:--:|-----------:|------------------------:|---------------------------------:|
  | D1 | ncol       | nCells                  | lat x lon                        |
  | D2 | ncol       | nEdges                  | (levgrnd or levdcmp) x lat x lon |
  | D3 | lev x ncol | nCells x nVertLevels    | levlak x lat x lon               |
  | D4 |            | nEdges x nVertLevels    | ltype x lat x lon                |
  | D5 |            | nVertices x nVertLevels | natpft x lat x lon               |
  | D6 |            | nCells x nVertLevelsP1  |                                  |

### Dimension sizes of decomposition maps available in this folder
  | map_f_case_16p | map_g_case_16p      | map_i_case_16p |
  |:---------------|:--------------------|:---------------|
  | ncol = 866     | nCells = 285        | lat = 96       |
  | lev = 72       | nEdges = 879        | lon = 144      |
  |                | nVertLevels = 100   | levgrnd = 15   |
  |                | nVertices = 593     | levdcmp = 15   |
  |                | nVertLevelsP1 = 101 | levlak = 10    |
  |                |                     | ltype = 9      |
  |                |                     | natpft = 17    |

### Decomposition Map's dimension sizes used in production runs (standard resolution)
  | F case 21600p | G case 9600p        | I case 1344p |
  |:--------------|:--------------------|:-------------|
  | ncol = 777602 | nCells = 3693225    | lat = 360    |
  | lev = 72      | nEdges = 11135652   | lon = 720    |
  |               | nVertLevels = 80    | levgrnd = 15 |
  |               | nVertices = 7441216 | levdcmp = 15 |
  |               | nVertLevelsP1 = 81  | levlak = 10  |
  |               |                     | ltype = 9    |
  |               |                     | natpft = 17  |

* For information about the statistics of decomposition maps and variables that
  are partitioned using the maps, readers are referred to
  [../docs/variables.md](../docs/variables.md).

