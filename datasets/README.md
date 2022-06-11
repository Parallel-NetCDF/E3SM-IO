
* Dimensions in Decomposition Maps
  |    | F case     | G case                  | I case                           |
  |:--:|-----------:|------------------------:|---------------------------------:|
  | D1 | ncol       | nCells                  | lat x lon                        |
  | D2 | ncol       | nEdges                  | (levgrnd or levdcmp) x lat x lon |
  | D3 | lev x ncol | nCells x nVertLevels    | levlak x lat x lon               |
  | D4 |            | nEdges x nVertLevels    | ltype x lat x lon                |
  | D5 |            | nVertices x nVertLevels | natpft x lat x lon               |
  | D6 |            | nCells x nVertLevelsP1  |                                  |

* Dimension sizes of decomposition maps available in this folder
  | f_case_866x72_16p | g_case_cmpaso_16p   | i_case_f19_g16_16p |
  |:------------------|:--------------------|:-------------------|
  | ncol = 866        | nCells = 28571      | lat = 96           |
  | lev = 72          | nEdges = 87980      | lon = 144          |
  |                   | nVertLevels = 100   | levgrnd = 15       |
  |                   | nVertices = 59329   | levdcmp = 15       |
  |                   | nVertLevelsP1 = 101 | levlak = 10        |
  |                   |                     | ltype = 9          |
  |                   |                     | natpft = 17        |

* Decomposition Map's dimension sizes used in production runs (standard resolution)
  | F case 21600p | G case 9600p        | I case 1344p |
  |:--------------|:--------------------|:-------------|
  | ncol = 48602  | nCells = 3693225    | lat = 360    |
  | lev = 72      | nEdges = 11135652   | lon = 720    |
  |               | nVertLevels = 80    | levgrnd = 15 |
  |               | nVertices = 7441216 | levdcmp = 15 |
  |               | nVertLevelsP1 = 81  | levlak = 10  |
  |               |                     | ltype = 9    |
  |               |                     | natpft = 17  |

