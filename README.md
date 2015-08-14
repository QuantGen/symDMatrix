## symMatrix

### Memory-mapped distributed symmetric matrix


Class: ```symMatrix``` 

Slots: 
      - names (character)
      - data (list) each element of the list is an ff object
      - centers (numeric) column-means used in the computation of the matrix
      - scales (numeric) column-standard deviations used to scale the matrix.

#### Examples:

```R
 source('~/GitHub/symDMatrix/definitions.r')
 library(BGLR)
 data(wheat)
 A=as.symDMatrix(x=wheat.A,nChunks=10,folder='wheatA')
 round(wheat.A[1:2,1:2],6)==round(A@data[[1]][1:2,1:2],6)

```
