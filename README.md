## symDMatrix

### Memory-mapped distributed symmetric matrix

**Contact**: Gustavo de los Campos (gdeloscampos@gmail.com), Paulino Perez-Rodriguez (perpdgo@gmail.com  )

**Class**: symMatrix 

**Slots**:

     * names (character)
     * data (list) each element of the list is an ff object
     * centers (numeric) column-means used in the computation of the matrix
     * scales (numeric) column-standard deviations used to scale the matrix.

#### Examples

``` R
 # loading genotypes from a mice data set
  library(BGLR)
  data(mice)

 # Computing a G-matrix (in ram)
  X=mice.X
  p=ncol(X);n=nrow(X)
  rownames(X)=paste0('ID_',1:nrow(X))
  G=tcrossprod(scale(X))/p

 # Converting G into a symDMatrix
   G2=as.symDMatrix(G,folder="mice",nChunks=5)

 # Basic operators applied to a matrix in RAM and to a symDMatrix
  # dimension operators
   dim(G)==dim(G2)
   nrow(G)==nrow(G2)
   ncol(G)==ncol(G2)
   plot(diag(G),diag(G2))

  # names operators
   all(dimnames(G)[[1]]==dimnames(G2)[[1]])
   all(dimnames(G)[[2]]==dimnames(G2)[[2]])
   all(rownames(G)==rownames(G2))
   all(colnames(G)==rownames(G2))

 # Chunk operators
   nChunks(G2)
   chunkSize(G2)
   chunks.symDMatrix(G2)

 # Indexing (can use booleans, integers or labels)
   G2[1:2,1:2]
   G2[c("ID_1","ID_2"),c("ID_1","ID_2")]
   tmp=c(T,T,rep(F,nrow(G2)-2))
   G2[tmp,tmp]


 # Checking indexing
   for(i in 1:100){
      n1=sample(1:50,size=1)
      n2=sample(1:50,size=1)
      rows=sample(1:n,size=n1)
      cols=sample(1:n,size=n2)

      TMP1=G[rows,cols]
      TMP2=G2[rows,cols]
      stopifnot(round(cor(as.vector(TMP1),as.vector(TMP2)),5)==1)
      print(i)
   }

  ## Computing a G matrix (returns a symDMatrix)
   G3=getG(X,scaleCol=T,centerCol=T,folder='tmp',chunkSize=300)
   class(G3)
   cor(diag(G),diag(G3))
   
```


#### Pending

     - Replacement: function to update content of the matrix
     - Add chunk:   function to add one chunk
     - chol:        A recursive method to compute a cholesky decomposition
     - updateChol:  A method for updating a cholesky when a chunk is added
     - Parallel: modify some mehtods (getG, as.symDMatrix) to run in parallel
     
     
