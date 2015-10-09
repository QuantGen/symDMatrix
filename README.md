## symDMatrix

**Contact**: Gustavo de los Campos (gdeloscampos@gmail.com)
**Class**: symMatrix 

### A memory-mapped-distributed symmetric matrix

A symmeryc matrix is partitioned into blocks as follows:

| G11 | G12 | G33 |
|:---:|:---:|:---:|
| G21 | G22 | G33 |
| G31 | G32 | G33 |

Because the matrix is assumed to be symmetric Gij=Gji; therefore, only the upper-triangular blocks are stored. Each block is stored in an ff object. The package defines the class and multiple methods that allow treating this memmory-mapped matrix as a standard RAM matrix.



**Slots**:

     * @names   (character)
     * @data    (list) each element of the list is an ff object
     * @centers (numeric) column-means used in the computation of the matrix
     * @scales  (numeric) column-standard deviations used to scale the matrix.

### (0) Creating a symmetric matrix in RAM

Before we start, let's create a symmetric matrix in RAM.

```R
 # loading genotypes from a mice data set
  library(BGLR)
  data(mice)
 
 # Computing a G-matrix (in ram)
  X=mice.X
  p=ncol(X);n=nrow(X)
  rownames(X)=paste0('ID_',1:nrow(X))
  G=tcrossprod(scale(X))
  G<-G/mean(diag(G))
```  

### (1) Converting a RAM matrix into a symDMatrix

In practice, if we can hold a matrix in RAM there is not much of a point to convert it to a symDMatrix, however, this will help us to get started.

```R
  G2=as.symDMatrix(G,folder="mice",nChunks=5,vmode='double') # can use single for lighter files.
```

### (2) Exploring operators
```R
 # Basic operators applied to a matrix in RAM and to a symDMatrix
  # dimension operators
   dim(G)==dim(G2)
   nrow(G)==nrow(G2)
   ncol(G)==ncol(G2)
   all.equal(diag(G),diag(G2))
     
  # names operators
   all(dimnames(G)[[1]]==dimnames(G2)[[1]])
   all(dimnames(G)[[2]]==dimnames(G2)[[2]])
   all(rownames(G)==rownames(G2))
   all(colnames(G)==rownames(G2))

 # Chunk operators
   nChunks(G2)
   chunkSize(G2)
   chunks(G2)
  
 # Indexing (can use booleans, integers or labels)
   G2[1:2,1:2]
   G2[c("ID_1","ID_2"),c("ID_1","ID_2")]
   tmp=c(T,T,rep(F,nrow(G2)-2))
   G2[tmp,tmp]
   head(G2[tmp,])
  
 # Exhaustive check of indexing
   for(i in 1:100){
   	  n1=sample(1:50,size=1)
   	  n2=sample(1:50,size=1)
   	  i1=sample(1:n,size=n1)
   	  i2=sample(1:n,size=n2)
   	  
   	  TMP1=G[i1,i2,drop=FALSE]
   	  TMP2=G2[i1,i2,drop=FALSE]
   	  stopifnot(all.equal(TMP1,TMP2))
   }
   
```
### (3) Creating a symDMatrix from genotypes

The function ```getG.symDMatrix``` computes G=XX' (with options for centering and scaling) without ever loading G in RAM, it creates the symDMatrix directly. In this example X is a matrix in RAM, for large genotype data sets X could be a mmemory-mapped matrix, ff object, or part of a BGData object.

```R
  G3=getG.symDMatrix(X,scaleCol=T,centerCol=T,folder='tmp',chunkSize=300,mc.cores=6,vmode='double')
  class(G3)
  all.equal(diag(G),diag(G3))
  
 for(i in 1:10){ 
  n1<-sample(1:25,size=1)
  i1<-sample(1:25,size=n1)
  for(j in 1:10){
    n2<-sample(1:nrow(G),size=1)
    i2<-sample(1:nrow(G),size=n2)
    tmp1=G[i1,i2]
    tmp2<-G3[i1,i2]
    stopifnot(all.equal(tmp1,tmp2))
  }
 }

```
### (4) Creating a symDMatrix from ff files containing the blocks.

For very large G-matrices, computation of the blocks of the symDMatrix can be done in parallel (e.g., in an HPC). The function getGij is similar to `getG()` (see [BGData](https://github.com/quantgen/bgdata) package) but accepts arguments i1 and i2 which define a block of G (i.e., rows of X).

```R
 nBlocks=3
 dir.create('ff_files')
 setwd('ff_files')

 n<-nrow(G)
 stepSize<-ceiling(n/nBlocks)
 
 ## This loop may be executed in parallel
 for(i in 1:nBlocks){
   i_ini=(i-1)*stepSize+1
   if(i_ini<=n){
     i_end<-min(n,i_ini+stepSize-1)
     for(j in i:nBlocks){
      j_ini<-(j-1)*stepSize+1
      if(j_ini<=n){
       j_end<-min(n,j_ini+stepSize-1)
       Gij=as.ff(G[i_ini:i_end,j_ini:j_end],file=paste0('data_',i,'_',j,'.bin'),vmode='double')
       save(Gij,file=paste0('data_',i,'_',j,'.ff'))
       print(paste(i_ini,i_end,' ; ', j_ini,j_end))
      }
     }
   }
 }
 ## end of blocks
 
 ## Now we create the object (centers, scales, etc can be also added)
 G5=symDMatrix(dataFiles=list.files(pattern='*.ff'),names=rownames(X))
 all.equal(diag(G5),diag(G))
```



#### Pending

     - Replacement: function to update content of the matrix.
     - Add chunk:   function to add one chunk.
     - chol:        A recursive method to compute a cholesky decomposition
     - updateChol:  A method for updating a cholesky when a chunk is added

     
     
