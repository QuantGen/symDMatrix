### (1) Converting a RAM matrix into a symDMatrix

```R

 dir.create("~/test")
 setwd("~/test")
 library(BGData)
 source('~/GitHub/symDMatrix/definitions.r')
 
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
```

### (2) Exploring operators
```R
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
   head(G2[tmp,])
  
 # Exhaustive check of indexing
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
```

### (3) Creating a symDMatrix from ff files containing the blocks

```R
 # 1st, let's create the blocks
 
 nBlocks=10
 
 n<-nrow(G)
 stepSize<-ceiling(n/nBlocks)
 
 for(i in 1:nBlocks){
   i_ini=(i-1)*stepSize+1
   if(i_ini<=n){
     i_end<-min(n,i_ini+stepSize-1)
     for(j in i:nBlocks){
      j_ini<-(j-1)*stepSize+1
      if(j_ini<=n){
       j_end<-min(n,j_ini+stepSize-1)
       Gij=as.ff(G[i_ini:i_end,j_ini:j_end],file=paste0('data_',i,'_',j,'.bin'))
       save(Gij,file=paste0('data_',i,'_',j,'.ff'))
       print(paste(i_ini,i_end,' ; ', j_ini,j_end))
      }
     }
   }
 }
 
 ## code for new
 
 

counter=0
dataList=list()
nChunks=(-1+sqrt(1+4*2*length(fileList)))/2
blockID=rep(1:Chunks,each=ceiling(n/nBlocks))[1:n]


for(i in 1:nChunks){
	nRow=sum(blockID==i)
    dataList[[i]]<-list()
    for(j in i:nChunks){
    	nCol=sum(blockID==j)
    	counter=counter+1
    	dataList[[i]][[j-i+1]]<-ff(dim=c(nRow,nCol),filename=fileList[counter],vmode='double')
       
	}
 }
 
```
