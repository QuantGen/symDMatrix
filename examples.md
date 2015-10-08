### (1) Converting a RAM matrix into a symDMatrix

```R

 dir.create("~/test")
 setwd("~/test")
 library(LinkedMatrix)
 library(BGData)
 source('~/GitHub/symDMatrix/definitions.r')
 
 # loading genotypes from a mice data set
  library(BGLR)
  data(mice)
 
 # Computing a G-matrix (in ram)
  X=mice.X
  p=ncol(X);n=nrow(X)
  rownames(X)=paste0('ID_',1:nrow(X))
  G=tcrossprod(scale(X))
  G<-G/mean(diag(G))
  

 # Converting G into a symDMatrix
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

### (3) Creating a symDMatrix from ff files containing the blocks

```R
 # 1st, let's create the blocks
 
 nBlocks=3
 dir.create('ff_files')
 setwd('ff_files')

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
       Gij=as.ff(G[i_ini:i_end,j_ini:j_end],file=paste0('data_',i,'_',j,'.bin'),vmode='double')
       save(Gij,file=paste0('data_',i,'_',j,'.ff'))
       print(paste(i_ini,i_end,' ; ', j_ini,j_end))
      }
     }
   }
 }
 
 
 ## code for new
 
 counter=1
 fileList=list.files(pattern='*.ff')
 dataList=list()
 n_chunks=(-1+sqrt(1+4*2*length(fileList)))/2
 for(i in 1:n_chunks){
    dataList[[i]]<-list()
    for(j in i:n_chunks){
        load(fileList[counter])
        counter=counter+1
        dataList[[i]][[j-i+1]]<-Gij
        rm(Gij)
    }
 }
 G4<-new('symDMatrix',names=rownames(G),centers=colMeans(X),
         scales=apply(X=X,MARGIN=2,FUN=sd)*sqrt((n-1)/n),data=dataList)
  all.equal(diag(G),diag(G4)) 
```
