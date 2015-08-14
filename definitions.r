
library(ff)

setClass('symDMatrix',slots=c(names='character',data='list') )

nChunks=function(A) (-1+sqrt(1+8*length(A@data)))/2
chunkSize=function(A) nrow(A@data[[1]])
rownames.symDMatrix<-function(x) x@names
colnames.symDMatrix=rownames.symDMatrix

 
as.symDMatrix<-function(x,nChunks=3,vmode='single',folder=randomString()){
	#nChunks=3;vmode='single';folder=randomString(5)
	n=nrow(x)
	if(ncol(x)!=n){ stop('x must by a square matrix') }

	tmpDir=getwd()
	dir.create(folder)
	setwd(folder)
	
	chunkSize=ceiling(n/nChunks)
	
	## Determining chunk size and subjects in each chunk
	TMP=matrix(nrow=nChunks,ncol=3)
	TMP[,1]<-1:nChunks
	TMP[1,2]<-1
	TMP[1,3]<-chunkSize
	if(nChunks>1){
		for(i in 2:nChunks){
			TMP[i,2]=TMP[(i-1),3]+1
			TMP[i,3]=min(TMP[i,2]+chunkSize-1, n)
		}
	}
	ini=1
	end=0
	DATA=list()
	
	counter=1
	for(i in 1:nChunks){
		rowIndex=eval(parse(text=paste0(TMP[i,2],":",TMP[i,3])))
		for(j in i:nChunks){
			colIndex=eval(parse(text=paste0(TMP[j,2],":",TMP[j,3])))
			DATA[[counter]]=ff(dim=c(length(rowIndex),length(colIndex)),
			                   vmode=vmode,initdata=as.vector(x[rowIndex,colIndex]),
			                   filename=paste0('data_',i,'_',j,'.bin')
			                   )
			counter=counter+1
		}
	}
	setwd(tmpDir)
	out=new('symDMatrix',names=rownames(x),data=DATA)
	return(out)
}

chunks.symDMatrix<-function(x){
    if(class(x)!='symDMatrix'){ stop(' the input must be a symDMatrix object.') }
    nFiles=length(x@data)
    n=sqrt(1+8*nFiles)/2 -1
    OUT=matrix(nrow=n,ncol=3)
    OUT[,1]<-1:n
    OUT[1,1]=1
    OUT[1,2]=nrow(x@data[[1]])
    
    colnames(OUT)<-c('chunk','ini','end')
    slotNum=1
    if(n>1){
	    for(i in 2:n){
    		slotNum=slotNum+(n-i+1)
    		OUT[i,1]=OUT[(i-1),2]+1
    		OUT[i,2]=OUT[i,1]+nrow(x@data[[slotNum]])-1
    	     }
	}
      return(OUT)
}


randomString<-function(n=10)    paste(sample(c(0:9,letters,LETTERS),size=n,replace=TRUE),collapse="")


if(FALSE){
  ## Code for subseting
i=sample(1:599,size=2)
j=sample(1:599,size=2)

j=ifelse(j<=i,i,j)

chunk.i=ceiling(i/chunkSize(A))
chunk.j=ceiling(j/chunkSize(A))

local.i=i-(chunk.i-1)*chunkSize(A)
local.j=j-(chunk.j-1)*chunkSize(A)

INDEX=cbind(i,j,chunk.i,chunk.j,local.i,local.j,out.i=1:length(i),out.j=1:length(j))


uniqueSets=unique(apply(X=INDEX[,c('chunk.i','chunk.j')],FUN=paste,MARGIN=1,collapse='-'))
uniqueSets=matrix(data=as.integer(unlist(strsplit(uniqueSets,split='-'))),ncol=2,byrow=T) 
uniqueSets=uniqueSets[order(uniqueSets[,1],uniqueSets[,2]),]
colnames(uniqueSets)<-c('rowChunk','colChunk')
 
OUT=matrix(nrow=length(i),ncol=length(j),NA)
colnames(OUT)<-colnames.symDMatrix(A)[j]
rownames(OUT)<-rownames.symDMatrix(A)[i]


## Have to fix this

for(k in 1:nrow(uniqueSets)){
	tmp=(INDEX[,'chunk.i']==uniqueSets[k,'rowChunk'])&(INDEX[,'chunk.j']==uniqueSets[k,'colChunk'])
	
	previousRow=uniqueSets[k,1]-1
	filesPreviousRow=previousRow*nChunks(A)-(previousRow)*(previousRow-1)/2
	filesInRow=uniqueSets[k,'colChunk']-(uniqueSets[k,'rowChunk']-1)
	fileNumber=filesPreviousRow+filesInRow
	
	print(c(filesPreviousRow,filesInRow,fileNumber))          
	
	OUT[INDEX[tmp,'out.i'],INDEX[tmp,'out.i']]=A@data[[fileNumber]][INDEX[tmp,'local.i'],INDEX[tmp,'local.j']]
}

}
