
library(ff)

setClass('symDMatrix',slots=c(names='character',data='list') )

chunks.symDMatrix<-function(x){
    if(class(x)!=){ stop(' the input must be a symDMatrix object.') }
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
    		OUT[i,2]=OUT[i,1]+nrow(x@data[[slotNum]]-1
    	}
	}
	return(OUT)
}

rownames.symDMatrix(x)<-function(x){
	return(x@names)
}

colnames.symDMatrix(x)=rownames.symDMatrix


as.symDMatrix(x,nChunks=3,vmode='double'){
	n=nrow(x)
	if(ncol(x)!=n){ stop('x must by a symmetric matrix') }
	chunkSize=ceiling(n/nChunks)
	
	## Determining chunk size
	TMP=matrix(nrow=nChunks,ncol=3)
	TMP[,1]<-1:nChunks
	TMP[1,1]<-1
	TMP[1,2]<-chunkSize
	if(nChunks>1){
		TMP[i,1]=TMP[(i-1),2]+1
		TMP[i,2]=min(TMP[i,1]+chunkSize-1, n)
	}
	
	ini=1
	end=0
	OUT=list()
	
	counter=1
	for(i in 1:nChunks){
		rowIndex=eval(parse(text=paste0(TMP[i,1],":",TMP[i,2])))
		for(j in i:nChunks){
			colIndex=eval(parse(text=paste0(TMP[j,1],":",TMP[j,2])))
			OUT[[counter]]=ff(dim=c(length(rowIndex),length(colIndex)),vmode=vmode,initdata=as.vector(x[rowIndex,colIndex]))
		}
	}
}

randomString<-function(){
    paste(sample(c(0:9,letters,LETTERS),size=5,replace=TRUE),collapse="")
}
