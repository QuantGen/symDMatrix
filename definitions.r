
library(ff)

setOldClass('ff_matrix')


`colnames<-.symDMatrix`<-function(x,value) {
    x@names=values
    return(x)
}

setClass('symDMatrix',slots=c(names='character',centers='numeric',
                              scales='numeric',data='list') )

nChunks=function(x) length(x@data[[1]])
chunkSize=function(x) nrow(x@data[[1]][[1]])

setMethod('rownames',signature='symDMatrix',definition=function(x) x@names)
setMethod('colnames',signature='symDMatrix',definition=function(x) x@names)

nrow.symDMatrix<-function(x) length(x@names)
setMethod('nrow',signature='symDMatrix',definition=nrow.symDMatrix)
setMethod('ncol',signature='symDMatrix',definition=nrow.symDMatrix)
setMethod('dim',signature='symDMatrix',definition=function(x) rep(nrow(x),2) )

diag.ff<-function(x){
	if(class(x)[1]!='ff_matrix'){ stop('x must be an ff_matrix object')}
	n=min(dim(x))
	out=rep(NA,n)
	for(i in 1:n){
		out[i]=x[i,i]
	}
	return(out)
}
setMethod('diag',signature='ff_matrix',definition=diag.ff)


diag.symDMatrix<-function(x){
	n=min(dim(x))
	out=rep(NA,n)
	
	nChunks=nChunks(x)
	end=0
	for(i in 1:nChunks){
		tmp=diag(x@data[[i]][[1]])
		ini<-end+1
		end<-ini+length(tmp)-1
		out[ini:end]<-tmp
	}
	names(out)=x@names
	return(out)
}
setMethod('diag',signature='symDMatrix',definition=diag.symDMatrix)

as.symDMatrix<-function(x,nChunks=3,vmode='single',folder=randomString(),saveRData=TRUE){
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
		DATA[[i]]=list()
		for(j in i:nChunks){
			colIndex=eval(parse(text=paste0(TMP[j,2],":",TMP[j,3])))
			k=j-i+1
			DATA[[i]][[k]]=ff( dim=c(length(rowIndex),length(colIndex)),
			                   vmode=vmode,initdata=as.vector(x[rowIndex,colIndex]),
			                   filename=paste0('data_',i,'_',j,'.bin')
			                 )
			colnames(DATA[[i]][[k]])<-colnames(x)[colIndex]
			rownames(DATA[[i]][[k]])<-rownames(x)[rowIndex]
		}
	}
	setwd(tmpDir)
	out=new('symDMatrix',names=rownames(x),data=DATA,centers=0,scales=0)
	if(saveRData){save(out,file='symDMatrix.RData') }
	return(out)
}

chunks.symDMatrix<-function(x){
    if(class(x)!='symDMatrix'){ stop(' the input must be a symDMatrix object.') }
    
    n=length(x@data)
    OUT=matrix(nrow=n,ncol=3)
    OUT[,1]<-1:n    
    colnames(OUT)<-c('chunk','ini','end')
    end=0
    for(i in 1:n){
    	ini=end+1
    	end=ini+nrow(x@data[[i]][[1]])-1
    	OUT[i,2]=ini
    	OUT[i,3]=end
    	ini=end+1
    }
    return(OUT)
}

setMethod('chunks',signature='symDMatrix',definition=chunks.symDMatrix)

subset.symDMatrix=function(x,i,j){
  #i=sample(1:nrow(x),size=5)
  #j=sample(1:ncol(x),size=3)
  tmpClass<- class(i)
   if(tmpClass=='factor')    i=as.character(i)
   if(tmpClass=='character') i=which(x@names%in%i)
   if(tmpClass=='logical')   i=which(i)
   if(tmpClass=='numeric')   i=as.integer(i)
   
 tmpClass<- class(j)
  if(tmpClass=='factor')    j=as.character(j)
  if(tmpClass=='character')	j=which(x@names%in%j)
  if(tmpClass=='logical')	j=which(j)
  if(tmpClass=='numeric')	j=as.integer(j)
  
 nChunks=nChunks(x)
 chunkSize=ncol(x@data[[1]][[1]])
 i0=i
 j0=j
  
 i=rep(i0,each=length(j0))
 j=rep(j0,times=length(i0))
 
 # switching indexes ( i must be smaller than j)
 tmp=i>j
 tmp2=i[tmp]
 i[tmp]=j[tmp]
 j[tmp]=tmp2

 out.i=rep(1:length(i0),each=length(j0))
 out.j=rep(1:length(j0),length(i0))
 
 row.chunk=ceiling(i/chunkSize)
 col.chunk=ceiling(j/chunkSize)
 local.i=i-(row.chunk-1)*chunkSize
 local.j=j-(col.chunk-1)*chunkSize
 
 cbind(i,j,row.chunk,col.chunk,local.i,local.j,out.i,out.j)



 OUT=matrix(nrow=length(i0),ncol=length(j0),NA)
 rownames(OUT)=x@names[i0]
 colnames(OUT)=x@names[j0]

 for(i in unique(row.chunk)){
	tmp=which(row.chunk==i)
	for(j in unique(col.chunk[tmp])){
	   k=which(row.chunk==i & col.chunk==j)
	   tmp.row.in=local.i[k]
	   tmp.col.in=local.j[k]
 	   tmp.in=(tmp.col.in-1)*nrow(x@data[[i]][[j-i+1]])+tmp.row.in
 
 	   tmp.row.out=out.i[k]
 	   tmp.col.out=out.j[k]
 	   tmp.out=(tmp.col.out-1)*nrow(OUT)+tmp.row.out
 	   OUT[tmp.out]=x@data[[i]][[j-i+1]][tmp.in]
	   
	}
 }
 return(OUT)
}

setMethod('[',signature='symDMatrix',definition=subset.symDMatrix)





getG=function(X,chunkSize=100,centers=NULL, scales=NULL,centerCol=T,scaleCol=T,folder=randomString(5),vmode='single',
				verbose=TRUE,saveRData=TRUE){
 
    timeIn=proc.time()[3]
	n<-nrow(X)
	p<-ncol(X)
	
	if( (centerCol|scaleCol)&(is.null(centers)|is.null(scales))){		
		if(is.null(centers)&is.null(scales)){
			centers=rep(NA,p); scales=rep(NA,p)	
			for(i in 1:p){	
				xi=X[,i]
				scales[i]=sd(xi,na.rm=TRUE)*sqrt((n-1)/n)*sqrt(p)
				centers[i]=mean(xi,na.rm=TRUE)
			}
		}
		if(is.null(centers)&(!is.null(scales))){
			scales=rep(NA,p)	
			for(i in 1:p){	
				xi=X[,i]
				scales[i]=sd(xi,na.rm=TRUE)*sqrt((n-1)/n)*sqrt(p)
			}
		}
		if((!is.null(centers))&is.null(scales)){
			centers=rep(NA,p)
			for(i in 1:p){	
				xi=X[,i]
				centers[i]=mean(xi,na.rm=TRUE)
			}
		}		
	}
	if(!centerCol) centers<-rep(0,p)
	if(!scaleCol)  scales<-rep(0,p)
	 
	chunkID=ceiling(1:n/chunkSize)
	nChunks=max(chunkID)
	nFiles=nChunks*(nChunks+1)/2
	DATA=list()
	counter=1

	tmpDir=getwd()
    dir.create(folder)
    setwd(folder)
     
	for(i in 1:nChunks){
		DATA[[i]]=list()
    	rowIndex_i=which(chunkID==i)
    	Xi=X[rowIndex_i,] 
    	
    	# centering/scaling
    	for(k in 1:p){
    		xik=Xi[,k]
    		xik=(xik-centers[k])/scales[k]
    		xik[is.na(xik)]=0
    		Xi[,k]=xik
    	}
        
        for(j in i:nChunks){
            rowIndex_j=which(chunkID==j)
    		Xj=X[rowIndex_j,] 
    		
			# centering/scaling
	    	for(k in 1:p){
    			xjk=Xj[,k]
    			xjk=(xjk-centers[k])/scales[k]
    			xjk[is.na(xjk)]=0
    			Xj[,k]=xjk
    		}            
            
            Gij=tcrossprod(Xi,Xj)
            
            DATA[[i]][[j-i+1]]=ff(dim=dim(Gij),
                                  vmode=vmode,initdata=as.vector(Gij),
                                  filename=paste0('data_',i,'_',j,'.bin')
                                 )
            colnames(DATA[[i]][[j-i+1]])<-colnames(X)[rowIndex_j]
            rownames(DATA[[i]][[j-i+1]])<-rownames(X)[rowIndex_i]
            counter=counter+1
            if(verbose){ cat(' Done with pair ', i,'-',j,' (',round(100*counter/(nChunks*(nChunks+1)/2)),'% ', round(proc.time()[3]-timeIn,3),' seconds).\n',sep='')}
        }
    }
    setwd(tmpDir)
    out=new('symDMatrix',names=rownames(X),data=DATA,centers=centers,scales=scales)
	if(saveRData){save(out,file='G.RData') }
	return(out)
}

randomString<-function(n=10)    paste(sample(c(0:9,letters,LETTERS),size=n,replace=TRUE),collapse="")

