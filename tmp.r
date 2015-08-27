getG.symDMatrix2=function(X,nChunks=5,chunkSize=NULL,centers=NULL, scales=NULL,centerCol=T,scaleCol=T,
						folder=randomString(5),vmode='single',verbose=TRUE,saveRData=TRUE,nCores=4){
 
 	if(nCores>1) library(snow) 
 	
 	
    timeIn=proc.time()[3]
	n<-nrow(X)
	p<-ncol(X)
	if(is.null(chunkSize)) nChunks=ceiling(n/nChunks)
	
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
	if(!scaleCol)  scales<-rep(1,p)
	 
	chunkID=ceiling(1:n/chunkSize)
	nChunks=max(chunkID)
	nFiles=nChunks*(nChunks+1)/2
	DATA=list()
	counter=1

	tmpDir=getwd()
    dir.create(folder)
    setwd(folder)
    
    ## Cluster applly here... 
    
    CL <- makeSOCKcluster(rep("localhost",nCores))
    clusterExport(CL, "X")

	TMP=clusterApply(CL, 1:((nChunks*(nChunks+1)/2)), getGij,X=X,scales=scales,centers=centers,nChunks=nChunks )

nChunks=10
system.time(TMP<-clusterApply(CL, 1:((nChunks*(nChunks+1)/2)), getGij,X=X,scales=scales,centers=centers,nChunks=nChunks ))
system.time(G<-tcrossprod(scale(X)))

