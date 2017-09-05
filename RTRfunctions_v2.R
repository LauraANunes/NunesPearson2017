
library('doParallel')
library('doRNG')

library('raster')
library('dismo')
library('rgeos') #for convexhull and centroid
library('maptools')

### Trimming occurences, remove occurences which represent outliers in environmental hyperspace (Farber and Kadmon, 2003)
trim.occ<-function(bios,occ, plot=TRUE){
  S=matrix(c(extract(bios,occ)),ncol=nlayers(bios),nrow=nrow(occ))
  C=cov(S) #covariance matrix
  M=apply(S,2,mean) #mean conditions of each climate indice
  
  D<-rep(NA,nrow(S))
  for(i in 1:nrow(S)){
    
    Tr=t(S[i,]-M) #transpose opeerator
    D[i]=Tr%*%((C)^-1)%*%(S[i,]-M) #malahanobis distance
  }
  
  D<-cbind(c(1:nrow(S)),D) #add locality index to D distances
  
  qD<-D[c(which(D[,2]>quantile(D,0.05))),] #trim 5%
  qD<-qD[c(which(qD[,2]<quantile(D,0.95))),] #trim 95%
  
  occ_trim<-occ[qD[,1],] #keep locality points within 5-95th percentile
  if(plot==TRUE){
    plot(D[,2], xlab='', ylab='Mahalanobis Distance')
    abline(h=quantile(D,0.95), col='red')
    abline(h=quantile(D,0.05), col='red')
    plot(bios[[1]], col='grey')
    points(occ)
    points(occ_trim,col='red')
  }
  return(occ_trim)
}
########## COMPARE NICHE
#### new metric - niche overlap function ###
mo.metric<-function(sp1,sp2,bios){ 
  
  e_sp1<-na.omit(extract(bios,sp1)) #extract climatic variables for species 1
  e_sp2<-na.omit(extract(bios,sp2)) #extract climatic variables for species 1
  
  e1_max<-apply(e_sp1,2,max) #max of species 1
  e1_min<-apply(e_sp1,2,min) #min of species 1
  
  e2_max<-apply(e_sp2,2,max) #max of species 2
  e2_min<-apply(e_sp2,2,min) #min of species 2
  
  overlap<-rep(NA,ncol(e_sp1)) #to store  overlap of enviromental varialble (Axis overlap, Fig. 1)
  for(n in 1:ncol(e_sp1)){    
    overlap[n]<-(min(e1_max[n],e2_max[n])-max(e1_min[n],e2_min[n]))/(max(e1_max[n],e2_max[n])-min(e1_min[n],e2_min[n]))
    if(overlap[n]<0){overlap[n]<-0} #there is no overlap between the axis if overlap is negative
  }
  
  sumoverlap<-sum(overlap)/ncol(e_sp1) # Cumulative sum of axis overlap, relative to potential overal
  return(list(overlap,sumoverlap)) 
}


################# end function ############

################ rtr function
rtr<-function(sp1,sp2,bios,rotation=TRUE,translation=TRUE){
  library('raster')
  library('maptools')
  library('dismo')
  library('rgeos') #for convexhull and centroid
  
  ### NEW ADDITION added post publication in Nunes & Pearson, 2016 :
  
  ##this makes sure that both sp1 and sp2 have same names so that they can be bound together
  
  colnam<-c("Longitude","Latitude")
  names(sp1)=c(colnam)
  names(sp2)=c(colnam)  
  
  ########
  sp3<- rbind(unique(sp1),unique(sp2)) # unify both clouds of points to maintain spatial configuration
  poly<-SpatialPoints(sp3) 
  poly<-gConvexHull(poly) #make minimum convex polygon
  cent<-gCentroid(poly) #find centroid of polygon
  
  x_centre<-as.vector(extent(cent))[1]  #find x- centre of polygon
  y_centre<-as.vector(extent(cent))[3]  #find y- centre of polygon
  
  x<-sp3[,1] #x-coordinates
  y<-sp3[,2] #y-coordinates
  
  cat(paste('',i,'-RUN  ',sep=''))
  repeat{
    ###### rotation
    if(rotation==TRUE){ ##option to rotate the points, no translation
      deg<-runif(1,0,360) #random angle in degrees
      
      rad<-deg*0.0174532925 #coverts degrees to radian
      
      R<-matrix(c(cos(rad),sin(rad),0,-sin(rad),cos(rad),0,0,0,1),byrow=T,ncol=3,nrow=3) #clockwise
      
      a<-matrix(c(1,0,0,0,1,0,x_centre,y_centre,1),nrow=3,ncol=3,byrow=T) #rotation at the centre of origin
      
      cn<-matrix(c(1,0,0,0,1,0,-x_centre,-y_centre,1),nrow=3,ncol=3,byrow=T) #rotation at the centre of origin
      
      R<-cn%*%R%*%a
      
      x<-sp3[,1] #x-coordinates
      y<-sp3[,2] #y-coordinates
      
      rot<-matrix(rep(NA),ncol=3,nrow=length(x))
      for(j in 1:length(x)){
        rot[j,]<-R%*%matrix(c(x[j],y[j],1),ncol=1,nrow=3)    #apply rotation to points    
      }
      
      #apply rotation matrix to the original matrix
      repmat<-matrix(rep(c(x_centre,y_centre)),ncol=2,nrow=nrow(sp3),byrow=T)
      
      s<-sp3-repmat
      s<-cbind(s,1)
      s<-matrix(unlist(s),ncol=3,byrow=F)
      so<-tcrossprod(s,round(R,5))
      
      newpool1<-so[,1:2]+repmat     
      
      x<-newpool1[,1]
      y<-newpool1[,2]
    }
    if(translation==TRUE){ #option to translate the points, no rotation
      #### translation 
      xmin=extent(bios)[1];xmax=extent(bios)[2];ymin=extent(bios)[3];ymax=extent(bios)[4]
      
      x_trans<-runif(1,-(xmax-xmin),(xmax-xmin)) #random long trans
      y_trans<-runif(1,-(ymax-ymin),(ymax-ymin)) #random lat trans
      
      x<-round(x+x_trans) #add translation vector
      
      y<-round(y+y_trans) #add translation vector
      
      
      newpool1<-cbind(x,y) #combine new coordinates
    }
    
    # plot(bios[[1]])
    #points(newpool1)
    
    climate_niche = extract(bios, newpool1) #if NA, then a point is in the ocean/outside study region - no variables
    if(isTRUE(length(which(is.na(climate_niche)==FALSE))==length(climate_niche))){  #if TRUE no NAs then break loop, keep RTR polygon
      break}}
  
  gridcells<-extract(bios, newpool1,cellnumbers=TRUE)[,1] #store gridcells to check for duplicates
  sim.niche<-mo.metric(newpool1[1:nrow(sp1),],newpool1[(nrow(sp1)+1):nrow(climate_niche),],bios)[[2]] #measure niche overlap of simulated polygon
  print(noquote(paste('Simulated Niche Overlap:',sim.niche))) #print measured niche overlap
  
  return(list(sim.niche,as.vector(gridcells)))} #return niche overlap of RTR and the gridcells occupied by the simulated points

####################### end functions ########################

################### PLOT SIGNIFICANCE
RTRsignificance<-function(sp1,sp2,bios,rtrs, tails=TRUE, divergence=TRUE){
  
  climate_sp1 = na.omit(extract(bios, sp1)) #extract niches
  climate_sp2 = na.omit(extract(bios, sp2)) #extract niches
  
  observed<-mo.metric(sp1,sp2,bios)[[2]] #observed niche overlap
  meaniche<-mean(rtrs)
  
  ### is PNC, one-tailed test
  
  nineperc<-as.numeric(quantile(rtrs,0.95)) #define critical value limits for up
  top<-observed>as.numeric(quantile(rtrs,0.95)) #niche conserved
  Fn<-ecdf(as.numeric(rtrs))
  location<-Fn(observed)  
  pvalue<-min((sum(rtrs <= observed)/length(rtrs)),(sum(rtrs >= observed)/length(rtrs))) #one tailed each time
  bottom<-NA
  fiveperc<-NA
  
  if(tails==TRUE){ ##two-tailed test
    fiveperc<-as.numeric(quantile(rtrs,0.025)) #define critical value limits for bottom
    nineperc<-as.numeric(quantile(rtrs,0.975)) #define critical value limits for up
    bottom<-observed<as.numeric(quantile(rtrs,0.025)) #niche diverge
    top<-observed>as.numeric(quantile(rtrs,0.975)) #niche conserved
    Fn<-ecdf(as.numeric(rtrs))
    location<-Fn(observed)
    pvalue<-min((sum(rtrs <= observed)/length(rtrs)),(sum(rtrs >= observed)/length(rtrs)))*2 #two tailed p-value
  }
  
  if(divergence==TRUE){ ## is PND, one tailed test
    fiveperc<-as.numeric(quantile(rtrs,0.05)) #define critical value limits for bottom
    bottom<-observed<as.numeric(quantile(rtrs,0.05)) #niche diverge
    Fn<-ecdf(as.numeric(rtrs))
    location<-Fn(observed)
    pvalue<-min((sum(rtrs <= observed)/length(rtrs)),(sum(rtrs >= observed)/length(rtrs))) #one tailed
    top<-NA
    nineperc<-NA
  }
  
  hist(rtrs, xlab='Niche Overlap Value', ylab='Frequency',main=NULL,xlim=c(0,1))
  legend('topright',legend=c('Observed Niche Overlap Value'), lty=c(1),lwd=c(2),col='red')
  abline(v=observed, col='red', lwd=2)
  return(data.frame(observed,bottom,top,meaniche,fiveperc,nineperc,pvalue,location))}

################## end function

#### function to store grid cells to test for spatial bias
rtr.bias<-function(sp1,sp2,bios,rotation=TRUE,translation=TRUE){
  library('raster')
  library('maptools')
  library('dismo')
  library('rgeos') #for convexhull and centroid
  
  ### NEW ADDITION added post publication in Nunes & Pearson, 2016 :
  
  ##this makes sure that both sp1 and sp2 have same names so that they can be bound together
  
  colnam<-c("Longitude","Latitude")
  names(sp1)=c(colnam)
  names(sp2)=c(colnam)  
  
  ########
  sp3<- rbind(unique(sp1),unique(sp2)) # unify both clouds of points to maintain spatial configuration
  poly<-SpatialPoints(sp3) 
  poly<-gConvexHull(poly) #make minimum convex polygon
  cent<-gCentroid(poly) #find centroid of polygon
  
  x_centre<-as.vector(extent(cent))[1]  #find x- centre of polygon
  y_centre<-as.vector(extent(cent))[3]  #find y- centre of polygon
  
  x<-sp3[,1] #x-coordinates
  y<-sp3[,2] #y-coordinates
  
  repeat{
    ###### rotation
    if(rotation==TRUE){ ##option to rotate the points, no translation
      deg<-runif(1,0,360) #random angle in degrees
      
      rad<-deg*0.0174532925 #coverts degrees to radian
      
      R<-matrix(c(cos(rad),sin(rad),0,-sin(rad),cos(rad),0,0,0,1),byrow=T,ncol=3,nrow=3) #clockwise
      
      a<-matrix(c(1,0,0,0,1,0,x_centre,y_centre,1),nrow=3,ncol=3,byrow=T) #rotation at the centre of origin
      
      cn<-matrix(c(1,0,0,0,1,0,-x_centre,-y_centre,1),nrow=3,ncol=3,byrow=T) #rotation at the centre of origin
      
      R<-cn%*%R%*%a
      
      x<-sp3[,1] #x-coordinates
      y<-sp3[,2] #y-coordinates
      
      rot<-matrix(rep(NA),ncol=3,nrow=length(x))
      for(j in 1:length(x)){
        rot[j,]<-R%*%matrix(c(x[j],y[j],1),ncol=1,nrow=3)    #apply rotation to points    
      }
      
      #apply rotation matrix to the original matrix
      repmat<-matrix(rep(c(x_centre,y_centre)),ncol=2,nrow=nrow(sp3),byrow=T)
      
      s<-sp3-repmat
      s<-cbind(s,1)
      s<-matrix(unlist(s),ncol=3,byrow=F)
      so<-tcrossprod(s,round(R,5))
      
      newpool1<-so[,1:2]+repmat     
      
      x<-newpool1[,1]
      y<-newpool1[,2]
    }
    if(translation==TRUE){ #option to translate the points, no rotation
      #### translation 
      xmin=extent(bios)[1];xmax=extent(bios)[2];ymin=extent(bios)[3];ymax=extent(bios)[4]
      
      x_trans<-runif(1,-(xmax-xmin),(xmax-xmin)) #random long trans
      y_trans<-runif(1,-(ymax-ymin),(ymax-ymin)) #random lat trans
      
      x<-round(x+x_trans) #add translation vector
      
      y<-round(y+y_trans) #add translation vector
      
      
      newpool1<-cbind(x,y) #combine new coordinates
    }
    
    # plot(bios[[1]])
    #points(newpool1)
    
    climate_niche = extract(bios, newpool1) #if NA, then a point is in the ocean/outside study region - no variables
    if(isTRUE(length(which(is.na(climate_niche)==FALSE))==length(climate_niche))){  #if TRUE no NAs then break loop, keep RTR polygon
      break}}
  
  gridcells<-extract(bios, newpool1,cellnumbers=TRUE)[,1] #store gridcells to check for duplicates or spatial bias
  
  return(as.vector(gridcells))} #return the gridcells occupied by the simulated points
#### end function ####

#### function to measure and visualise spatial bias ###
rtr.bias.map<-function(sp1,sp2,bios,iter,plot=TRUE){
  
  r<-raster(ncol=ncol(bios[[1]]),nrow=nrow(bios[[1]]))
  extent(r)<-extent(bios[[1]])
  res(r)<-res(bios[[1]])
  projection(r)<-projection(bios[[1]])
  values(r)<-0
  r[which(is.na(values(bios[[1]])))]<-NA
  r[rtr.bias(sp1,sp2,bios[[1]])]<-1
  r3<-r
  for(i in 1:iter){
    cat(paste('',i,'-RUN  ',sep=''))
    values(r)<-0
    r[rtr_bias(sp1,sp2,bios[[1]])]<-1
    r2<-r
    r3<-mosaic(r3,r2,fun=sum)
    removeTmpFiles(1) #remove temp files from 1 hour ago - to save memory
  }
  r3[which(is.na(values(bios[[1]])))]<-NA
  r4<-r3
  if(plot==TRUE){
    plot(r3)
  }
  return(r4)}
#### end function


############ RTR script for plotting percentile and stability througout simulations
library('doParallel')
library('doRNG')

cluster<-3 #number of CPUs to run simulation -advice: always leave one CPU free to not overlaod the computer
cl<-makeCluster(cluster) # add number of CPUs
registerDoParallel(cl)
seeded=10 #set seed for Mesenne-Twister, same seed=repeatability; NA = always different seeds,non-repeatable

trials=10  #number of RTRs per batch, (e.g do 100 RTRs each time, to measure percentile, stability, duplicates etc) 
cutoff=999  # stop running when reaches this number+1
sample=5 #define how many percentiles necessary to measure stability

start.time<-Sys.time() # record start time

observed<-mo.metric(sp1,sp2,bios)[[2]] #observed niche overlap

set.seed(seeded, kind="Mersenne-Twister") #set seed

RTRs<-foreach(i=1:trials, .combine='cbind', .errorhandling='remove') %dopar% {
  rtr(sp1,sp2,bios,rotation=TRUE, translation=TRUE)}

## check for duplicates
if(length(which(duplicated(RTRs[2,])==TRUE))>0){ #is number of duplicates > 0 ?
  non.dup<-RTRs[,-(which(duplicated(RTRs[2,])==TRUE))]### if YES then remove duplicates to keep unique NOs (niche overlap)
}else{non.dup<-RTRs} #if no duplicates then no removal


null.no<-c(non.dup[1,]) ##store NO values only for building null library


Fn<-ecdf(as.numeric(null.no)) #Empirical Cumulative Distribution Function
perc<-Fn(observed) #Note: this is not p-value, it is the location of the observed niche overlap in the null distribution (the percentile)

plot(1, type="n", xlab="Number of Iterations", ylab="Percentile", xlim=c(0,cutoff+1), ylim=c(0, 1))
points(length(null.no),perc, pch=19)
abline(0.05,0,col='red')
abline(0.95,0,col='red')
abline(0,0,col='green')
legend(cutoff/2,0.8, c('Percentile','Stability','Significance Threshold','Stability Threshold'), pch=c(19,19,NA,NA),lty = c(0, 0, 1,1),lwd = c(0, 0, 1,1),bty = "n",col = c("black","blue","red","green"))


#_prev ==stands for previous value

null.no_prev<-length(null.no) #number of unique NOs of previous step

perc_prev<-perc #store percentile of previous step

stability_prev<-c(1) #start of stability string

#checkpoints
print(noquote(paste('Number of Unique RTRs:','',length(null.no))))#print number of Iterations (number of unique NOs)

end.time<-Sys.time() #Register end of loop
time.taken<-end.time-start.time ##Register time it took to run the batch
print(noquote(paste('Duration of Analysis:',time.taken))) # print duration of batch so can estimate how long it will take to finish all the repeats

repeat{
  RTRs<-foreach(i=1:trials, .combine='cbind', .errorhandling='remove') %dopar% {
    rtr(sp1,sp2,bios,rotation=TRUE, translation=TRUE)}
  
  RTRs<-cbind(non.dup,RTRs) ##bind new null NO with previous batch of unique NOs
  
  ## check for duplicates
  if(length(which(duplicated(RTRs[2,])==TRUE))>0){ #is number of duplicates > 0 ?
    non.dup<-RTRs[,-(which(duplicated(RTRs[2,])==TRUE))]### if YES then remove duplicates to keep unique NOs
  }else{non.dup<-RTRs} #if no duplicates then no removal
  
  
  null.no<-c(non.dup[1,]) ##store NO values only for building null library
  
  
  Fn<-ecdf(as.numeric(null.no)) #Empirical Cumulative Distribution Function
  perc<-Fn(observed) #Note: this is not p-value, it is the location of the observed niche overlap in the null distribution (the percentile)
  
  points(length(null.no),perc,pch=19) # plot percentile
  
  null.no_prev<-c(null.no_prev,length(null.no)) # #number of combined unique null nicheo verlap values 
  
  perc_prev<-c(perc_prev,perc) #store all percentile locations
  
  
  df<-data.frame(null.no_prev,perc_prev) #table showing numble of RTRs and percentile of observed niche overlap
  colnames(df)<-c('Iterations','Percentile')
  
  
  lines(df[,1],df[,2]) #plot line showing number of NOs and percentile
  if(isTRUE(length(null.no_prev)>sample)){
    points(length(null.no),round(sd(df[(nrow(df)-sample):nrow(df),2]),2),col='blue',pch=19) #plot stability of percentile
    
    stability<-c(round(sd(df[(nrow(df)-sample):nrow(df),2]),2)) #calculate stability of percentile => standard deviation of percentile of the last *sample* batches
    stability_prev<-c(stability_prev,stability) #store stability measure
  }
  #checkpoint for stability
  print(noquote(paste('Stable:',isTRUE(abs(df[nrow(df),2]-df[(nrow(df)-1),2])<0.01)))) #measure of stability
  #checkpoints
  print(noquote(paste('Number of Unique RTRs:','',length(null.no))))#print number of Iterations (number of unique NOs)
  
  end.time<-Sys.time() #Register end of loop
  time.taken<-end.time-start.time ##Register time it took to run the batch
  print(noquote(paste('Duration of Analysis:',time.taken))) # print duration of batch so can estimate how long it will take to finish all the repeats
  
  if(isTRUE(length(null.no)>cutoff)) break } #if number of NOs reaches cutoff, stop the simulation
stopCluster(cl)  

stability<-stability_prev[-1] #print measure of stability at the end of the analysis

final.table<-cbind(df,c(rep('NA',sample),stability)) #table showing number of NOs, percentile and stability of percentile
colnames(final.table)<-c('Iterations','Percentile','Stability')
final.table #print table# NunesPearson2016
