options(java.parameters = "-Xmx4g")

library(EBImage)
library(RBioFormats)
library(MiXR)

library(sp)
library(doParallel)
library(pbapply)

#Necessary function==================================================================

invisible(sapply(list.files('Functions',full.names=T),source))

#====================================================================================
# Structure data

vsi = list.files('Acquisitions',pattern='.vsi',full.names=T);vsi=vsi[!grepl('Temoin',vsi)]
TumID = unique(sapply(strsplit(gsub('.vsi','',basename(vsi)),'_'),function(x)x[[1]]))
vsi.t = data.frame(TumID,do.call('rbind',lapply(TumID,function(x)grep(x,vsi,value=T))))

Marks = c('CD11c','CD3','F480','HES')
colnames(vsi.t)=c('TumID',Marks)  

#For testing
#vsi.t = head(vsi.t,5)

#====================================================================================
#PART1 : Gates definition ##
#====================================================================================

full.g=list()  #list of all gates
full.spt=list() #list of all images splits
SptF=5 # nrow for splitting

Pl=list() #List of paths for each Tumor
Ml=list() #Lits of metadata for each Tumor

for(Ti in vsi.t$TumID){ #Studied TumorID
  
  P = vsi.t[which(vsi.t$TumID==Ti),Marks]
  P = sapply(P, as.character)
  #
  nCores=length(Marks)-1
  cl=makeCluster(nCores)
  registerDoParallel(cl,nCores)
  M = pblapply(P,function(x){
    library(RBioFormats)
    read.metadata(x)
    },cl=cl)
  stopCluster(cl);rm(cl)
  #
  Pl = c(Pl,list(P))
  Ml = c(Ml,list(M))
  
  #--------------------------------
  #VSI#
  #Block 1 to 6 = low zoom sets
  #Block 7 to 13 = high zoom sets
  #Block 14 = thumb image
  
  Block = list(S1=1:6,S2=7:13,S3=14)
  #--------------------------------
  
  BigDim = lapply(M,function(x)x[[Block[[2]][1]]]$coreMetadata[c('sizeX','sizeY')]) # Get dimensions of our best resolved image
  SptDim = lapply(BigDim,function(x)lapply(x,function(y)floor(seq(1,y,length.out=SptF+1)))) #Blocks definition for analysis
  full.spt=c(full.spt,list(SptDim))
  
  SmallDim = lapply(M,function(x)x[[Block[[2]][4]]]$coreMetadata[c('sizeX','sizeY')]) # Get dimensions of our low res image (used for gating)
  scaling = lapply(seq_along(BigDim),function(i)c(ax=BigDim[[i]]$sizeX/SmallDim[[i]]$sizeX,ay=BigDim[[i]]$sizeY/SmallDim[[i]]$sizeY))#Scaling from high to low Res
  names(scaling)=names(BigDim)
  
  #--------------------------------
  
  #Define ROIs 
  HE=read.image(as.character(P['HES']),series=2,res=4)
  h=10;w=h*(nrow(HE)/ncol(HE))
  windows(width=w,height=h)
  display(HE,method='raster')
  
  g.marks = list();i=1
  for(m in subset(Marks,!grepl('HES',Marks))){
    I.thb = read.image(as.character(P[m]),series=2,res=4)
    h=10;w=h*(nrow(I.thb)/ncol(I.thb))
    windows(width=w,height=h)
    display(I.thb,method='raster')
    title(main=m,col='red')
    g = locator(n=512,type='o',lwd=2,col='red',xpd=T)
    #
    g.marks=c(g.marks,list(g))
    i=i+1
  }
  
  sapply(dev.list(),dev.off)
  rm(list=c('i','I.thb'));gc()
  names(g.marks)=subset(Marks,!grepl('HES',Marks))
  
  #Transforms gate coordinates to high res image coordinates
  g.hmarks = lapply(names(g.marks),function(p)list(x=g.marks[[p]]$x*scaling[[p]]['ax'],y=g.marks[[p]]$y*scaling[[p]]['ay']))
  names(g.hmarks)=names(g.marks)
  
  #Add info to the final list
  full.g=c(full.g,list(g.hmarks))
}

rm(list=setdiff(ls(),c('vsi.t','Marks','full.g','full.spt','SptF','Pl','Ml','VSI.celldetect')))
names(Ml)=names(Pl)=names(full.spt)=names(full.g)=vsi.t$TumID

#====================================================================================
#PART2 :  Image analysis#
#====================================================================================

full.dat=list()
#dir.create('Segs')

for(Ti in vsi.t$TumID){ #Studied TumorID
  
  print(paste0('Analyzing Tumor ',Ti,'...'))
  full=list()
  P = Pl[[Ti]]
  
  for(m in subset(Marks,!grepl('HES',Marks))){
    
    cl=makeCluster(13)
    registerDoParallel(cl,13)
    clusterExport(cl,'VSI.celldetect')
    invisible(clusterEvalQ(cl,{library('EBImage');library('MiXR');library('sp')}))
    
    dat=pblapply(1:(SptF^2),function(i,Ipath,SptF,SptD,Gate,Stain,ID){
      options(java.parameters = "-Xmx4g")
      library(RBioFormats)
      expath=paste0('Segs/',ID,'_',Stain,'_I',i,'.tif')
      print(SptD$sizeX)
      exp=F;if(i==median(1:SptF^2)){exp=T}
      VSI.celldetect(Ipath,iSpt=i,SptF,SptD,Gate,Stain,ID,export=exp,expath=expath)
    },SptF=SptF,SptD=full.spt[[Ti]][[m]],Gate=full.g[[Ti]][[m]],Ipath=as.character(P[m]),Stain=m,ID=Ti,cl=cl)
    
    stopCluster(cl);rm(cl);gc()
    dat=data.frame(do.call('rbind',dat))
    
    #dat check (for testing)------------------------------
    # dat$col='black'
    # dat$col[which(dat$ROI==1)]='red'
    # color.gradient <- function(x, colors=c("red","yellow","green"), colsteps=100) {
    #   return(colorRampPalette(colors)(colsteps)[findInterval(x, seq(min(x),max(x),length.out=colsteps))])
    # }
    # dat$col=color.gradient(dat$Cytoplasm.Brown.b.mean,
    #                        colors=c('white','cornsilk','red','orange','yellow'),
    #                        1000)
    # windows()
    # plot(dat$Centroid.x,max(dat$Centroid.y)-dat$Centroid.y,
    #      pch=19,cex=0.4,col=dat$col,axes=F,xlab='',ylab='',
    #      main=paste(dat$Stain[1],dat$ID[1]))
    #-----------------------------------------------------
    full=c(full,list(dat))
    
  }
  names(full)=subset(Marks,!grepl('HES',Marks))
  full.dat=c(full.dat,list(full))
  
}

rm(list=c('dat','full','m','Ti'))
names(full.dat)=vsi.t$TumID

