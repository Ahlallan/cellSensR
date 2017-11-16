options(java.parameters = "-Xmx36g")

library(EBImage)
library(RBioFormats)
library(MiXR)

library(doParallel)
library(foreach)
library(sp)


#Define filters =====================================================================

Lap = matrix(c(-1,-1,-1,-1,8,-1,-1,-1,-1),nrow=3)

#====================================================================================

#setwd('')

#====================================================================================
# Structure data

vsi = list.files('Acquisitions',pattern='.vsi',full.names=T);vsi=vsi[!grepl('Temoin',vsi)]
TumID = unique(sapply(strsplit(gsub('.vsi','',basename(vsi)),'_'),function(x)x[[1]]))
vsi.t = data.frame(TumID,do.call('rbind',lapply(TumID,function(x)grep(x,vsi,value=T))))

Marks = c('CD11c','CD3','F480','HES')
colnames(vsi.t)=c('TumID',Marks)  

#====================================================================================
#Play with data
full.dat=list()

for(Ti in vsi.t$TumID[1:2]){ #Studied TumorID
  
  P = vsi.t[which(vsi.t$TumID==Ti),Marks]
  M = lapply(P,function(x)read.metadata(as.character(x)))
  
  #--------------------------------
  #VSI#
  #Block 1 to 6 = low zoom sets
  #Block 7 to 13 = high zoom sets
  #Block 14 = thumb image
  
  Block = list(S1=1:6,S2=7:13,S3=14)
  #--------------------------------
  
  BigDim = lapply(M,function(x)x[[Block[[2]][1]]]$coreMetadata[c('sizeX','sizeY')]) # Get dimensions of our best resolved image
  SptF=4;SptDim = lapply(BigDim,function(x)lapply(x,function(y)floor(seq(1,y,length.out=SptF+1)))) #Blocks definition for analysis
  
  SmallDim = lapply(M,function(x)x[[Block[[2]][4]]]$coreMetadata[c('sizeX','sizeY')]) # Get dimensions of our low res image (used for gating)
  scaling = lapply(seq_along(BigDim),function(i)c(ax=BigDim[[i]]$sizeX/SmallDim[[i]]$sizeX,ay=BigDim[[i]]$sizeY/SmallDim[[i]]$sizeY))#Scaling from high to low Res
  names(scaling)=names(BigDim)
  
  #--------------------------------
  #--------------------------------
  
  #Define ROI and calculate movement between stains
  I.thb = read.image(as.character(P[,'HES']),series=2,res=4) #Reads HES image with low res
  windows(width=nrow(I.thb),height=ncol(I.thb))
  display(I.thb,method='raster')
  g = locator(n=512,type='o',lwd=2,col='red')
  #dev.off()
  
  #HES centroid
  otsu.adj=0.95
  I.det = channel(I.thb,'grey');q = tail(quantile(I.det,probs=seq(0,1,0.1)),n=2)[1];I.det[which(I.det>q)]=q #Correct addition of white pixels
  I.det = NormalizeIm(1-medianFilter(I.det,5),autoRange=T,verbose=F);q=quantile(I.det,probs=seq(0,1,0.01));q=min(subset(q,q>0))
  I.det=NormalizeIm(log10(I.det+q),autoRange=T,verbose=F)
  I.det = I.det>otsu(I.det)*otsu.adj;I.det=opening(fillHull(I.det),makeBrush(25,'disc'))
  HE.geo = computeFeatures.moment(I.det,I.thb)
  #display(paintObjects(I.det,I.thb,thick=T),method='raster')
  #polygon(g$x,g$y)
  
  #...And the others ...
  O.geo = list();g.marks = list();i=1
  for(m in subset(Marks,!grepl('HES',Marks))){
    I.thb = read.image(as.character(P[,m]),series=2,res=4)
    I.det = channel(I.thb,'grey'); q = tail(quantile(I.det,probs=seq(0,1,0.1)),n=2)[1];I.det[which(I.det>q)]=q 
    I.det = NormalizeIm(1-medianFilter(I.det,5),autoRange=T,verbose=F);q=quantile(I.det,probs=seq(0,1,0.01));q=min(subset(q,q>0))
    I.det = NormalizeIm(log10(I.det+q),autoRange=T,verbose=F)
    I.det = I.det>otsu(I.det)*otsu.adj;I.det=opening(fillHull(I.det),makeBrush(25,'disc'))
    #
    O.geo=c(O.geo,list(computeFeatures.moment(I.det,I.thb)))
    g.marks=c(g.marks,list(g))
    #
    newx = (g.marks[[i]]$x-(HE.geo[1]-O.geo[[i]][1]))
    newy = (g.marks[[i]]$y-(HE.geo[2]-O.geo[[i]][2]))
    
    #theta = O.geo[[i]][5]-HE.geo[5]
    theta = 0
    g.marks[[i]]$x = newx*cos(theta)-newy*sin(theta)
    g.marks[[i]]$y = newx*sin(theta)+newy*cos(theta)
    
    #
    #windows(width=nrow(I.thb),height=ncol(I.thb))
    #display(paintObjects(I.det,I.thb,thick=T),method='raster')
    #polygon(g.marks[[i]]$x,g.marks[[i]]$y)
    #
    i=i+1
  }
  
  rm(list=c('i','I.thb','I.det','newx','newy','theta'));gc()
  names(O.geo)=names(g.marks)=subset(Marks,!grepl('HES',Marks))
  Sys.sleep(1)
  dev.off()
  
  #Transforms gate coordinates to high res image coordinates
  g.hmarks = lapply(names(g.marks),function(p)list(x=g.marks[[p]]$x*scaling[[p]]['ax'],y=g.marks[[p]]$y*scaling[[p]]['ay']))
  names(g.hmarks)=names(g.marks)
  
  #--------------------------------
  #--------------------------------
  
  full=list()
  for(m in subset(Marks,!grepl('HES',Marks))){
    #cl=makeCluster(2)
    #registerDoParallel(cl,2)
    dat=foreach(i=1:(SptF^2),.packages=c('EBImage','RBioFormats'),.combine='rbind')%do%{
      
      Yp= ceiling(i/SptF)
      Xp = i-(SptF*(Yp-1))
      #==
      I.mat = read.image(as.character(P[,m]),series=2,res=1,subset=list(X=SptDim[[m]]$sizeX[[Xp]]:SptDim[[m]]$sizeX[[Xp+1]],
                                                                        Y=SptDim[[m]]$sizeX[[Yp]]:SptDim[[m]]$sizeX[[Yp+1]]))
      
      #==Nucleus detection==#
      I1 = channel(1-I.mat[,,1],'grey') #Contains nucleus Info
      I2 = channel(1-I.mat[,,2],'grey') #Contains mix info
      I3 = channel(1-I.mat[,,3],'grey') #Contains Cyto Info
      
      #Positive cells: weak nuclear stain
      Hole = NormalizeIm(I3-opening(I3,makeBrush(9,'box'))) #roughly brown stain
      Hmask = closing(thresh(Hole,w=25,h=25,offset=0.1),makeBrush(7,'box')) #The mask of cyto + holes
      lowN = opening(fillHull(Hmask)-Hmask,makeBrush(5,'disc'))
      
      #Other Nuclei
      N = NormalizeIm(medianFilter((NormalizeIm(I1)+NormalizeIm(I2)-NormalizeIm(I3)),3),autoRange=T,verbose=F) 
      N=filter2(gblur(N,sigma=5),Lap);N[which(N<0)]=0;N=NormalizeIm(N,autoRange=T,verbose=F)
      Nmask = opening(fillHull(closing(thresh(N,w=21,h=21,offset=0.1),makeBrush(5,'box'))),makeBrush(7,'disc'))
      Nmask = Nmask|lowN
      Nmask = propagate(x=Nmask,seeds=bwlabel(erode(Nmask,makeBrush(5,'disc'))),mask=Nmask)
      #display(paintObjects(Nmask,I.mat,col=c('purple',NA)))
      
      #==Cytoplasm detection==#
      Ring = dilate(Nmask, makeBrush(7,'disc'))-Nmask
      Cyto = opening(fillHull(Hmask),makeBrush(3,'disc'))
      
      Cell = (Ring|Cyto)|Nmask;Cell=closing(Cell,makeBrush(5,'box'))
      Cell.mask = propagate(Cell,Nmask,mask=Cell)
      Cyto.mask = Cell.mask-Nmask
      #display(paintObjects(x=Cyto.mask,tgt=paintObjects(x=Nmask,tgt=I.mat,col=c('red',NA)),col=c('limegreen',NA)))
      
      #==Features extraction==#
      Feat = data.frame(Stain=m,ID=Ti,
                        computeFeatures(Nmask,I1,xname='Nucleus',refnames='Blue',
                                        methods.noref=c("computeFeatures.moment", "computeFeatures.shape"),
                                        methods.ref=c("computeFeatures.basic"),basic.quantiles=0.25),
                        computeFeatures(Nmask,I3,xname='Cytoplasm',refnames='Brown',
                                        methods.noref=c("computeFeatures.shape"),
                                        methods.ref=c("computeFeatures.basic"),basic.quantiles=0.25)
      )
      Feat$Centroid.x=Feat$Nucleus.0.m.cx+SptDim[[m]]$sizeX[[Xp]]-1
      Feat$Centroid.y=Feat$Nucleus.0.m.cy+SptDim[[m]]$sizeY[[Yp]]-1
      Feat$ROI=point.in.polygon(Feat$Centroid.x,Feat$Centroid.y,g.hmarks[[m]]$x,g.hmarks[[m]]$y)
    }
    #stopCluster(cl);rm(cl)
    full=c(full,list(dat))
  }
  names(full)=subset(Marks,!grepl('HES',Marks))
  full.dat=c(full.dat,list(full))
  
}

names(full.dat)=vsi.t$TumID[1:2]

