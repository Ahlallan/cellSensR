options(java.parameters = "-Xmx8000m")

library(EBImage)
library(RBioFormats)
library(MiXR)

library(doParallel)
library(foreach)

#Define filters =====================================================================

Lap = matrix(c(-1,-1,-1,-1,8,-1,-1,-1,-1),nrow=3)

#====================================================================================

setwd('')

#====================================================================================
# Structure data

vsi = list.files('Acquisitions',pattern='.vsi',full.names=T)
TumID = unique(sapply(strsplit(gsub('.vsi','',basename(vsi)),'_'),function(x)x[[1]]))
vsi.t = data.frame(TumID,do.call('rbind',lapply(TumID,function(x)grep(x,vsi,value=T))))
colnames(vsi.t)=c('TumID','CD3','F480','HES')  

#====================================================================================
#Play with data

#for blabla {}
Ti = 35165 #Studied TumorID
Marks = c('F480','CD3')

P = vsi.t[which(vsi.t$TumID==Ti),c('HES','CD3','F480')]
M = lapply(P,function(x)read.metadata(as.character(x)))

#--------------------------------
#VSI#
#Block 1 to 6 = low zoom sets
#Block 7 to 13 = high zoom sets
#Block 14 = thumb image

Block = list(S1=1:6,S2=7:13,S3=14)
#--------------------------------
BigDim = lapply(M,function(x)x[[Block[[2]][1]]]$coreMetadata[c('sizeX','sizeY')]) # Get dimensions of our best resolved image
SptF=4;SptDim = lapply(BigDim,function(x)lapply(x,function(y)floor(seq(1,y,length.out=SptF+1))))

#Define ROI and calculate movement between stains
I.thb = read.image(as.character(P[,'HES']),series=2,res=4) #Reads HES image with low res
windows(width=ncol(I.thb),height=nrow(I.thb))
display(I.thb,method='raster')
g = locator(n=512,type='o',lwd=2,col='red')
dev.off()

#HES centroid
I.det = NormalizeIm(1-medianFilter(channel(I.thb,'grey'),5),autoRange=T,verbose=F)
I.det = I.det>otsu(I.det)*0.8;I.det=opening(fillHull(I.det),makeBrush(25,'disc'))
HE.geo = computeFeatures.moment(I.det,I.thb)
#...And the others ...
for(m in Marks){
  I.det = channel(read.image(as.character(P[,m]),series=2,res=4),'grey')
  ...
  ...
  I.det = NormalizeIm(1-medianFilter(I.det,5),autoRange=T,verbose=F)
  I.det = I.det>otsu(I.det)*0.8;I.det=opening(fillHull(I.det),makeBrush(25,'disc'))
}

#--------------------------------

for(m in Marks){
  for(i in 1:(SptF^2)){
    
    Yp= ceiling(i/SptF)
    Xp = i-(SptF*(Yp-1))
    #==
    I.mat = read.image(as.character(P[,'F480']),series=2,res=1,subset=list(X=SptDim[[m]]$sizeX[[Xp]]:SptDim[[m]]$sizeX[[Xp+1]],
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
    Feat$Centroid.x=Feat$Nucleus.0.m.cx+SptDim$F480$sizeX[[Xp]]-1
    Feat$Centroid.y=Feat$Nucleus.0.m.cy+SptDim$F480$sizeY[[Yp]]-1
    
  }
}


