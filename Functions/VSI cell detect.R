VSI.celldetect=function(Ipath,SptF,iSpt,SptD,Gate,Stain,ID,Th=0.04){
  
  options(java.parameters = "-Xmx4g")
  library(RBioFormats)
  library(EBImage)
  library(MiXR)
  library(sp)
  
  Lap = matrix(c(-1,-1,-1,-1,8,-1,-1,-1,-1),nrow=3)
  #==
  Yp = ceiling(iSpt/SptF)
  Xp = iSpt-(SptF*(Yp-1))
  #==
  I.mat = read.image(Ipath,series=2,res=1,subset=list(X=SptD$sizeX[[Xp]]:SptD$sizeX[[Xp+1]],
                                                      Y=SptD$sizeY[[Yp]]:SptD$sizeY[[Yp+1]]))
  
  #==Nucleus detection==#
  I1 = channel(1-I.mat[,,1],'grey') #Contains nucleus Info
  I2 = channel(1-I.mat[,,2],'grey') #Contains mix info
  I3 = channel(1-I.mat[,,3],'grey') #Contains Cyto Info
  
  #Positive cells: weak nuclear stain
  Hole = NormalizeIm(I3-opening(I3,makeBrush(9,'box'))) #roughly brown stain
  Hmask = closing(thresh(Hole,w=25,h=25,offset=0.1),makeBrush(7,'box')) #The mask of cyto + holes
  lowN = opening(fillHull(Hmask)-Hmask,makeBrush(5,'disc'))
  
  #Other Nuclei
  N=NormalizeIm(medianFilter((NormalizeIm(I1)+NormalizeIm(I2)-NormalizeIm(I3)),3),autoRange=F,verbose=F)
  
  #Go or NoGo for LoG filter
  if(sd(N)>Th){
    N=filter2(gblur(NormalizeIm(N,autoRange=T,verbose=F),sigma=5),Lap);N[which(N<0)]=0;N=NormalizeIm(N,autoRange=T,verbose=F)
  }else{
    warning(paste0('Warning! LoG filter not applied for image #',iSpt,' from ',basename(Ipath)))
  }
  ##
  Nmask = opening(fillHull(closing(thresh(N,w=21,h=21,offset=0.1),makeBrush(5,'box'))),makeBrush(7,'disc'))
  Nmask = Nmask|lowN
  Nmask = propagate(x=Nmask,seeds=bwlabel(erode(Nmask,makeBrush(5,'disc'))),mask=Nmask)
  
  #==Cytoplasm detection==#
  Ring = dilate(Nmask, makeBrush(7,'disc'))-Nmask
  Cyto = opening(fillHull(Hmask),makeBrush(3,'disc'))
  
  Cell = (Ring|Cyto)|Nmask;Cell=closing(Cell,makeBrush(5,'box'))
  Cell.mask = propagate(Cell,Nmask,mask=Cell)
  Cyto.mask = Cell.mask-Nmask
  
  #==Export segmented Image==#
  #display(paintObjects(x=Cyto.mask,tgt=paintObjects(x=Nmask,tgt=I.mat,col=c('red',NA)),col=c('limegreen',NA)))
  
  
  #==Features extraction==#
  Feat = data.frame(computeFeatures(Nmask,I1,xname='Nucleus',refnames='Blue',
                                    methods.noref=c("computeFeatures.moment", "computeFeatures.shape"),
                                    methods.ref=c("computeFeatures.basic"),basic.quantiles=0.25),
                    computeFeatures(Nmask,I3,xname='Cytoplasm',refnames='Brown',
                                    methods.noref=c("computeFeatures.shape"),
                                    methods.ref=c("computeFeatures.basic"),basic.quantiles=0.25)
  )
  Feat$Centroid.x=Feat$Nucleus.0.m.cx+SptD$sizeX[[Xp]]-1
  Feat$Centroid.y=Feat$Nucleus.0.m.cy+SptD$sizeY[[Yp]]-1
  Feat$ROI=point.in.polygon(Feat$Centroid.x,Feat$Centroid.y,Gate$x,Gate$y)
  
  if(!missing(Stain) & !missing(ID)){
    Feat=data.frame(Stain,ID,Feat)
  }
  return(Feat)
}

