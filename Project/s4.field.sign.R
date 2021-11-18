setwd("C:/Users/BANZH/Downloads/2019spring/AOS 219/AOS project")
###Step1: nmRead in TCWV and LAI data as 3D matrix
library(ncdf4)
LAInc=nc_open("1982-2015.grow.LAI.nc")
LAI=ncvar_get(LAInc,"Band1")
lon=ncvar_get(LAInc,"longitude")
lat=ncvar_get(LAInc,"latitude")
nc_close(LAInc)
TCWVnc=nc_open("1982-2015.grow.TCWV.nc")
TCWV=ncvar_get(TCWVnc,"tcwv")
nc_close(TCWVnc)
###set <0 LAI to be 0
ind=which(LAI<0, arr.ind=TRUE)
LAI[ind]=0
LAI=LAI/100
dim=dim(TCWV)
laiclim=array(0,dim=c(dim[1],dim[2],204))
tcwvclim=laiclim
LAIlatcor=array(0,dim=dim)
TCWVlatcor=array(0,dim=dim)
for(i in 1:length(lat)){
  LAIlatcor[,i,]=LAI[,i,]*abs(cos(lat[i]/180*pi)) ##latitude correction
  TCWVlatcor[,i,]=TCWV[,i,]*abs(cos(lat[i]/180*pi))
}
#for(i in 1:6){
#  laiclim[,,seq(i,198+i,by=6)]=rowMeans(LAIlatcor[,,seq(i,198+i,by=6)],dims=2)
#  tcwvclim[,,seq(i,198+i,by=6)]=rowMeans(TCWVlatcor[,,seq(i,198+i,by=6)],dims=2)
#}
#laianom=LAIlatcor-laiclim
#tcwvanom=TCWVlatcor-tcwvclim
Cll=array(NA,dim=dim)
Ctt=array(NA,dim=dim)
for(i in 1:dim[1]){
  for(j in 1:dim[2]){
    if(!is.na(LAIlatcor[i,j,1])){
      if(sum(LAIlatcor[i,j,])!=0){
      tmp1=acf(LAIlatcor[i,j,],lag.max=204,plot=FALSE)
      tmp2=acf(TCWVlatcor[i,j,],lag.max=204,plot=FALSE)
      Cll[i,j,]=tmp1$acf
      Ctt[i,j,]=tmp2$acf
      }
    }
  }
}
#laiind=which(!is.na(laianom[,,1]),arr.ind=TRUE)
#laianomuse=array(0,dim=c(length(laiind[,1]),dim[3]))
#tcwvind=which(!is.na(tcwvanom[,,1]),arr.ind=TRUE)
#tcwvanomuse=array(0,dim=c(length(laiind[,1]),dim[3]))
#for(i in 1:dim[3]){
#  laianomuse[,i]=laianom[,,i][laiind]
#  tcwvanomuse[,i]=tcwvanom[,,i][laiind]
#}
#Cll=array(0,dim=c(length(laianomuse[,1]),204))
#Ctt=array(0,dim=c(length(laianomuse[,1]),204))
#for(i in 1:length(laianomuse[,1])){
#  tmp1=acf(laianomuse[i,],lag.max=204)
#  tmp2=acf(tcwvanomuse[i,],lag.max=204)
#  Cll[i,]=as.numeric(tmp1$acf)
#  Ctt[i,]=as.numeric(tmp2$acf)
#}
tau=array(NA,dim=c(dim[1],dim[2]))
for(i in 1:dim[1]){
  for(j in 1:dim[2]){
    if(!is.na(sum(Cll[i,j,]))){
      tau[i,j]=1+2*sum(Cll[i,j,]*Ctt[i,j,])
    }
  }
}
Neff=array(NA,dim=c(dim[1],dim[2]))
Neff=204/tau  ##Neff
getp <- function(ts1,ts2,Neff){
  r=cor(ts1,ts2)
  t=r*sqrt(Neff-2)/sqrt(1-r*r)
  p=pnorm(t)
  vec=c(p,r)
  return(vec)
}
testpr=array(NA,dim=c(dim[1],dim[2],2))
for(i in 1:dim[1]){
  for(j in 1:dim[2]){
     if(!is.na(sum(Cll[i,j,]))){
#    testp[(i-1)*31+j]=getp(ninodjf,hgtdjfanom[i,j,],Neff[i,j])
     testpr[i,j,]=getp(TCWVlatcor[i,j,],LAIlatcor[i,j,],Neff[i,j])
     }
  }
}
laiind=which(!is.na(LAIlatcor[,,1]),arr.ind=TRUE)
L=length(laiind[,1])
##Chen's method
check <- function(ts1,ts2){
  test=cor.test(ts1,ts2)
  if(((test$p.value)<=0.025) || (test$p.value>=0.975)){
    return(1)
  }else{
    return(0)
  }
}
check1 <- function(ts1,ts2,Neff){
  r=cor(ts1,ts2)
  t=r*sqrt(Neff-2)/sqrt(1-r*r)
  p=pnorm(t)
  if((p<=0.025) || (p>=0.975)){
    return(1)
  }else{
    return(0)
  }
}
test=array(0,dim=c(dim[1],dim[2]))
testran=test
for(i in 1:dim[1]){
  for(j in 1:dim[2]){
    if(!is.na(sum(Cll[i,j,]))){
    gau=rnorm(204,mean=0,sd=1)
    test[i,j]=check1(LAIlatcor[i,j,],TCWVlatcor[i,j,],Neff[i,j])
    testran[i,j]=check(gau,TCWVlatcor[i,j,])
    }
  }
}
sum=sum(test)
sumran=sum(testran)
ratio=sum/L
ratioran=sumran/L
print(ratio)
print(ratioran)
library(lattice)
library(raster)
library(rasterVis)
library(maps)
library(mapdata)
library(maptools)
rotate <- function(x) t(apply(x, 2, rev))
plotmat1 <- function(LAIp,x){
  r <-raster(rotate(LAIp[720:1,]),  xmn=-180, xmx=180, ymn=-90, ymx=90, 
             crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=180,0,0"))
  ext <- as.vector(extent(r))
  boundaries <- map('worldHires', fill=TRUE,
                    xlim=ext[1:2], ylim=ext[3:4],
                    plot=FALSE)
  IDs <- sapply(strsplit(boundaries$names, ":"), function(x) x[1])
  bPols <- map2SpatialPolygons(boundaries, IDs=IDs,
                               proj4string=CRS(projection(r)))
  #rgb.palette <- colorRampPalette(c("red", "yellow","green"))
  levelplot(r,margin=FALSE,scales=list(cex=3),at=seq(-1,1,by=0.001),col.regions=colorRampPalette(c("red", "yellow", "blue"))(2000),
            xlab=list(cex=3),ylab=list(cex=3),main=list(x,cex=3),
            colorkey=list(labels=list(cex=3))) + layer(sp.polygons(bPols))
}
plotmat2 <- function(LAIp,x){
  r <-raster(rotate(LAIp[720:1,]),  xmn=-180, xmx=180, ymn=-90, ymx=90, 
             crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=180,0,0"))
  ext <- as.vector(extent(r))
  boundaries <- map('worldHires', fill=TRUE,
                    xlim=ext[1:2], ylim=ext[3:4],
                    plot=FALSE)
  IDs <- sapply(strsplit(boundaries$names, ":"), function(x) x[1])
  bPols <- map2SpatialPolygons(boundaries, IDs=IDs,
                               proj4string=CRS(projection(r)))
  #rgb.palette <- colorRampPalette(c("red", "yellow","green"))
  levelplot(r,margin=FALSE,scales=list(cex=3),at=seq(-1,1,by=0.001),col.regions=colorRampPalette(c("red", "yellow", "blue"))(2000),
            xlab=list(cex=3),ylab=list(cex=3),main=list(x,cex=3),
            colorkey=list(labels=list(cex=3))) + layer(sp.polygons(bPols))
}
testprsig=array(NA,dim=c(dim[1],dim[2],2))
for(i in 1:dim[1]){
  for(j in 1:dim[2]){
    if(!is.na(testpr[i,j,1])){
      if((testpr[i,j,1]<=0.025) || (testpr[i,j,1]>=0.975)){
      testprsig[i,j,]=testpr[i,j,]
      }
    }
  }
}
testprcor=testpr[c(361:720,1:360),,2]
testprp=testpr[c(361:720,1:360),,1]
testprsigcor=testprsig[c(361:720,1:360),,2]
jpeg("Figures/correlation.jpg",width=1200,height=600)
plotmat2(testprcor,"correlation value between TCWV and LAI")
dev.off()
jpeg("Figures/p_corr.jpg",width=1200,height=600)
plotmat1(testprp,"p value for correlation")
dev.off()
jpeg("Figures/corr_sign.jpg",width=1200,height=600)
plotmat2(testprsig[,,2],"correlation value at only significant locations")
dev.off()
