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
#laianomcs=array(0,dim=c(360,180,204))
#tcwvanomcs=array(0,dim=c(360,180,204))
laics=array(0,dim=c(360,180,204))
tcwvcs=array(0,dim=c(360,180,204))
##Spatial average to 1 deg##
lo=seq(1,719,by=2)
la=seq(1,359,by=2)
for(i in 1:length(lo)){
  for(j in 1:length(la)){
    for(t in 1:204){
      laics[i,j,t]=mean(LAIlatcor[((i-1)*2+1):(i*2),((j-1)*2+1):(j*2),t])
      tcwvcs[i,j,t]=mean(TCWVlatcor[((i-1)*2+1):(i*2),((j-1)*2+1):(j*2),t])
    }
  }
}
##Spatial average to 2 deg##
lo=seq(1,359,by=2)
la=seq(1,179,by=2)
laics1=array(0,dim=c(180,90,204))
tcwvcs1=array(0,dim=c(180,90,204))
for(i in 1:length(lo)){
  for(j in 1:length(la)){
    for(t in 1:204){
      laics1[i,j,t]=mean(laics[((i-1)*2+1):(i*2),((j-1)*2+1):(j*2),t])
      tcwvcs1[i,j,t]=mean(tcwvcs[((i-1)*2+1):(i*2),((j-1)*2+1):(j*2),t])      
    }
  }
}
tcwvcs2=array(0,dim=c(90,45,204))
lo=seq(1,179,by=2)
la=seq(1,89,by=2)
for(i in 1:length(lo)){
  for(j in 1:length(la)){
    for(t in 1:204){
      tcwvcs2[i,j,t]=mean(tcwvcs1[((i-1)*2+1):(i*2),((j-1)*2+1):(j*2),t])      
    }
  }
}
cnt=1
for(i in 1:180){
  for(j in 1:90){
    if(!is.na(laics1[i,j,1]) & (sum(laics1[i,j,])!=0)){
      tmp=c(i,j)
      if(cnt==1){
        laiind=tmp
      }else{
        laiind=rbind(laiind,tmp)
      }
      cnt=cnt+1
    }
  }
}
laiuse=array(0,dim=c(length(laiind[,1]),dim[3]))
tcwvind=which(!is.na(tcwvcs2[,,1]),arr.ind=TRUE)
tcwvuse=array(0,dim=c(length(tcwvind[,1]),dim[3]))
for(i in 1:dim[3]){
  #laiuse[,i]=laics1[,,i][laiind]
  tcwvuse[,i]=tcwvcs2[,,i][tcwvind]
}
memory.limit(size = 56000)
#Clai=laiuse %*% t(laiuse)/204
Ctcwv=tcwvuse %*% t(tcwvuse)/204
#eiglai<- eigen(Clai)
eigtcwv<- eigen(Ctcwv)
#eigveclai <- eiglai$vectors[,1:3]
eigvectcwv <- eigtcwv$vectors[,1:3]
#laiPCs <- t(eigveclai)%*%laiuse
tcwvPCs <- t(eigvectcwv)%*%tcwvuse
dim=dim(laics1)
#laiEOF1 <- array(NA,dim=c(dim[1],dim[2]))
#laiEOF2 <- laiEOF1
#laiEOF3 <- laiEOF1
tcwvEOF1 <- array(NA,dim=c(90,45))
tcwvEOF2 <- array(NA,dim=c(90,45))
tcwvEOF3 <- array(NA,dim=c(90,45))
for(i in 1:length(tcwvind[,1])){
  tcwvEOF1[tcwvind[i,1],tcwvind[i,2]]=eigvectcwv[i,1]
  tcwvEOF2[tcwvind[i,1],tcwvind[i,2]]=eigvectcwv[i,2]
  tcwvEOF3[tcwvind[i,1],tcwvind[i,2]]=eigvectcwv[i,3]
}
for(i in 1:length(laiind[,1])){
  laiEOF1[laiind[i,1],laiind[i,2]]=eigveclai[i,1]
  laiEOF2[laiind[i,1],laiind[i,2]]=eigveclai[i,2]
  laiEOF3[laiind[i,1],laiind[i,2]]=eigveclai[i,3]
}
for(j in 1:45){
  for(i in 1:90){
    tcwvEOF1[i,j]=eigvectcwv[(j-1)*90+i,1]
    tcwvEOF2[i,j]=eigvectcwv[(j-1)*90+i,2]
    tcwvEOF3[i,j]=eigvectcwv[(j-1)*90+i,3]
  }
}
write.table(laiEOF1,"laiEOF1.txt",col.names=FALSE,row.names=FALSE)
write.table(laiEOF2,"laiEOF2.txt",col.names=FALSE,row.names=FALSE)
write.table(laiEOF3,"laiEOF3.txt",col.names=FALSE,row.names=FALSE)
write.table(tcwvEOF1,"tcwvEOF1.txt",col.names=FALSE,row.names=FALSE)
write.table(tcwvEOF2,"tcwvEOF2.txt",col.names=FALSE,row.names=FALSE)
write.table(tcwvEOF3,"tcwvEOF3.txt",col.names=FALSE,row.names=FALSE)
write.table(laiPCs,"laiPCs.txt",col.names=FALSE,row.names=FALSE)
write.table(tcwvPCs,"tcwvPCs.txt",col.names=FALSE,row.names=FALSE)
######Start to plot the figures#######
library(lattice)
library(raster)
library(rasterVis)
library(maps)
library(mapdata)
library(maptools)
rotate <- function(x) t(apply(x, 2, rev))
plotmat2 <- function(LAIEOF,x){
  r <-raster(rotate(LAIEOF[180:1,]),  xmn=-180, xmx=180, ymn=-90, ymx=90, 
             crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  ext <- as.vector(extent(r))
  boundaries <- map('worldHires', fill=TRUE,
                    xlim=ext[1:2], ylim=ext[3:4],
                    plot=FALSE)
  IDs <- sapply(strsplit(boundaries$names, ":"), function(x) x[1])
  bPols <- map2SpatialPolygons(boundaries, IDs=IDs,
                               proj4string=CRS(projection(r)))
  #rgb.palette <- colorRampPalette(c("red", "yellow","green"))
  levelplot(r,margin=FALSE,scales=list(cex=3),at=seq(-0.08,0.08,by=0.0002),col.regions=colorRampPalette(c("red", "yellow", "green"))(1000),
            xlab=list(cex=3),ylab=list(cex=3),main=list(x,cex=3),
            colorkey=list(labels=list(cex=3))) + layer(sp.polygons(bPols))
}
plotmat3 <- function(TCWVEOF,x){
  r <-raster(rotate(TCWVEOF[90:1,]),  xmn=-180, xmx=180, ymn=-90, ymx=90, 
             crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  ext <- as.vector(extent(r))
  boundaries <- map('worldHires', fill=TRUE,
                    xlim=ext[1:2], ylim=ext[3:4],
                    plot=FALSE)
  IDs <- sapply(strsplit(boundaries$names, ":"), function(x) x[1])
  bPols <- map2SpatialPolygons(boundaries, IDs=IDs,
                               proj4string=CRS(projection(r)))
  #rgb.palette <- colorRampPalette(c("red", "yellow","green"))
  levelplot(r,margin=FALSE,scales=list(cex=3),at=seq(-0.08,0.08,by=0.0002),col.regions=colorRampPalette(c("red", "yellow", "blue"))(1000),
            xlab=list(cex=3),ylab=list(cex=3),main=list(x,cex=3),
            colorkey=list(labels=list(cex=3))) + layer(sp.polygons(bPols))
}
LAI1=laiEOF1[c(91:180,1:90),]
LAI2=laiEOF2[c(91:180,1:90),]
LAI3=laiEOF3[c(91:180,1:90),]
jpeg("Figures/LAIEOF1.jpg",width=1200,height=600)
plotmat2(LAI1,"LAI EOF1 96.45%")
dev.off()
jpeg("Figures/LAIEOF2.jpg",width=1200,height=600)
plotmat2(LAI2,"LAI EOF2 2.12%")
dev.off()
jpeg("Figures/LAIEOF3.jpg",width=1200,height=600)
plotmat2(LAI3,"LAI EOF3 0.54%")
dev.off()
TCWV1=tcwvEOF1[c(46:90,1:45),]
TCWV2=tcwvEOF2[c(46:90,1:45),]
TCWV3=tcwvEOF3[c(46:90,1:45),]
jpeg("Figures/TCWVEOF1.jpg",width=1200,height=600)
plotmat3(TCWV1,"TCWV EOF1 ")
dev.off()
jpeg("Figures/TCWVEOF2.jpg",width=1200,height=600)
plotmat3(TCWV2,"TCWV EOF2 ")
dev.off()
jpeg("Figures/TCWVEOF3.jpg",width=1200,height=600)
plotmat3(TCWV3,"TCWV EOF3 ")
dev.off()
laiPCs=as.matrix(read.table("laiPCs.txt",header=FALSE))
jpeg("Figures/LAIPC1.jpg",width=1200,height=600)
plot(laiPCs[1,],type='l',lwd=2,xlab="Grow season months since May 1982",ylab="LAI PC1",main="LAI PC1")#,cex.lab=3,cex.axis=3)
abline(v=114,col="red",lwd=2)
mtext("yr2001",side=3,at=114)
dev.off()
jpeg("Figures/LAIPC2.jpg",width=1200,height=600)
plot(laiPCs[2,],type='l',lwd=2,xlab="Grow season months since May 1982",ylab="LAI PC2",main="LAI PC2")#,cex.lab=3,cex.axis=3)
abline(v=114,col="red",lwd=2)
mtext("yr2001",side=3,at=114)
dev.off()
jpeg("Figures/LAIPC3.jpg",width=1200,height=600)
plot(laiPCs[3,],type='l',lwd=2,xlab="Grow season months since May 1982",ylab="LAI PC3",main="LAI PC3")#,cex.lab=3,cex.axis=3)
abline(v=114,col="red",lwd=2)
mtext("yr2001",side=3,at=114)
dev.off()
jpeg("Figures/tcwvPC1.jpg",width=1200,height=600)
plot(tcwvPCs[1,],type='l',lwd=2,xlab="Grow season months since May 1982",ylab="TCWV PC1",main="TCWV PC1")#,cex.lab=3,cex.axis=3)
abline(v=114,col="red",lwd=2)
mtext("yr2001",side=3,at=114)
dev.off()
jpeg("Figures/tcwvPC2.jpg",width=1200,height=600)
plot(tcwvPCs[2,],type='l',lwd=2,xlab="Grow season months since May 1982",ylab="TCWV PC2",main="TCWV PC2")#,cex.lab=3,cex.axis=3)
abline(v=114,col="red",lwd=2)
mtext("yr2001",side=3,at=114)
dev.off()
jpeg("Figures/tcwvPC3.jpg",width=1200,height=600)
plot(tcwvPCs[3,],type='l',lwd=2,xlab="Grow season months since May 1982",ylab="TCWV PC3",main="TCWV PC3")#,cex.lab=3,cex.axis=3)
abline(v=114,col="red",lwd=2)
mtext("yr2001",side=3,at=114)
dev.off()

