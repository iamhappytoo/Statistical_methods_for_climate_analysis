setwd("C:/Users/BANZH/Downloads/2019spring/AOS 219/AOS project")
LAIp=as.matrix(read.table("LAIp.txt",header=FALSE))
LAIslope=as.matrix(read.table("LAIslope.txt",header=FALSE))
LAIslopesig=as.matrix(read.table("LAIslopesig.txt",header=FALSE))
TCWVp=as.matrix(read.table("TCWVp.txt",header=FALSE))
TCWVslope=as.matrix(read.table("TCWVslope.txt",header=FALSE))
TCWVslopesig=as.matrix(read.table("TCWVslopesig.txt",header=FALSE))
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
  rgb.palette <- colorRampPalette(c("green", "yellow"))
  levelplot(r,col.regions=rgb.palette(120),margin=FALSE,scales=list(cex=3),xlab=list(cex=3),ylab=list(cex=3),main=list(x,cex=3),colorkey=list(labels=list(cex=3))) + layer(sp.polygons(bPols))
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
  levelplot(r,margin=FALSE,scales=list(cex=3),at=seq(-0.003,0.003,by=0.0002),col.regions=colorRampPalette(c("red", "yellow", "green"))(100),
            xlab=list("lon",cex=3),ylab=list("lat",cex=3),main=list(x,cex=3),
            colorkey=list(labels=list(cex=3))) + layer(sp.polygons(bPols))
}
plotmat3 <- function(LAIp,x){
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
  levelplot(r,margin=FALSE,scales=list(cex=3),at=seq(-0.02,0.02,by=0.0005),col.regions=colorRampPalette(c("red", "yellow", "blue"))(100),
            xlab=list(cex=3),ylab=list(cex=3),main=list(x,cex=3),
            colorkey=list(labels=list(cex=3))) + layer(sp.polygons(bPols))
}
jpeg("Figures/LAIp.jpg",width=1200,height=600)
LAI1=LAIp[c(360:720,1:359),]
plotmat2(LAI1,"LAI p value of trend")
dev.off()
jpeg("Figures/LAIslope.jpg",width=1200,height=600)
LAI2=LAIslope[c(360:720,1:359),]
plotmat2(LAI2, "LAI trend")
dev.off()
jpeg("Figures/LAIslopesig.jpg",width=1200,height=600)
LAI3=LAIslopesig[c(360:720,1:359),]
plotmat2(LAI3,"LAI trend at signficant location only")
dev.off()
jpeg("Figures/TCWVp.jpg",width=1200,height=600)
TCWV1=TCWVp[c(360:720,1:359),]
plotmat1(TCWV1,"TCWV p value of trend")
dev.off()
jpeg("Figures/TCWVslope.jpg",width=1200,height=600)
TCWV2=TCWVslope[c(360:720,1:359),]
plotmat3(TCWV2,"TCWV trend")
dev.off()
jpeg("Figures/TCWVslopesig.jpg",width=1200,height=600)
TCWV3=TCWVslopesig[c(360:720,1:359),]
plotmat3(TCWV3,"TCWV trend at significant location only")
dev.off()
