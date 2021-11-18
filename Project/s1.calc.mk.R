setwd("/home/zban/Downloads/2019spring/AOS219/project/")
###Step1: Read in TCWV and LAI data as 3D matrix
library(ncdf4)
LAInc=nc_open("1982-2015.grow.LAI.nc")
LAI=ncvar_get(LAInc,"Band1")
lon=ncvar_get(LAInc,"longitude")
lat=ncvar_get(LAInc,"latitude")
nc_close(LAInc)
TCWVnc=nc_open("1982-2015.grow.TCWV.nc")
TCWV=ncvar_get(TCWVnc,"tcwv")
###set <0 LAI to be 0
ind=which(LAI<0, arr.ind=TRUE)
LAI[ind]=0
LAI=LAI/100
###Step2: start modified mk test
library(modifiedmk)
dim=dim(TCWV)
LAIp=array(NA,dim=c(dim[1],dim[2]))
LAIslope=array(NA,dim=c(dim[1],dim[2]))
LAIslopesig=LAIslope
TCWVp=LAIp
TCWVslope=LAIslope
TCWVslopesig=TCWVslope
for(i in 1:dim[1]){
  for(j in 1:dim[2]){
    if((sum(LAI[i,j,])!=0) & (!is.na(sum(LAI[i,j,])))){
      tmp=mmky(LAI[i,j,])
      if((tmp[2]<=0.05) | (tmp[2]>=0.95)){
        LAIslopesig[i,j]=tmp[7]
      }
      LAIp[i,j]=tmp[2]
      LAIslope[i,j]=tmp[5]
    }
    if(!is.na(sum(TCWV[i,j,]))){
      tmp=mmky(TCWV[i,j,])
      if((tmp[2]<=0.05)|(tmp[2]>=0.95)){
        TCWVslopesig[i,j]=tmp[7]
      }
      TCWVp[i,j]=tmp[2]
      TCWVslope[i,j]=tmp[5]
    }
  }
}
write.table(LAIp,"LAIp.txt",col.names=FALSE,row.names=FALSE)
write.table(LAIslope,"LAIslope.txt",col.names=FALSE,row.names=FALSE)
write.table(LAIslopesig,"LAIslopesig.txt",col.names=FALSE,row.names=FALSE)
write.table(TCWVp,"TCWVp.txt",col.names=FALSE,row.names=FALSE)
write.table(TCWVslope,"TCWVslope.txt",col.names=FALSE,row.names=FALSE)
write.table(TCWVslopesig,"TCWVslopesig.txt",col.names=FALSE,row.names=FALSE)
#library(lattice)
#library(raster)
#library(rasterVis)
#library(maps)
#library(mapdata)
#library(maptools)
#rotate <- function(x) t(apply(x, 2, rev))
#plot1 <- levelplot(LAIp, )
#plotmat <- function(MAT1,x){
#  r <-raster(rotate(MAT1[180:1,]),  xmn=0, xmx=360, ymn=-90, ymx=90, 
#             crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
#  ext <- as.vector(extent(r))
#  boundaries <- map('worldHires', fill=TRUE,
#                    xlim=ext[1:2], ylim=ext[3:4],
#                    plot=FALSE)
#  IDs <- sapply(strsplit(boundaries$names, ":"), function(x) x[1])
#  bPols <- map2SpatialPolygons(boundaries, IDs=IDs,
#                               proj4string=CRS(projection(r)))
#  levelplot(r,margin=FALSE,main=paste0(x)) + layer(sp.polygons(bPols))
#}
