setwd("C:/Users/BANZH/Downloads/2019spring/AOS 219")
nc=nc_open("sst.mnmean.v4.nc")
sst=ncvar_get(nc,"sst")
sstuse=sst[,,565:1824]
dim=dim(sstuse)
JAmean=array(0,dim=c(dim[1],dim[2],105))
for(i in 1:dim[1]){
      for(j in 1:dim[2]){
        for(z in 1:105){
          JAmean[i,j,z]=mean(sstuse[i,j,((z-1)*12+6):((z-1)*12+8)])
        }
      }
}
JAind=which(!is.na(JAmean[,,1]),arr.ind=TRUE)
JAuse=array(0,dim=c(length(JAind[,1]),105))
JAuse=JAmean[JAind,]