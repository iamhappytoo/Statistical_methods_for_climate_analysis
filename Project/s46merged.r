setwd("/home/zban/Downloads/Master_thesis/Colosubbasin_experiment/data/awco/")
awcodata=array(0,dim=c(15,4))
for(i in 1:15){
  data=as.numeric(as.matrix(read.table(paste0(i,"/ROawco/annual_awco"),header=FALSE)))
  awcodata[i,]=data
}
rf=awcodata[,4]
total=sum(rf)
ratio=rf/total
ratio=as.data.frame(ratio)
ind=seq(1:15)
ratio$ind=paste(rep("basin",15),seq(1:15))
###plot a column figure.
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
png("/home/zban/Downloads/Master_thesis/Colosubbasin_experiment/plot/Figures/s4.colo.contribution.png")
plot4 <- ratio %>%
  mutate(ind = reorder(ind, ratio)) %>%
  ggplot( aes(x=ind, y=ratio, fill=ind)) +
  geom_bar(stat="identity") +
  xlab("basin index") +
  ylab("RF contribution ratio")+
  coord_flip()
dev.off()
setwd("/home/zban/Downloads/Master_thesis/Colosubbasin_experiment/data/awco")
awcodata=array(0,dim=c(15,3,4))
for(i in 1:15){
  data=as.numeric(as.matrix(read.table(paste0(i,"/ROawco/annual_awco"),header=FALSE)))
  awcodata[i,1,]=data
  data=as.numeric(as.matrix(read.table(paste0(i,"/ROawco/warm_awco"),header=FALSE)))
  awcodata[i,2,]=data
  data=as.numeric(as.matrix(read.table(paste0(i,"/ROawco/cold_awco"),header=FALSE)))
  awcodata[i,3,]=data
}
annualratio=array(0,dim=c(15,3))
warmratio=array(0,dim=c(15,3))
coldratio=array(0,dim=c(15,3))
for(i in 1:15){
  tmpbase=awcodata[i,1,4]  ##for basin i, annual_awco, origin RO
  annualratio[i,1]=(awcodata[i,1,1]-awcodata[i,1,4])/tmpbase ##annual V rel change under annual warming
  annualratio[i,2]=(awcodata[i,1,2]-awcodata[i,1,4])/tmpbase ##annual V rel change under warm warming
  annualratio[i,3]=(awcodata[i,1,3]-awcodata[i,1,4])/tmpbase ##annual V rel change under cold warming
  warmratio[i,1]=(awcodata[i,2,1]-awcodata[i,2,4])/tmpbase   ##warm seas V rel change under annual warming
  warmratio[i,2]=(awcodata[i,2,2]-awcodata[i,2,4])/tmpbase   ##warm seas V rel change under warm warming
  warmratio[i,3]=(awcodata[i,2,3]-awcodata[i,2,4])/tmpbase   ##warm seas V rel change under cold warming
  coldratio[i,1]=(awcodata[i,3,1]-awcodata[i,3,4])/tmpbase
  coldratio[i,2]=(awcodata[i,3,2]-awcodata[i,3,4])/tmpbase
  coldratio[i,3]=(awcodata[i,3,3]-awcodata[i,3,4])/tmpbase
}
colnames(annualratio)=c("annual3d","warm3d","cold3d")
rownames(annualratio)=paste( rep("basin",15) , c(1:15) , sep=" ")
annualratio=as.data.frame(annualratio)
annualratio$ratio=ratio$ratio
annualratio=annualratio[order(annualratio$ratio,decreasing = TRUE),]
colnames(warmratio)=c("annual3d","warm3d","cold3d")
rownames(warmratio)=paste( rep("basin",15) , c(1:15) , sep=" ")
warmratio=as.data.frame(warmratio)
warmratio$ratio=ratio$ratio
warmratio=warmratio[order(warmratio$ratio,decreasing = TRUE),]
colnames(coldratio)=c("annual3d","warm3d","cold3d")
rownames(coldratio)=paste( rep("basin",15) , c(1:15) , sep=" ")
coldratio=as.data.frame(coldratio)
coldratio$ratio=ratio$ratio
coldratio=coldratio[order(coldratio$ratio,decreasing=TRUE),]

library(lattice)
jpeg("/home/zban/Downloads/Master_thesis/Colosubbasin_experiment/plot/Figures/s6.seaschange.jpg",width=1500,height=600)
plot1 <- levelplot(t(as.matrix(annualratio[c(15:1) ,1:3])),xlab="warming scenarios",ylab="basin index",main = list('Annual streamflow change',side=1,line=0.5),
                   at=seq(-0.25,0.25,by=0.001), col.regions = colorRampPalette(c("red","white","blue")))
plot2 <- levelplot(t(as.matrix(warmratio[c(15:1) , 1:3])),xlab="warming scenarios",ylab="",main = list('warm season streamflow change',side=1,line=0.5),
                   at=seq(-0.25,0.25,by=0.001), col.regions = colorRampPalette(c("red","white","blue")))
plot3 <- levelplot(t(as.matrix(coldratio[c(15:1) , 1:3])),xlab="warming scenarios",ylab="",main = list('cold season streamflow change',side=1,line=0.5),
                   at=seq(-0.07,0.07,by=0.001), col.regions = colorRampPalette(c("red","white","blue")))
library(gridExtra)
library(grid)
grid.arrange(arrangeGrob(plot1,plot2,plot3, nrow=1, ncol=3),
             arrangeGrob(plot4, nrow=1, ncol=1), widths=c(1.5,1)) 
dev.off()
