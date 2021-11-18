setwd("C:/Users/BANZH/Downloads/2019spring/AOS 219")
nino=read.table("NINO34.txt")
nino=nino[1:70,]
library(ncdf4)
nc=nc_open("outhgt.nc")
hgt=ncvar_get(nc,"hgt")
hgt=hgt[,,24:842]  #Cut the hgt data from 194912-291802
hgtdjf=array(0,dim=c(41,31,69))
for(i in 1:69){
  hgtdjf[,,i]=rowMeans(hgt[,,((i-1)*12+1):((i-1)*12+3)],dim=2)
}
hgtdjfmean=rowMeans(hgtdjf,dim=2)
hgtdjfanom=array(0,dim=c(41,31,69))
for(i in 1:69){
  hgtdjfanom[,,i]=hgtdjf[,,i]-hgtdjfmean
}
ninots=array(0,dim=c(69,3))
ninodjf=rep(0,69)
for(i in 1:69)
{
  ninots[i,1]=nino[i,13]
  ninots[i,2]=nino[i+1,2]
  ninots[i,3]=nino[i+1,3]
  ninodjf[i]=mean(ninots[i,])
}
acfnino=acf(ninodjf,lag.max=69)
CHH=as.numeric(acfnino$acf)
CSS=array(0,dim=c(41,31,69))
for(i in 1:41){
  for(j in 1:31){
    acftmp=acf(hgtdjfanom[i,j,],lag.max=69,plot=FALSE)
    CSS[i,j,]=as.numeric(acftmp$acf)
  }
}
tau=array(0,dim=c(41,31))
for(i in 1:41){
  for(j in 1:31){
    tau[i,j]=1+2*sum(CHH*CSS[i,j,])
  }
}
Neff=array(0,dim=c(41,31))
Neff=69/tau  ##Neff
getp <- function(ts1,ts2,Neff){
  r=cor(ts1,ts2)
  t=r*sqrt(Neff-2)/sqrt(1-r*r)
  p=pnorm(t)
  return(p)
}
testp=rep(0,41*31)
for(i in 1:41){
  for(j in 1:31){
    testp[(i-1)*31+j]=getp(ninodjf,hgtdjfanom[i,j,],Neff[i,j])
  }
}
threshold=rep(0,41*31)
N=41*31
alpha=0.05
for(i in 1:(41*31)){
  threshold[i]=alpha*i/N
}
sortp=sort(testp)
index=seq(1,1271,by=1)
plot(x=index[1:10],y=sortp[1:10],col="red",xlab="Rank",ylab="p value")
lines(threshold[1:10],col="black",lwd=2,lty=2)
legend("bottomright",c("sort p","FDR threshold"),lwd=c(NA,2),pch=c(1,NA),col=c("red","black"),lty=c(NA,2),cex=1.5)

##Chen's method
check <- function(ts1,ts2){
  test=cor.test(ts1,ts2)
  if(((test$p.value)<0.025) || (test$p.value>0.975)){
    return(1)
  }else{
    return(0)
  }
}
check1 <- function(ts1,ts2,Neff){
  r=cor(ts1,ts2)
  t=r*sqrt(Neff-2)/sqrt(1-r*r)
  p=pnorm(t)
  if((p<0.025) || (p>0.975)){
    return(1)
  }else{
    return(0)
  }
}
gau=rnorm(69,mean=0,sd=1)
test=array(0,dim=c(41,31))
testran=test
for(i in 1:41){
  for(j in 1:31){
    test[i,j]=check1(ninodjf,hgtdjfanom[i,j,],Neff[i,j])
    testran[i,j]=check(gau,hgtdjfanom[i,j,])
  }
}
sum=sum(test)
sumran=sum(testran)
ratio=sum/N
ratioran=sumran/N
nc_close(nc)
