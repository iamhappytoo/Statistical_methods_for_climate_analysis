setwd("C:/Users/BANZH/Downloads/2019spring/AOS 219")
data=read.table("LA.txt",header = FALSE)
len=2018-1981+1
rf=rep(0,len)
for(i in 1981:2018){
  y1=i-1
  y=i
  ind1=which(data[,1]==y1)
  ind=which(data[,1]==y)
  rf[i-1980]=sum(c(data[ind1[3],3],data[ind[1:2],3]))
}
##Question 1##
##Fit expected Gamma distribution
rf1=rf[which(rf!=0)]
meanrf1=mean(rf1)
D=log(meanrf1)-1/length(rf1)*sum(log(rf1))
alpha=(1+sqrt(1+4*D/3))/(4*D)
beta=meanrf1/alpha
xbar=alpha*beta
s2=alpha*beta*beta
sd=sqrt(s2/length(rf1))
##There are 38 years, for the chi-square method to be valid, at least each bin has 5 frequency, so #bin=38/5~8
binwd=sd
lwd=xbar-4*sd
uwd=xbar+4*sd
x=seq(lwd,uwd,sd)
N=len
No=rep(0,8)
Negau=rep(0,8)
Negam=rep(0,8)
for(j in 1:8){
  No[j]=length(which((rf>x[j])&(rf<=x[j+1])))
  Negau[j]=(pnorm((x[j+1]-xbar)/sd)-pnorm((x[j]-xbar)/sd))*N
  Negam[j]=(pgamma(x[j+1],shape=alpha,scale=beta)-pgamma(x[j],shape=alpha,scale=beta))*N
}
chisgam=sum((No-Negam)^2/Negam)
chisgau=sum((No-Negau)^2/Negau)
df=8-2-1=5 ##Degree of freedom for both Gaussian and Gamma
##At 95% confidence##
print(paste0("At 99% confidence, null hypothesis for gaussian distribution is rejected, with X2=",chisgau,"> X2_99%=15.086"))
print(paste0("At 99% confidence, null hypothesis for gamma distribution is accepted, with X2=",chisgam,"< X2_99%=15.086"))
##Therefore, at 99% confidence, the rainfall data can be accepted as following Gamma distribution
###Question 2###
lanina=c(1985,  1989,   1999,   2008,   2009,   2011)
select <- function(data,l){
  y1=l-1
  y=l
  ind1=which(data[,1]==y1)
  ind=which(data[,1]==y)
  tmp=sum(c(data[ind1[3],3],data[ind[1:2],3]))
  return(tmp)
}
data2=rep(0,length(lanina))
for(i in 1:length(lanina)){
  data2[i]=select(data,lanina[i])
}
d=unique(data[,1])
nonlanina=d[!(d %in% lanina)]
lenn=length(nonlanina)
nonlanina=nonlanina[2:lenn] ##delete the year 1980
data3=rep(0,length(nonlanina))
for(i in 1:length(nonlanina)){
  data3[i]=select(data,nonlanina[i])
}
N1=length(data2)
N2=length(data3)
#Use t test
x1bar=mean(data2)
x2bar=mean(data3)
s1sq=var(data2)
s2sq=var(data3)
sigma1=sqrt((N1*s1sq+N2*s2sq)/(N1+N2-2))
t=(x1bar-x2bar)/(sigma1*sqrt(1/N1+1/N2))
print(paste0("t=",t))
N=N1+N2-2
print(paste0("N=",N))
##For 95% confidence level and 37 samples, t0.025=1.960>abs(t), therefore, cannot reject null hypothesis
##Not significantly different from other years.
###Smirnov two-sample test
srt_l=sort(data2)
srt_nl=sort(data3)
ecdfl=ecdf(srt_l)
ecdfnl=ecdf(srt_nl)
rf_l=ecdfl(rf)
rf_nl=ecdfnl(rf)
plot(x=sort(rf),y=sort(rf_l),type='l',lwd=2,xlab="Sorted or ranked Rainfall value (inches)",ylab="CDF",pch="+",cex=2)
points(x=sort(rf),y=sort(rf_l),lwd=2,pch="+",cex=2)
lines(x=sort(rf),y=sort(rf_nl),lwd=2,col="blue",type="o",cex=2)
legend("bottomright",lty=1,pch=c("+","o"),col=c("black","blue"),lwd=2,legend=c("lanina","non-lanina"),cex=c(2,2))
Ds=max(abs(sort(rf_l)-sort(rf_nl)))
Dsc=sqrt(-0.5*(1/N1+1/N2)*log(0.05/2))
print(paste0("Because Ds=",Ds,"< Dsc=",Dsc,", cannot reject the null hypothesis at 0.95 confidence level."))
##Question 3##
srt_rf=sort(rf)
rankl=rep(0,N1)
ranknl=rep(0,N2)
for(i in 1:N1){
  rankl[i]=mean(which(srt_rf==srt_l[i]))
}
for(i in 1:N2){
  ranknl[i]=mean(which(srt_rf==srt_nl[i]))
}
R1=sum(rankl)
R2=sum(ranknl)
U1=N1*N2+N1*(N1+1)/2-R1
U2=N1*N2+N2*(N2+1)/2-R2
U=min(U1,U2)
mu=0.5*N1*N2
t=2 ##tie 
Tc=t^3-t
N=N1+N2
sigma=sqrt(N1*N2*(N^3-N-Tc)/(12*N*(N-1)))
z=1.96 #for a 0.05 significance level
Uc=mu-z*sigma-0.5
print(paste0("U=",U,">Uc=",Uc,"cannot reject the null hypothesis that the lanina and non-lanina are from the same distribution"))
###Question 4###
S=0
for(i in 1:(N-1)){
  S=S+length(which(rf[(i+1):N]>rf[i]))-length(which(rf[(i+1):N]<rf[i]))
}
varS=(N*(N-1)*(2*N-5)-2*(2-1)*(2*2-5)*2)/18
z=(S+1)/sqrt(varS)
##For p<0.05 two tailed test, Z0.025=-1.96
print(paste0("|z|<=Z_1-p/2 when p=0.05, because |z|=",abs(z),", and Z_1-p/2=",qnorm(1-0.05/2),"cannot reject the null hypothesis"))
print("Therefore, the trend is not significant under 95% confidence.")
