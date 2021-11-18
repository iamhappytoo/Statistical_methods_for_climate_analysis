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
install.packages(c("modifiedmk"))
library(modifiedmk)
mmky(rf)
plot(rf)
#   This procedure computes the nonparametric estimator
#   of the long run variance using Andrews (1991) AR(1)
#   plug in optimal truncation lag.  v is the vector of
#   residuals and M is the trucation lag.  M<0 uses the
#   automatic bandwidth.  Kernel is an integer representing
#   the kernel that is used
#   1 = bartlett, 2 = parzen, 3 = Tukey-Hanning, 4 = quadratic spectral

#   if prewhite == 1 then prewhitening using an AR(1) model is used
#   if prewhite == 0 no prewhitening is used.

#   8/13/96


sig2np <- function(v, M, my.kernel, prewhite){
  
  #    v = uhat;
  #    prewhite = 0;
  #    M = -1
  T = nrow(v);
  
  if(exists("my.kernel")) {
    if(my.kernel < 0 | my.kernel > 5){
      my.kernel=4;
      cat("using the Quadratic spectral kernel\n");
    }
  }
  if(!exists("my.kernel")) {
    my.kernel = 4;
    cat("using the Quadratic spectral kernel\n");
  }
  
  R = matrix(0,nrow = (T-1), ncol = 1);
  R0 = t(v[1:T-1]) %*% v[1:T-1];
  R[1]=t(v[2:T]) %*% v[1:T-1];
  rho=R[1]/R0;
  
  if(prewhite == 1){
    rho2 = rho;
    v = cbind(c(0,(v[2:T]-rho*v[1:T-1])))
    R0 = t(v[1:T-1]) %*% v[1:T-1];
    R[1]=t(v[2:T]) %*% v[1:T-1];
    rho=R[1]/R0;
  }
  
  if(M < 0.0){
    a1 = 4.0*rho*rho/(1.0-rho*rho)^2;
    a2 = 4.0*rho*rho/(1.0-rho)^4;
    if(my.kernel == 1) { ST=1.1447*(T*a1)^(1.0/3.0) }
    if(my.kernel == 2) { ST=2.6614*(T*a2)^(0.2) }
    if(my.kernel == 3) { ST=1.7462*(T*a2)^(0.2) }
    if(my.kernel == 4) { ST=1.3221*(T*a2)^(0.2) }
    if(my.kernel == 5) { ST=1 }
  }
  if(M>=0.0){
    ST=M;
  }
  R0 = t(v) %*% v / T;
  
  
  i = 1;
  for(i in 1:(T-1)){
    R[i,1] = t(v[(i+1):T]) %*% v[1:(T-i)] / T
  }
  
  sigma2 = R0 + 2.0 * t(K(seq(1,T-1)/ST, my.kernel)) %*% R;
  
  if(prewhite != 1) { return(sigma2) }
  if(prewhite == 1) { return(sigma2/(1.0-rho2)^2) }
}



# Kernels
K = function(x, my.kernel){
  
  x = as.matrix(x)
  x = abs(x);
  xx = ((x > 0) * x) + (x < 0.0000001)
  
  if(my.kernel == 1) { # bartlett
    y= (x <= 1.0) * (1.0 - x);
  }
  if(my.kernel == 2) { # parzen
    y = (x >= 0.0) * (x <= 0.5) * (1.0 - 6.0 * x * x * (1.0 - x));
    y = y + (x > 0.5) * (x <= 1.0) * 2.0 * (1.0 - x)^3;
  }
  if(my.kernel == 3) { # Tukey-hanning
    y = (x <= 1.0) * 0.5 * (1.0 + cos(pi * x));
  }
  if(my.kernel == 4) { # quadratic spectral
    y = (x > 0) * (25.0 / (12.0 * pi^2 * x * x)) * ((sin(1.2 * pi * x) / (1.2 * pi * x)) - cos(1.2 * pi * x));
    y = y + (x == 0.0);
  }
  if(my.kernel == 5) { # Daniells
    y = ((x > 0) * sin(pi * xx) / (pi * xx)) + (abs(x) < 0.00000001);
  }
  return(y);
}

v.trend <- function(y) {
  prewhite = 0; # set to 1 to use prewhitening for the t-dan test
  
  # critical values for right tailed 5% 2.5% and 1% test and one-sided
  #   confidence intervals
  cvtps =  c(1.331, 1.720, 2.152, 2.647);
  btps =   c(0.494, 0.716, 0.995, 1.501);
  #   cvtdan = c(1.710, 2.052, 2.462);
  #   btdan =  c(1.32227, 1.79541, 2.46582);
  #   Jcv =    c(0.488, 0.678, 0.908);
  
  z = cumsum(y);
  T = nrow(y);
  t1 = 1:nrow(y)
  X1 = cbind(rep(1, nrow(y)), t1)
  
  # Compute the partial sums of X1
  nr = nrow(X1);
  nc = ncol(X1);                                           
  X2 = apply(X1, 2, cumsum);
  
  # Compute inverse matrices
  
  X1inv = as.matrix(solve(pdMat(t(X1) %*% X1)))
  X2inv = as.matrix(solve(pdMat(t(X2) %*% X2)))
  
  # Compute orthonormal trends for the J statistics
  Xj = matrix(1, nrow = nr) / sqrt(nr)
  i = 1;
  for(i in 1:9){
    tt = t1^i
    ehat = tt - Xj %*% as.matrix(solve(pdMat(t(Xj) %*% Xj))) %*% (t(Xj) %*% tt);
    ehat = ehat / as.vector(sqrt(t(ehat) %*% ehat));
    Xj = cbind(Xj, ehat);
  }
  Xjinv = as.matrix(solve(pdMat(t(Xj) %*% Xj)))
  
  # Compute stats in standard regressions
  bhat = X1inv %*% t(X1) %*% as.matrix(y);
  uhat = y - X1 %*% bhat;
  rss1 = t(uhat) %*% uhat;
  
  s2dan = sig2np(v = uhat, M = max(trunc(0.02*T),2), my.kernel = 1, prewhite = 0);
  #Kernel 1 is chosen based on Formby and Vogelsang (2002) which uses the Bartlett kernel
  
  
  btild = X2inv %*% t(X2) %*% z;
  jbeta = Xjinv %*% t(Xj) %*% y;
  rssj = t((y - Xj %*% jbeta)) %*% (y - Xj %*% jbeta);
  s2z = t((z - X2 %*% btild)) %*% (z - X2 %*% btild) / (nrow(X2) - ncol(X2));
  J = (rss1-rssj) / rssj;
  
  #   tdan = (bhat[2] / sqrt(s2dan %*% X1inv[2,2])) %*% exp(-btdan * J);
  #   tps = (btild[2] / sqrt(T*s2z*X2inv[2,2])) %*% exp(-btps * J);
  tps = (btild[2]/sqrt(T*s2z*X2inv[2,2]))*exp(-btps*J);
  
  # Write bhat, sig level for tdan and tps,
  tps.100 <- ifelse(tps[1] > cvtps[1] | tps[1] < -1*cvtps[1], 1, 0)
  tps.050 <- ifelse(tps[2] > cvtps[2] | tps[1] < -1*cvtps[2], 1, 0)
  tps.025 <- ifelse(tps[3] > cvtps[3] | tps[2] < -1*cvtps[3], 1, 0)
  tps.010 <- ifelse(tps[4] > cvtps[4] | tps[3] < -1*cvtps[4], 1, 0)
  data2export <- c(bhat[2], tps.100, tps.050, tps.025, tps.010)
  return(data2export)
}

#############################################################
#############################################################
# Start Here
# Do Vogelsang time series analysis
#download nlme package
install.packages(c("nlme"))
library(nlme)
#inputTs = 'read your input time series here';
inputTs<-read.delim("/Users/rfu/Desktop/Rongsfiles/AOS-219-statistics-lecture/HW assignments/AOS219-HW1/DFJLArainfall.txt").  # you need to replace this path by your own path to the data
dataMat = as.matrix(rf,nrow=length(inputTs),ncol=1); # reformat the inputTs as Matrix

results <- v.trend(dataMat);# run the vogelsang test
trend <- results[1]; # get the trend value (slope)

# It should be noted that there are only two possible values for the output pValue: 0 and 1. 
# pValue equals to 0 means the trend is not significant at the corresponding confidence level, i.e. cannot pass the test
# pValue equals to 1 means the trend is significant at the corresponding cofindence level, i.e. pass the test
pValue.100 <- results[2]; # get the pValue at 0.1 confidence level
pValue.050 <- results[3]; # get the pValue at 0.05 confidence level
pValue.025 <- results[4]; # get the pValue at 0.025 confidence level
pValue.010 <- results[5]; # get the pValue at 0.01 confidence level