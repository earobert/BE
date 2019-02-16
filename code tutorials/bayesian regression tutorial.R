#bayesian regression tutorial

## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
?gl # generate factor levels
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lmfit <- lm(weight ~ group)
summary(lmfit)
plot(weight ~ group)

plot(ctl, trt)

QR<-lmfit$qr
df.residual<-lmfit$df.residual
R<-qr.R(QR) ## R component
coef<-lmfit$coef
Vb<-chol2inv(R) ## variance(unscaled)
s2<-(t(lmfit$residuals)%*%lmfit$residuals)
s2<-s2[1,1]/df.residual

## function to compute the bayesian analog of the lmfit
## using non-informative priors and Monte Carlo scheme
## based on N samples

bayesfit<-function(lmfit,N){
  QR<-lmfit$qr
  df.residual<-lmfit$df.residual
  R<-qr.R(QR) ## R component
  coef<-lmfit$coef
  Vb<-chol2inv(R) ## variance(unscaled)
  s2<-(t(lmfit$residuals)%*%lmfit$residuals)
  s2<-s2[1,1]/df.residual
  
  ## now to sample residual variance
  sigma<-df.residual*s2/rchisq(N,df.residual)
  coef.sim<-sapply(sigma,function(x) mvrnorm(1,coef,Vb*x))
  ret<-data.frame(t(coef.sim))
  names(ret)<-names(lmfit$coef)
  ret$sigma<-sqrt(sigma)
  ret
}

Bayes.sum<-function(x)
{
  c("mean"=mean(x),
    "se"=sd(x),
    "t"=mean(x)/sd(x),
    "median"=median(x),
    "CrI"=quantile(x,prob=0.025),
    "CrI"=quantile(x,prob=0.975)
  )
}

set.seed(1234)  ## reproducible sim
lmfit <- lm(weight ~ group)
bf<-bayesfit(lmfit,10000)
?t
t(apply(bf,2,Bayes.sum))
# mean        se         t     median    Cr.2.5%  Cr.97.5%
#   (Intercept)  5.0332172 0.2336049 21.545857  5.0332222  4.5651327 5.4902380
# groupTrt    -0.3720698 0.3335826 -1.115375 -0.3707408 -1.0385601 0.2895787
# sigma        0.7262434 0.1275949  5.691789  0.7086832  0.5277051 1.0274083
summary(lmfit)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.0710 -0.4938  0.0685  0.2462  1.3690 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   5.0320     0.2202  22.850 9.55e-15 ***
#   groupTrt     -0.3710     0.3114  -1.191    0.249  
