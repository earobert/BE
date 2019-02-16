## ============================================================================
## FISH 454 Ecological Modeling - Winter 2014
## Lab #3 - Monte Carlo simulations
## January 2014
## ============================================================================

setwd(" NEED TO SET THIS") # set your working directory to the location you 
                                # are storing your files
source("Lab3_2function.R") 
par(ask=F)

MC.func <- function(nsims=500) {
  output <- rep(NA, nsims)
  a.vec <- runif(nsims, -1.648, -0.412)
  b.vec <- runif(nsims, -0.688, -0.172)
  alpha.vec <- runif(nsims, 0.38, 1.5)
  ar.vec <- runif(nsims, 10, 40)
  M.vec <- runif(nsims, 0.05, 0.15)
  lambda.vec <- exp(a.vec+b.vec*1.35)
    
  for(i in 1:nsims) {
    temp <- snake.func(N.start,  N.target, Yp, a=a.vec[i], b=b.vec[i], alpha=alpha.vec[i], ar=ar.vec[i], br, N.max, N.min, damage, M=M.vec[i])
    output[i] = temp[[2]]
  }   

  output.mat <- cbind(output, a.vec, b.vec, lambda.vec, alpha.vec, ar.vec, M.vec)
  colnames(output.mat) <- c("Annual cost", "a", "b", "lambda", "alpha", "ar", "M")
  
  plot(output~a.vec, main = "Parameter=a", xlab = "Parameter value", ylab = "Annual Cost (USD Million)")
  plot(output~b.vec, main = "Parameter=b", xlab = "Parameter value", ylab = "Annual Cost (USD Million)")
  plot(output~lambda.vec, main = "Parameter=lambda", xlab = "Parameter value", ylab = "Annual Cost (USD Million)")
  plot(output~alpha.vec, main = "Parameter=alpha", xlab = "Parameter value", ylab = "Annual Cost (USD Million)")
  plot(output~ar.vec, main = "Parameter=ar", xlab = "Parameter value", ylab = "Annual Cost (USD Million)")
  plot(output~M.vec, main = "Parameter=M", xlab = "Parameter value", ylab = "Annual Cost (USD Million)")
  return(output.mat)
  
}



runMC <- MC.func(nsims = 1000)
# run this and type to view the plots of each (press 'Enter')

# to look at each plot individually:
#par(ask = T)
#plot.a <- plot(runMC[,1]~runMC[,2], main = "Parameter=a", xlab = "Parameter value", ylab = "Annual Cost (USD Million)")
#plot.b <- plot(runMC[,1]~runMC[,3], main = "Parameter=b", xlab = "Parameter value", ylab = "Annual Cost (USD Million)")
#plot.lambda <- plot(runMC[,1]~runMC[,4], main = "Parameter=lambda", xlab = "Parameter value", ylab = "Annual Cost (USD Million)")
#plot.alpha <- plot(runMC[,1]~runMC[,5], main = "Parameter=alpha", xlab = "Parameter value", ylab = "Annual Cost (USD Million)")
#plot.ar <- plot(runMC[,1]~runMC[,6], main = "Parameter=ar", xlab = "Parameter value", ylab = "Annual Cost (USD Million)")
#plot.M <- plot(runMC[,1]~runMC[,7], main = "Parameter=M", xlab = "Parameter value", ylab = "Annual Cost (USD Million)")




par(mfrow=c(4,4),oma=c(4,4,0,0),mar=c(4,2,1,1),las=1,cex.lab=1.5)
# make plots

max.matrix=matrix(apply(abs(runMC),2,max),nrow=nrow(runMC),ncol=7,byrow=TRUE)
std.pars<-runMC/max.matrix

for (i in 3:6){
  ylab.txt="Total Cost"
  main.txt.vector=c('a','b','lambda','alpha','ar','M')
  for (j in 3:6){
    if (i==j){
      plot(c(0,1),c(0,1),type='n',ylab='',xlab='',main='',axes=FALSE)
      text(0.33,0.5,labels=main.txt.vector[i],cex=3)
      
    } else {
      radius=std.pars[,j+1]^2
      symbols(runMC[,i+1],runMC[,1],circles=radius,inches=0.05,fg='blue',bg='blue',xlab=main.txt.vector[i],ylab='',xaxs="i",yaxs='i')
    }
    }
  }
par(las=0)  
mtext(outer=TRUE,side=2,text="Annual Cost",cex=2)

# END