# eivtools used for the Autumn dataset
# here we are estimating the SE using error in the variables through a total least squares estimation (EIV)
# Note there is a single paper recommending this approach, so it is worth checking out other methods
# Not sure this is the best way to do this... it's adding a normal random component to the 45 data points
# I wonder if estimating without a would be better
# Mod2 at the end isn't working - it calculates 
# JAGS = just another gibbs sampler, MCMC
# estimate X's and B's


rm(list=ls())

# devtools::install_github("jrlockwood/eivtools")
# install.packages("R2jags")
# to get this to work I installed rjags from the sourceforge website AND R2jags!?
# it's odd b/c you install R2jags and maybe rjags, though rjags might not be necesary but then don't call the library R2jags
rm(list=ls())

library(eivtools)
#library(rjags)

p.Aut <- data.frame(
  names = c("a","b","d","e","C","del"),
  val = c(14,10.3,0.69,1,0.182,.333),
  se = c(7,4.8,.01,.0001,.01,.006),
  season = rep("Autumn", times = length(names))
)
p.Aut
p.Spr <- data.frame(
  names = c("a","b","d","e","C","del"),
  val = c(45,34.4,0.69,1,0.182,.295),
  se = c(14,9.2,.01,.0001,.01,.004),
  season = rep("Spring", times = length(names))
)
p.Spr
p.Aut <- data.frame(
  names = c("a","b","d","e","C","del"),
  val = c(14,10.3,0.69,1,0.182,.333),
  se = c(1,1,.01,.0001,.001,.001),
  season = rep("Autumn", times = length(names))
)
p.Aut

p <- p.Aut

setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="practice_data.csv")
# season = "Autumn"
# df <- df[df$season==season,]

a <- p$val[1]
b <- p$val[2]
d <- p$val[3]
e <- p$val[4]
C <- p$val[5]
del <- p$val[6]

a.se <- p$se[1]
b.se <- p$se[2]
d.se <- p$se[3]
e.se <- p$se[4]
C.se <- p$se[5]
del.se <- p$se[6]

set.seed(1001)
## simulate data with covariates x1, x2 and z.
.n    <- length(mass_init)
# .d    <- data.frame(x1 = rnorm(.n))
# .d$x2 <- sqrt(0.5)*.d$x1 + rnorm(.n, sd=sqrt(0.5))
# .d$z  <- as.numeric(.d$x1 + .d$x2 > 0)

#rates
growth_obs <- del^3*((df$len_final_QC/10)^3-(df$len_init_QC/10)^3)/df$expt_length
mass_init <- del^3*(df$len_init_QC/10)^3
Nth <- df$thread_count_QC/df$expt_length

.d <- data.frame(
  x1 = C * a * mass_init ^ d, #f is this coeff
  x2 = - C * b * mass_init ^ 1, # this coeff should be 1
  x3 = - C * Nth #h is this coeff
)

.d$y <- growth_obs


se.w1 <- sqrt((C.se/C)^2+(a.se/a)^2) #need to remember how to incorporate an exponent
se.w2 <- sqrt((C.se/C)^2+(b.se/b)^2)
se.w3 <- C.se
Sigma_error <- diag(c(se.w1,se.w2,se.w3))
dimnames(Sigma_error) <- list(c("w1","w2","w3"), c("w1","w2","w3"))
#   w1        w2   w3
# w1 0.5030099 0.0000000 0.00
# w2 0.0000000 0.4692473 0.00
# w3 0.0000000 0.0000000 0.01

.d$w1 <- .d$x1 + rnorm(.n, sd = sqrt(se.w1)) #Is this right? sd? vs. se? 
.d$w2 <- .d$x2 + rnorm(.n, sd = sqrt(se.w2))
.d$w3 <- .d$x3 + rnorm(.n, sd = sqrt(se.w3))


## generate outcome
## true regression parameters are c(2,1,1,-1)
# .d$y  <- 2.0 + 1.0*.d$x1 + 1.0*.d$x2 - 1.0*.d$z + rnorm(.n)

## generate error-prone covariates w1 and w2
# Sigma_error <- diag(c(0.20, 0.30))
# dimnames(Sigma_error) <- list(c("w1","w2"), c("w1","w2"))
# .d$w1 <- .d$x1 + rnorm(.n, sd = sqrt(Sigma_error["w1","w1"]))
# .d$w2 <- .d$x2 + rnorm(.n, sd = sqrt(Sigma_error["w2","w2"]))


# This is odd because there is a z
## fit EIV regression specifying known measurement error covariance matrix
.mod1 <- eivreg(y ~ w1 + w2 + w3, data = .d, Sigma_error = Sigma_error)
print(class(.mod1))
.tmp <- summary(.mod1)
print(class(.tmp))
print(.tmp)

#issue negatives in place of positives - need starting value? 

# Low variance version
# Call:
#   eivreg(formula = y ~ w1 + w2 + w3, data = .d, Sigma_error = Sigma_error)
# Error Covariance Matrix
# w1      w2    w3
# w1 0.07164 0.00000 0.000
# w2 0.00000 0.09724 0.000
# w3 0.00000 0.00000 0.001
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0091121 -0.0033861  0.0000926  0.0035039  0.0112590 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.022387   0.016552   1.352    0.184
# w1          -0.014587   0.013978  -1.044    0.303
# w2          -0.008938   0.009025  -0.990    0.328
# w3           0.001684   0.001143   1.474    0.149
# Number of observations used: 43 
# Latent residual standard deviation: 0.001906 
# Latent R-squared: 0.608, (df-adjusted: 0.5678)
# 
# EIV-Adjusted vs Unadjusted Coefficients:
#   Adjusted Unadjusted
# (Intercept)  0.022387   0.008773
# w1          -0.014587  -0.002506
# w2          -0.008938  -0.001077
# w3           0.001684   0.001020

## fit EIV regression specifying known reliabilities.  Note that
## point estimator is slightly different from .mod1 because
## the correction matrix S must be estimated when the reliability
## is known.
.lambda <- c(1,1,1) / (c(1,1,1) + diag(Sigma_error))
.mod2 <- eivreg(y ~ w1 + w2 + w3, data = .d, reliability = .lambda)
print(summary(.mod2))

# Call:
#   eivreg(formula = y ~ w1 + w2 + w3, data = .d, reliability = .lambda)
# 
# Reliability:
#   w1     w2     w3 
# 0.9331 0.9114 0.9990 
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0044774 -0.0018029 -0.0002702  0.0015590  0.0083741 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.0091612  0.0020298   4.513 5.74e-05 ***
#   w1          -0.0028628  0.0012344  -2.319  0.02571 *  
#   w2          -0.0013204  0.0008961  -1.474  0.14864    
# w3           0.0010392  0.0003403   3.054  0.00406 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Number of observations used: 43 
# Latent residual standard deviation: 0.002727 
# Latent R-squared: 0.1978, (df-adjusted: 0.1155)
# 
# EIV-Adjusted vs Unadjusted Coefficients:
#   Adjusted Unadjusted
# (Intercept)  0.009161   0.008773
# w1          -0.002863  -0.002506
# w2          -0.001320  -0.001077
# w3           0.001039   0.001020

