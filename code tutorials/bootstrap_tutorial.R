#bootstrapping tutorial
#https://www.statmethods.net/advstats/bootstrapping.html


# Bootstrap 95% CI for R-Squared
library(boot)
# function to obtain R-Squared from the data 
rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
} 
# bootstrapping with 1000 replications 
results <- boot(data=mtcars, statistic=rsq, 
                R=1000, formula=mpg~wt+disp)

model <- lm(mpg~wt+disp,data = mtcars)
summary(model)
summary(model)$r.square
plot(mpg~wt+disp,data = mtcars)
Call:
  lm(formula = mpg ~ wt + disp, data = mtcars)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.4087 -2.3243 -0.7683  1.7721  6.3484 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 34.96055    2.16454  16.151 4.91e-16 ***
#   wt          -3.35082    1.16413  -2.878  0.00743 ** 
#   disp        -0.01773    0.00919  -1.929  0.06362 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.917 on 29 degrees of freedom
# Multiple R-squared:  0.7809,	Adjusted R-squared:  0.7658 
# F-statistic: 51.69 on 2 and 29 DF,  p-value: 2.744e-10


# view results
results 
plot(results)

# get 95% confidence interval 
boot.ci(results, type="bca")
# CALL : 
#   boot.ci(boot.out = results, type = "bca")
# 
# Intervals : 
#   Level       BCa          
# 95%   ( 0.6413,  0.8522 )  
# Calculations and Intervals on Original Scale
#Some BCa intervals may be unstable
#-------------

# Bootstrap 95% CI for regression coefficients 
library(boot)
# function to obtain regression weights 
bs <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- lm(formula, data=d)
  return(coef(fit)) 
} 
# bootstrapping with 1000 replications 
results <- boot(data=mtcars, statistic=bs, 
                R=1000, formula=mpg~wt+disp)

# view results
results
plot(results, index=1) # intercept 
plot(results, index=2) # wt 
plot(results, index=3) # disp 

# get 95% confidence intervals 
boot.ci(results, type="bca", index=1) # intercept 
boot.ci(results, type="bca", index=2) # wt 
boot.ci(results, type="bca", index=3) # disp


