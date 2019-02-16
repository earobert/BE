set.seed(194812)
n <- 100
x <- rnorm(n)
residual_sd <- exp(x)
y <- 2*x + residual_sd*rnorm(n)

plot(x,y)

mod <- lm(y~x)
summary(mod)

library(sandwich)
vcovHC(mod, type = "HC")

sandwich_se <- diag(vcovHC(mod, type = "HC"))^.5
sandwich_se

coef(mod)-1.96*sandwich_se
coef(mod)+1.96*sandwich_se
