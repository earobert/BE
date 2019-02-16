# linear model
# B1 = -thread_num
# B2 = a'*mass^d
# B3 = -b*mass^e
# y = B1 * h + B2 * food + B3

rm(list=ls())

# Model parameters ####
# conversion factors
conversion_gWW_per_gDW <- 3.918 # From Summer 2015 collection, shape coeff estimation
shape_coeff <- .304 # From Summer 2015 collection, shape coeff estimation 
mass_ww_per_dw <- 3.9 #coverts from mass_DW (dry weight) to mass_WW (wet weight)

# exponents ####
d <- 0.67 # intake; the model is very sensitive to this. (hope this is an exponent for g WW not mg WW)
e <- 1 # the model is very sensitive to this (hope this is an exponent for g WW not mg WW)

# calculation of b
respiration_reference_J_per_day <- 0.07*4.75*4.184*24 # Units: J / (day) from 0.07mlO2/hr, Fly and Hilbish 2013
respiration_reference_gDW <- 0.24 # Fly and Hilbish 2013
respiration_reference_gWW <- respiration_reference_gDW * 3.9 # Using wet:dry conversion
b_permass <- respiration_reference_J_per_day / (respiration_reference_gWW)^e # Note that e is just 1
b_J_per_g_WW_per_day <- b_permass # Units: J / (gWW * day)

# Temp response, other  ####
Tmult_cost <- 1
Tmult_int <- 1
reprod_multiplier <- 1 # reproduction multiplier, brachi was 1, but if k is about half and half then this should be more like .5
opt_size_reducer <- 1 # lets reproduction equal surplus early

# other parameters ####
en_density_g_p_J <- .002 #energy_density
size_maturity <- 10 # mg, size at maturity, 3000mg from tross del_M calc worksheet, somatic<full tissue weight

# calculated scalers
b <- b_J_per_g_WW_per_day*Tmult_cost #Wait the byssus multiplier is earlier here than I thought

Wopt_measured_gDW <- 0.8 #gDW from sample collection
Wopt_measured_gWW <- Wopt_measured_gDW * conversion_gWW_per_gDW
a_fromWopt <- (b*e)/((Wopt_measured_gWW)^(d-e)*d) # backwards calculation of a
a_J_per_day_per_food_scalar <- a_fromWopt
a_prime <- a_J_per_day_per_food_scalar*Tmult_int

# import data ####
setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="practice_data.csv", stringsAsFactors = FALSE)
head(df)
df$treatment <- as.factor(df$treat)
#df$season <- as.factor(df$season)
#df$treatment <- ordered(df$treatment, levels = c("never", "weekly","daily"))
#df <- df[df$season==season,]

# all data
#df <- df[!is.na(df$len_init_QC)&!is.na(df$len_final_QC),]
#len_init <- df[df$season==season,]$len_init_QC/10 #converted from mm to cm
#len_final <- df[df$season==season,]$len_final_QC/10# converted from mm to cm
thread_num <- df$thread_num
mass_ww_init <- df$mass_ww_init
#mass_ww_final <- (len_final*.304)^3
mass_ww <- mass_ww_init
#growth_shell <- len_final-len_init
growth_tissue <- df$model.growth

practice_data_params<-c(3,.1)
cost_per_thread <-practice_data_params[2]
food_scalar <- practice_data_params[1]

B2 <- en_density_g_p_J*a_prime*(mass_ww^d) #good, this is already in mass_ww... hope this is g not mg
B3 <- -en_density_g_p_J*b*(mass_ww^e) #*(1-baseline_byssus_multiplier)
# byssus_baseline <- b*(mass_ww^e)*baseline_byssus_multiplier
#byssus_induced <- thread_num*cost_per_thread
B1 <- -en_density_g_p_J*thread_num
reproduction <- 0
pred_growth <- cost_per_thread*B1+food_scalar*B2+B3
compare <- data.frame(growth_tissue,
                      pred_growth)
plot(growth_tissue, pred_growth,
     ylim = c(.02,.4),
     xlim = c(0.02, .4)
     )
x <- seq(from = -.4, to=.4, by=.1)
lines(x,x, lty = 2)
#This shows that if I use the costs used for the practice data, I get a 1:1 relationship!  

dat_lm <- data.frame(
  growth_tissue = growth_tissue,
  B1 = B1,
  B2 = B2,
  B3 = B3
)

y <- lm(growth_tissue ~ B1 + B2 + offset(B3), data = dat_lm)
length(B1)
length(growth_tissue)
sum_y <- summary(y)
pred_growth_2 <-sum_y$coefficients[2]*B1+sum_y$coefficients[3]*B2+B3
plot(growth_tissue~pred_growth_2)
x <- seq(from = -.4, to=.4, by=.1)
lines(x,x, lty = 2)

# this is really great - It's giving .1 for the cost and 3 for the food availability, which are my practive values! Yay, this is working! :)
# I'm also noticing the standard error is really small... kind of not in a believable way.
# > sum_y
# 
# Call:
#   lm(formula = growth_tissue ~ B1 + B2 + offset(B3), data = dat_lm)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -5.102e-16 -1.571e-16 -1.634e-17  1.782e-16  5.369e-16 
# 
# Coefficients:
#               Estimate Std. Error   t value Pr(>|t|)    
#   (Intercept) 3.312e-16  2.802e-16 1.182e+00    0.249    
#   B1          1.000e-01  1.521e-16 6.575e+14   <2e-16 ***
#   B2          3.000e+00  2.846e-15 1.054e+15   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.796e-16 on 24 degrees of freedom
# Multiple R-squared:      1,	Adjusted R-squared:      1 
# F-statistic: 6.371e+29 on 2 and 24 DF,  p-value: < 2.2e-16


# This is probably more correct. Setting the intercept to 0.
# Both values are significant. 
# Why isn't this working with the real data?

y <- lm(growth_tissue ~ B1 + B2 + offset(B3)+0, data = dat_lm)
length(B1)
length(growth_tissue)
sum_y <- summary(y)
pred_growth_2 <-sum_y$coefficients[1]*B1+sum_y$coefficients[2]*B2+B3
plot(growth_tissue~pred_growth_2)
x <- seq(from = -.4, to=.4, by=.1)
lines(x,x, lty = 2)

# > sum_y
# 
# Call:
#   lm(formula = growth_tissue ~ B1 + B2 + offset(B3) + 0, data = dat_lm)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -5.432e-16 -1.461e-16  3.926e-17  1.783e-16  5.751e-16 
# 
# Coefficients:
#   Estimate Std. Error   t value Pr(>|t|)    
#   B1 1.000e-01  1.409e-16 7.096e+14   <2e-16 ***
#   B2 3.000e+00  9.443e-16 3.177e+15   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.751e-16 on 25 degrees of freedom
# Multiple R-squared:      1,	Adjusted R-squared:      1 
# F-statistic: 7.035e+30 on 2 and 25 DF,  p-value: < 2.2e-16


B1 <- -thread_num
B2 <- a_prime*(mass_ww^d) #good, this is already in mass_ww... hope this is g not mg
B3 <- -b*(mass_ww^e) #*(1-baseline_byssus_multiplier)
growth_tissue <- df$model.growth

dat_lm <- data.frame(
  growth_tissue = growth_tissue,
  B1 = B1,
  B2 = B2,
  B3 = B3
)

y <- lm(growth_tissue ~ B1 + B2 + B3+0, data = dat_lm)
sum_y <- summary(y)
growth_efficiency <-sum_y$coefficients[3] #.002
cost_of_thread <- sum_y$coefficients[1]/growth_efficiency #0.1
food_scalar <- sum_y$coefficients[2]/growth_efficiency #3
pred_growth_2 <-sum_y$coefficients[1]*B1+sum_y$coefficients[2]*B2+sum_y$coefficients[3]*B3
plot(growth_tissue~pred_growth_2)
x <- seq(from = -.4, to=.4, by=.1)
lines(x,x, lty = 2)

# > sum_y
# 
# Call:
#   lm(formula = growth_tissue ~ B1 + B2 + B3 + 0, data = dat_lm)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
#   -0.118658 -0.031204 -0.009152  0.021569  0.182601 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
#   B1 1.232e-04  5.531e-05   2.226  0.03168 * 
#   B2 9.200e-03  2.752e-03   3.343  0.00181 **
#   B3 1.846e-02  7.037e-03   2.623  0.01228 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.06224 on 40 degrees of freedom
# Multiple R-squared:  0.6835,	Adjusted R-squared:  0.6597 
# F-statistic: 28.79 on 3 and 40 DF,  p-value: 4.382e-10

