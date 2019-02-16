# This is not just linear model code (current working code.... more or less matches model 5), 
# but also I looked at two different estimations of growth, 
# but the one that is now here is 0.35 and matches the data, 
# not 0.304, which matches the Penn Cove mussels in Summer 2015 or so. 
# I use my original energy density estimate of 0.002, and then use an estimate. 



# linear model now using real data
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
season="Autumn"
setwd("~/BE/BE/Datasets")
df <- read.csv(file="Spring_Fall.csv", stringsAsFactors = FALSE)
df$treatment <- as.factor(df$treatment)
df$season <- as.factor(df$season)
df$treatment <- ordered(df$treatment, levels = c("never", "weekly","daily"))
df <- df[df$season==season,]
df.never <- df[df$treatment=="never",]

# all data
df <- df[!is.na(df$len_init_QC)&!is.na(df$len_final_QC),]
len_init <- df[df$season==season,]$len_init_QC/10 #converted from mm to cm
len_final <- df[df$season==season,]$len_final_QC/10# converted from mm to cm
thread_num <- df[df$season==season,]$thread_count_QC
mass_ww_init <- (len_init*.35)^3
mass_ww_final <- (len_final*.35)^3

plot(df$total_wt_dry,df$total_wt_wet)
relWW <- lm(df$total_wt_wet~df$total_wt_dry+0)
summary(relWW)
# 
# Call:
#   lm(formula = df$total_wt_wet ~ df$total_wt_dry + 0)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.16095 -0.02694  0.03677  0.06975  0.21294 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# df$total_wt_dry  3.98420    0.06567   60.67   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.08795 on 42 degrees of freedom
# Multiple R-squared:  0.9887,	Adjusted R-squared:  0.9884 
# F-statistic:  3681 on 1 and 42 DF,  p-value: < 2.2e-16


head(df)

time <- as.Date(df$date_final, "%d-%b-%y")-as.Date(df$date_init,"%d-%b-%y")
en_density_g_p_J <- .002*as.numeric(time) #energy_density

par(mfrow = c(1,1))
#Mass is actual measured mass???####
mass_ww_final_meas <- df[df$season==season,]$total_wt_wet
plot(mass_ww_final,mass_ww_final_meas,
     xlim = c(0.3,1.5),
     ylim = c(0.3,1.5))
x <- seq(from = .3, to=1.5, by=.1)
lines(x,x, lty = 2)

plot(len_final, mass_ww_final_meas)
x <- seq(from = 2, to=3.2, by=.1)
y <- (x*.34)^3
lines(x,y, lty = 2)
y <- (x*.304)^3
lines(x,y, lty = 1)

# import data here and make sure this is the same relationshipo
#plot(len_initial, )
x <- seq(from = 2, to=3.2, by=.1)
y <- (x*.34)^3
lines(x,y, lty = 2)
y <- (x*.304)^3
lines(x,y, lty = 1)

#I'm trying to estimate del with SE, but this is tricky because it is not linear
# y = (len * del)^3
# ln(y) = 3ln(len*del)
# ln(y) = 3ln(len)+3ln(del)
y = mx + b

nonlin_mod=nls(y~a*exp(b*x),start=list(a=13,b=0.1))
nonlin_mod=nls(mass_ww_final_meas~(del*len_final)^3,start=list(del=.35))
summary(nonlin_mod)
# Formula: mass_ww_final_meas ~ (del * len_final)^3
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# del 0.340354   0.002743   124.1   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1282 on 42 degrees of freedom
# 
# Number of iterations to convergence: 3 
# Achieved convergence tolerance: 1.3e-09

#====#
mass_ww <- mass_ww_init
growth_shell <- len_final-len_init
growth_tissue <- mass_ww_final-mass_ww_init


B2 <- en_density_g_p_J*a_prime*(mass_ww^d) #good, this is already in mass_ww... hope this is g not mg
B3 <- -en_density_g_p_J*b*(mass_ww^e) #*(1-baseline_byssus_multiplier)
# byssus_baseline <- b*(mass_ww^e)*baseline_byssus_multiplier
#byssus_induced <- thread_num*cost_per_thread
B1 <- -en_density_g_p_J*thread_num
reproduction <- 0

dat_lm <- data.frame(
  growth_tissue = growth_tissue,
  B1 = B1,
  B2 = B2,
  B3 = B3
)

y <- lm(growth_tissue ~ B1 + B2 + offset(B3) + 0, data = dat_lm)
length(B1)
length(growth_tissue)
sum_y <- summary(y)
pred_growth_2 <-sum_y$coefficients[1]*B1+sum_y$coefficients[2]*B2+B3
par(mfrow = c(3,1), mar = c(4, 4, 0, 0) + 0.1)
plot(pred_growth_2~growth_tissue,
      ylim = c(-0.1,.5),
      xlim = c(-0.1,.5),
     col=df$treatment
)
x <- seq(from = -.4, to=.4, by=.1)
lines(x,x, lty = 2)

plot(pred_growth_2~thread_num,col=df$treatment
)

plot(pred_growth_2~mass_ww, col=df$treatment)



y1 <-y
# # Hmm... here if we make the intercept 0, then B1 is plus or minus 0.02... yikes. 
# # I'm not sure if this means that we should estimate an additional variable. 
# 
# # > old sum_y
# # 
# # Call:
# #   lm(formula = growth_tissue ~ B1 + B2 + offset(B3) + 0, data = dat_lm)
# # 
# # Residuals:
# #   Min       1Q   Median       3Q      Max 
# # -0.12043 -0.04059 -0.01417  0.03571  0.18014 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # B1  0.02862    0.02506   1.142     0.26    
# # B2  1.40318    0.16635   8.435 1.69e-10 ***
# #   ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# # 
# # Residual standard error: 0.06554 on 41 degrees of freedom
# # Multiple R-squared:  0.6364,	Adjusted R-squared:  0.6186 
# # F-statistic: 35.87 on 2 and 41 DF,  p-value: 9.862e-10
# Call:
#   lm(formula = growth_tissue ~ B1 + B2 + offset(B3) + 0, data = dat_lm)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.18378 -0.06195 -0.02163  0.05450  0.27491 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# B1  0.04368    0.03825   1.142     0.26    
# B2  1.61323    0.19125   8.435 1.69e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1 on 41 degrees of freedom
# Multiple R-squared:  0.6364,	Adjusted R-squared:  0.6186 
# F-statistic: 35.87 on 2 and 41 DF,  p-value: 9.862e-10

# This is to test later on...
# This is to test later on...

# practice_data_params<-c(.19,.061)
# cost_per_thread <-practice_data_params[2]
# food_scalar <- practice_data_params[1]
# 
# pred_growth <- cost_per_thread*B1+food_scalar*B2+B3
# compare <- data.frame(growth_tissue,
#                       pred_growth)
# plot(growth_tissue, pred_growth,
#      ylim = c(0,.6),
#      xlim = c(0, .6)
# )
# x <- seq(from = -.4, to=.4, by=.1)
# lines(x,x, lty = 2)



# practice_data_params<-c(1.5,.06)
# cost_per_thread <-practice_data_params[2]
# food_scalar <- practice_data_params[1]
# 
# pred_growth <- cost_per_thread*B1+food_scalar*B2+B3
# compare <- data.frame(growth_tissue,
# pred_growth)
# plot(growth_tissue, pred_growth,
#      ylim = c(0,.6),
#      xlim = c(0, .6)
# )
# x <- seq(from = -.4, to=.4, by=.1)
# lines(x,x, lty = 2)

B1 <- -thread_num
B2 <- a_prime*(mass_ww^d) #good, this is already in mass_ww... hope this is g not mg
B3 <- -b*(mass_ww^e) #*(1-baseline_byssus_multiplier)
growth_tissue <- mass_ww_final-mass_ww_init

dat_lm <- data.frame(
  growth_tissue = growth_tissue,
  B1 = B1,
  B2 = B2,
  B3 = B3
)

y <- lm(growth_tissue ~ B1 + B2 + B3 + 0, data = dat_lm)
sum_y <- summary(y)
growth_efficiency_time <- sum_y$coefficients[3]
growth_efficiency <-sum_y$coefficients[3] / as.numeric(time[1])
growth_efficiency_SE <-growth_efficiency*sqrt((sum_y$coefficients[3,2]/sum_y$coefficients[3])^2+(.1/29)^2) # I think that time doesn't change the SE though it maybe should... ?
#growth_efficiency_SE <- sum_y$coefficients[3,2] 
cB1 <- sum_y$coefficients[1]
cB2 <- sum_y$coefficients[2]
cB1_SE <- sum_y$coefficients[1,2]
cB2_SE <- sum_y$coefficients[2,2]
cost_of_thread <- sum_y$coefficients[1]/(growth_efficiency_time)
cost_of_thread_SE <- cost_of_thread * sqrt((cB1_SE/cB1)^2+(growth_efficiency_SE/growth_efficiency_time)^2)
food_scalar <- sum_y$coefficients[2]/growth_efficiency_time
food_scalar_SE <- food_scalar * sqrt((cB2_SE/cB2)^2+(growth_efficiency_SE/growth_efficiency_time)^2)

# > sum_y
# 
# Call:
#   lm(formula = growth_tissue ~ B1 + B2 + B3 + 0, data = dat_lm)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.18108 -0.04762 -0.01397  0.03292  0.27867 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
#   B1 1.879e-04  8.442e-05   2.226  0.03168 * 
#   B2 1.058e-02  3.164e-03   3.343  0.00181 **
#   B3 1.846e-02  7.037e-03   2.623  0.01228 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.09498 on 40 degrees of freedom
# Multiple R-squared:  0.6835,	Adjusted R-squared:  0.6597 
# F-statistic: 28.79 on 3 and 40 DF,  p-value: 4.382e-10


stats <- data.frame(
  names = c("growth_efficiency","cost_of_thread","food_scalar"),
  means = c(growth_efficiency,cost_of_thread,food_scalar),
  SE = c(growth_efficiency_SE,cost_of_thread_SE,food_scalar_SE)
)

stats
# # > stats old
# # names       means          SE
# # 1 growth_efficiency 0.018457318 0.007037088
# # 2    cost_of_thread 0.006672387 0.003931023
# # 3       food_scalar 0.498433799 0.241544571
# > stats new
# names      means          SE
# 1 growth_efficiency 0.01845732 0.007037088
# 2    cost_of_thread 0.01018274 0.005999139
# 3       food_scalar 0.57304668 0.277702505

par(mfrow = c(3,1))

pred_growth_2 <-sum_y$coefficients[1]*B1+sum_y$coefficients[2]*B2+sum_y$coefficients[3]*B3
max(pred_growth_2)
max(growth_tissue)
plot(pred_growth_2~growth_tissue,
      ylim = c(-.01,.5),
      xlim = c(-.01,.5),
     col=df$treatment
     )
x <- seq(from = -.4, to=.4, by=.1)
lines(x,x, lty = 2)

plot(pred_growth_2~thread_num,col=df$treatment
)

plot(pred_growth_2~mass_ww, col=df$treatment)

par(mfrow = c(1,1))
plot(growth_tissue~mass_ww,col=df$treatment)

#Why is there so much variability???

y2 <-y
anova(y1, y2)
aov(y1)
aov(y2)

# Wow, so there is a significantly better estimation 
# when lm is allowed to estimate the growth efficiency.
# I wonder if I should allow there to also be an error estimate??? 
# But also is it weird that there are more DFs for model 1 when this model just has the offset? 
# # > anova(y1, y2)
# # Analysis of Variance Table
# # 
# # Model 1: growth_tissue ~ B1 + B2 + offset(B3) + 0
# # Model 2: growth_tissue ~ B1 + B2 + B3 + 0
# # Res.Df     RSS Df Sum of Sq      F  Pr(>F)  
# # 1     41 0.17611                              
# # 2     40 0.15493  1  0.021184 5.4693 0.02444 *
# #   ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# > anova(y1, y2)
# Analysis of Variance Table
# 
# Model 1: growth_tissue ~ B1 + B2 + offset(B3) + 0
# Model 2: growth_tissue ~ B1 + B2 + B3 + 0
# Res.Df     RSS Df Sum of Sq      F  Pr(>F)  
# 1     41 0.41017                              
# 2     40 0.36083  1  0.049338 5.4693 0.02444 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

