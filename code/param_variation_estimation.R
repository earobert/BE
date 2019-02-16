# Need to change these this practice dataset
# Found issue that intake and cost are very much correlated. 
# What does this mean? Should I use partial least squares regression???
# http://blog.minitab.com/blog/statistics-and-quality-data-analysis/giving-thanks-for-the-regression-menu-v2

library(eivtools)
library(lavaan)

# calculate error in b relative to wet weight ####
# b [J/day/gWW]

W_inf <- c(0.72,0.06) #mean + SD
R <- c(0.037, 0.017) #mean + SE
ratioWWpDW <- c(3.98,0.07) #mean + SE
R_gDW <- .429 # mean weight of tissue
J_per_cal <- 4.184
cal_per_mlO <- 4.75
hr_per_day <- 24

#high estimate
R_upper <- R[1]+R[2]
ratioWWpDW_lower <- ratioWWpDW[1]-ratioWWpDW[2]

#low estimate
R_lower <- R[1]-R[2]
ratioWWpDW_upper <- ratioWWpDW[1]+ratioWWpDW[2]

#b = R / (R_gDW * ratioWWpDW) * J

(b <- (R[1] * J_per_cal *cal_per_mlO *hr_per_day) / (R_gDW[1] * ratioWWpDW[1]))
(b_upper <- (R_upper * J_per_cal *cal_per_mlO *hr_per_day) / (R_gDW[1] * ratioWWpDW_lower))
(b_lower <- (R_lower * J_per_cal *cal_per_mlO *hr_per_day) / (R_gDW[1] * ratioWWpDW_upper))

#Range 5.4 to 15.3 This is in part because of the dry to wet weight ratio. 
b_upper-b

5.0/10.3
#49% error

p.practice <- data.frame( 
  names = c("a","b","d","e","C","del"),
  val = c(14,10.3,0.65,1,0.000182,.333),
  se = c(7,5.0,.01,0,.00001,.006),
  season = rep("Autumn", times = length(names))
)
p.practice


# calculate error in b, using dry weight only ####
# b [J/day/gDW]

# W_inf <- c(0.72,0.06) #mean + SD
# R <- c(0.037, 0.017) #mean + SE
# R_gDW <- .429 # mean weight of tissue
# J_per_cal <- 4.184
# cal_per_mlO <- 4.75
# hr_per_day <- 24
# 
# #high estimate
# R_upper <- R[1]+R[2]
# 
# #low estimate
# R_lower <- R[1]-R[2]
# 
# #b = R / (R_gDW * ratioWWpDW) * J
# 
# (b <- (R[1] * J_per_cal *cal_per_mlO *hr_per_day) / (R_gDW[1]))
# (b_upper <- (R_upper * J_per_cal *cal_per_mlO *hr_per_day) / (R_gDW[1]))
# (b_lower <- (R_lower * J_per_cal *cal_per_mlO *hr_per_day) / (R_gDW[1]))
# 
# #Range 5.4 to 15.3 This is in part because of the dry to wet weight ratio. 
# b_upper-b
# b-b_lower
# 
# 18.9/41.1
# # 46% error
# 
# p.practice <- data.frame( 
#   names = c("a","b","d","e","C","del"),
#   val = c(14,41.1,0.69,1,0.182,.333),
#   se = c(7,18.9,.01,.0001,.01,.006),
#   season = rep("Autumn", times = length(names))
# )
# p.practice

p <- p.practice
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

# Import dataset

# lm 1
season="Autumn"
setwd("~/BE/BE/Datasets")
df <- read.csv(file="Spring_Fall.csv", stringsAsFactors = FALSE)
df$treatment <- as.factor(df$treatment)
df$season <- as.factor(df$season)
df$treatment <- ordered(df$treatment, levels = c("never", "weekly","daily"))
df <- df[df$season==season,]
df <- df[(!is.na(df$len_final_QC))&(!is.na(df$len_init_QC)),]
df <- df[df$len_init_QC<30,]
df$len_final_QC
plot(df$len_init_QC,df$len_final_QC)
#df.subset <- df[df$treatment=="never"|df$treatment=="daily",]
#df <- df.subset

growth_obs <- ((del*(df$len_final_QC/10))^3-(del*(df$len_init_QC/10))^3)/df$expt_length
growth_obs_upper <- ((del*(df$len_final_QC/10))^3-(del*(df$len_init_QC/10))^3)/df$expt_length

mass_init <- (del*(df$len_init_QC/10))^3
mass_init_upper <- ((del+del.se)*(df$len_init_QC/10))^3
mass_init_lower <- ((del-del.se)*(df$len_init_QC/10))^3
mass_final <- (del*(df$len_final_QC/10))^3
mass_final_upper <- ((del+del.se)*(df$len_final_QC/10))^3
mass_final_lower <- ((del-del.se)*(df$len_final_QC/10))^3

growth_obs <- mass_final - mass_init
growth_obs_upper <- mass_final_upper - mass_init_lower
growth_obs_lower <- mass_final_lower - mass_init_upper

mass_final_real <- df$total_wt_dry*3.98
Nth <- df$thread_count_QC

length(growth_obs)
length(mass_init)
length(Nth)

C #g/J
C_upper <- C + C.se
C_lower <- C - C.se
#C_mg <- C/1000

#growth_tissue <- growth_obs*1000


# just using normal C [g/J]
# non_byss (negative)
non_byss <- b*C*(mass_init)^e
non_byss_upper <- b_upper*C_upper*(mass_init_upper)
non_byss_lower <- b_lower*C_lower*(mass_init_lower)

plot(mass_init,non_byss,
     ylim = c(.0003,.0028),
     col = df$treatment)
points(mass_init,non_byss_lower, col = df$treatment)
points(mass_init,non_byss_upper, col = df$treatment)

# byss (negative)
byss <- C*Nth
byss_upper <- C_upper*Nth
byss_lower <- C_lower*Nth

plot(byss,byss, col = df$treatment)
points(byss,byss_upper, col = df$treatment)
points(byss,byss_lower, col = df$treatment)

# intake (positive)
intake <- C*(mass_init)^d 
plot(intake,intake, col = df$treatment)
intake_upper <- C_upper*(mass_init_upper)^(d)
intake_lower <- C_lower*(mass_init_lower)^(d)
points(intake,intake_upper, col = df$treatment)
points(intake,intake_lower, col = df$treatment)
intake_upper <- C_upper*(mass_init_upper)^(d+d.se)
intake_lower <- C_lower*(mass_init_lower)^(d-d.se)
points(intake,intake_upper, col = 'blue')
points(intake,intake_lower, col = 'blue')


# What about the relationship between mass and intake-non_byss? Which plays out?
a<- 14
SFBG <- a*intake*1000 - non_byss*1000
SFBG_lower<-(a*intake_lower - non_byss_upper)*1000
SFBG_upper<-(a*intake_upper - non_byss_lower)*1000

plot(mass_init,SFBG,
     ylim = c(-.6,4), col = df$treatment)
points(mass_init,SFBG_upper)
points(mass_init,SFBG_lower)
# looks like neither dominates, the line seems pretty flat, 
# with increasing varience with increasing mass. The slope of this line also depends on the intake coeff.
# the low end estimates negative growth at most mass values. 

plot(mass_init,growth_obs, col = df$treatment)
legend(x = "topleft", legend = unique(df$treatment), 
       col = as.numeric(unique(df$treatment)), pch = 1)

# Just checking how this works, visually. 
# Assuming the energy efficiency is not known. 
SFG <- 20*intake*1000 - non_byss*1000 - .01*byss*1000
plot(mass_init,SFG, col = df$treatment)

#checking how this would be different if using mg rather than grams 
# because the exponent would have a more linear effect
# input <- C*(mass_init*1000)^d 
# plot(input,input, col = df$treatment)
# input_upper <- C_upper*(mass_init_upper*1000)^(d)
# input_lower <- C_lower*(mass_init_lower*1000)^(d)
# points(input,input_upper, col = df$treatment)
# points(input,input_lower, col = df$treatment)
# input_upper <- C_upper*(mass_init_upper*1000)^(d+d.se)
# input_lower <- C_lower*(mass_init_lower*1000)^(d-d.se)
# points(input,input_upper, col = 'blue')
# points(input,input_lower, col = 'blue')
# turns out the exponent has a larger effect... but doesn't really apply
# because if you check out Jones et al., 1997, they use g of tissue and 0.7

cov(x=data.frame(input, byss, non_byss, growth_obs))
cor(x=data.frame(input, byss, non_byss, growth_obs))
# input       byss    non_byss  growth_obs
# input       1.00000000 -0.3553659  0.99920619  0.04640339
# byss       -0.35536592  1.0000000 -0.35610221 -0.31467507
# non_byss    0.99920619 -0.3561022  1.00000000  0.04037684
# growth_obs  0.04640339 -0.3146751  0.04037684  1.00000000
# Note that non_byss and input are 99.9% correlated
# Note these are only 4% correlated with growth observed
# Note that byssus is 30% correlated with growth observed
mod1 <- lm(growth_obs ~ input + byss + non_byss)
summary(mod1)
vcov(mod1)
# > summary(mod1)
# 
# Call:
#   lm(formula = growth_obs ~ input + byss + non_byss + 0)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.15491 -0.04195 -0.01388  0.02952  0.24107 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
#   input    4128.6916  1381.9776   2.988  0.00484 **
#   byss       -0.9062     0.4146  -2.186  0.03491 * 
#   non_byss -336.5907   146.7327  -2.294  0.02727 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.08282 on 39 degrees of freedom
# Multiple R-squared:  0.6814,	Adjusted R-squared:  0.6569 
# F-statistic: 27.81 on 3 and 39 DF,  p-value: 8.705e-10
# 
# > vcov(mod1)
# input         byss      non_byss
# input    1909862.0723 -346.2114457 -201561.08304
# byss        -346.2114    0.1719028      32.61977
# non_byss -201561.0830   32.6197673   21530.47162

lavaan(mod1)
# The Holzinger and Swineford (1939) example
HS.model <- ' visual  =~ x1 + x2 + x3,
              textual =~ x4 + x5 + x6,
              speed =~ x7 + x8 + x9'
head(HolzingerSwineford1939)
fit <- lavaan(HS.model, data=HolzingerSwineford1939,
              auto.var=TRUE, auto.fix.first=TRUE,
              auto.cov.lv.x=TRUE)
summary(fit, fit.measures=TRUE)

cost_model <- '
x1 =~ C*mass_init^d
x2 =~ C*Nth
x3 =~ 1*b*C
y1 =~ x1-x2-x3
# intercept
y1 ~ 0
# covariances
x1 ~~ x2
x1 ~~ x3
x2 ~~ x3
'


data(Demo.growth)
describe(Demo.growth)
pairs.panels(Demo.growth)

