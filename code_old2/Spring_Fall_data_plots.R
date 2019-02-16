# These are going to be bar graphs + stats of:
# 1. the number of threads produced
# currently stuck on trying to get the graph to plot
# 2. shell growth
# 3. estimate tissue growth, and plot

rm(list=ls())

#------------------------#
# Thread production graph####
#------------------------#

setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is
season <- "Autumn"
df.Aut <- df[df$season==season,]
season <- "Spring"
df.Spr <- df[df$season==season,]


setwd("~/BE/BE_2018_12_3/figs")
png(filename = "thread_boxplots.png", width = 2*480, height = 480, units = "px", pointsize = 18)
par(mfrow = c(1,2), oma = c(0,0,0,0), mar = c(3, 4, 2, 0)+.1)

df <- df.Aut
length(df$treatment)
length(df$thread_count_QC)
df$treatment <- ordered(df$treatment, levels = c("never", "weekly","daily"))
boxplot(df$thread_count_QC~df$treatment,
        ylim = c(0,600),
        ylab = "Thread production (#)")

df <- df.Spr
length(df$treatment)
length(df$thread_count_QC)
df$treatment <- ordered(df$treatment, levels = c("never", "weekly","daily"))
boxplot(df$thread_count_QC~df$treatment, ylim = c(0,600), yaxt = "n")

dev.off()

# 2-way ANOVA
setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is
df$treatment <- as.factor(df$treatment)
df$treatment <- ordered(df$treatment, levels = c("never", "weekly","daily"))
shell_growth <- (df$len_final_QC-df$len_init_QC)/df$expt_length
model <- lm(df$thread_count_QC~df$treatment*df$season)
anova(model)

# Analysis of Variance Table
# 
# Response: df$thread_count_QC
# Df  Sum Sq Mean Sq  F value    Pr(>F)    
# df$treatment            2 1592171  796085 150.7573 < 2.2e-16 ***
#   df$season               1  357130  357130  67.6309 3.192e-12 ***
#   df$treatment:df$season  2   51414   25707   4.8682   0.01015 *  
#   Residuals              79  417166    5281                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#--------------------------#
#### Shell growth graph ####
#--------------------------#

setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is
season <- "Autumn"
df.Aut <- df[df$season==season,]
season <- "Spring"
df.Spr <- df[df$season==season,]


setwd("~/BE/BE_2018_12_3/figs")
png(filename = "shell_growth_boxplots.png", width = 2*480, height = 480, units = "px", pointsize = 18)
par(mfrow = c(1,2), oma = c(0,0,0,0), mar = c(3, 4, 2, 0)+.1)

df <- df.Aut
shell_growth <- (df$len_final_QC-df$len_init_QC)/df$expt_length
df$treatment <- ordered(df$treatment, levels = c("never", "weekly","daily"))
boxplot(shell_growth~df$treatment,
        ylim = c(0,0.15),
        ylab = "Shell growth (mm/day)")

df <- df.Spr
shell_growth <- (df$len_final_QC-df$len_init_QC)/df$expt_length
df$treatment <- ordered(df$treatment, levels = c("never", "weekly","daily"))
boxplot(shell_growth~df$treatment, 
        ylim = c(0,.15), 
        yaxt = "n"
        )

dev.off()

# 2-way ANOVA
setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is
df$treatment <- ordered(df$treatment, levels = c("never", "weekly","daily"))
shell_growth <- (df$len_final_QC-df$len_init_QC)/df$expt_length
model <- lm(shell_growth~df$treatment*df$season)
anova(model)

# Response: shell_growth
# Df   Sum Sq   Mean Sq F value   Pr(>F)   
# df$treatment            2 0.014911 0.0074555  7.5027 0.001079 **
#   df$season               1 0.006923 0.0069226  6.9664 0.010122 * 
#   df$treatment:df$season  2 0.002548 0.0012742  1.2823 0.283494   
# Residuals              74 0.073535 0.0009937                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#--------------------------------------#
#### Calculate conversions DW to WW ####
#--------------------------------------#
rm(list=ls())

# Autumn only, no data in Spring
season <- "Autumn"

setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is

head(df,1)
df <- df[df$season==season,]
#ratio<- df$total_wt_wet/df$total_wt_dry
#conversion_gWW_per_gDW <- mean(ratio, na.rm = TRUE)


# Model parameters ####
# conversion factors
conversion_gWW_per_gDW <- 3.918 # From Summer 2015 collection, shape coeff estimation
#shape_coeff <- .304 # From Summer 2015 collection, shape coeff estimation 
#mass_ww_per_dw <- 3.9 #coverts from mass_DW (dry weight) to mass_WW (wet weight)

# import data ####

df$treatment <- as.factor(df$treatment)
df$season <- as.factor(df$season)
df$treatment <- ordered(df$treatment, levels = c("never", "weekly","daily"))
df <- df[df$season==season,]
#df.never <- df[df$treatment=="never",]

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

#--------------------------------------#
#### All: del ####
#--------------------------------------#

#--------------------------------------#
#### Calculate del and SE using final data ####
#--------------------------------------#
season <- "Autumn"
setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is

df <- df[df$season==season,]
df <- df[!is.na(df$len_init_QC)&!is.na(df$len_final_QC),]
len_init <- df[df$season==season,]$len_init_QC/10 #converted from mm to cm
len_final <- df[df$season==season,]$len_final_QC/10# converted from mm to cm
mass_dry_final <- df$total_wt_dry

#Estimate del with SE, but this is tricky because it is not linear
# y = (len * del)^3
nonlin_mod=nls((mass_dry_final*conversion_gWW_per_gDW)~(del*len_final)^3,start=list(del=.35))
summary(nonlin_mod)

# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# del 0.335981   0.004185   80.28   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1906 on 42 degrees of freedom

plot(len_final, (mass_dry_final*conversion_gWW_per_gDW),
     ylab = "Total tissue mass (g WW)",
     xlab = "Length")
x <- seq(from = 2, to=3.2, by=.1)
y <- (x*.34)^3
lines(x,y, lty = 2)
y <- (x*.332)^3
lines(x,y, lty = 2)

#---
season <- "Spring"
setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is

df <- df[df$season==season,]
df <- df[!is.na(df$len_init_QC)&!is.na(df$len_final_QC),]
len_init <- df[df$season==season,]$len_init_QC/10 #converted from mm to cm
len_final <- df[df$season==season,]$len_final_QC/10# converted from mm to cm
mass_dry_final <- df$total_wt_dry

#Estimate del with SE, but this is tricky because it is not linear
# y = (len * del)^3
nonlin_mod=nls((mass_dry_final*conversion_gWW_per_gDW)~(del*len_final)^3,start=list(del=.35))
summary(nonlin_mod)

# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# del 0.297987   0.002356   126.5   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.08139 on 36 degrees of freedom



points(len_final, (mass_dry_final*conversion_gWW_per_gDW))
points(len_final, (mass_dry_final*conversion_gWW_per_gDW), pch = 20)
x <- seq(from = 2, to=3.2, by=.1)
y <- (x*.300)^3
lines(x,y, lty = 1)
y <- (x*.296)^3
lines(x,y, lty = 1)


#--------------------------------------#
#### Calculate del and SE using never data ####
#--------------------------------------#
season <- "Autumn"
setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is

df <- df[df$season==season,]
df <- df[df$treatment=="never",]
df <- df[!is.na(df$len_init_QC)&!is.na(df$len_final_QC),]
len_init <- df[df$season==season,]$len_init_QC/10 #converted from mm to cm
len_final <- df[df$season==season,]$len_final_QC/10# converted from mm to cm
mass_dry_final <- df$total_wt_dry

#Estimate del with SE, but this is tricky because it is not linear
# y = (len * del)^3
nonlin_mod=nls((mass_dry_final*conversion_gWW_per_gDW)~(del*len_final)^3,start=list(del=.35))
summary(nonlin_mod)

# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# del 0.333256   0.006388   52.17   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1861 on 14 degrees of freedom

points(len_final, (mass_dry_final*conversion_gWW_per_gDW),
     ylab = "Total tissue mass (g WW)",
     xlab = "Length",
     col = "red")
x <- seq(from = 2, to=3.5, by=.1)
y <- (x*.339)^3
lines(x,y, lty = 2, col = "red")
y <- (x*.327)^3
lines(x,y, lty = 2, col = "red")



season <- "Spring"
setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is

df <- df[df$season==season,]
df <- df[df$treatment=="never",]
df <- df[!is.na(df$len_init_QC)&!is.na(df$len_final_QC),]
len_init <- df[df$season==season,]$len_init_QC/10 #converted from mm to cm
len_final <- df[df$season==season,]$len_final_QC/10# converted from mm to cm
mass_dry_final <- df$total_wt_dry

#Estimate del with SE, but this is tricky because it is not linear
# y = (len * del)^3
nonlin_mod=nls((mass_dry_final*conversion_gWW_per_gDW)~(del*len_final)^3,start=list(del=.35))
summary(nonlin_mod)

# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# del 0.294994   0.003719   79.32   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07869 on 12 degrees of freedom

points(len_final, (mass_dry_final*conversion_gWW_per_gDW), col = "red")
points(len_final, (mass_dry_final*conversion_gWW_per_gDW), pch = 20, col = "red")
x <- seq(from = 2, to=3.5, by=.1)
y <- (x*.298)^3
lines(x,y, lty = 1, col = "red")
y <- (x*.292)^3
lines(x,y, lty = 1, col = "red")


#--------------------------------------#
#### Calculate del and SE using initial reference data ####
#--------------------------------------#

# Spring has no initial tissue weight measurements
season <- "Autumn"
setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="Spring_Fall_reference.csv")  #or wherever this file is
head(df)

df <- df[df$season==season,]
head(df,5)
#df <- df[df$treatment=="never",]
#df <- df[!is.na(df$len_init_QC)&!is.na(df$len_final_QC),]
len_init <- df[df$season==season,]$Length.initial #converted from mm to cm
#len_final <- df[df$season==season,]$len_final_QC/10# converted from mm to cm
mass_dry_init <- df$total_wt_dry

#Estimate del with SE, but this is tricky because it is not linear
# y = (len * del)^3
nonlin_mod=nls((mass_dry_init*conversion_gWW_per_gDW)~(del*len_init)^3,start=list(del=.35))
summary(nonlin_mod)

# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# del 0.303766   0.007193   42.23 3.66e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2355 on 14 degrees of freedom

points(len_init, (mass_dry_init*conversion_gWW_per_gDW),
       ylab = "Total tissue mass (g WW)",
       xlab = "Length",
       col = "blue",
       pch = 20)

x <- seq(from = 2, to=3.5, by=.1)
y <- (x*.311)^3
lines(x,y, lty = 2, col = "blue")
y <- (x*.297)^3
lines(x,y, lty = 2, col = "blue")


#--------------------------------------#
#### Autumn: del ####
#--------------------------------------#
par(mfrow = c(1,2))
#--------------------------------------#
#### Calculate del and SE using final data ####
#--------------------------------------#
season <- "Autumn"
setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is

df <- df[df$season==season,]
df <- df[!is.na(df$len_init_QC)&!is.na(df$len_final_QC),]
len_init <- df[df$season==season,]$len_init_QC/10 #converted from mm to cm
len_final <- df[df$season==season,]$len_final_QC/10# converted from mm to cm
mass_dry_final <- df$total_wt_dry

#Estimate del with SE, but this is tricky because it is not linear
# y = (len * del)^3
nonlin_mod=nls((mass_dry_final*conversion_gWW_per_gDW)~(del*len_final)^3,start=list(del=.35))
summary(nonlin_mod)

# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# del 0.335981   0.004185   80.28   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1906 on 42 degrees of freedom

plot(len_final, (mass_dry_final*conversion_gWW_per_gDW),
     ylab = "Total tissue mass (g WW)",
     xlab = "Length",
     pch = 16,
     xlim = c(2.0,3.5))
x <- seq(from = 2, to=3.5, by=.1)
y <- (x*.34)^3
lines(x,y, lty = 2)
y <- (x*.332)^3
lines(x,y, lty = 2)




#--------------------------------------#
#### Calculate del and SE using never data ####
#--------------------------------------#
season <- "Autumn"
setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is

df <- df[df$season==season,]
df <- df[df$treatment=="never",]
df <- df[!is.na(df$len_init_QC)&!is.na(df$len_final_QC),]
len_init <- df[df$season==season,]$len_init_QC/10 #converted from mm to cm
len_final <- df[df$season==season,]$len_final_QC/10# converted from mm to cm
mass_dry_final <- df$total_wt_dry

#Estimate del with SE, but this is tricky because it is not linear
# y = (len * del)^3
nonlin_mod=nls((mass_dry_final*conversion_gWW_per_gDW)~(del*len_final)^3,start=list(del=.35))
summary(nonlin_mod)

# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# del 0.333256   0.006388   52.17   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1861 on 14 degrees of freedom

points(len_final, (mass_dry_final*conversion_gWW_per_gDW),
       ylab = "Total tissue mass (g WW)",
       xlab = "Length",
       col = "red",
       pch = 16)
x <- seq(from = 2, to=3.5, by=.1)
y <- (x*.339)^3
lines(x,y, lty = 1, col = "red")
y <- (x*.327)^3
lines(x,y, lty = 1, col = "red")






#--------------------------------------#
#### Calculate del and SE using initial reference data ####
#--------------------------------------#

# Spring has no initial tissue weight measurements
season <- "Autumn"
setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="Spring_Fall_reference.csv")  #or wherever this file is
head(df)

df <- df[df$season==season,]
head(df,5)

len_init <- df[df$season==season,]$Length.initial #converted from mm to cm
mass_dry_init <- df$total_wt_dry

#Estimate del with SE, but this is tricky because it is not linear
# y = (len * del)^3
nonlin_mod=nls((mass_dry_init*conversion_gWW_per_gDW)~(del*len_init)^3,start=list(del=.35))
summary(nonlin_mod)

# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# del 0.303766   0.007193   42.23 3.66e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2355 on 14 degrees of freedom

points(len_init, (mass_dry_init*conversion_gWW_per_gDW),
       ylab = "Total tissue mass (g WW)",
       xlab = "Length",
       col = "blue",
       pch = 20)

x <- seq(from = 2, to=3.5, by=.1)
y <- (x*.311)^3
lines(x,y, lty = 2, col = "blue")
y <- (x*.297)^3
lines(x,y, lty = 2, col = "blue")

#--------------------------------------#
#### Spring: del ####
#--------------------------------------#

#--------------------------------------#
#### Calculate del and SE using final data ####
#--------------------------------------#


#---
season <- "Spring"
setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is

df <- df[df$season==season,]
df <- df[!is.na(df$len_init_QC)&!is.na(df$len_final_QC),]
len_init <- df[df$season==season,]$len_init_QC/10 #converted from mm to cm
len_final <- df[df$season==season,]$len_final_QC/10# converted from mm to cm
mass_dry_final <- df$total_wt_dry

#Estimate del with SE, but this is tricky because it is not linear
# y = (len * del)^3
nonlin_mod=nls((mass_dry_final*conversion_gWW_per_gDW)~(del*len_final)^3,start=list(del=.35))
summary(nonlin_mod)

# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# del 0.297987   0.002356   126.5   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.08139 on 36 degrees of freedom



plot(len_final, (mass_dry_final*conversion_gWW_per_gDW), 
       pch = 16,
       ylim = c(0.3,1.4),
       xlim = c(2, 3.5),
     ylab = "Total tissue mass (g WW)",
     xlab = "Length",)
x <- seq(from = 2, to=3.5, by=.1)
y <- (x*.300)^3
lines(x,y, lty = 1)
y <- (x*.296)^3
lines(x,y, lty = 1)



#--------------------------------------#
#### Calculate del and SE using never data ####
#--------------------------------------#


season <- "Spring"
setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is

df <- df[df$season==season,]
df <- df[df$treatment=="never",]
df <- df[!is.na(df$len_init_QC)&!is.na(df$len_final_QC),]
len_init <- df[df$season==season,]$len_init_QC/10 #converted from mm to cm
len_final <- df[df$season==season,]$len_final_QC/10# converted from mm to cm
mass_dry_final <- df$total_wt_dry

#Estimate del with SE, but this is tricky because it is not linear
# y = (len * del)^3
nonlin_mod=nls((mass_dry_final*conversion_gWW_per_gDW)~(del*len_final)^3,start=list(del=.35))
summary(nonlin_mod)

# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# del 0.294994   0.003719   79.32   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07869 on 12 degrees of freedom

points(len_final, (mass_dry_final*conversion_gWW_per_gDW), col = "red")
points(len_final, (mass_dry_final*conversion_gWW_per_gDW), pch = 20, col = "red")
x <- seq(from = 2, to=3.5, by=.1)
y <- (x*.298)^3
lines(x,y, lty = 1, col = "red")
y <- (x*.292)^3
lines(x,y, lty = 1, col = "red")






#--------------------------#
#### Tissue growth graph ####
#--------------------------#

season <- "Autumn"
setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is
season <- "Autumn"
df.Aut <- df[df$season==season,]
season <- "Spring"
df.Spr <- df[df$season==season,]

sigma_Aut <- 0.333
sigma_Spr <- 0.295
conversion_gWW_per_gDW <- 3.918

setwd("~/BE/BE_2018_12_3/figs")
png(filename = "tissue_growth_boxplots.png", width = 2*480, height = 480, units = "px", pointsize = 18)
par(mfrow = c(1,2), oma = c(0,0,0,0), mar = c(3, 4, 2, 0)+.1)

sigma <- sigma_Spr
df <- df.Aut

final_gWW <- df$total_wt_dry*conversion_gWW_per_gDW
init_gWW <- (df$len_init_QC/10*sigma)^3 # length is in mm, convert to cm. 
tissue_growth <- (final_gWW*1000-init_gWW*1000)/df$expt_length
df.growth.Aut <- data.frame(tissue_growth=tissue_growth,
                            treatment=df$treatment,
                            season=df$season
                            )

df$treatment <- ordered(df$treatment, levels = c("never", "weekly","daily"))
boxplot(tissue_growth~df$treatment,
        #ylim = c(0,0.15),
        ylab = "Tissue growth (mgWW / day)")


sigma <- sigma_Spr
df <- df.Spr

final_gWW <- df$total_wt_dry*conversion_gWW_per_gDW
init_gWW <- (df$len_init_QC/10*sigma)^3 # length is in mm, convert to cm. 
tissue_growth <- 10000*(final_gWW-init_gWW)/df$expt_length
df.growth.Spr <- data.frame(tissue_growth=tissue_growth,
                            treatment=df$treatment,
                            season=df$season
                            )

df$treatment <- ordered(df$treatment, levels = c("never", "weekly","daily"))
boxplot(tissue_growth~df$treatment
        #, 
        #ylim = c(0,.15), 
        #yaxt = "n"
)

dev.off()


# 2-way ANOVA
setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is
df$treatment <- ordered(df$treatment, levels = c("never", "weekly","daily"))

new <- rbind(df.growth.Spr, df.growth.Aut)
model <- lm(new$tissue_growth~new$treatment*new$season)
anova(model)
# Analysis of Variance Table

# Response: new$tissue_growth
# Df  Sum Sq Mean Sq F value    Pr(>F)    
# new$treatment             2  1173.9   586.9  1.5605    0.2168    
# new$season                1 23010.0 23010.0 61.1760 2.607e-11 ***
#   new$treatment:new$season  2   839.9   419.9  1.1165    0.3328    
# Residuals                75 28209.6   376.1     


