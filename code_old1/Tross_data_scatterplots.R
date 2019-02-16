# Scater plot of threads vs. growth for both seasons. 
# Note that baseline thread production ("never") was generally lower in the Spring.  
# Mussels are coming out of a period of low food availability 
# and may be allocating more energy towards reproduction at this time of the year. 

season <- "Spring"
setwd("~/BE/BE/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is
head(df,1)
df <- df[df$season==season,]
len_init <- df[df$season==season,]$len_init_QC
len_final <- df[df$season==season,]$len_final_QC
thread_num <- df[df$season==season,]$thread_count_QC
mass_ww_init <- len_init^3*.304
mass_ww_final <- len_final^3*.304 #using the shape coefficient and DW to WW conversion to g
mass_ww <- mass_ww_init
growth_shell <- len_final-len_init
growth_tissue <- mass_ww_final-mass_ww_init
gonad_wt_dry <- df[df$season==season,]$gonad_wt_dry
total_wt_dry <- df[df$season==season,]$total_wt_dry
gonad_proportion <- gonad_wt_dry / total_wt_dry

par(mfrow = c(1,1))
plot(len_init, thread_num, 
     col = as.numeric(df$treatment))
rel <- lm(thread_num~ len_init)
summary(rel)
abline(rel)

plot(log(df$thread_count_QC), log(growth_tissue), col =df$treatment, ylim = c(5,11), xlim = c(2,7))

season <- "Spring"
setwd("~/BE/BE/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is
head(df,1)
df <- df[df$season==season,]
len_init <- df[df$season==season,]$len_init_QC
len_final <- df[df$season==season,]$len_final_QC
thread_num <- df[df$season==season,]$thread_count_QC
mass_ww_init <- len_init^3*.304*3.9
mass_ww_final <- len_final^3*.304*3.9 #using the shape coefficient and DW to WW conversion to g
mass_ww <- mass_ww_init
growth_shell <- len_final-len_init
growth_tissue <- mass_ww_final-mass_ww_init
gonad_wt_dry <- df[df$season==season,]$gonad_wt_dry
total_wt_dry <- df[df$season==season,]$total_wt_dry
gonad_proportion <- gonad_wt_dry / total_wt_dry

points(log(df$thread_count_QC), log(growth_tissue), col =df$treatment, pch = 15)

# gonad proportion ####
season <- "Autumn"
setwd("~/BE/BE/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is
head(df,1)
df <- df[df$season==season,]
len_init <- df[df$season==season,]$len_init_QC
len_final <- df[df$season==season,]$len_final_QC
thread_num <- df[df$season==season,]$thread_count_QC
mass_ww_init <- len_init^3*.304*3.9
mass_ww_final <- len_final^3*.304*3.9 #using the shape coefficient and DW to WW conversion to g
mass_ww <- mass_ww_init
growth_shell <- len_final-len_init
growth_tissue <- mass_ww_final-mass_ww_init
gonad_wt_dry <- df[df$season==season,]$gonad_wt_dry
total_wt_dry <- df[df$season==season,]$total_wt_dry
gonad_proportion <- gonad_wt_dry / total_wt_dry




plot(log(gonad_proportion+1), log(growth_tissue), col =df$treatment,
     xlim = c(.01,.30),
     ylim = c(5,10))

#season <- "Spring"
season <- "Autumn"

setwd("~/BE/BE/Datasets")
df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is
head(df,1)
df <- df[df$season==season,]
len_init <- df[df$season==season,]$len_init_QC
len_final <- df[df$season==season,]$len_final_QC
thread_num <- df[df$season==season,]$thread_count_QC
mass_ww_init <- len_init^3*.304*3.9
mass_ww_final <- len_final^3*.304*3.9 #using the shape coefficient and DW to WW conversion to g
mass_ww <- mass_ww_init
growth_shell <- len_final-len_init
growth_tissue <- mass_ww_final-mass_ww_init
gonad_wt_dry <- df[df$season==season,]$gonad_wt_dry
total_wt_dry <- df[df$season==season,]$total_wt_dry
gonad_proportion <- gonad_wt_dry / total_wt_dry

points(log(gonad_proportion+1), log(growth_tissue), 
       col =df$treatment, pch = 15)

plot(len_init, gonad_proportion,
     col =df$treatment)
rel <- lm(gonad_proportion~len_init*df$treatment)
abline(rel)
anova(rel)
summary(rel)
