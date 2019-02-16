#_____________________________________________________________________----
# This is the main working model right now. Now I'm adding time as an aspect in... 
#===============================================#
# Model 5 - there is no gamma, and no cost of threads
# MLE estimation of cost of byssus given growth and thread production####
# With separate estimates of a and the cost of thread production 
# BUT including a baseline value
# -use Ken's suggestion of having their be a baseline... Tried that here. Doesn't change anything.
# I could assume that 0 thread cutting was for that baseline value. 
#===============================================#



mod.1 <- function(params, season) {
  food_scalar <- params[1] 
  cost_per_thread <- params[2]
  #baseline_byssus_multiplier <- 0.08
  sigma <- params[3]
  season <- season
  
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
  a <- a_J_per_day_per_food_scalar*food_scalar*Tmult_int
  
  # import data ####
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
  mass_ww_init <- (len_init*.304)^3
  mass_ww_final <- (len_final*.304)^3
  mass_ww <- mass_ww_init
  growth_shell <- len_final-len_init
  growth_tissue <- mass_ww_final-mass_ww_init
  gonad_wt_dry <- df[df$season==season,]$gonad_wt_dry
  total_wt_dry <- df[df$season==season,]$total_wt_dry
  gonad_proportion <- gonad_wt_dry / total_wt_dry
  
  # model####
  intake <- a*(mass_ww^d) #good, this is already in mass_ww... hope this is g not mg
  cost <- b*(mass_ww^e) #*(1-baseline_byssus_multiplier)
  # byssus_baseline <- b*(mass_ww^e)*baseline_byssus_multiplier
  byssus_induced <- thread_num*cost_per_thread
  reproduction <- 0
  model.predG_J <- intake-cost-byssus_induced-reproduction #predicts mussel growth in J
  model.predG_g <- model.predG_J*en_density_g_p_J  #predicts mussel growth in g DW or WW???
  out <- data.frame(
    growth_tissue = growth_tissue,
    model.predG_g = model.predG_g
  )
  return(out)
}

#For troubleshooting:
params <- c(food = 3, cost_induced_byssus = 0.01, sigma = 100)
season <- "Autumn"
out <- mod.1(params,season)
growth_tissue <- out$growth_tissue
model.predG_g <- out$model.predG_g
#plot(growth_tissue,model.predG_g)

#
a.est.NLL1 <- function(params, season) {
  out <- mod.1(params,season)
  growth_tissue <- out$growth_tissue
  model.predG_g <- out$model.predG_g
  sigma <- params[3]
  NLL <- -sum(dnorm(x=growth_tissue, mean=model.predG_g, sd=sigma, log=TRUE))
  return(NLL)
}

Autumn.est <- optim(fn=a.est.NLL1, par=c(food = 3, cost_induced_byssus = 0.01, sigma = 100), season = "Autumn") # par are the starting values
Spring.est <- optim(fn=a.est.NLL1, par=c(food = 3, cost_induced_byssus = 0.01, sigma = 100), season = "Spring") # par are the starting values


#=====================
#plot####
#=====================

a.plot.NLL1 <- function(params, season) {
  out <- mod.1(params,season)
  growth_tissue <- out$growth_tissue
  model.predG_g <- out$model.predG_g
  
  # import data ####
  setwd("~/BE/BE/Datasets")
  df <- read.csv(file="Spring_Fall.csv")
  df <- df[df$season==season,]
  df <- df[!is.na(df$len_init_QC)&!is.na(df$len_final_QC),]
  
  p.1 <- lm(growth_tissue~ df$thread_count_QC)
  plot(x = df$thread_count_QC, y = growth_tissue, 
       col = as.numeric(df$treatment)+26,
       pch = 19,
       xlab = "Thread production (#)",
       ylab = "Observed growth (gWW)"
  )
  abline(p.1)
  
  p.2 <- lm(model.predG_g~growth_tissue)
  plot(x = growth_tissue, y = model.predG_g, 
       col = as.numeric(df$treatment)+26,
       pch = 19,
       ylim = c(0,.3), xlim = c(0,.3),
       xlab = "Observed growth (gWW)",
       ylab = "Predicted growth (gWW)"
  )
  abline(p.2)
  x <- seq(from = -.4, to=.4, by=.1)
  lines(x,x, lty = 2)
}

dev.off()
par(mfrow = c(2,2),mar = c(5,4,4,2)+0.1)
a.plot.NLL1(par=c(food = Autumn.est$par[1], cost_per_byssus = Autumn.est$par[2], sigma = Autumn.est$par[3]), season = "Autumn")
a.plot.NLL1(par=c(food = Spring.est$par[1], cost_per_byssus = Spring.est$par[2], sigma = Spring.est$par[3]), season = "Spring")


# plot cost using estimated costs. 
#Spring####
season <- "Spring"

# import data ####
setwd("~/BE/BE/Datasets")
df <- read.csv(file="Spring_Fall.csv", stringsAsFactors = FALSE)
df$treatment <- as.factor(df$treatment)
df$season <- as.factor(df$season)
df$treatment <- ordered(df$treatment, levels = c("never", "weekly","daily"))
df <- df[df$season==season,]
df <- df[!is.na(df$len_init_QC)&!is.na(df$len_final_QC),]
df.never <- df[df$treatment=="never",]


plot(df$treatment, df$thread_count_QC*Spring.est$par[2], ylab = "Cost (J)")

params <- Spring.est$par

food_scalar <- params[1] 
cost_per_thread <- params[2]
#baseline_byssus_multiplier <- 0.08
sigma <- params[3]
season <- season

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
en_density_g_p_J <- .002 *29 #energy_density
size_maturity <- 10 # mg, size at maturity, 3000mg from tross del_M calc worksheet, somatic<full tissue weight

# calculated scalers
b <- b_J_per_g_WW_per_day*Tmult_cost #Wait the byssus multiplier is earlier here than I thought

Wopt_measured_gDW <- 0.8 #gDW from sample collection
Wopt_measured_gWW <- Wopt_measured_gDW * conversion_gWW_per_gDW
a_fromWopt <- (b*e)/((Wopt_measured_gWW)^(d-e)*d) # backwards calculation of a
a_J_per_day_per_food_scalar <- a_fromWopt
a <- a_J_per_day_per_food_scalar*food_scalar*Tmult_int

# import data ####
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
mass_ww_init <- (len_init*.304)^3
mass_ww_final <- (len_final*.304)^3
mass_ww <- mass_ww_init
growth_shell <- len_final-len_init
growth_tissue <- mass_ww_final-mass_ww_init
gonad_wt_dry <- df[df$season==season,]$gonad_wt_dry
total_wt_dry <- df[df$season==season,]$total_wt_dry
gonad_proportion <- gonad_wt_dry / total_wt_dry

# model####
intake <- a*(mass_ww^d) #good, this is already in mass_ww... hope this is g not mg
cost <- b*(mass_ww^e) #*(1-baseline_byssus_multiplier)
# byssus_baseline <- b*(mass_ww^e)*baseline_byssus_multiplier
byssus_induced <- thread_num*cost_per_thread
reproduction <- 0
model.predG_J <- intake-cost-byssus_induced-reproduction #predicts mussel growth in J
model.predG_g <- model.predG_J*en_density_g_p_J  #predicts mussel growth in g DW or WW???

length(df$treatment)
length(model.predG_J)


new <- data.frame(treat = df$treatment, pred_J =model.predG_J)
new$treat <- as.factor(new$treat)

dev.off()
par(mfrow = c(2,2))

plot(new$treat, intake, ylim = c(0,120), ylab = "intake (J)")

plot(new$treat,byssus_induced, ylim = c(0,120), ylab = "cost byssus (J)")

plot(new$treat,cost, ylim = c(0,120), ylab = "cost non-byssus (J)")

plot(new$treat,new$pred_J, ylim = c(0,120), ylab = "surplus (J)")


dev.off()
par(mfrow = c(2,2))

prop_byss <- byssus_induced / (intake)
plot(new$treat, prop_byss*100, ylim = c(0,40),  ylab = "% energy to byssus")

prop_predG_J <- new$pred_J/ intake
plot(new$treat,prop_predG_J*100, ylim = c(0,80), ylab = "% energy to growth")

prop_byss <- cost / (intake)
plot(new$treat, prop_byss*100, ylim = c(0,100), ylab = "% energy to non-byssus costs")

plot(df$treatment, gonad_proportion*100*prop_byss, ylim = c(0,100), ylab = "% energy to reproduction")

# plot cost using estimated costs. 
#Autumn####
season <- "Autumn"

# import data ####
setwd("~/BE/BE/Datasets")
df <- read.csv(file="Spring_Fall.csv", stringsAsFactors = FALSE)
df$treatment <- as.factor(df$treatment)
df$season <- as.factor(df$season)
df$treatment <- ordered(df$treatment, levels = c("never", "weekly","daily"))
df <- df[df$season==season,]
df <- df[!is.na(df$len_init_QC)&!is.na(df$len_final_QC),]
df.never <- df[df$treatment=="never",]


plot(df$treatment, df$thread_count_QC*Autumn.est$par[2], ylab = "Cost (J)")

params <- Autumn.est$par

food_scalar <- params[1] 
cost_per_thread <- params[2]
#baseline_byssus_multiplier <- 0.08
sigma <- params[3]
season <- season

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
a <- a_J_per_day_per_food_scalar*food_scalar*Tmult_int

# import data ####
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
mass_ww_init <- (len_init*.304)^3
mass_ww_final <- (len_final*.304)^3
mass_ww <- mass_ww_init
growth_shell <- len_final-len_init
growth_tissue <- mass_ww_final-mass_ww_init
gonad_wt_dry <- df[df$season==season,]$gonad_wt_dry
total_wt_dry <- df[df$season==season,]$total_wt_dry
gonad_proportion <- gonad_wt_dry / total_wt_dry

# model####
intake <- a*(mass_ww^d) #good, this is already in mass_ww... hope this is g not mg
cost <- b*(mass_ww^e) #*(1-baseline_byssus_multiplier)
# byssus_baseline <- b*(mass_ww^e)*baseline_byssus_multiplier
byssus_induced <- thread_num*cost_per_thread
reproduction <- 0
model.predG_J <- intake-cost-byssus_induced-reproduction #predicts mussel growth in J
model.predG_g <- model.predG_J*en_density_g_p_J  #predicts mussel growth in g DW or WW???

length(df$treatment)
length(model.predG_J)


new <- data.frame(treat = df$treatment, pred_J =model.predG_J)
new$treat <- as.factor(new$treat)

dev.off()
par(mfrow = c(2,2))

plot(new$treat, intake, ylim = c(0,120), ylab = "intake (J)")

plot(new$treat,byssus_induced, ylim = c(0,120), ylab = "cost byssus (J)")

plot(new$treat,cost, ylim = c(0,120), ylab = "cost non-byssus (J)")

plot(new$treat,new$pred_J, ylim = c(0,120), ylab = "surplus (J)")


dev.off()
par(mfrow = c(2,2))

prop_byss <- byssus_induced / (intake)
plot(new$treat, prop_byss*100, ylim = c(0,40),  ylab = "% energy to byssus")

prop_predG_J <- new$pred_J/ intake
plot(new$treat,prop_predG_J*100, ylim = c(0,80), ylab = "% energy to growth")

prop_byss <- cost / (intake)
plot(new$treat, prop_byss*100, ylim = c(0,100), ylab = "% energy to non-byssus costs")

plot(df$treatment, gonad_proportion*100*prop_byss, ylim = c(0,100), ylab = "% energy to reproduction")




# Some things to try... 
# -use frequency of manipulation rather than number of threads (.3per week, 1per week, 7per week) ... 
# -use Ken's suggestion of having their be a baseline... Tried that here.
# -calculate with same cost of threads
# -perform a monte carlo on a and threads



# Some things to try... 
# -use frequency of manipulation rather than number of threads (.3per week, 1per week, 7per week) ... 
# -use Ken's suggestion of having their be a baseline... Tried that here.
# -calculate with same cost of threads
# -perform a monte carlo on a and threads
