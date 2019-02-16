# Oct 31, 2018
# Import the two datasets
# Create two dataframes, one for each season, with (init_len, est_init_mass_, est_change_mass, est_change_len)
# Maybe set up an MLE estimation of average "a" from ultimate size?
# Set up an MLE estimation of food level given growth, either of length or mass. (I'm not sure which?) 
#### Later
# Also look at differences btw. season and thread production
# Also look at differences btw. thread production and growth... would be modifying to be more than 0.08.


# Import the two datasets ####
setwd("~/BE/BE/Datasets")
df <- read.csv("Spring_Fall.csv", stringsAsFactors = FALSE)
head(df,1)
#df <- df[!is.na(df$len_init_QC)|!is.na(df$len_final_QC),]
df$season <- as.factor(df$season)
df$treatment <- as.factor(df$treatment)

shell_growth <- df$len_final_QC-df$len_init_QC
df <- cbind(df,shell_growth)

plot(df$len_init_QC,shell_growth, col = 19+as.numeric(df$season), pch = 1+as.numeric(df$treatment))

plot(shell_growth~season*treatment, data=df, notch=FALSE, 
        col=(c("gold","darkgreen","blue")),
        main="Shell growth", xlab="Treat")

plot(thread_count_QC~season*treatment, data=df, notch=FALSE, 
     col=(c("gold","darkgreen","blue")),
     main="Shell growth", xlab="Treat")

Spring<-df[df$season=="Spring",]
Fall<-df[df$season=="Autumn",]


# Example MLE ####
VB.NLL4 <- function(params, gender, ...) {
  Linfinity <- params[1]
  K <- params[2]
  sigma <- params[3]
  gender <- gender
  LA <- read.csv(file="LengthAge.csv")  #or wherever this file is
  ages <- LA[LA$Gender==gender,]$Ages
  lengths <- LA[LA$Gender==gender,]$Lengths
  model.predL <- Linfinity*(1-exp(-K*ages))
  ndata <- length(ages)
  NLL <- -sum(dnorm(x=lengths, mean=model.predL, sd=sigma, log=TRUE))
  return(NLL)
}

optim(fn=VB.NLL4, # Specify function
      par=c(100, 0.2, 10), gender = "female" ) # PAR are the starting values





#===============================================
# MLE estimation of food level given growth ####
#===============================================
rm(list = ls())
a.est.NLL1 <- function(params, season, reproduction_toggle, ...) {
  # parameters to optimize ###
  food_scalar <- params[1] # intake multiplier
  #cost_per_thread <- params[2]
  sigma <- params[2]
  season <- season
  cost_per_thread=.001 # Will estimate this later
  
  # conversion factors
  conversion_gWW_per_gDW <- 3.918 # From Summer 2015 collection, shape coeff estimation
  shape_coeff <- .304 # From Summer 2015 collection, shape coeff estimation 
  
  # scalars ####
 # a_J_per_mgfood <- 2*4.2 *(1/1000)^(2/3) #2cal/g I X 4.2J/cal X (1g/1000 mg)^2/3 #0.25, brachi, will use M. edulis, Widdows and Bayne 1971
  a_J_per_gfood <- 2*4.2 #2cal/g I X 4.2J/cal X (1g/1000 mg)^2/3 #0.25, brachi, will use M. edulis, Widdows and Bayne 1971
  a_J_per_gfood_per_day <- 2*4.2*24
  #b_per_mg_DW <- 1.7*4.2*1/1000 #.0145, brachi, will use M. edulis, was about .007 per mg dry weight I think , I calculated 0.006 for M. trossulus :)
  #b_J_per_g_WW <-1.7*4.2*1*3.9 
  
  # exponents ####
  d <- 0.67 #0.67 # intake #the model is very sensitive to this. (hope this is an exponent for g WW not mg WW)
  e <- 1 # cost # the model is very sensitive to this (hope this is an exponent for g WW not mg WW)
  
  # calculation of b
  respiration_reference_J_per_day <- 0.07*4.75*4.184*24 # Units: J / (day) from 0.07mlO2/hr, Fly and Hilbish 2013
  respiration_reference_gDW <- 0.24 # Fly and Hilbish 2013
  respiration_reference_gWW <- respiration_reference_gDW * 3.9 # Using wet:dry conversion
  b_permass <- respiration_reference_J_per_day / (respiration_reference_gWW)^e # Note that e is just 1
  b_J_per_g_WW_per_day <- b_permass # Units: J / (gWW * day)
  
  
  # calculation of a
  #eggs <- 0.11 # eggs, brachi
  #AE <- 0.75 # assimilation, Fly & Hilbish - AE vs. assimilation?
  baseline_byssus_multiplier <- 0.08
  

  
  
   
  # Temp response ####
  Tmult_cost <- 1
  Tmult_int <- 1
  
  # other parameters ####
  #en_density_mg_p_J <- 2 #500j/g =.5j/mg, so 2mg/j; #.38 # mg/J, tissue per energy unit, brachi - probably conserved but how did they get this number?
  en_density_g_p_J <- .002 #energy_density
  size_maturity <- 10 # mg, size at maturity, 3000mg from tross del_M calc worksheet, somatic<full tissue weight
  #byssus_multiplier <- 0.10 # fraction of respiration to byssus
  
  
  #stress <- 1 # cost multiplier,  brachi
  reprod_multiplier <- 1 # reproduction multiplier, brachi was 1, but if k is about half and half then this should be more like .5
  opt_size_reducer <- 1 # lets reproduction equal surplus early
  max_time <- 356 # days, How many days to use in the model
  initial_mass_ww <- 1 # mg, What is the initial mass 
  
  shape_coeff <- .304 #converts from L^3 to mass_DW, and is from different dataset
  mass_ww_per_dw <- 3.9 #coverts from mass_DW (dry weight) to mass_WW (wet weight)
  

  
  # calculated scalers
  #a <- a_J_per_food*food*Tmult_int 
  b <- b_J_per_g_WW_per_day*(1+baseline_byssus_multiplier)*Tmult_cost #Wait the byssus multiplier is earlier here than I thought
  #Wopt <- ((b*e)/(a*d))^(1/(d-e)) #Everything is in gWW, so Wopt should be in gWW. 
  Wopt_measured_gDW <- 0.8 #gDW from sample collection
  Wopt_measured_gWW <- Wopt_measured_gDW * conversion_gWW_per_gDW
  
  a_fromWopt <- (b*e)/((Wopt_measured_gWW)^(d-e)*d) # backwards calculation of a
  AE <- 0.75
  I_fromWopt <- a_fromWopt/AE
  a_J_per_day_per_food_scalar <- a_fromWopt
  a <- a_J_per_day_per_food_scalar*food_scalar*Tmult_int
  
  
  
  #season <- "Spring"
  #season <- "Autumn"
  setwd("~/BE/BE/Datasets")
  df <- read.csv(file="Spring_Fall.csv")  #or wherever this file is
  head(df,1)
  df <- df[df$season==season,]
  df <- df[!is.na(df$len_init_QC)&!is.na(df$len_final_QC),]
  len_init <- df[df$season==season,]$len_init_QC/10 #converted from mm to cm
  len_final <- df[df$season==season,]$len_final_QC/10# converted from mm to cm
  thread_num <- df[df$season==season,]$thread_count_QC
  mass_ww_init <- (len_init*.304)^3#*3.9
  mass_ww_final <- (len_final*.304)^3#*3.9 #using the shape coefficient and DW to WW conversion to g
  mass_ww <- mass_ww_init
  growth_shell <- len_final-len_init
  growth_tissue <- mass_ww_final-mass_ww_init
  gonad_wt_dry <- df[df$season==season,]$gonad_wt_dry
  total_wt_dry <- df[df$season==season,]$total_wt_dry
  gonad_proportion <- gonad_wt_dry / total_wt_dry
  
  #plot(as.numeric(df$treatment), gonad_proportion, col = df$len_init_QC)
  # model####
  intake <- a*(mass_ww^d) #good, this is already in mass_ww... hope this is g not mg
  cost <- b*(mass_ww^e)*(1-baseline_byssus_multiplier)
  byssus <- thread_num*cost_per_thread*b*(mass_ww)^e
  # These mussels are maybe reproductive... could look at gonad index. 
  # how to deal with reproduction? Could either toggle this on in the spring and off in the autumn,
  # or could say if gonad proportion is > than 0.1, say then calculate this. Could check out Emily's paper.
  # or could just have it "on" regardless of the conditions
   reproduction <- 0
  # reproduction_mult <- 1
  # if(reproduction_toggle==1){
  #   reproduction <- (intake-cost)*(1-reproduction_mult*(Wopt-mass_ww)/Wopt)
  # }


  
  model.predG_J <- intake-cost-byssus-reproduction #predicts mussel growth in J
  model.predG_g <- model.predG_J*en_density_g_p_J  #predicts mussel growth in g DW or WW???
  ndata <- length(len_init)
  NLL <- -sum(dnorm(x=growth_tissue, mean=model.predG_g, sd=sigma, log=TRUE))
  return(NLL)
}




optim(fn=a.est.NLL1, par=c(food = 1, sigma = 100), season = "Spring") # PAR are the starting values
optim(fn=a.est.NLL1, par=c(food = 1, sigma = 100), season = "Autumn")
# params <- c(food = 1, cost_per_thread=.001, sigma = 100)
# season <- "Autumn"
#optim(fn=a.est.NLL1, par=c(food = 1, byssus_per_thread=.001, sigma = 100), season = "Spring", reproduction_toggle=1) # PAR are the starting values

#==============
# Results
#==============
#==========
# Spring
#==========
# > optim(fn=a.est.NLL1, par=c(food = 1, sigma = 100), season = "Spring") # PAR are the starting values
# $par
# food     sigma 
# 1.5183766 0.0554419 
# 
# $value
# [1] -54.51853
# 
# $counts
# function gradient 
# 109       NA 
# 
# $convergence
# [1] 0
# 
# $message
# NULL

#==========
# Autumn
#==========
# > optim(fn=a.est.NLL1, par=c(food = 1, sigma = 100), season = "Autumn")
# $par
# food      sigma 
# 1.24675183 0.06381785 
# 
# $value
# [1] -57.31035
# 
# $counts
# function gradient 
# 129       NA 
# 
# $convergence
# [1] 0
# 
# $message
# NULL

a.season <- data.frame(
  season = c("Spring","Autumn"),
  a_est = c(1.52,1.25)
)



