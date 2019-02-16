# NEW_but older version of M. trossulus MODEL FROM Sebens et al ####
# This time I'm doing an optimization to estimate reproduction 
# Question - if I change food, then that changes my estimate of reproduction, right? 
# Food should just be 1... weird that it is 1.5 in the paper... I guess I'll leave it at 

# age at repro
# setting reproduction
# check labels on graph - they are mislabeled
# 


# scalars ####
a_perfood <- 0.25 #0.25, brachi
b_permass <- 0.0145 #brachi
eggs <- 0.11 # eggs, brachi
AE <- 0.75 # assimilation, Fly & Hilbish - AE vs. assimilation?
mortality <- 0.896 # mortality, brachi

# exponents ####
d <- 0.67 # intake, brachi, makes sense with surface area
e <- 1 # cost, brachi, makes sense with volume

# temp parameters #### 
# guestimated. 
# Right I'm now using 10deg as a ref temp, a factor of 1.
temp <- 10 # temp (was in "other parameters category")
Aarhenius_T <- 6842.7 # Used to change cost with T # From Fly & Hilbish
ref_T_cost <- 283 #283K=10C # Used to change cost with T # At this temp, coeff is 1
ambient_T_cost <- temp # ambient temp
critical_T_intake <- 20 # this is the optimal intake temp
max_T_intake <- 30 # max temp for intake, shut down
const_Aarhenius <- 2
# Could fit these parameters for mussels

# Checked fcn here, used this to see if it is reasonable:
# temp<- seq(from = 5, to = 25, by = 1)
# 
# Tmult_cost <- const_Aarhenius^(Aarhenius_T/ref_T_cost-Aarhenius_T/(temp+273))
# for(i in 1:length(temp)){
# if(temp[i]<critical_T_intake){
#   Tmult_int[i] <- Tmult_cost[i]
# } else {
#   Tmult_int[i] <- Tmult_cost[i] - Tmult_cost[i]*(temp[i]-critical_T_intake)/(max_T_intake-critical_T_intake)
# }
# }
# plot(temp, Tmult_int, ylim = c(0,2))

Tmult_cost <- const_Aarhenius^(Aarhenius_T/ref_T_cost-Aarhenius_T/(temp+273))
if(temp<critical_T_intake){
  Tmult_int <- Tmult_cost
} else {
  Tmult_int <- Tmult_cost - Tmult_cost*(temp-critical_T_intake)/(max_T_intake-critical_T_intake)
}


# other parameters ####
food <- 1.5 # intake multiplier, brachi
en_density <- 0.38 # mg/J, tissue per energy unit, brachi
stress <- 1 # cost multiplier,  brachi
size_maturity <- 3000 # mg, size at maturity, from tross del_M calc worksheet, somatic<full tissue weight
reprod_multiplier <- .5 # reproduction multiplier, brachi was 1, but if k is about half and half then this should be more like .5
opt_size_reducer <- 1 # lets reproduction equal surplus early
byssus_multiplier <- 0.10 # fraction of respiration to byssus
max_time <- 1700 # days, How many days to use in the model
initial_mass_ww <- 1 # mg, What is the initial mass 

# calculated scalers
a <- a_perfood*food*Tmult_int
b <- b_permass*stress*(1+byssus_multiplier)*Tmult_cost
(Wopt <- ((b*e)/(a*d))^(1/(d-e)))

# set up vectors ####
time <- seq(from = 1, to = max_time, by = 1)

# mass mg wet weight
mass_ww <- rep(x = -99, len = max_time)

# assimilated intake with temp effect, 
# minus basal metabolic cost, 
# mult by temp effect and byssus effect
E_surplus <- rep(x = -99, len = max_time)

intake <- rep(x = -99, len = max_time)
cost <- rep(x = -99, len = max_time)
byssus <- rep(x = -99, len = max_time)
reproduction <- rep(x = -99, len = max_time)

# run model ####
for(i in 1:(max_time)){
  # Calculate mass
  if(i == 1){
    mass_ww[i] <- initial_mass_ww
  } else {
    mass_ww[i] <- mass_ww[i-1]+en_density*(intake[i-1]-cost[i-1]-reproduction[i-1]) #in mg
  }
  
  intake[i] <- a*(mass_ww[i]^d)
  cost[i] <- b*(mass_ww[i]^e)
  #E_surplus[i] <- a*(mass_ww[i]^d)-b*(mass_ww[i]^e) #in J, Is this necessary?
  E_surplus[i] <- intake[i]-cost[i]
  byssus[i] <- cost[i]-cost[i]/(1+byssus_multiplier) # Here byssus is a fraction f the cost
  
}

E_surplus_max <- max(E_surplus)
reprod_multiplier <- E_surplus_max

for(i in 1:(max_time)){
# Calculate reproduction
if(mass_ww[i]<size_maturity){
  reproduction[i] <- 0
} else {
  reproduction[i] <- E_surplus[i]*(1-reprod_multiplier*((Wopt-opt_size_reducer*mass_ww[i])/Wopt))
}
}



# calculate relationships
E_growth <- rep(x = -99, len = max_time)
E_growth[2:max_time] <- mass_ww[2:max_time] - mass_ww[1:(max_time-1)]
E_growth[1:1] <- 0  # initial growth at very first timepoint is 0

# Sebens plots ####  
length(time)
length(mass_ww)
plot(time,E_growth, type = "l", 
     lty = 3,
     lwd = 5,
     ylim = c(0,40),
     ylab = "Energy Units (J)",
     xlab = "Time")
lines(time,byssus, lwd = 5, col = "grey")
lines(time,E_surplus, lwd = 2)
lines(time,reproduction, lwd = 1, lty = 3, col = "black")

# Von-Bertlanffy growth
plot(time,mass_ww, type = "l")

# EOS within realistic range of survival
plot(mass_ww, cost, lty = 2, type = "l", ylim = c(0,140))
lines(mass_ww,intake)
lines(mass_ww, E_surplus, lwd = 5)