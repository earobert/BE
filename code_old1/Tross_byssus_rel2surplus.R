# NEW M. trossulus MODEL FROM Sebens et al_new_byssus_allocation ####
# scalars ####
a_perfood <- 0.25
b_permass <- 0.0145
eggs <- 0.11 # eggs
AE <- 0.55 # assimilation
mortality <- 0.896 # mortality

# exponents ####
d <- 0.67 # intake
e <- 1 # cost

# temp parameters ####
temp <- 16 # temp (was in "other parameters category")
Aarhenius_T <- 6842.7 # Used to change cost with T # From Fly & Hil
ref_T_cost <- 293 #285K=12C # Used to change cost with T # How to calc?
ambient_T_cost <- temp # ambient temp
critical_T_intake <- 10 # this is the optimal intake temp
max_T_intake <- 25 # max temp for intake, shut down
const_Aarhenius <- 0.14
# Could fit these parameters for mussels
Tmult_cost <- const_Aarhenius^(Aarhenius_T/ref_T_cost-Aarhenius_T/(temp+273))
if(temp<critical_T_intake){
  Tmult_int <- Tmult_cost
} else {
  Tmult_int <- Tmult_cost - Tmult_cost*(temp-critical_T_intake)/(max_T_intake-critical_T_intake)
}

# other parameters ####
food <- 1.5 # intake multiplier
en_density <- 0.38 # mg/J, tissue per energy unit 
stress <- 1 # cost multiplier
size_maturity <- 300 # size at maturity
reprod_multiplier <- 1 # reproduction multiplier
opt_size_reducer <- 1 # lets reproduction equal surplus early
byssus_multiplier <- 0.25 # fraction of respiration to byssus
max_time <- 1500 # days, How many days to use in the model
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
    mass_ww[i] <- mass_ww[i-1]+en_density*(intake[i-1]-cost[i-1]-reproduction[i-1]-byssus[i-1]) #in mg
  }
  
  intake[i] <- a*(mass_ww[i]^d)
  cost[i] <- b*(mass_ww[i]^e)
  E_surplus[i] <- intake[i]-cost[i]
  
  # Calculate reproduction
  if(mass_ww[i]<size_maturity){
    reproduction[i] <- 0
   } else {
     reproduction[i] <- E_surplus[i]*(1-reprod_multiplier*((Wopt-opt_size_reducer*mass_ww[i])/Wopt))
   }
  byssus[i] <- (E_surplus[i]-reproduction[i])*0.1 # Here byssus is upregulated when there is more energy available
  
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
     ylim = c(0,10),
     ylab = "Energy Units (J)",
     xlab = "Time")
lines(time,byssus, lwd = 5, col = "grey")
lines(time,E_surplus, lwd = 2)
lines(time,reproduction, lwd = 1, lty = 3, col = "black")

plot(time,byssus, type = "l")

# Von-Bertlanffy growth
plot(time,mass_ww, type = "l")

# EOS within realistic range of survival
plot(mass_ww, cost, lty = 2, type = "l", ylim = c(0,100))
lines(mass_ww,intake)
lines(mass_ww, E_surplus, lwd = 5)

# Maybe this relates to the number of threads but not the type of threads...
plot(E_surplus,byssus)

