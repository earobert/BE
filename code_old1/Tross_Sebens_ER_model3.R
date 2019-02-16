# NEW M. trossulus MODEL FROM Sebens et al_new_byssus_allocation ####
# scalars ####
# a per food
AE <- 0.8 # ranges from 0.75 to 0.85 in Fly and Hilbish
CR_max_ml_pmin <- 5.003 #ml/min from Fly and Hilbish
CR_max_ml_pday <- CR_max_ml_pmin * 60 * 24 #ml/day
CR_max_L_pday <- CR_max_ml_pday/1000 #L/day
#food <- 1.6 #mg per L of SW in Fly and Hilbish
food <- 1.3*.1 #mg/L Estimate for experiment... probably different
#maybe try converting mg to J

a_perfood <- CR_max_L_pday
  #references: Fly and Hilbish, (Bayne, Arifin and Bendell-Young)

# b per mass
# No reference for tross, maybe bayne
b_permass <- 0.0145 # convert cost of tissue from ___ into J

# exponent d
# No reference for tross? Check Bayne, Denny
d <- 0.67 # intake
e <- 1 # cost, this makes sense, cost relates to volume


# energy density
# No reference for tross, probably bayne
en_density <- 0.38 #mg/J

# reproduction
# No reference for tross, probably bayne
eggs <- 11
reprod_multiplier <- 1
opt_size_reducer <- 1

# stress
stress <- 1

# byssus
byssus_multiplier <- 0.1 # Fraction of respiration to byssus


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

# calculated scalers
a <- a_perfood*food*Tmult_int
b <- b_permass*stress*(1+byssus_multiplier)*Tmult_cost
(Wopt <- ((b*e)/(a*d))^(1/(d-e)))

# model scope
max_time <- 1500 # days, How many days to use in the model
initial_mass_ww <- 1 # mg, What is the initial mass 

# set up vectors ####
time <- seq(from = 1, to = max_time, by = 1)
mass_ww <- rep(x = -99, len = max_time)
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
plot(time,E_growth*1000, type = "l", 
     lty = 3,
     lwd = 5,
     ylim = c(0,10),
     ylab = "Energy Units (J)",
     xlab = "Time")
lines(time,byssus*1000, lwd = 5, col = "grey")
lines(time,E_surplus*1000, lwd = 2)
lines(time,reproduction*1000, lwd = 1, lty = 3, col = "black")


