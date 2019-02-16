# ORIGINAL Brachi MODEL FROM Sebens et al ####

# Modified the growth calculation 
# to be in J not weight. (after transferring over)


# scalars ####
a_perfood <- 0.25 #the model is relatively sensitive, don't see AE anywhere
b_permass <- 0.0145 #the model is also relatively sensitive to this
eggs <- 0.11 # eggs
AE <- 0.75 # assimilation, was 55, but in the byssus paper it's 75, this has no affect on this graph...
mortality <- 0.896 # mortality #this has no affect on the EOS graph

# exponents ####
d <- 0.67 #0.67 # intake #the model is very sensitive to this. 
e <- 1 # cost # the model is very sensitive to this

# temp parameters ####
temp <- 16 # temp (was in "other parameters category")
Aarhenius_T <- 8000 # Used to change cost with T
ref_T_cost <- 285 #285K=12C # Used to change cost with T
ambient_T_cost <- temp # ambient temp
critical_T_intake <- 16 # this is the optimal intake temp
max_T_intake <- 25 # max temp for intake, shut down
const_Aarhenius <- 2.71816
# Could fit these parameters for mussels
Tmult_cost <- const_Aarhenius^(Aarhenius_T/ref_T_cost-Aarhenius_T/(temp+273))
if(temp<critical_T_intake){
  Tmult_int <- Tmult_cost
  } else {
    Tmult_int <- Tmult_cost - Tmult_cost*(temp-critical_T_intake)/(max_T_intake-critical_T_intake)
  }
#
#Tmult_cost <- 1

# other parameters ####
food <- 1.5 # intake multiplier
en_density <- 0.38 # mg/J, tissue per energy unit 
stress <- 1 # cost multiplier
size_maturity <- 300 #300 # size at maturity
reprod_multiplier <- 1 # reproduction multiplier, could be advantageous for species who had high mortality and were optimized for reproduction
opt_size_reducer <- 1 # lets reproduction equal surplus early, would reduce total size which could be advantagous for some species
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
      reproduction[i] <- 1*E_surplus[i]*(1-reprod_multiplier*((Wopt-opt_size_reducer*mass_ww[i])/Wopt))
      #reproduction[i] <- 0

            }
}


# calculate relationships
W_growth <- rep(x = -99, len = max_time)
W_growth[2:max_time] <- mass_ww[2:max_time] - mass_ww[1:(max_time-1)]
W_growth[1:1] <- 0  # initial growth at very first timepoint is 0
E_growth <- W_growth / en_density
  
# Sebens plots ####  
length(time)
length(mass_ww)
plot(time,E_growth, type = "l", 
     lty = 3,
     lwd = 5,
     ylim = c(0,80),
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
lines(mass_ww, E_growth, lwd = 3)
lines(mass_ww, E_surplus, lwd = 5)

# MODIFICATION_Past ultimate size ####
# Longer time for optimal EOS ####
max_time <- 1700

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

# run model #####
for(i in 1:(max_time)){
  # Calculate mass
  if(i == 1){
    mass_ww[i] <- initial_mass_ww
  } else {
    mass_ww[i] <- mass_ww[i-1]+50
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
# plot EOS ####
plot(mass_ww, cost, lty = 2, 
     type = "l", xlim = c(0,11000), 
     ylim = c(0,350),
     yaxs = "i",
     xaxs = "i")
lines(mass_ww,intake)
lines(mass_ww, E_surplus, lwd = 5)

# MODIFICATION_Temperatures ####
# temp parameters ####
temp <- 16 # temp (was in "other parameters category")
Aarhenius_T <- 8000 # Used to change cost with T
ref_T_cost <- 285 #285K=12C # Used to change cost with T 
ambient_T_cost <- temp # ambient temp
critical_T_intake <- 16 # this is the optimal intake temp
max_T_intake <- 25 # max temp for intake, shut down
const_Aarhenius <- 2.71816
# Could fit these parameters for mussels

temp <- seq(from = 5, to = 20, length.out = 20)
# set up vectors ####

# assimilated intake with temp effect, 
# minus basal metabolic cost, 
# mult by temp effect and byssus effect
E_surplus <- seq(from = 5, to = 30, length.out = 20)

intake <- seq(from = 5, to = 30, length.out = 20)
cost <- seq(from = 5, to = 30, length.out = 20)
byssus <- seq(from = 5, to = 30, length.out = 20)
reproduction <- seq(from = 5, to = 30, length.out = 20)

mass_ww <- 300
for(i in 1:length(temp)){
  Tmult_cost[i] <- const_Aarhenius^(Aarhenius_T/ref_T_cost-Aarhenius_T/(temp[i]+273))
  if(temp[i]<critical_T_intake){
    Tmult_int[i] <- Tmult_cost[i]
  } else {
    Tmult_int[i] <- Tmult_cost[i] - Tmult_cost[i]*(temp[i]-critical_T_intake)/(max_T_intake-critical_T_intake)
  }
  a[i] <- a_perfood*food*Tmult_int[i]
  b[i] <- b_permass*stress*(1+byssus_multiplier)*Tmult_cost[i]
  Wopt[i] <- ((b[i]*e)/(a[i]*d))^(1/(d-e))
  
  intake[i] <- a[i]*(mass_ww^d)
  cost[i] <- b[i]*(mass_ww^e)
  E_surplus[i] <- intake[i]-cost[i]
  byssus[i] <- cost[i]-cost[i]/(1+byssus_multiplier) # Here byssus is a fraction f the cost
  
  # Calculate reproduction
  if(mass_ww<size_maturity){
    reproduction[i] <- 0
  } else {
    reproduction[i] <- E_surplus[i]*(1-reprod_multiplier*
             ((Wopt[i]-opt_size_reducer*mass_ww[i])/Wopt[i]))
  }
  
}


plot(temp,Tmult_cost, type = "l", ylim = c(0,3), xlim = c(5,20))
lines(temp,Tmult_int)

plot(temp,intake, type = "l", ylim = c(0,30))
lines(temp,cost)
lines(temp,E_surplus)

plot(temp,byssus, type = "l")

plot(E_surplus,byssus,type = 'l', xlim = c(0,18))










# Plot ####
# Plot the two temperature correction functions
plot(temp, Tmult_cost, type = "l")
lines(temp, Tmult_int)

# MODIFICATION_Lifetime Temperatures ####
# run model ####
# set up vectors ####
max_time = 1500
time <- seq(from = 1, to = max_time, by = 1)
temp <- seq(from = 5, to = 20, length.out = 20)
a <- rep(x = 0, length.out = length(temp))
b <- rep(x = 0, length.out = length(temp))
Wopt <- rep(x = 0, length.out = length(temp))
# mass mg wet weight


# assimilated intake with temp effect, 
# minus basal metabolic cost, 
# mult by temp effect and byssus effect
mass_ww <- matrix(data=NA, nrow = length(time), ncol = length(temp))
intake <- matrix(data=NA, nrow = length(time), ncol = length(temp))
cost <- matrix(data=NA, nrow = length(time), ncol = length(temp))
E_surplus <- matrix(data=NA, nrow = length(time), ncol = length(temp))
byssus <- matrix(data=NA, nrow = length(time), ncol = length(temp))
reproduction <- matrix(data=NA, nrow = length(time), ncol = length(temp))
str(mass_ww)

# run model ####
for(j in 1:length(temp)){
  temp[j]
  Tmult_cost[j] <- const_Aarhenius^(Aarhenius_T/ref_T_cost-Aarhenius_T/(temp[j]+273))
  if(temp[j]<critical_T_intake){
    Tmult_int[j] <- Tmult_cost[j]
  } else {
    Tmult_int[j] <- Tmult_cost[j] - Tmult_cost[j]*(temp[j]-critical_T_intake)/(max_T_intake-critical_T_intake)
  }
  a[j] <- a_perfood*food*Tmult_int[j]
  b[j] <- b_permass*stress*(1+byssus_multiplier)*Tmult_cost[j]
  Wopt[j] <- ((b[j]*e)/(a[j]*d))^(1/(d-e))
for(i in 1:(max_time)){
  # Calculate mass
  if(i == 1){
    mass_ww[i,j] <- initial_mass_ww
  } else {
    mass_ww[i,j] <- mass_ww[i-1,j]+en_density*(intake[i-1,j]-cost[i-1,j]-reproduction[i-1,j]) #in mg
  }
  intake[i,j] <- a[j]*(mass_ww[i,j]^d)
  cost[i,j] <- b[j]*(mass_ww[i,j]^e)
  E_surplus[i,j] <- intake[i,j]-cost[i,j]
  byssus[i,j] <- cost[i,j]-cost[i,j]/(1+byssus_multiplier) # Here byssus is a fraction f the cost
  
  # Calculate reproduction
  if(mass_ww[i,j]<size_maturity){
    reproduction[i,j] <- 0
  } else {
    reproduction[i,j] <- E_surplus[i,j]*(1-reprod_multiplier*((Wopt[j]-opt_size_reducer*mass_ww[i,j])/Wopt[j]))
  }
}
}


# Plot ####
# If mussels held at one temperature:
plot(temp,Wopt, type = "l", ylim = c(0,3000), lty = 2, ylab = "Yearly mass (ww)")
lines(temp, mass_ww[1350,], lwd = 4)
lines(temp, mass_ww[1000,], lwd = 3)
lines(temp, mass_ww[700,], lwd = 2)
lines(temp, mass_ww[350,], lwd = 1)



