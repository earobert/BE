# NEW M. trossulus MODEL FROM Sebens et al ####
# I'M NOT getting stuff that makes sense. Why does b_permass have to be so small, is it per mg not per g?
# 
# age at repro
# setting reproduction
# check labels on graph - they are mislabeled
# 


# scalars ####
a_perfood <- 2*4.2 *(1/1000)^(2/3) #2cal/g I X 4.2J/cal X (1g/1000 mg)^2/3 #0.25, brachi, will use M. edulis, Widdows and Bayne 1971
b_permass <- 1.7*4.2*1/1000 #.0145, brachi, will use M. edulis, was about .007, I calculated 0.006 for M. trossulus :)
eggs <- 0.11 # eggs, brachi
AE <- 0.75 # assimilation, Fly & Hilbish - AE vs. assimilation?
mortality <- 0.896 # mortality, brachi

# exponents ####
d <- 0.67 # intake, brachi, makes sense with surface area, tross might be 0.70
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
food <- 1 # intake multiplier, brachi was 1.5
en_density <- 2 #500j/g =.5j/mg, so 2mg/j; #.38 # mg/J, tissue per energy unit, brachi - probably conserved but how did they get this number?
stress <- 1 # cost multiplier,  brachi
size_maturity <- 10 # mg, size at maturity, 3000mg from tross del_M calc worksheet, somatic<full tissue weight
reprod_multiplier <- 1 # reproduction multiplier, brachi was 1, but if k is about half and half then this should be more like .5
opt_size_reducer <- 1 # lets reproduction equal surplus early
byssus_multiplier <- 0.10 # fraction of respiration to byssus
max_time <- 356 # days, How many days to use in the model
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
plot(time,E_surplus, type = "l", 
     lty = 1,
     lwd = 2,
     #ylim = c(0,40),
     ylab = "Energy Units (J)",
     xlab = "Time")
lines(time,byssus, lwd = 5, col = "grey")
lines(time,E_growth, lwd = 5, lty = 3)
lines(time,reproduction, lwd = 1, lty = 3, col = "black")

# Von-Bertlanffy growth
plot(time,mass_ww, type = "l")

# EOS within realistic range of survival
plot(mass_ww, intake, lty = 1, type = "l")
lines(mass_ww,cost, lty = 2)
lines(mass_ww, E_surplus, lwd = 5)
lines(mass_ww, reproduction, lwd = 1, lty = 4)
lines(mass_ww, byssus, lwd = 1, lty = 3)

# MODIFICATION_Temperatures ####
# temp parameters ####
Aarhenius_T <- 6842.7 # Used to change cost with T
ref_T_cost <- 285 #285K=12C # Used to change cost with T 
ambient_T_cost <- temp # ambient temp
critical_T_intake <- 10 # this is the optimal intake temp
max_T_intake <- 25 # max temp for intake, shut down
const_Aarhenius <- 2.71816
# Could fit these parameters for mussels

temp <- seq(from = 12, to = 21, length.out = 20)
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


plot(temp,Tmult_cost, type = "l", ylim = c(0,3), xlim = c(10,25))
lines(temp,Tmult_int)

plot(temp,intake, type = "l", ylim = c(0,30))
lines(temp,cost)
lines(temp,E_surplus)

plot(temp,byssus, type = "l")

plot(E_surplus,byssus,type = 'l', xlim = c(0,18))

# MODIFICATION_Temperatures_different_a ####
# temp parameters ####
Aarhenius_T <- 6842.7 # Used to change cost with T
ref_T_cost <- 285 #285K=12C # Used to change cost with T 
ambient_T_cost <- temp # ambient temp
critical_T_intake <- 10 # this is the optimal intake temp
max_T_intake <- 25 # max temp for intake, shut down
const_Aarhenius <- 2.71816
# Could fit these parameters for mussels

temp <- seq(from = 12, to = 21, length.out = 20)
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

plot(temp,Tmult_cost, type = "l", ylim = c(0,3), xlim = c(10,25))
lines(temp,Tmult_int)

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

intake2 <- intake/2
E_surplus2 <- intake2-cost

lines(temp,Tmult_cost)
lines(temp,Tmult_int)
lines(temp,Tmult_int2)

plot(temp,intake, type = "l", ylim = c(-30,30))
lines(temp,cost)
plot(temp,E_surplus, type = "l")
lines(temp,E_surplus2)

plot(temp,byssus, type = "l")

plot(E_surplus,byssus,type = 'l')
lines(E_surplus2,byssus,type = 'l')
# This gives us the OPPOSITE prediction... 

plot(temp,byssus,type = 'l')
lines(temp,byssus,type = 'l')

