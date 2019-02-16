## ============================================================================
## FISH 454 Ecological Modeling - Winter 2014
## Lab #3
## January 2014
## ============================================================================

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
# Policy Variables
Yp <- 1.35      # annual expenditures for prevention (millions USD)
N.target <- 2   # target level after removal


# Main Parameters
a <- -0.824     # Ln(annual average invasion rate)
b <- -0.344     # effect of each million USD spent on invasion rate
alpha <- 0.75   # maximum reproductive rate (year)
ar <- 10        # cost to remove last snake
br <- 0.83      # power function, describing how population density effects 
                # cost of removal
N.max <- 7500000 # carrying capacity
N.min <- 2      # cost per year for each snake
damage <- 122.31 # cost per year for each snake
M <- 0.1        # annual mortality rate

# Other Parameters
N.start <- 0    # initial population size
lambda <- exp(a + b*Yp)           # annual invation rate
beta <- (alpha - M) / (M * N.max) # density dependence in reproduction

# -----------------------------------------------------------------------------
# Calculating annual values
# -----------------------------------------------------------------------------
# this model is set up very similarly to that in the excel file for Lab3
# it is not the most efficient way to write the code - but hopefully will be
# more intuitive when looking at all the outputs

n.years <- 100 # we will run this model forward 100 years
pop.matrix <- matrix(ncol = 7, nrow = n.years)
colnames(pop.matrix) <- c("N(t)", "Invaded", "Births", "N(t+1)", "Removals",
                          "Removal cost", "Damages")

# now we need to fill in the top row of the matrix, and that will be used to
# calculate all the subsequent rows
pop.matrix[1,1] <- N.start # starting population at N.start
pop.matrix[, 2] <- lambda   # every year has the same lambda, so the column can 
                           # be filled in with this value
births.y1 <- max(0, alpha * (pop.matrix[1,1] - N.min) / 
                   (1 + beta * (pop.matrix[1,1] - N.min))) # births the first year
pop.matrix[1, 3] <- births.y1
pop.matrix[1, 4] <- N.start*(1-M) + lambda + births.y1  # the new population size
pop.matrix[1, 5] <- max(0, (pop.matrix[1, 4] - N.target)) # removals
pop.matrix[1, 6] <- ar / N.target^br*pop.matrix[1,5] # removal cost 
pop.matrix[1, 7] <- damage * pop.matrix[1,1]/1000000

head(pop.matrix,3) # you can look to see what the first few rows of the matrix 
                   # look like


# now in order to fill in all the rows of the population matrix, we need to 
# make a for loop to calculate each subsequent year

for (i in 1:(n.years-1)) {
  # first we calculate N(t+1) so the new population size based on 
  # N(t) + invasions + births - deaths
  pop.matrix[i+1, 1] <- pop.matrix[i, 4] - pop.matrix[i,5]
  pop.matrix[i+1, 3] <- max(0, alpha * (pop.matrix[i+1,1] - N.min) / 
                              (1 + beta * (pop.matrix[i+1,1] - N.min)))
  pop.matrix[i+1, 4] <- pop.matrix[i+1,1] * (1-M) + pop.matrix[i+1,2] + pop.matrix[i+1,3]
  pop.matrix[i+1, 5] <- max(0, (pop.matrix[i+1, 4] - N.target)) 
  pop.matrix[i+1, 6] <- ar / N.target^br*pop.matrix[i+1,5] 
  pop.matrix[i+1, 7] <- damage * pop.matrix[i+1,1]/1000000
  
}

# -----------------------------------------------------------------------------
# Calculating costs
# -----------------------------------------------------------------------------
total.prevention <- Yp*100
total.removal <- sum(pop.matrix[, 6])
total.damages <- sum(pop.matrix[, 7])

total.cost <- sum(total.prevention, total.removal, total.damages)/100
# print output
print(paste('Total Cost=$',round(total.cost,3)," Million USD",sep=""))

