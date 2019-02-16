#_____________________________________________________________________----
#===============================================#
# Model 1, in the format of Model 2. Gamma is constrained, but cancels out
# MLE estimation of cost of byssus given growth and thread production####
# With separate estimates of a and the cost of thread production 
# BUT including a baseline value
# -use Ken's suggestion of having their be a baseline... Tried that here. Doesn't change anything.
# I could assume that 0 thread cutting was for that baseline value. 
#===============================================#
rm(list = ls())

#====



filled.legend <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                                                         length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes, ...) 
  {
    # modification of filled.contour by Carey McGilliard and Bridget Ferris
    # designed to just plot the legend
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
      stop("increasing 'x' and 'y' values expected")
    #  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    #  on.exit(par(par.orig))
    #  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    #layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
    #  par(las = las)
    #  mar <- mar.orig
    #  mar[4L] <- mar[2L]
    #  mar[2L] <- 1
    #  par(mar = mar)
    # plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
                yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
    if (missing(key.axes)) {
      if (axes) 
        axis(4)
    }
    else key.axes
    box()
  }









filled.contour3 <- function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...) 
  {
    # modification by Ian Taylor of the filled.contour function
    # to remove the key and facilitate overplotting with contour()
    # further modified by Carey McGilliard and Bridget Ferris
    # to allow multiple plots on one page
    
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
      stop("increasing 'x' and 'y' values expected")
    # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    # on.exit(par(par.orig))
    # w <- (3 + mar.orig[2]) * par("csi") * 2.54
    # par(las = las)
    # mar <- mar.orig
    plot.new()
    # par(mar=mar)
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
      stop("no proper 'z' matrix specified")
    if (!is.double(z)) 
      storage.mode(z) <- "double"
    .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                    col = col)
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(x, side = 1)
        Axis(y, side = 2)
      }
    }
    else plot.axes
    if (frame.plot) 
      box()
    if (missing(plot.title)) 
      title(...)
    else plot.title
    invisible()
}





#=====


monte_carlo_mod.1 <- function(shape_coeff_min, shape_coeff_max, d_min, d_max,
                              respiration_reference_J_per_day_min, respiration_reference_J_per_day_max, 
                              en_density_g_p_J_min, en_density_g_p_J_max, 
                              mass_ww_min, mass_ww_max, 
                              thread_num_min, thread_num_max,
                              params, season) {
  nsims<-10000
  food_scalar <- params[1] 
  cost_per_thread <- params[2]
  sigma <- params[3]
  season <- season
  output<-matrix(NA,nrow=nsims, ncol=6)
  
  baseline_byssus_multiplier <- 0.08
  sims <- 1
  for(sims in 1:nsims){
    shape_coeff <- (shape_coeff_max-shape_coeff_min)*runif(1)+shape_coeff_min
    d <- (d_max-d_min)*runif(1)+d_min
    respiration_reference_J_per_day <- (respiration_reference_J_per_day_max-respiration_reference_J_per_day_min)*runif(1)+respiration_reference_J_per_day_min
    en_density_g_p_J <- (en_density_g_p_J_max-en_density_g_p_J_min)*runif(1)+en_density_g_p_J_min
    mass_ww <- (mass_ww_max-mass_ww_min)*runif(1)+mass_ww_min
    thread_num <- (thread_num_max-thread_num_min)*runif(1)+thread_num_min
  
  # conversion factors
  mass_gWW_per_gDW <- 3.918 # From Summer 2015 collection, shape coeff estimation
  #shape_coeff <- .304 # From Summer 2015 collection, shape coeff estimation 

  # exponents ####
  #d <- 0.67 # intake; the model is very sensitive to this. (hope this is an exponent for g WW not mg WW)
  e <- 1 # the model is very sensitive to this (hope this is an exponent for g WW not mg WW)
  
  # calculation of b
  #respiration_reference_J_per_day <- 0.07*4.75*4.184*24 # Units: J / (day) from 0.07mlO2/hr, Fly and Hilbish 2013
  respiration_reference_gDW <- 0.24 # Fly and Hilbish 2013
  respiration_reference_gWW <- respiration_reference_gDW * mass_gWW_per_gDW # Using wet:dry conversion
  b_permass <- respiration_reference_J_per_day / (respiration_reference_gWW)^e # Note that e is just 1
  b_J_per_g_WW_per_day <- b_permass # Units: J / (gWW * day)
  
  # Temp response, other  ####
  Tmult_cost <- 1
  Tmult_int <- 1
  reprod_multiplier <- 1 # reproduction multiplier, brachi was 1, but if k is about half and half then this should be more like .5
  opt_size_reducer <- 1 # lets reproduction equal surplus early
  
  # other parameters ####
  #en_density_g_p_J <- .002 #energy_density
  size_maturity <- 10 # mg, size at maturity, 3000mg from tross del_M calc worksheet, somatic<full tissue weight
  
  # calculated scalers
  b <- b_J_per_g_WW_per_day*Tmult_cost #Wait the byssus multiplier is earlier here than I thought
  
  Wopt_measured_gDW <- 0.8 #gDW from sample collection
  Wopt_measured_gWW <- Wopt_measured_gDW * mass_gWW_per_gDW
  a_fromWopt <- (b*e)/((Wopt_measured_gWW)^(d-e)*d) # backwards calculation of a
  a_J_per_day_per_food_scalar <- a_fromWopt
  a <- a_J_per_day_per_food_scalar*food_scalar*Tmult_int
  

  # model####
  intake <- a*(mass_ww^d) #good, this is already in mass_ww... hope this is g not mg
  cost <- b*(mass_ww^e)*(1-baseline_byssus_multiplier)
  byssus_baseline <- b*(mass_ww^e)*baseline_byssus_multiplier
  byssus_induced <- b*(mass_ww^e)*thread_num*cost_per_thread
  reproduction <- 0
  model.predG_J <- intake-cost-byssus_baseline-byssus_induced-reproduction #predicts mussel growth in J
  model.predG_g <- model.predG_J*en_density_g_p_J  #predicts mussel growth in g DW or WW???
  
  output[sims,] <-c(shape_coeff, d, respiration_reference_J_per_day, en_density_g_p_J, mass_ww, model.predG_g)

  
  }
  output.df <- as.data.frame(output)
  colnames(output.df) <- c(
    "shape_coeff", 
    "d", 
    "respiration_reference_J_per_day", 
    "en_density_g_p_J", 
    "mass_ww", 
    "model.predG_g"
    )
  return(output.df)
}  

#Set up parameters ####

shape_coeff_min <- .28
shape_coeff_max <- .31
d_min <- .6
d_max <- .75
respiration_reference_J_per_day_min <- 20
respiration_reference_J_per_day_max <- 40
en_density_g_p_J_min <- 0.001
en_density_g_p_J_max <- 0.050


# Autumn ####
season <- "Autumn"

# import data ####
setwd("~/BE/BE/Datasets")
df <- read.csv(file="Spring_Fall.csv")
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
  
params <- c(1,.002,100)

mass_ww_min <- min(mass_ww)
mass_ww_max <- max(mass_ww)
thread_num_min <- min(thread_num)
thread_num_max <- max(thread_num)

output.A <- monte_carlo_mod.1(shape_coeff_min, shape_coeff_max, d_min, d_max,
                  respiration_reference_J_per_day_min, respiration_reference_J_per_day_max, 
                  en_density_g_p_J_min, en_density_g_p_J_max, mass_ww_min, mass_ww_max, thread_num_min, thread_num_max,
                  params, season)

# Spring ####
season <- "Spring"

# import data ####
setwd("~/BE/BE/Datasets")
df <- read.csv(file="Spring_Fall.csv")
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

params <- c(1,.002,100)

mass_ww_min <- min(mass_ww)
mass_ww_max <- max(mass_ww)
thread_num_min <- min(thread_num)
thread_num_max <- max(thread_num)

output.S <- monte_carlo_mod.1(shape_coeff_min, shape_coeff_max, d_min, d_max,
                              respiration_reference_J_per_day_min, respiration_reference_J_per_day_max, 
                              en_density_g_p_J_min, en_density_g_p_J_max, mass_ww_min, mass_ww_max, thread_num_min, thread_num_max,
                              params, season)

head(output.A)
# par(mfrow = c(4,2),mai = c(.5, .5, 0.2, 0.1), oma = c(.1,.1,1,.1))
# plot(output.A$model.predG_g~output.A[,"d"], pch = ".", ylab = "Growth", xlab = "d", ylim = c(0,4))
# title(main = "Autumn")
# plot(output.S$model.predG_g~output.S[,"d"], pch = ".", ylab = "Growth", xlab = "d", ylim = c(0,4))
# title(main = "Spring")
# plot(output.A$model.predG_g~output.A[,"respiration_reference_J_per_day"], pch = ".", ylab = "Growth", xlab = "Respiration (J/d)", ylim = c(0,4))
# plot(output.S$model.predG_g~output.S[,"respiration_reference_J_per_day"], pch = ".", ylab = "Growth", xlab = "Respiration (J/d)", ylim = c(0,4))
# plot(output.A$model.predG_g~output.A[,"en_density_g_p_J"], pch = ".", ylab = "Growth", xlab = "Energy density (g per J)", ylim = c(0,4))
# plot(output.S$model.predG_g~output.S[,"en_density_g_p_J"], pch = ".", ylab = "Growth", xlab = "Energy density (g per J)", ylim = c(0,4))
# plot(output.A$model.predG_g~output.A[,"mass_ww"], pch = ".", ylab = "Growth", xlab = "mass (g WW)", ylim = c(0,4))
# plot(output.S$model.predG_g~output.S[,"mass_ww"], pch = ".", ylab = "Growth", xlab = "mass (g WW)", ylim = c(0,4))

dev.off()

library(MASS)
par(mfrow = c(4,2),mai = c(.5, .5, 0.2, 0.1), oma = c(.1,.1,1,.1))
a <- output.A[,"mass_ww"]
b <- output.A$model.predG_g
fA1 <- kde2d(output.A[,"d"], output.A$model.predG_g, n = 100)
fA2 <- kde2d(output.A[,"respiration_reference_J_per_day"], output.A$model.predG_g, n = 100)
fA3 <- kde2d(output.A[,"en_density_g_p_J"], output.A$model.predG_g, n=100)
fA4 <- kde2d(output.A[,"mass_ww"], output.A$model.predG_g, n = 100)
fS1 <- kde2d(output.S[,"d"], output.S$model.predG_g, n = 100)
fS2 <- kde2d(output.S[,"respiration_reference_J_per_day"], output.S$model.predG_g, n = 100)
fS3 <- kde2d(output.S[,"en_density_g_p_J"], output.S$model.predG_g, n = 100)
fS4 <- kde2d(output.S[,"mass_ww"], output.S$model.predG_g, n = 100)
par(mfrow = c(4,2))
#Source the following functions (change the paths as necessary)
filled.contour3(fA1, ylim = c(0,1), xlab = "d")
filled.contour3(fS1, ylim = c(0,1), xlab = "d")
filled.contour3(fA2, ylim = c(0,1), xlab = "respiration (J/day)")
filled.contour3(fS2, ylim = c(0,1), xlab = "respiration (J/day)")
filled.contour3(fA3, ylim = c(0,1), xlab = "energy density (g/J)")
filled.contour3(fS3, ylim = c(0,1), xlab = "energy density (g/J)")
filled.contour3(fA4, ylim = c(0,1), xlab = "mass (gWW)")
filled.contour3(fS4, ylim = c(0,1), xlab = "mass (gWW)")

