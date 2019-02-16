###################################################################################
#Scope for Growth model estimation of a and cost_byssus
#Created by Emily A. Roberts earobert@uw.edu
#Starting Fall 2018
#Latest modification: 3 November 2018 for the cost of byssus manuscript
###################################################################################

setwd("~/BE/BE")

set.seed(1)

source("stackpolyTB2.r")
source("random.status.pauly.2007.R")
source("lognorm.catch.R")
source("st.autocorrel.catch.R")
source("simulated.plot2.R")

pdf("Fig 1 v1.pdf",width=8,height=6)
simulated.plot2(autocorrel=0.5, nsims=2000)
dev.off()
