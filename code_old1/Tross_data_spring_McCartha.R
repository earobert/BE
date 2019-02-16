# Fitting the tross model to data, and estimating the cost of byssal thread production
# Checks the statistical significance of the length to tissue weight relationship (p=.01)




rm(list=ls())
setwd("~/BE/BE_2018_12_3/Datasets")
df <- read.csv("Spring_mussel_byssus_energetics_simple.csv", stringsAsFactors = FALSE)
df$treatment <- as.factor(df$treatment)

par(mfrow = c(1,2))
# Estimate initial mass ####
plot(df$length_final_QC^3,df$all_tissue_wt, col = df$treatment)
plot(df$length_final_QC,df$all_tissue_wt, col = df$treatment)
plot(df$length_initial,df$length_final, col = df$treatment, ylim = c(0,31), xlim = c(0,31))
lines(x=c(22,31),y=c(22,31))
points(df$length_initial_QC,df$length_final_QC, col = df$treatment, pch = 20)

legend(x = 'topleft', legend = unique(df$treatment), pch =1, col=unique(df$treatment))

subset <- df[df$treatment=="never",]
plot(subset$length_final_QC,subset$all_tissue_wt)
plot(subset$length_final_QC^3,subset$all_tissue_wt)
legend(x = 'topleft', legend = unique(subset$treatment), pch =1, col=unique(subset$treatment))

v.l <- subset$length_final^3/1000

# one model (#2) uses length cubed, 
# and the other (#1) just looks at length and mass
mod.1 <- lm(subset$all_tissue_wt~subset$length_final)
mod.2 <- lm(subset$all_tissue_wt~v.l)
anova(mod.1)
anova(mod.2)

summary(mod.1)
summary(mod.2)

coef(mod.1)
coef(mod.2)

# Plot the two fits
plot(subset$length_final_QC,subset$all_tissue_wt)
abline(mod.1)

plot(df$length_final^3/1000,df$all_tissue_wt, col = df$treatment, xlim = c(10,30), ylim = c(.05,0.2))
points(subset$length_final_QC^3/1000,subset$all_tissue_wt, pch = 20, col = "blue")
abline(mod.2)

# Estimate initial mass and plot
est_initial_mass <- (df$length_initial^3/1000)*coef(mod.2)[2]+coef(mod.2)[1]
points(df$length_initial^3/1000,est_initial_mass, col = "blue")
points(df$length_final^3/1000,df$all_tissue_wt, col = df$treatment)

df <- cbind(df,est_initial_mass)

plot(df$est_initial_mass,df$all_tissue_wt)
plot(df$length_initial,df$est_initial_mass)

est_change_mass <- df$all_tissue_wt - df$est_initial_mass


df <- cbind(df,est_change_mass)

plot(df$length_initial,df$est_change_mass,col = df$treatment)

df_QC<- df[df$length_initial_QC>0,]

# Check effect of treatment on change in mass
boxplot(df_QC$est_change_mass~df_QC$treatment)
fit <- aov(df_QC$est_change_mass~df_QC$treatment)
summary(fit) # not significant p=0.3
TukeyHSD(fit)

# Check ANCOVA that takes into account initial length
fit <- aov(df_QC$est_change_mass~df_QC$treatment*df_QC$length_initial)
summary(fit) # Interaction term is 0.3, so will pool
fit <- aov(df_QC$est_change_mass~df_QC$treatment+df_QC$length_initial)
summary(fit) # Still not significant, p=0.263 for treatment effect with sig effect of length

# Change in length ####

change_length <- df$length_final-df$length_initial
df <- cbind(df,change_length)
plot(df$length_initial_QC,df$change_length, col = df$treatment)
legend(x = 'topright', legend = unique(df$treatment), pch =1, col=unique(df$treatment))

df_QC <- df[df$length_initial_QC>0,]
df_QC <- df_QC[df_QC$length_initial<31,]

boxplot(df_QC$change_length~df_QC$treatment, ylim = c(0,5))
fit <- aov(df_QC$change_length~df_QC$treatment)
summary(fit)
TukeyHSD(fit)

fit <- aov(df_QC$change_length~df_QC$treatment*df_QC$length_initial)
summary(fit)



# Change in length is not right.... 

# Should I be using condition index

# Using % change in length
perc_change_length <- change_length/df$length_initial
df <- cbind(df,perc_change_length)
boxplot(df$perc_change_length~df$treatment)
fit <- aov(df$perc_change_length~df$treatment)
summary(fit)
TukeyHSD(fit)

ss.manip <- df[df$treatment!="never",]
ss.manip$treatment <- as.character(ss.manip$treatment)
ss.manip$treatment <- as.factor(ss.manip$treatment)
boxplot(ss.manip$perc_change_length~ss.manip$treatment)
fit <- aov(ss.manip$perc_change_length~ss.manip$treatment)
summary(fit)
TukeyHSD(fit)

# Relative change in length ####
# What about relative change at that length compared to the never treatment
# Start with length
plot(df$length_initial_QC,df$change_length, col = df$treatment)
legend(x = 'topright', legend = unique(df$treatment), pch =1, col=unique(df$treatment))
subset <- df[df$treatment=="never",]
mod.3 <- lm(subset$change_length~subset$length_initial)
summary(mod.3)
abline(mod.3)
subset <- df[df$treatment=="weekly",]
mod.4 <- lm(subset$change_length~subset$length_initial)
abline(mod.4)
summary(mod.4)
subset <- df[df$treatment=="daily",]
mod.5 <- lm(subset$change_length~subset$length_initial)
abline(mod.5)
summary(mod.5)
abline(mod.3)
mod.3
# This is actually problematic... 
# not really sure if this is a strong direction to go, 
# seems to be directed by a few points not a hypothesis

# Next steps ####
# Next step will be to make the tross model and change the cost value. 

