# Regression comparison
# Start with one linear model: y = f*x1 + h * x2 + 1 * x3
# lm 1
# eiv 1
# SIMEX

# Set parameters




# Need to change these this practice dataset
p.practice <- data.frame( 
  names = c("a","b","d","e","C","del"),
  val = c(14,10.3,0.69,1,0.000182,.333),
  se = c(7,4.8,.01,.0001,.01,.006),
  season = rep("Autumn", times = length(names))
)
p.practice

p.Aut <- data.frame(
  names = c("a","b","d","e","C","del"),
  val = c(14,10.3,0.69,1,0.182,.333),
  se = c(7,4.8,.01,.0001,.01,.006),
  season = rep("Autumn", times = length(names))
)
p.Aut
p.Spr <- data.frame(
  names = c("a","b","d","e","C","del"),
  val = c(45,34.4,0.69,1,0.182,.295),
  se = c(14,9.2,.01,.0001,.01,.004),
  season = rep("Spring", times = length(names))
)
p.Spr
p.Aut <- data.frame(
  names = c("a","b","d","e","C","del"),
  val = c(14,10.3,0.69,1,0.182,.333),
  se = c(1,1,.01,.0001,.001,.001),
  season = rep("Autumn", times = length(names))
)
p.Aut

p <- p.practice
a <- p$val[1]
b <- p$val[2]
d <- p$val[3]
e <- p$val[4]
C <- p$val[5]
del <- p$val[6]

a.se <- p$se[1]
b.se <- p$se[2]
d.se <- p$se[3]
e.se <- p$se[4]
C.se <- p$se[5]
del.se <- p$se[6]

# Import dataset

# lm 1
season="Autumn"
setwd("~/BE/BE/Datasets")
df <- read.csv(file="Spring_Fall.csv", stringsAsFactors = FALSE)
df$treatment <- as.factor(df$treatment)
df$season <- as.factor(df$season)
df$treatment <- ordered(df$treatment, levels = c("never", "weekly","daily"))
df <- df[df$season==season,]
df <- df[!is.na(df$len_final_QC)&&!is.na(df$len_init_QC),]
df.subset <- df[df$treatment=="never"|df$treatment=="daily",]
df <- df.subset

growth_obs <- ((del*(df$len_final_QC/10))^3-(del*(df$len_init_QC/10))^3)/df$expt_length
mass_init <- (del*(df$len_init_QC/10))^3
mass_final_real <- df$total_wt_dry*3.98
Nth <- df$thread_count_QC

length(growth_obs)
length(mass_init)
length(Nth)

C #g/J
C_mg <- C/1000

growth_tissue <- growth_obs*1000
input <- C_mg*a*(mass_init*1000)^d #to make sense that ^2/3 of the mass is less than ^1 of the mass, this has to be in mg
byss <- -C_mg*Nth
non_byss <- -b*C_mg*(mass_init*1000)^e



growth_tissue <- growth_obs
input <- C*a*(mass_init)^d #to make sense that ^2/3 of the mass is less than ^1 of the mass, this has to be in mg
byss <- -C*(Nth/100)
non_byss <- -b*C*(mass_init)^e

dat <- data.frame(
  input = input,
  byss = byss,
  non_byss = non_byss
)

input.sd <- 1
byss.sd <- 1
non_byss.sd <- 1


mod2 <- lm(growth_tissue~ input + byss + offset(non_byss) + 0)
summary(mod2)
sum_mod2 <- summary(mod2)

plot(input, growth_tissue)
abline()
plot(-h*byss, growth_tissue, col = log(mass_init*10+2), pch = 20)
plot(-h*byss, growth_tissue, col = df$cage, pch = 20)

plot(growth_tissue, -h*byss, col = df$treatment, pch = 20,
     xlab = "Tissue growth",
     ylab = "Energy to byssus")

plot(-non_byss, f*input-h*byss, col = df$treatment, pch = 20,
     xlab = "Energy needed for maintenance",
     ylab = "Energy available after byssal threads",
     xlim = c(.0005,.01),
     ylim = c(.0005,.01)
     )

plot(growth_tissue, -h*byss, col = df$cage, pch = 20,
     xlab = "Tissue growth",
     ylab = "Energy to byssus")

plot(f*input-h*byss, -non_byss, col = df$cage, pch = 20,
     xlab = "Energy available after byssal threads",
     ylab = "Energy needed for maintenance")

plot(mass_init, f*input-non_byss+h*byss, col = df$treatment, pch = 20)

plot(mass_init, growth_tissue, col = df$treatment, pch = 20)


plot(mass_final_real, growth_tissue, col = df$treatment, pch = 20)

plot(df$len_init_QC, growth_tissue, col = df$treatment, pch = 20)


plot(df$gonad_wt_dry/df$somatic_wt_dry, growth_tissue, col = df$treatment, pch = 20)

plot(df$gonad_wt_dry, growth_tissue, col = df$treatment, pch = 20)

plot(df$len_init_QC, df$gonad_wt_dry, col = df$treatment, pch = 20)
plot(df$len_init_QC, df$somatic_wt_dry, col = df$treatment, pch = 20)
plot(df$somatic_wt_dry)

glm(mass_init)

plot(mass_init, growth_tissue, col = df$cage, pch = 20)



plot(f*input, -non_byss, col = df$treatment, pch = 20)


plot(mass_init, mass_init, col = log(mass_init*10))
plot(byss + non_byss, growth_tissue)
plot(input + byss + non_byss, growth_tissue)
plot(input + byss + non_byss, growth_tissue)

f <- sum_mod2$coefficients[1]
h <- sum_mod2$coefficients[2]

plot(input*f + byss*h + non_byss, growth_tissue,
     ylim = c(.0,.01),
     xlim = c(.0,.01)
     )



mod.ancova <- glm(growth_tissue~df$treatment*byss) #but then cage matters too

mod.ancova <- glm(growth_tissue~df$treatment*byss*mass_init) #but then cage matters too
summary(mod.ancova)

plot(growth_tissue~byss*1000, col = df$treatment)

mod1 <- lm(growth_tissue ~ input + byss + offset(non_byss), x=T, y=T, data = dat)

plot(input, growth_tissue)
plot(byss, growth_tissue)
plot(non_byss, growth_tissue)

summary(mod1)
sum_mod1 <- summary(mod1)
# pred_growth_2 <-sum_mod1$coefficients[1]+sum_mod1$coefficients[2]*input+sum_mod1$coefficients[3]*byss+
#   sum_mod1$coefficients[4]*non_byss

pred_growth_2 <-sum_mod1$coefficients[1]+sum_mod1$coefficients[2]*input+sum_mod1$coefficients[3]*byss + non_byss
plot(pred_growth_2~growth_tissue,
     ylim = c(.0,.05),
     xlim = c(.0,.05),
     col=df$treatment
)
x <- seq(from = -.4, to=.4, by=.1)
lines(x,x, lty = 2)

#SIMEX
naive <- mod1
mod.sim <- simex(naive,
                 measurement.error = c(input.sd),
                 SIMEXvariable = c("input"))
mod.sim
n<-nrow(dat)
mat <- data.frame(input = rep(input.sd,times = n), 
                  byss = rep(input.sd,times = n), 
                  non_byss= rep(input.sd,times = n))
nrow(mat)
mod.sim <- simex(naive,
                 measurement.error = c(1,1,1), # or mat?
                 SIMEXvariable = c("input", "byss", "non-byss"))
mod.sim

nrow(matrix(c(input.sd, byss.sd, non_byss.sd)))
length(c("input", "byss", "non-byss"))
nrow(dat)

