# Deming - this works if just one x 

library(deming)

fit <- deming(aes ~ aas, data=arsenate, xstd=se.aas, ystd=se.aes)
print(fit)

model <- lm(mpg~wt+disp,data = mtcars)
summary(model)
summary(model)$r.square
plot(mpg~wt+disp,data = mtcars)

 fit <- deming(mpg ~ wt + disp, data=mtcars, xstd=se.aas, ystd=se.aes)
 print(fit)
 
 
 data(ferritin)
 temp <- ferritin[ferritin$period <4,]
 plot(temp$old.lot, temp$new.lot, type='n', log='xy',
      xlab="Old lot", ylab="New Lot")
 text(temp$old.lot, temp$new.lot, temp$period,
      col=temp$period)
 