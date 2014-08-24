
library(minimal)
data(AuJC)

gold <- predict.material(AuJC, range=c(0.35, 1.2), n=300)
head(gold)
raw.gold <- subset(normalise.material(AuJC), wavelength < 1 & wavelength > 0.4)

library(ggplot2)

compare.material()
