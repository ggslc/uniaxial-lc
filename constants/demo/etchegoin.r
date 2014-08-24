library(constants)
data(AuJC)

gold.raw <- predict.material(AuJC, range=c(0.3, 1)) # original values for fit
gold <- predict.material(AuJC, range=c(0.3, 1), n=100) # original values in wider range
 
model <- function(p, lambda=gold$wavelength){
  epsilon <- drude(lambda=lambda*1e-6,
                   interband=TRUE, basic=TRUE, surf=FALSE,
                   omega.p = p[1]*1e+16, 
                   gamma.p = p[2]*1e+14,
                   eps.inf=p[3], 
                   A.i = p[4:5], 
                   phi.i =c(-pi/4, -pi/4), 
                   lambda.i = p[6:7]*  1e-09, 
                   gamma.i = p[8:9]* 1e-09)

  data.frame(wavelength=lambda, epsilon=epsilon)
}

objective <- function(p){
  
	model <- model(p, lambda=gold.raw$wavelength)$epsilon
	sos.real <- sum( (Re(model) - Re(gold.raw$epsilon))^2 ) / sum(Re(gold.raw$epsilon)^2)
	sos.imag <- sum( (Im(model) - Im(gold.raw$epsilon))^2 ) / sum(Im(gold.raw$epsilon)^2)
	sos <-  sos.real + sos.imag
	sos
}


p0 <- c(1.32, 1.3, 1.54, 1.27, 1.1, 470, 325, 1900, 1060 )

guess <- model(p0)

res <- optim(p0, objective)
fit <- model(res$par)

splitComplex <- function(d){

  transform(d, real=Re(epsilon), imag=Im(epsilon))
}

comparison <- lapply(list(raw=gold.raw,
                          data=gold,
                          guess=guess,
                          fit=fit), splitComplex)

library(ggplot2)

m <- melt(comparison[-1], meas=c("real", "imag"))

str(m)

ggplot(m, aes(wavelength, value)) +
  facet_grid(variable~., scales="free_y")+
  geom_path(aes(colour=L1)) +
  geom_point(data=melt(comparison[1], meas=c("real", "imag"))) +
  ylab(expression(epsilon)) +
  xlab(expression(wavelength*" / "*mu*m)) +
  scale_colour_discrete("") +
  theme_bw()

