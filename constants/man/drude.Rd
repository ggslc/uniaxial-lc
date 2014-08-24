\name{drude}
\alias{drude}
\title{drude}
\usage{drude(lambda, omega, omega.p=1.3172e+16, gamma.p=1.2991e+14,
    basic=FALSE, surf=FALSE, interband=FALSE, eps.inf=1.54, A=1,
    nu.f=1.4 * 1e+06, leff=1, A.i=c(1.27, 1.1), phi.i=c(-pi/4, -pi/4),
    lambda.i=c(470, 325) * 1e-09, gamma.i=c(1900, 1060) * 1e-09)}
\value{...}
\arguments{\item{lambda}{}
\item{omega}{}
\item{omega.p}{}
\item{gamma.p}{}
\item{basic}{}
\item{surf}{}
\item{interband}{}
\item{eps.inf}{ionbackground}
\item{nu.f}{}
\item{leff}{surface scattering}}
\examples{library(constants)
data(AuJC)

gold.raw <- predict.material(AuJC, range=c(0.7, 1.2)) # original values for fit
gold <- predict.material(AuJC, range=c(0.3, 1.2), n=500) # interpolated values in wider range

model <- function(p, lambda=gold$wavelength){
epsilon <- drude(lambda=lambda*1e-6, # in meters
interband=FALSE, basic=TRUE, surf=FALSE,
omega.p = p[1]*1e+16, 
gamma.p = p[2]*1e+14,
eps.inf=p[3])

data.frame(wavelength=lambda, epsilon=epsilon)
}



objective <- function(p){

model <- model(p, lambda=gold.raw$wavelength)$epsilon
sos.real <- sum( (Re(model) - Re(gold.raw$epsilon))^2 ) / sum(Re(gold.raw$epsilon)^2)
sos.imag <- sum( (Im(model) - Im(gold.raw$epsilon))^2 ) / sum(Im(gold.raw$epsilon)^2)
sos <-  sos.real + sos.imag
sos
}


p0 <- c(1, 1, 1)
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
\dontrun{
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
}}
