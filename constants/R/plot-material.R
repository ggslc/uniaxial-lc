##' @alias fortify.material, ggplot.material
##' @param model material data.frame
##' @param data irrelevant here
##' @param ... additional ggplot2 commands
##' @return ggplot object
##' @example 
##' library(constants)
##' data(AuJC)
##' material.predict(AuJC)
##' ggplot(material.predict(AuJC)) +
##'   geom_path( aes(wavelength, value, colour=variable) )

fortify.material <- function (model, data, ...) 
{
  require(plyr)
  dwide <- with(model, data.frame(wavelength, epsilon.real=Re(epsilon),
                                  epsilon.imag=Im(epsilon)))
  melt(dwide, id="wavelength")
}

ggplot.material <- function (data = NULL, ...) {
  require(ggplot2)
  ggplot(fortify(data), ...) 
}




compare.material <- function(material=AuJC, range=c(0.4, 1), ...){

ggplot(subset(normalise.material(material), wavelength < range[2] & wavelength > range[1])) +
  geom_path(data=material.predict(material, range=range), aes(wavelength, value, colour=variable))+
  geom_point(aes(wavelength, value, colour=variable))

  }

