##' @alias normalise.material, spline.material, predict.material, permittivity
##' @param material
##' @param ... control parameters to spline
##' @param sp smoothing spline
##' @param range range of the predicted values, in um
##' @param n length of the predicted values
##' @return normalise and predict return a data.frame; spline returns a list of two smooth.spline.
##' @example 
##' library(constants)
##' data(AuJC)
##' material.predict(AuJC)

normalise.material <- function(material=AuJC){
  structure(with(material, data.frame(wavelength = lambda*1e6,
                                     epsilon = complex(real=epsr, imag=epsi))),
            class=c("normalised", "material", "data.frame"))
}

spline.material <- function(material, ...){
  if(!"normalised" %in% class(material))
    material <- normalise.material(material)
  
  list(real=smooth.spline(material$wavelength, Re(material$epsilon), ...),
       imag=smooth.spline(material$wavelength, Im(material$epsilon), ...))
  
}

predict.material <- function(material, sp=spline.material(material),
                             range=c(0.4, 0.9), n=NULL, wavelength=NULL, ...){

  if(is.null(range))
    range <- range(sp$real$data$x)

  if(is.null(wavelength)){ 
    wavelength <- if(is.null(n)) { # use the original points
      sp$real$data$x[sp$real$data$x > range[1] & sp$real$data$x < range[2]]
      
    } else { 
      seq(range[1], range[2], length=n)
    }
  }
  
  structure(with(sp,
                 data.frame(wavelength=wavelength,
                            epsilon=complex(real=predict(real, wavelength, ...)$y,
                              imag=predict(imag, wavelength, ...)$y))),
            class=c("smoothed", "normalised", "material", "data.frame"))
  
}

permittivity <- function(material=AuJC, wavelength=0.6328, ...){
  material.predict(material, range=rep(wavelength, 2), n=1)
}

# dumb naming scheme...
material.predict <- predict.material
