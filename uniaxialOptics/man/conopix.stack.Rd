\name{conopix.stack}
\alias{conopix.stack}
\title{Conoscopy of uniaxial multilayer stacks: pixel calculation}
\description{Computes m pixels in the conoscopy (interference) figure
  of a multilayer stack. Multiple (l) wavelengths can be computed
  simultaneously for the same stack.

  A conoscopy figure is formed when a cone of linearly polarized
  light passes through a sample (which should be uniform across
  the whole beam spot) and then, following collimation,
  through a second, crossed polarizer onto a screen.
  Provided the beam spot is some distance from the focus,
  the cone can be built up from rays, each of which strikes
  a differnt point on the screen. The trajectory of these
  rays can be computed by matrix methods for inifinite
  plane waves, such as the Jones or Berreman methods.
  A coordinate sytem x,y is defined such that 
  \code{r = x^2 + y^2  = sin(angle of incidence)},
  \code{phi = tan(y/x) = }
  angle between plane of incidence and stack y-axis.

  A list l including l\$crossed, a (m) x (l) matrix
  of intensities, is returned.
}
\usage{conopix.stack(lambda, depth, tilt, twist, eo, ee,
x, y, az0 = pi/4,method = "ko")}
\arguments{
\item{lambda}{vector of wavelengths, in metres, length = l}
\item{depth}{vector of layer depths, length n }
\item{tilt}{vector of tilt angles between optic axis and z-axis, length n }
\item{twist}{vector of twist angles between optic axis and y-axis, length n }
\item{eo}{ordinary relative permittivity. Either a complev vector, length = n, or an n x l complex matrix if dispersion is to be included}
\item{ee}{extraordinary relative permittivity. See eo.}
\item{x}{vector, length m, of pixel x-coordinates to compute}
\item{y}{vector, length m, of pixel y-coordinates to compute}
\item{az0}{scalar, angle between the polarizers transmission axis and
  the x-axis}
\item{method}{which underlying method to use. Currently "ko" for Berreman 4x4 treatment, and "lien" for Jones treatment}
}
\seealso{optics.stack, conorect.stack, conopix.lc, conorect.lc}
\examples{

#compute a section along the x-axis of the conoscopy figure
#shows the difference bewteen Ko/Berreman and Lien/Jones
#methods. Ko/Berreman includes thin film inteference.
#simple stack, one birefringent, near homeotropic layer between
#ito coated glass

depth <- c(0,.03,15,.03,0) * 1e-6
glass <- 1.52^2 + 0.0i
ito <- 3.8 + 0.08i
eo <- c(glass,ito,glass,ito,glass)
ee <- c(glass,ito,(sqrt(Re(glass)) + 0.2)^2 + 0i,ito,glass)
tilt <- c(0,0,pi/20,0,0)
twist <- rep(0,length(tilt))

angle <- seq(0,pi/3,l=100)
azimuth <- rep(pi/4,length(angle))
lambda <- c(632.8e-9,400e-9)

numap <- .55/1.52 # numerical aperture, in glass, of ~50X lens
x <- seq(-numap,numap,l=100)
y <- rep(0.1,length(x))

cko <- conopix.stack(lambda,depth,tilt,twist,eo,ee, x, y,  method="ko")
#clien <- conopix.stack(lambda,depth,tilt,twist,eo,ee, x, y, method="lien")
#lien not working for now

#define matrix plotting function for convenience
mp <- function(x,y,add=FALSE,lty=1,...){
  col <- colorRampPalette(c("red","blue"),space="Lab")(ncol(y))
  if (add){
    matlines(x,y,type='l',lty=lty,col=col,...)
  } else {
    matplot(x,y,type='l',lty=lty,col=col,...)
  }
}

mp(x,cko$crossed,xlab='x', ylab=expression(Irradiance))


}
\keyword{}
