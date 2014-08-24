\name{conorect.stack}
\alias{conorect.stack}
\title{Conoscopy of uniaxial multilayer stacks: image calculation}
\description{Computes a complete rectangular conoscopy image
  for given minumum and maximum angles of incindence

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

  A list \code{l} is returned, including an \code{nx x ny x nl} array
  \code{l$crossed}, of which the element \code{l$crossed[,,l]}
  is a matrix if pixel values. Pixels outside. Pixels outside
  the annulus defined by \code{pmin,pmax} are set to 0.

  Be warned: it is easy to ask for a lot of pixels, in which
  case the function will take a long time to evaluate. When
  trying out new code using it, keep the resolution down.
  
  
}
\usage{conorect.stack(lambda, depth, tilt, twist, eo, ee,
nx, ny, pmin, pmax, az0 = pi/4,method = "ko")}
\arguments{
\item{lambda}{vector of wavelengths, in metres, length = l}
\item{depth}{vector of layer depths, length n }
\item{tilt}{vector of tilt angles between optic axis and z-axis, length n }
\item{twist}{vector of twist angles between optic axis and y-axis, length n }
\item{eo}{ordinary relative permittivity. Either a complev vector, length = n, or an n x l complex matrix if dispersion is to be included}
\item{ee}{extraordinary relative permittivity. See eo.}
\item{nx}{number of horizontal pixels}
\item{ny}{number of vertical pixels}
\item{pmin}{minuimum angle of incidence}
\item{pmax}{maximum angle of incidence}
\item{az0}{scalar, angle between the polarizers transmission axis and
  the x-axis}
\item{method}{which underlying method to use. Currently "ko" for
  Berreman 4x4 treatment, and "lien" for Jones treatment}
}

\seealso{conopix.stack, conorect.lc}
\examples{

#compute a conoscopy figure of a uniaxial, homeotropic
#layer between ito coated glass

depth <- c(0,.03,15,.03,0) * 1e-6
glass <- 1.52^2 + 0.0i
ito <- 3.8 + 0.08i
eo <- c(glass,ito,glass,ito,glass)
ee <- c(glass,ito,(sqrt(Re(glass)) + 0.2)^2 + 0i,ito,glass)
tilt <- c(0,0,pi/20,0,0)
twist <- rep(0,length(tilt))

angle <- seq(0,pi/3,l=100)
azimuth <- rep(pi/4,length(angle))
lambda <- 633e-9

numap <- .55/1.52 # numerical aperture, in glass, of ~50X lens
n <- 96 #horizontal and vertical resolution

rko <- conorect.stack(lambda,depth,tilt,twist,eo,ee,n,n, 0, numap, method="ko")


#define image plotting function 
imp <- function(r,l,bright,...){
  image(t(r$crossed[,,l]),
        col=colorRampPalette(c("black",bright),bias=1/2.2)(64),
        axes=FALSE,...)
}
imp(rko,1,bright="red",main=paste(lambda[1]*1e+9," nm, ko"))






}
\keyword{}
