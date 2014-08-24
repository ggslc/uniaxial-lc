\name{conopix.lc}
\alias{conopix.lc}
\title{Conoscopy of liquid crystal cell: pixel calculation}
\description{Computes m pixels in the conoscopy (interference) figure
  of a liquid crystal cell. Multiple (l) wavelengths can be computed
  simultaneously for the same stack. See conopix.stack for more details
}
\usage{conopix.lc(lambda, lc, eolc, eelc, db,eb,dt,et,x, y,
az0 = pi/4, i = 1, method = "ko")}
\arguments{
\item{lambda}{vector of wavelengths, in metres, length = l}
\item{lc}{a list containing at least depthz (sublayer depths) ,a vector of length n, and
  two n x m matrices, tilt and twist, which define m director profiles. The function
  lc\_switching\_sequence() in the EricksenLeslie package provdes just
  such a list. It is assumed that tilt[1] and tilt[m] are defined at the
  edge of the LC film, ie at the bottom and top of the end sublayers,
  whereas every other tilt/twist is defined in the middle of its sublayer.}
\item{eolc}{ordinary relative permittivity of the LC. Either a complex
  scalar, or a  complex vector of length l if dispersion is to be included}
\item{eelc}{extraordinary relative permittivity of the LC. See eo.}
\item{db}{vector, length nt, of depths for the top isotropic layers}
\item{eb}{relative permittivity of the top layers. Either a complex
  vector of length nt or a complex nt x l matrix if dispersion is to be
  included}
\item{dt}{vector, length nb, of depths for the bottom isotropic layers}
\item{et}{relative permittivity of the top layers. See et}
\item{x}{vector, length m, of pixel x-coordinates to compute}
\item{y}{vector, length m, of pixel y-coordinates to compute}
\item{az0}{scalar, angle between the polarizers transmission axis and
  the x-axis}
\item{i}{index of tilt,twist matrices to extract - ie director profile
  at seqeunce index i.Defaults to 1}
\item{method}{which underlying method to use. Currently "ko" for Berreman 4x4 treatment, and "lien" for Jones treatment}
}
\seealso{optics.lc, conopix.stack, conorect.lc}
\examples{


glass <- 1.52^2 + 0.0i
ito <- 3.8 + 0.08i
nz <- 48
d <- 5e-6
eolc <- 1.52^2 + 0i
eelc <- 1.68^2 + 0i
di <- c(0,30e-9)
ei <- c(glass,ito)

lc <- list(depthz=c(0.5,rep(1,nz-2),0.5) * d/nz,
           tilt=matrix(seq(0,pi/2,l=nz),ncol=1),
           twist=matrix(seq(0,0,l=nz),ncol=1))

lambda <- c(632.8e-9,500e-9,400e-9)

numap <- .55/1.52 # numerical aperture, in glass, of ~50X lens
x <- seq(-numap,numap,l=100)
y <- x


cko <- conopix.lc(lambda,lc,eolc,eelc, di,ei,rev(di),rev(ei), x, y,  method="ko")

clien <- conopix.lc(lambda,lc,eolc,eelc, di,ei,rev(di),rev(ei), x, y,  method="lien")

#define matrix plotting function for convenience
mp <- function(x,y,add=FALSE,lty=1,...){
  col <- rainbow(ncol(y))
  if (add){
    matlines(x,y,type='l',lty=lty,col=col,...)
  } else {
    matplot(x,y,type='l',lty=lty,col=col,...)
  }
}

mp(x,cko$crossed,xlab='x', ylab=expression(Irradiance),ylim=c(0,1),
   main ="x=y section through HAN cell conoscopy figure")

mp(x,clien$crossed,xlab='x', ylab=expression(Irradiance),add=TRUE,lty=2)

legend(x="bottomright",leg=lambda*1e+9,title=expression(lambda),lwd=1,
       col= rainbow(3),
       bg="#eeeeee")

legend(x="bottomleft",leg=c("ko","lien"),title="method",lty=c(1,2),
       bg="#eeeeee")

}
\keyword{}