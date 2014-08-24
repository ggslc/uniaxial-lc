\name{optics.lc}
\alias{optics.lc}
\title{Optics of liquid crystal cells}
\description{Calculate the complex and intensity reflection and transmission coefficients
  through a liquid crystal cell. The cell comprises a liquid crytsal (LC) film
sandwiched between two isotropic multi-layers, which include the
substrates, electrodes, aligners, and so on}
\usage{optics.lc(lambda,azimuth,angle,lc,eolc,eelc,db,eb,dt,et, method="ko" )}
\arguments{
\item{lambda}{vector of wavelengths, in metres, length = l}
\item{azimuth}{vector of angles between the y axis of the system and the plane of incidence, length = m }
\item{angle}{vector of angles of incidence, length = m }
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
\item{method}{which underlying method to use. Currently "ko"
  for Berreman 4x4 treatment, and "lien" for Jones treatment}
}
\seealso{optics.stack}
\examples{

#shows the difference bewteen Ko/Berreman and Lien/Jones
#methods. Ko/Berreman includes thin film inteference.

#5CB hybrid-aligned LC film between ITO coated glass
glass <- 1.52^2 + 0.0i
ito <- 3.8 + 0.08i
nz <- 48
d <- 3e-6
eolc <- 1.52^2 + 0i
eelc <- 1.68^2 + 0i
di <- c(0,30e-9)
ei <- c(glass,ito)

lc <- list(depthz=c(0.5,rep(1,nz-2),0.5) * d/nz,
           tilt=matrix(seq(0,pi/2,l=nz),ncol=1),
           twist=matrix(seq(0,0,l=nz),ncol=1))

angle <- seq(0,pi/3,l=100)
azimuth <- rep(pi/4,length(angle))

rko <- optics.lc(lambda=632.8e-9,azimuth,angle,lc,eolc,eelc,
                 di,ei,rev(di),rev(ei),method="ko")
rlien <- optics.lc(lambda=632.8e-9,azimuth,angle,lc,eolc,eelc,
                 di,ei,rev(di),rev(ei),method="lien")

plot(rko$angle * 180/pi, rko$Tpp, type='l',
xlab='Angle of incidence (degrees)', ylab=expression(T[pp]),ylim=c(0,1))
lines(rlien$angle * 180/pi, rlien$Tpp,lty=2)

legend(x="bottomleft",lty=c(1,2),leg=c("ko","lien"))

}
\keyword{}
