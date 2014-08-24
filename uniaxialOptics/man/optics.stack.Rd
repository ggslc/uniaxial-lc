\name{optics.stack}
\alias{optics.stack}
\title{Optics of multilayer, anisotropic systems}
\description{Calculate the complex and intensity reflection and transmission coefficients
through a multilayer, anisotropic system.}
\usage{optics.stack(lambda,azimuth,angle,depth,tilt,twist,eo,ee, method="ko" )}
\arguments{
\item{lambda}{vector of wavelengths, in metres, length = l}
\item{azimuth}{vector of angles between the y axis of the system and the plane of incidence, length = m }
\item{angle}{vector of angles of incidence, length = m }
\item{depth}{vector of layer depths, length = n. First and last entries should be zero, and are interpreted
as semi-inifinite}
\item{tilt}{vector of angles between optic axis and z axis,  length = n}
\item{twist}{vector of angles between optic axis and y axis, length = n}
\item{eo}{ordinary relative permittivity. Either a complev vector, length = n, or an n x l complex matrix if dispersion is to be included}
\item{ee}{extraordinary relative permittivity. See eo.}
\item{method}{which underlying method to use. Currently "ko" for Berreman 4x4 treatment, and "lien" for Jones treatment}
}
\seealso{optics.lc}
\examples{

#shows the difference bewteen Ko/Berreman and Lien/Jones
#methods. Ko/Berreman includes thin film inteference.

#simple stack, one birefringent layer between ito coated glass
depth <- c(0,.03,3,.03,0) * 1e-6

glass <- 1.52^2 + 0.0i
ito <- 3.8 + 0.08i

eo <- c(glass,ito,glass,ito,glass)
ee <- c(glass,ito,glass + 0.04 + 0i,ito,glass)
tilt <- c(0,0,pi/2,0,0)
twist <- rep(0,length(tilt))

angle <- seq(0,pi/3,l=100)
azimuth <- rep(pi/4,length(angle))

rko <- optics.stack(lambda=632.8e-9,azimuth,angle,depth,tilt,twist,eo,ee, method="ko")
rlien <- optics.stack(lambda=632.8e-9,azimuth,angle,depth,tilt,twist,eo,ee, method="lien")

plot(rko$angle * 180/pi, rko$Tps, type='l',
xlab='Angle of incidence (degrees)', ylab=expression(T[ps]))
lines(rlien$angle * 180/pi, rlien$Tps,lty=2)

legend(x="bottomleft",lty=c(1,2),leg=c("ko","lien"))


}
\keyword{}
