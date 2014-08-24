\name{fineEps}
\alias{fineEps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Smoothing of permittivity values }
\description{
  Interpolates permittivity against wavelength values using cubic splines
}
\usage{
fineEps(lambda, eps, size = 100, deb = 3e-07, fin = 9e-07)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{ wavelength [m] }
  \item{eps}{ epsilon values (real vector) }
  \item{size}{ number of interpolated points }
  \item{deb}{ starting wavelength [m] }
  \item{fin}{ ending wavelength [m] }
}
\details{

}
\value{
list
  \item{lambdaFine }{}
  \item{epsFine }{}
}
\references{ }
\author{ baptiste }
\note{ 
}
\seealso{\code{\link{smooth.spline}} }
\examples{

# data(AuJChristy)
# attach(AuJChristy)
# 
# ls()
# 
# epsrf<-fineEps(lambda,epsr)
# epsif<-fineEps(lambda,epsi)
# par(mfrow=c(2,1))
# plot(lambda,epsr,type="p",pch=1)
# lines(epsrf$lambda.fine,epsrf$eps.fine)
# plot(lambda,epsi,type="p",pch=1)
# lines(epsif$lambda.fine,epsif$eps.fine)
}
