\name{conorect.lc}
\alias{conorect.lc}
\title{Conoscopy of liquid crystal cell: image calculation}
\description{Computes a complete rectangular conoscopy image
  for given minumum and maximum angles of incindence.Multiple (l)
  wavelengths can be computed simultaneously for the same stack.
  See conorect.stack for more details
  
  
}
\usage{conorect.lc(lambda, lc, eolc, eelc, db,eb,dt,et,
nx,ny,pmin, pmax, az0 = pi/4, i = 1, method = "ko")}
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
\item{nx}{number of horizontal pixels}
\item{ny}{number of vertical pixels}
\item{pmin}{minuimum angle of incidence}
\item{pmax}{maximum angle of incidence}
\item{az0}{scalar, angle between the polarizers transmission axis and
  the x-axis}
\item{i}{index of tilt,twist matrices to extract - ie director profile
  at seqeunce index i.Defaults to 1}
\item{method}{which underlying method to use. Currently "ko" for Berreman 4x4 treatment, and "lien" for Jones treatment}
}
\seealso{optics.lc, conopix.lc, conorect.stack}
\keyword{}